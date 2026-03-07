#' Generate PRS-CS compatible LD panels from 33KG GAUSS reference files
#'
#' Builds per-chromosome LD blocks directly from GAUSS 33KG reference files
#' (`33kg_geno.gz`, `33kg_index.gz`) without converting to PLINK/BGEN/VCF.
#' Genotypes are extracted window-by-window, signed LD correlations are computed
#' within each window, optional shrinkage / eigenvalue flooring is applied to
#' enforce positive-definiteness, and output is written to HDF5 in a block layout.
#'
#' @param reference_geno_file Path to `33kg_geno.gz` (or gzipped file with same layout).
#' @param reference_index_file Path to `33kg_index.gz` with at least columns
#'   `rsid`, `chr`, `bp`, `a1`, `a2`.
#' @param chromosome Integer chromosome to process.
#' @param output_h5 Path to output HDF5 file for one chromosome.
#' @param window_size Window size in base pairs.
#' @param step_size Step size in base pairs for sliding windows.
#' @param n_reference_populations Number of population genotype strings per line
#'   in `reference_geno_file`.
#' @param population_indices Integer indices of populations to use. Defaults to all.
#' @param shrinkage Numeric shrinkage in `[0, 1)`, applied as
#'   `R <- (1 - shrinkage) * R + shrinkage * I`.
#' @param min_eigenvalue Eigenvalue floor for positive-definite correction.
#' @param maf_min Minimum minor allele frequency filter (after selected populations).
#' @param missing_max Maximum allowed missing genotype fraction per SNP.
#' @param max_windows Optional cap on number of windows (useful for examples/tests).
#' @param dry_run Logical; if `TRUE`, only returns per-window SNP counts and does
#'   not read genotypes or write HDF5.
#' @param panel_name Character identifier stored in HDF5 file attributes.
#' @param dosage_allele Which allele genotype dosage corresponds to (`"A1"` or `"A2"`).
#' @param verbose Logical; print progress messages.
#'
#' @return If `dry_run = FALSE`, invisibly returns `output_h5`.
#'   If `dry_run = TRUE`, returns a data frame with block boundaries and SNP counts.
#' @export
#'
#' @examples
#' \dontrun{
#' generate_prscs_ld_panel(
#'   reference_geno_file = "33kg_geno.gz",
#'   reference_index_file = "33kg_index.gz",
#'   chromosome = 22,
#'   output_h5 = "chr22_prscs_ld.h5",
#'   window_size = 2e6,
#'   step_size = 1e6,
#'   n_reference_populations = 29,
#'   max_windows = 3,
#'   shrinkage = 0.01
#' )
#' }
generate_prscs_ld_panel <- function(reference_geno_file,
                                    reference_index_file,
                                    chromosome,
                                    output_h5,
                                    window_size = 3e6,
                                    step_size = 1e6,
                                    n_reference_populations = 29,
                                    population_indices = NULL,
                                    shrinkage = 0,
                                    min_eigenvalue = 1e-8,
                                    maf_min = 0,
                                    missing_max = 1,
                                    max_windows = NULL,
                                    dry_run = FALSE,
                                    panel_name = "GAUSS_33KG",
                                    dosage_allele = c("A1", "A2"),
                                    verbose = TRUE) {
  dosage_allele <- match.arg(dosage_allele)

  .validate_prscs_args(
    reference_geno_file = reference_geno_file,
    reference_index_file = reference_index_file,
    chromosome = chromosome,
    output_h5 = output_h5,
    window_size = window_size,
    step_size = step_size,
    n_reference_populations = n_reference_populations,
    population_indices = population_indices,
    shrinkage = shrinkage,
    min_eigenvalue = min_eigenvalue,
    maf_min = maf_min,
    missing_max = missing_max,
    max_windows = max_windows,
    dry_run = dry_run,
    panel_name = panel_name,
    verbose = verbose
  )

  if (!dry_run && !requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required. Please install it with install.packages('hdf5r').", call. = FALSE)
  }

  if (is.null(population_indices)) {
    population_indices <- seq_len(n_reference_populations)
  }

  idx <- .read_gauss_index(reference_index_file)
  idx$row_id <- seq_len(nrow(idx))

  idx_chr <- idx[idx$chr == chromosome, , drop = FALSE]
  if (nrow(idx_chr) == 0) {
    stop(sprintf("No SNPs found for chromosome %s in index file: %s", chromosome, reference_index_file), call. = FALSE)
  }
  idx_chr <- idx_chr[order(idx_chr$bp), , drop = FALSE]

  windows <- .build_sliding_windows(idx_chr$bp, window_size = window_size, step_size = step_size)
  if (!is.null(max_windows)) {
    windows <- windows[seq_len(min(nrow(windows), max_windows)), , drop = FALSE]
  }

  window_counts <- .window_counts(idx_chr$bp, windows)
  if (dry_run) {
    if (verbose) {
      message("Dry run complete. Returning block boundaries and SNP counts.")
    }
    return(window_counts)
  }

  if (file.exists(output_h5)) {
    file.remove(output_h5)
  }

  h5 <- hdf5r::H5File$new(output_h5, mode = "w")
  on.exit(h5$close_all(), add = TRUE)

  h5$attr_open("panel_name")$write(panel_name)
  h5$attr_open("chr")$write(as.integer(chromosome))
  h5$attr_open("window_bp")$write(as.integer(window_size))
  h5$attr_open("step_bp")$write(as.integer(step_size))
  h5$attr_open("shrinkage_lambda")$write(as.numeric(shrinkage))
  h5$attr_open("eig_floor")$write(as.numeric(min_eigenvalue))
  h5$attr_open("maf_min")$write(as.numeric(maf_min))
  h5$attr_open("missing_max")$write(as.numeric(missing_max))
  h5$attr_open("n_blocks")$write(as.integer(nrow(windows)))
  h5$attr_open("n_samples")$write(as.integer(length(population_indices) * .infer_samples_per_population(reference_geno_file, n_reference_populations)))
  h5$attr_open("n_selected_pops")$write(as.integer(length(population_indices)))
  h5$attr_open("selected_pops")$write(as.integer(population_indices))
  h5$attr_open("dosage_allele")$write(dosage_allele)

  blocks_group <- h5$create_group("blocks")

  for (i in seq_len(nrow(windows))) {
    wstart <- windows$start_bp[i]
    wend <- windows$end_bp[i]
    snp_rows <- idx_chr$bp >= wstart & idx_chr$bp <= wend
    block_idx <- idx_chr[snp_rows, , drop = FALSE]

    if (verbose) {
      message(sprintf("Processing chr%s block %s/%s: [%s, %s], SNPs(before filter)=%s",
        chromosome, i, nrow(windows), wstart, wend, nrow(block_idx)
      ))
    }

    block_group <- blocks_group$create_group(sprintf("block_%s", i))
    block_group$attr_open("start_bp")$write(as.integer(wstart))
    block_group$attr_open("end_bp")$write(as.integer(wend))

    if (nrow(block_idx) == 0) {
      block_group$attr_open("m_snps")$write(as.integer(0))
      block_group[["ld"]] <- matrix(numeric(0), nrow = 0, ncol = 0)
      block_group[["rsid"]] <- character(0)
      block_group[["bp"]] <- integer(0)
      block_group[["a1"]] <- character(0)
      block_group[["a2"]] <- character(0)
      next
    }

    geno_mat <- .extract_window_genotypes(
      reference_geno_file = reference_geno_file,
      target_row_ids = block_idx$row_id,
      n_reference_populations = n_reference_populations,
      population_indices = population_indices
    )

    filt <- .filter_snps_by_maf_missing(geno_mat, maf_min = maf_min, missing_max = missing_max)
    geno_mat <- geno_mat[, filt$keep, drop = FALSE]
    block_idx <- block_idx[filt$keep, , drop = FALSE]

    if (ncol(geno_mat) == 0) {
      block_group$attr_open("m_snps")$write(as.integer(0))
      block_group[["ld"]] <- matrix(numeric(0), nrow = 0, ncol = 0)
      block_group[["rsid"]] <- character(0)
      block_group[["bp"]] <- integer(0)
      block_group[["a1"]] <- character(0)
      block_group[["a2"]] <- character(0)
      next
    }

    ld <- stats::cor(geno_mat, use = "pairwise.complete.obs")
    ld[is.na(ld)] <- 0
    ld <- .shrink_and_make_pd(ld, shrinkage = shrinkage, min_eigenvalue = min_eigenvalue)

    block_group$attr_open("m_snps")$write(as.integer(ncol(geno_mat)))
    block_group[["ld"]] <- ld
    block_group[["rsid"]] <- as.character(block_idx$rsid)
    block_group[["bp"]] <- as.integer(block_idx$bp)
    block_group[["a1"]] <- as.character(block_idx$a1)
    block_group[["a2"]] <- as.character(block_idx$a2)
  }

  invisible(output_h5)
}

.validate_prscs_args <- function(reference_geno_file,
                                 reference_index_file,
                                 chromosome,
                                 output_h5,
                                 window_size,
                                 step_size,
                                 n_reference_populations,
                                 population_indices,
                                 shrinkage,
                                 min_eigenvalue,
                                 maf_min,
                                 missing_max,
                                 max_windows,
                                 dry_run,
                                 panel_name,
                                 verbose) {
  if (!is.character(reference_geno_file) || length(reference_geno_file) != 1 || !nzchar(reference_geno_file)) {
    stop("'reference_geno_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_geno_file)) {
    stop(sprintf("Reference genotype file does not exist: %s", reference_geno_file), call. = FALSE)
  }

  if (!is.character(reference_index_file) || length(reference_index_file) != 1 || !nzchar(reference_index_file)) {
    stop("'reference_index_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_index_file)) {
    stop(sprintf("Reference index file does not exist: %s", reference_index_file), call. = FALSE)
  }

  if (!is.numeric(chromosome) || length(chromosome) != 1 || is.na(chromosome) || chromosome < 1) {
    stop("'chromosome' must be a positive integer.", call. = FALSE)
  }
  if (!is.character(output_h5) || length(output_h5) != 1 || !nzchar(output_h5)) {
    stop("'output_h5' must be a non-empty output path.", call. = FALSE)
  }

  if (!is.numeric(window_size) || length(window_size) != 1 || is.na(window_size) || window_size <= 0) {
    stop("'window_size' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(step_size) || length(step_size) != 1 || is.na(step_size) || step_size <= 0) {
    stop("'step_size' must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(n_reference_populations) || length(n_reference_populations) != 1 || is.na(n_reference_populations) || n_reference_populations < 1) {
    stop("'n_reference_populations' must be a positive integer.", call. = FALSE)
  }

  if (!is.null(population_indices)) {
    if (!is.numeric(population_indices) || any(is.na(population_indices)) || length(population_indices) == 0) {
      stop("'population_indices' must be a non-empty integer vector when provided.", call. = FALSE)
    }
    if (any(population_indices < 1 | population_indices > n_reference_populations)) {
      stop("'population_indices' contains indices outside valid range [1, n_reference_populations].", call. = FALSE)
    }
  }

  if (!is.numeric(shrinkage) || length(shrinkage) != 1 || is.na(shrinkage) || shrinkage < 0 || shrinkage >= 1) {
    stop("'shrinkage' must be in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(min_eigenvalue) || length(min_eigenvalue) != 1 || is.na(min_eigenvalue) || min_eigenvalue <= 0) {
    stop("'min_eigenvalue' must be positive.", call. = FALSE)
  }
  if (!is.numeric(maf_min) || length(maf_min) != 1 || is.na(maf_min) || maf_min < 0 || maf_min >= 0.5) {
    stop("'maf_min' must be in [0, 0.5).", call. = FALSE)
  }
  if (!is.numeric(missing_max) || length(missing_max) != 1 || is.na(missing_max) || missing_max < 0 || missing_max > 1) {
    stop("'missing_max' must be in [0, 1].", call. = FALSE)
  }

  if (!is.null(max_windows) && (!is.numeric(max_windows) || length(max_windows) != 1 || is.na(max_windows) || max_windows < 1)) {
    stop("'max_windows' must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.logical(dry_run) || length(dry_run) != 1 || is.na(dry_run)) {
    stop("'dry_run' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(panel_name) || length(panel_name) != 1 || !nzchar(panel_name)) {
    stop("'panel_name' must be a non-empty string.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
}

.read_gauss_index <- function(reference_index_file) {
  idx <- utils::read.table(
    gzfile(reference_index_file, open = "rt"),
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (ncol(idx) < 7) {
    stop("Reference index file must have at least 7 whitespace-delimited columns.", call. = FALSE)
  }
  idx <- idx[, 1:7, drop = FALSE]
  names(idx) <- c("rsid", "chr", "bp", "a1", "a2", "af1ref", "fpos")
  idx$chr <- as.integer(idx$chr)
  idx$bp <- as.integer(idx$bp)
  idx
}

.build_sliding_windows <- function(bp, window_size, step_size) {
  starts <- seq(from = min(bp), to = max(bp), by = step_size)
  data.frame(
    start_bp = as.integer(starts),
    end_bp = as.integer(starts + window_size - 1),
    stringsAsFactors = FALSE
  )
}

.window_counts <- function(bp, windows) {
  n <- nrow(windows)
  counts <- integer(n)
  for (i in seq_len(n)) {
    counts[i] <- sum(bp >= windows$start_bp[i] & bp <= windows$end_bp[i])
  }
  data.frame(block = seq_len(n), windows, n_snps = counts, stringsAsFactors = FALSE)
}

.extract_window_genotypes <- function(reference_geno_file,
                                      target_row_ids,
                                      n_reference_populations,
                                      population_indices) {
  target_row_ids <- sort(unique(as.integer(target_row_ids)))
  row_set <- new.env(parent = emptyenv())
  for (rid in target_row_ids) assign(as.character(rid), TRUE, envir = row_set)

  con <- gzfile(reference_geno_file, open = "rt")
  on.exit(close(con), add = TRUE)

  collected <- vector("list", length(target_row_ids))
  names(collected) <- as.character(target_row_ids)

  current_row <- 0L
  while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
    current_row <- current_row + 1L
    key <- as.character(current_row)
    if (!exists(key, envir = row_set, inherits = FALSE)) next

    tokens <- strsplit(line, "[[:space:]]+")[[1]]
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens) < n_reference_populations) {
      stop(sprintf("Genotype line %s has fewer than %s population genotype strings.", current_row, n_reference_populations), call. = FALSE)
    }

    geno_strings <- tokens[population_indices]
    values <- unlist(strsplit(paste0(geno_strings, collapse = ""), ""), use.names = FALSE)
    values[values == "9"] <- NA_character_
    collected[[key]] <- as.numeric(values)
  }

  missing <- vapply(collected, is.null, logical(1))
  if (any(missing)) {
    stop(sprintf("Failed to retrieve genotype rows from reference_geno_file. Missing row ids include: %s",
      paste(names(collected)[missing][seq_len(min(5, sum(missing)))], collapse = ", ")), call. = FALSE)
  }

  mat <- do.call(cbind, collected)
  storage.mode(mat) <- "double"
  mat
}

.filter_snps_by_maf_missing <- function(geno_mat, maf_min, missing_max) {
  if (ncol(geno_mat) == 0) {
    return(list(keep = logical(0)))
  }
  missing_frac <- colMeans(is.na(geno_mat))
  af <- colMeans(geno_mat, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  keep <- (missing_frac <= missing_max) & is.finite(maf) & (maf >= maf_min)
  list(keep = keep)
}

.infer_samples_per_population <- function(reference_geno_file, n_reference_populations) {
  con <- gzfile(reference_geno_file, open = "rt")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  if (length(first) == 0) {
    stop("Reference genotype file is empty.", call. = FALSE)
  }
  tokens <- strsplit(first, "[[:space:]]+")[[1]]
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) < n_reference_populations) {
    stop("Unable to infer sample size from genotype file: fewer tokens than n_reference_populations.", call. = FALSE)
  }
  nchar(tokens[[1]])
}

.shrink_and_make_pd <- function(ld, shrinkage, min_eigenvalue) {
  p <- ncol(ld)
  if (p == 0) return(ld)

  ld <- (ld + t(ld)) / 2
  if (shrinkage > 0) ld <- (1 - shrinkage) * ld + shrinkage * diag(p)

  eig <- eigen(ld, symmetric = TRUE)
  eig$values[eig$values < min_eigenvalue] <- min_eigenvalue
  ld_pd <- eig$vectors %*% diag(eig$values, nrow = p) %*% t(eig$vectors)

  ld_pd <- (ld_pd + t(ld_pd)) / 2
  d <- sqrt(diag(ld_pd))
  d[d == 0] <- 1
  ld_pd <- sweep(sweep(ld_pd, 1, d, "/"), 2, d, "/")
  diag(ld_pd) <- 1
  ld_pd <- (ld_pd + t(ld_pd)) / 2
  ld_pd
}
