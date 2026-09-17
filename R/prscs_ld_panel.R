#' Write one chromosome's PRS-CS LD blocks for a superpopulation
#'
#' Partitions a chromosome into contiguous, non-overlapping blocks and computes
#' a signed, ancestry-weighted LD correlation matrix for each block via
#' \code{\link{computeSuperPopLDWindow}}. Output matches the on-disk format read
#' by PRS-CS's \code{parse_ldblk}/\code{parse_ldblk_mix}: an HDF5 file whose root
#' contains groups \code{blk_1, blk_2, ..., blk_N} (contiguous, 1-indexed, no
#' wrapping group), each holding an \code{ldblk} dataset (the LD matrix) and a
#' \code{snplist} dataset (rsids).
#'
#' Blocks must be non-overlapping so each SNP belongs to exactly one block --
#' PRS-CS's mixing logic (\code{parse_ldblk_mix}) aligns two populations' panels
#' by block index and tracks a running cumulative SNP position across blocks,
#' which breaks if a SNP appears in more than one block. Block boundaries are
#' derived from the full, population-agnostic SNP `bp` range in
#' `reference_index_file` (not from any one panel's QC'd SNP set), so panels
#' generated with the same `reference_index_file` and `block_size` share
#' identical block boundaries regardless of which superpopulation they're for --
#' required for `blk_i` to mean the same genomic region across panels meant to
#' be mixed together.
#'
#' @param chr Integer chromosome to process.
#' @param superpopulation_label Superpopulation label (e.g. AFR/AMR/ASN/EUR/SAS).
#' @param panel_dir Path to the panel's output directory (e.g. `ldblk_1kg_eur33kg`).
#'   Must already exist.
#' @param reference_index_file Reference panel index file.
#' @param reference_data_file Reference panel genotype data file.
#' @param reference_pop_desc_file Reference panel population description file.
#' @param block_size Block size in base pairs. Blocks are contiguous and
#'   non-overlapping: `[start, start + block_size - 1]`.
#' @param maf_cutoff Minimum minor allele frequency filter. Defaults to 0.01.
#' @param missing_cutoff Maximum allowed missing genotype fraction per SNP. Defaults to 1.0.
#' @param shrinkage Numeric shrinkage in `[0, 1)`, applied as
#'   `R <- (1 - shrinkage) * R + shrinkage * I`.
#' @param min_eigenvalue Eigenvalue floor for positive-definite correction.
#' @param panel_name Character identifier stored in HDF5 file attributes.
#'   Defaults to `sprintf("GAUSS_33KG_%s", toupper(superpopulation_label))`.
#' @param max_blocks Optional cap on number of blocks (useful for examples/tests).
#' @param write_snpinfo_sidecar Logical; if `TRUE` (default), also writes a
#'   per-chromosome snpinfo sidecar file (`.snpinfo_part_chr{chr}.tsv`) in
#'   `panel_dir`, for later assembly by \code{\link{finalize_prscs_snpinfo}}.
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly, a data frame with columns `chr, rsid, bp, a1, a2, maf`
#'   covering every SNP written across this chromosome's blocks.
#' @export
write_prscs_ld_chr <- function(chr,
                               superpopulation_label,
                               panel_dir,
                               reference_index_file,
                               reference_data_file,
                               reference_pop_desc_file,
                               block_size = 250000,
                               maf_cutoff = 0.01,
                               missing_cutoff = 1.0,
                               shrinkage = 0,
                               min_eigenvalue = 1e-8,
                               panel_name = NULL,
                               max_blocks = NULL,
                               write_snpinfo_sidecar = TRUE,
                               verbose = TRUE) {
  if (is.null(panel_name)) {
    panel_name <- sprintf("GAUSS_33KG_%s", toupper(superpopulation_label))
  }

  .validate_prscs_chr_args(
    chr = chr,
    superpopulation_label = superpopulation_label,
    panel_dir = panel_dir,
    reference_index_file = reference_index_file,
    reference_data_file = reference_data_file,
    reference_pop_desc_file = reference_pop_desc_file,
    block_size = block_size,
    maf_cutoff = maf_cutoff,
    missing_cutoff = missing_cutoff,
    shrinkage = shrinkage,
    min_eigenvalue = min_eigenvalue,
    panel_name = panel_name,
    max_blocks = max_blocks,
    write_snpinfo_sidecar = write_snpinfo_sidecar,
    verbose = verbose
  )

  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required. Please install it with install.packages('hdf5r').", call. = FALSE)
  }

  idx <- .read_gauss_index(reference_index_file)
  idx_chr <- idx[idx$chr == chr, , drop = FALSE]
  if (nrow(idx_chr) == 0) {
    stop(sprintf("No SNPs found for chromosome %s in index file: %s", chr, reference_index_file), call. = FALSE)
  }
  idx_chr <- idx_chr[order(idx_chr$bp), , drop = FALSE]

  blocks <- .build_contiguous_blocks(idx_chr$bp, block_size = block_size)
  if (!is.null(max_blocks)) {
    blocks <- blocks[seq_len(min(nrow(blocks), max_blocks)), , drop = FALSE]
  }

  output_h5 <- file.path(panel_dir, sprintf("ldblk_1kg_chr%d.hdf5", chr))
  if (file.exists(output_h5)) {
    file.remove(output_h5)
  }

  h5 <- hdf5r::H5File$new(output_h5, mode = "w")
  on.exit(h5$close_all(), add = TRUE)

  .h5_write_attr(h5, "panel_name", panel_name)
  .h5_write_attr(h5, "chr", as.integer(chr))
  .h5_write_attr(h5, "block_bp", as.integer(block_size))
  .h5_write_attr(h5, "maf_cutoff", as.numeric(maf_cutoff))
  .h5_write_attr(h5, "missing_cutoff", as.numeric(missing_cutoff))
  .h5_write_attr(h5, "shrinkage_lambda", as.numeric(shrinkage))
  .h5_write_attr(h5, "eig_floor", as.numeric(min_eigenvalue))
  .h5_write_attr(h5, "superpopulation_label", toupper(superpopulation_label))
  .h5_write_attr(h5, "n_blocks", as.integer(nrow(blocks)))

  snpinfo_parts <- vector("list", nrow(blocks))

  for (i in seq_len(nrow(blocks))) {
    start_bp <- blocks$start_bp[i]
    end_bp <- blocks$end_bp[i]

    if (verbose) {
      message(sprintf("Processing chr%s block %s/%s: [%s, %s]",
        chr, i, nrow(blocks), start_bp, end_bp
      ))
    }

    blk_group <- h5$create_group(sprintf("blk_%d", i))
    ld_window_fn <- function(chr, start_bp, end_bp) {
      computeSuperPopLDWindow(
        chr = chr,
        start_bp = start_bp,
        end_bp = end_bp,
        superpopulation_label = toupper(superpopulation_label),
        reference_index_file = reference_index_file,
        reference_data_file = reference_data_file,
        reference_pop_desc_file = reference_pop_desc_file,
        maf_cutoff = maf_cutoff,
        missing_cutoff = missing_cutoff
      )
    }

    snpinfo_parts[[i]] <- .compute_and_write_ld_block(
      blk_group = blk_group,
      chr = chr,
      start_bp = start_bp,
      end_bp = end_bp,
      ld_window_fn = ld_window_fn,
      shrinkage = shrinkage,
      min_eigenvalue = min_eigenvalue
    )
  }

  chr_snpinfo <- do.call(rbind, snpinfo_parts)
  if (is.null(chr_snpinfo)) {
    chr_snpinfo <- .empty_snpinfo()
  }

  if (write_snpinfo_sidecar) {
    utils::write.table(
      chr_snpinfo,
      file = file.path(panel_dir, sprintf(".snpinfo_part_chr%d.tsv", chr)),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  invisible(chr_snpinfo)
}

#' Assemble the combined PRS-CS snpinfo file for a panel
#'
#' Combines per-chromosome SNP info (either supplied directly in `snpinfo_list`,
#' or read from on-disk sidecar files written by
#' \code{\link{write_prscs_ld_chr}}) into the single `snpinfo_1kg_hm3` file
#' PRS-CS reads for the whole panel (all chromosomes together).
#'
#' @param panel_dir Path to the panel's output directory.
#' @param chromosomes Integer vector of chromosomes to include, in order.
#' @param snpinfo_list Optional list of per-chromosome data frames (as returned
#'   by \code{\link{write_prscs_ld_chr}}), in the same order as `chromosomes`.
#'   If `NULL` (default), per-chromosome sidecar files are read from `panel_dir`.
#' @param cleanup_sidecars Logical; if `TRUE` and sidecar files were read from
#'   disk, delete them after a successful write.
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly, the path to the written `snpinfo_1kg_hm3` file.
#' @export
finalize_prscs_snpinfo <- function(panel_dir,
                                   chromosomes = 1:22,
                                   snpinfo_list = NULL,
                                   cleanup_sidecars = TRUE,
                                   verbose = TRUE) {
  if (!is.character(panel_dir) || length(panel_dir) != 1 || !nzchar(panel_dir)) {
    stop("'panel_dir' must be a non-empty directory path.", call. = FALSE)
  }
  if (!is.numeric(chromosomes) || length(chromosomes) == 0 || anyNA(chromosomes)) {
    stop("'chromosomes' must be a non-empty integer vector.", call. = FALSE)
  }

  used_sidecars <- is.null(snpinfo_list)

  if (used_sidecars) {
    sidecar_paths <- file.path(panel_dir, sprintf(".snpinfo_part_chr%d.tsv", chromosomes))
    missing <- !file.exists(sidecar_paths)
    if (any(missing)) {
      stop(sprintf("Missing snpinfo sidecar file(s) for chromosome(s): %s",
        paste(chromosomes[missing], collapse = ", ")), call. = FALSE)
    }
    snpinfo_list <- lapply(sidecar_paths, function(p) {
      utils::read.table(p, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    })
  } else {
    if (!is.list(snpinfo_list) || length(snpinfo_list) != length(chromosomes)) {
      stop("'snpinfo_list' must be a list with the same length as 'chromosomes'.", call. = FALSE)
    }
  }

  combined <- do.call(rbind, snpinfo_list)
  if (is.null(combined)) {
    combined <- .empty_snpinfo()
  }

  out <- data.frame(
    CHR = combined$chr,
    SNP = combined$rsid,
    BP = combined$bp,
    A1 = combined$a1,
    A2 = combined$a2,
    MAF = combined$maf,
    stringsAsFactors = FALSE
  )

  out_path <- file.path(panel_dir, "snpinfo_1kg_hm3")
  utils::write.table(out, file = out_path, row.names = FALSE, quote = FALSE)

  if (used_sidecars && cleanup_sidecars) {
    file.remove(sidecar_paths)
  }

  if (verbose) {
    message(sprintf("Wrote %d SNPs to %s", nrow(out), out_path))
  }

  invisible(out_path)
}

#' Generate a complete PRS-CS LD reference panel for one superpopulation
#'
#' Builds a PRS-CS-compatible LD reference panel directory (matching the format
#' read by `parse_genet.py`/`v2parse_genet_mix.py`) for one superpopulation,
#' looping \code{\link{write_prscs_ld_chr}} over the requested chromosomes and
#' assembling the combined snpinfo file via \code{\link{finalize_prscs_snpinfo}}.
#'
#' @param reference_index_file Reference panel index file.
#' @param reference_data_file Reference panel genotype data file.
#' @param reference_pop_desc_file Reference panel population description file.
#' @param superpopulation_label Superpopulation label (e.g. AFR/AMR/ASN/EUR/SAS).
#' @param output_dir Parent directory the panel directory will be created under.
#'   The panel directory itself is
#'   `file.path(output_dir, sprintf("ldblk_1kg_%s33kg", tolower(superpopulation_label)))`.
#' @param chromosomes Integer vector of chromosomes to process. Defaults to 1:22.
#' @param block_size Block size in base pairs (see \code{\link{write_prscs_ld_chr}}).
#' @param maf_cutoff Minimum minor allele frequency filter. Defaults to 0.01.
#' @param missing_cutoff Maximum allowed missing genotype fraction per SNP. Defaults to 1.0.
#' @param shrinkage Numeric shrinkage in `[0, 1)` applied before writing each block.
#' @param min_eigenvalue Eigenvalue floor for positive-definite correction.
#' @param panel_name Character identifier stored in HDF5 file attributes.
#' @param max_blocks_per_chr Optional cap on number of blocks per chromosome
#'   (useful for examples/tests).
#' @param dry_run Logical; if `TRUE`, only returns block boundaries and raw SNP
#'   counts per chromosome, without reading genotypes or creating any directory.
#' @param verbose Logical; print progress messages.
#'
#' @return If `dry_run = FALSE`, invisibly returns the panel directory path.
#'   If `dry_run = TRUE`, returns a data frame with block boundaries and SNP
#'   counts across all requested chromosomes.
#' @export
#'
#' @examples
#' \dontrun{
#' generate_prscs_ld_panel(
#'   reference_index_file = "33kg_index.gz",
#'   reference_data_file = "33kg_geno.gz",
#'   reference_pop_desc_file = "33kg_pop_desc.txt",
#'   superpopulation_label = "EUR",
#'   output_dir = "prscs_ref",
#'   chromosomes = 22,
#'   block_size = 250000,
#'   shrinkage = 0.01
#' )
#' }
generate_prscs_ld_panel <- function(reference_index_file,
                                    reference_data_file,
                                    reference_pop_desc_file,
                                    superpopulation_label,
                                    output_dir,
                                    chromosomes = 1:22,
                                    block_size = 250000,
                                    maf_cutoff = 0.01,
                                    missing_cutoff = 1.0,
                                    shrinkage = 0,
                                    min_eigenvalue = 1e-8,
                                    panel_name = NULL,
                                    max_blocks_per_chr = NULL,
                                    dry_run = FALSE,
                                    verbose = TRUE) {
  .validate_prscs_panel_args(
    reference_index_file = reference_index_file,
    superpopulation_label = superpopulation_label,
    output_dir = output_dir,
    chromosomes = chromosomes,
    block_size = block_size,
    dry_run = dry_run,
    verbose = verbose
  )

  if (dry_run) {
    idx <- .read_gauss_index(reference_index_file)
    parts <- lapply(chromosomes, function(chr) {
      idx_chr <- idx[idx$chr == chr, , drop = FALSE]
      if (nrow(idx_chr) == 0) {
        return(data.frame(chr = integer(0), block = integer(0), start_bp = integer(0),
          end_bp = integer(0), n_snps_raw = integer(0)))
      }
      idx_chr <- idx_chr[order(idx_chr$bp), , drop = FALSE]
      blocks <- .build_contiguous_blocks(idx_chr$bp, block_size = block_size)
      if (!is.null(max_blocks_per_chr)) {
        blocks <- blocks[seq_len(min(nrow(blocks), max_blocks_per_chr)), , drop = FALSE]
      }
      counts <- .block_counts(idx_chr$bp, blocks)
      data.frame(chr = chr, counts, stringsAsFactors = FALSE)
    })
    if (verbose) {
      message("Dry run complete. Returning block boundaries and SNP counts.")
    }
    return(do.call(rbind, parts))
  }

  panel_dir <- file.path(output_dir, sprintf("ldblk_1kg_%s33kg", tolower(superpopulation_label)))
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

  snpinfo_list <- lapply(chromosomes, function(chr) {
    write_prscs_ld_chr(
      chr = chr,
      superpopulation_label = superpopulation_label,
      panel_dir = panel_dir,
      reference_index_file = reference_index_file,
      reference_data_file = reference_data_file,
      reference_pop_desc_file = reference_pop_desc_file,
      block_size = block_size,
      maf_cutoff = maf_cutoff,
      missing_cutoff = missing_cutoff,
      shrinkage = shrinkage,
      min_eigenvalue = min_eigenvalue,
      panel_name = panel_name,
      max_blocks = max_blocks_per_chr,
      write_snpinfo_sidecar = FALSE,
      verbose = verbose
    )
  })

  finalize_prscs_snpinfo(
    panel_dir = panel_dir,
    chromosomes = chromosomes,
    snpinfo_list = snpinfo_list,
    cleanup_sidecars = FALSE,
    verbose = verbose
  )

  invisible(panel_dir)
}

#' Generate aligned PRS-CS LD reference panels for multiple superpopulations
#'
#' Like \code{\link{generate_prscs_ld_panel}}, but builds panels for two or
#' more superpopulations together and reconciles them so every population's
#' panel retains the *same* SNPs in each block.
#'
#' This matters because \code{\link{computeSuperPopLDWindow}}'s MAF filter is
#' applied to each population's own allele frequency: the same nominal
#' \code{maf_cutoff} can retain a SNP in one superpopulation's panel while
#' dropping it from another's, since allele frequencies differ by ancestry.
#' Left unreconciled, PRS-CS's mixing logic (\code{parse_ldblk_mix}), which
#' aligns panels by block index and a running cumulative SNP position across
#' blocks, silently misaligns every SNP after the first such mismatch in a
#' block -- pairing one population's SNP with a different SNP's row in
#' another population's panel, which is exactly the kind of bug that can
#' surface downstream as spurious zero/inflated betas.
#'
#' For each block, every requested superpopulation's raw SNP set is computed
#' independently (same MAF/missingness filtering as
#' \code{\link{generate_prscs_ld_panel}}), then only the SNPs common to
#' *every* requested superpopulation are kept -- in identical bp order -- in
#' every population's output. Retained correlation submatrices are re-shrunk
#' and re-PD'd (via the same shrinkage/eigenvalue floor) after subsetting. A
#' block where fewer than 2 SNPs survive reconciliation is written as an
#' empty stub block in *every* population's panel, so block indices stay
#' aligned across all of them.
#'
#' @param reference_index_file Reference panel index file.
#' @param reference_data_file Reference panel genotype data file.
#' @param reference_pop_desc_file Reference panel population description file.
#' @param superpopulation_labels Character vector of two or more
#'   superpopulation labels (e.g. \code{c("EUR", "AFR")}) to build aligned
#'   panels for.
#' @param output_dir Parent directory each panel directory will be created
#'   under.
#' @param chromosomes Integer vector of chromosomes to process. Defaults to 1:22.
#' @param block_size Block size in base pairs (see \code{\link{write_prscs_ld_chr}}).
#' @param maf_cutoff Minimum minor allele frequency filter, applied per
#'   population before reconciliation. Defaults to 0.01.
#' @param missing_cutoff Maximum allowed missing genotype fraction per SNP.
#'   Defaults to 1.0.
#' @param shrinkage Numeric shrinkage in `[0, 1)` applied to each reconciled
#'   block's correlation matrix.
#' @param min_eigenvalue Eigenvalue floor for positive-definite correction.
#' @param panel_names Optional named character vector of panel names keyed by
#'   (uppercased) superpopulation label. Defaults to `GAUSS_33KG_<LABEL>` for
#'   each requested label.
#' @param max_blocks_per_chr Optional cap on number of blocks per chromosome
#'   (useful for examples/tests).
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly, a named list of panel directory paths, one per
#'   (uppercased) superpopulation label.
#' @export
generate_prscs_ld_panels <- function(reference_index_file,
                                     reference_data_file,
                                     reference_pop_desc_file,
                                     superpopulation_labels,
                                     output_dir,
                                     chromosomes = 1:22,
                                     block_size = 250000,
                                     maf_cutoff = 0.01,
                                     missing_cutoff = 1.0,
                                     shrinkage = 0,
                                     min_eigenvalue = 1e-8,
                                     panel_names = NULL,
                                     max_blocks_per_chr = NULL,
                                     verbose = TRUE) {
  .validate_prscs_panels_args(
    reference_index_file = reference_index_file,
    reference_data_file = reference_data_file,
    reference_pop_desc_file = reference_pop_desc_file,
    superpopulation_labels = superpopulation_labels,
    output_dir = output_dir,
    chromosomes = chromosomes,
    block_size = block_size,
    maf_cutoff = maf_cutoff,
    missing_cutoff = missing_cutoff,
    shrinkage = shrinkage,
    min_eigenvalue = min_eigenvalue,
    max_blocks_per_chr = max_blocks_per_chr,
    verbose = verbose
  )

  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required. Please install it with install.packages('hdf5r').", call. = FALSE)
  }

  labels <- toupper(superpopulation_labels)

  if (is.null(panel_names)) {
    panel_names <- stats::setNames(sprintf("GAUSS_33KG_%s", labels), labels)
  } else if (!all(labels %in% names(panel_names))) {
    stop("'panel_names' must have an entry for every superpopulation label.", call. = FALSE)
  }

  panel_dirs <- stats::setNames(
    file.path(output_dir, sprintf("ldblk_1kg_%s33kg", tolower(labels))),
    labels
  )
  for (d in panel_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  idx <- .read_gauss_index(reference_index_file)

  snpinfo_by_pop <- stats::setNames(vector("list", length(labels)), labels)
  for (label in labels) snpinfo_by_pop[[label]] <- vector("list", length(chromosomes))

  for (k in seq_along(chromosomes)) {
    chr_result <- .write_prscs_ld_chr_aligned(
      chr = chromosomes[k],
      labels = labels,
      panel_dirs = panel_dirs,
      panel_names = panel_names,
      reference_index_file = reference_index_file,
      reference_data_file = reference_data_file,
      reference_pop_desc_file = reference_pop_desc_file,
      block_size = block_size,
      maf_cutoff = maf_cutoff,
      missing_cutoff = missing_cutoff,
      shrinkage = shrinkage,
      min_eigenvalue = min_eigenvalue,
      max_blocks = max_blocks_per_chr,
      idx = idx,
      verbose = verbose
    )
    for (label in labels) {
      snpinfo_by_pop[[label]][[k]] <- chr_result[[label]]
    }
  }

  for (label in labels) {
    finalize_prscs_snpinfo(
      panel_dir = panel_dirs[[label]],
      chromosomes = chromosomes,
      snpinfo_list = snpinfo_by_pop[[label]],
      cleanup_sidecars = FALSE,
      verbose = verbose
    )
  }

  invisible(as.list(panel_dirs))
}

.validate_prscs_panel_args <- function(reference_index_file,
                                       superpopulation_label,
                                       output_dir,
                                       chromosomes,
                                       block_size,
                                       dry_run,
                                       verbose) {
  if (!is.character(reference_index_file) || length(reference_index_file) != 1 || !nzchar(reference_index_file)) {
    stop("'reference_index_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_index_file)) {
    stop(sprintf("Reference index file does not exist: %s", reference_index_file), call. = FALSE)
  }
  if (!is.character(superpopulation_label) || length(superpopulation_label) != 1 || !nzchar(superpopulation_label)) {
    stop("'superpopulation_label' must be a non-empty string.", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1 || !nzchar(output_dir)) {
    stop("'output_dir' must be a non-empty directory path.", call. = FALSE)
  }
  if (!is.numeric(chromosomes) || length(chromosomes) == 0 || anyNA(chromosomes)) {
    stop("'chromosomes' must be a non-empty integer vector.", call. = FALSE)
  }
  if (!is.numeric(block_size) || length(block_size) != 1 || is.na(block_size) || block_size <= 0) {
    stop("'block_size' must be a positive number.", call. = FALSE)
  }
  if (!is.logical(dry_run) || length(dry_run) != 1 || is.na(dry_run)) {
    stop("'dry_run' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
}

.validate_prscs_chr_args <- function(chr,
                                     superpopulation_label,
                                     panel_dir,
                                     reference_index_file,
                                     reference_data_file,
                                     reference_pop_desc_file,
                                     block_size,
                                     maf_cutoff,
                                     missing_cutoff,
                                     shrinkage,
                                     min_eigenvalue,
                                     panel_name,
                                     max_blocks,
                                     write_snpinfo_sidecar,
                                     verbose) {
  if (!is.numeric(chr) || length(chr) != 1 || is.na(chr) || chr < 1) {
    stop("'chr' must be a positive integer.", call. = FALSE)
  }
  if (!is.character(superpopulation_label) || length(superpopulation_label) != 1 || !nzchar(superpopulation_label)) {
    stop("'superpopulation_label' must be a non-empty string.", call. = FALSE)
  }
  if (!is.character(panel_dir) || length(panel_dir) != 1 || !nzchar(panel_dir)) {
    stop("'panel_dir' must be a non-empty directory path.", call. = FALSE)
  }
  if (!dir.exists(panel_dir)) {
    stop(sprintf("'panel_dir' does not exist: %s", panel_dir), call. = FALSE)
  }

  if (!is.character(reference_index_file) || length(reference_index_file) != 1 || !nzchar(reference_index_file)) {
    stop("'reference_index_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_index_file)) {
    stop(sprintf("Reference index file does not exist: %s", reference_index_file), call. = FALSE)
  }

  if (!is.character(reference_data_file) || length(reference_data_file) != 1 || !nzchar(reference_data_file)) {
    stop("'reference_data_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_data_file)) {
    stop(sprintf("Reference genotype file does not exist: %s", reference_data_file), call. = FALSE)
  }

  if (!is.character(reference_pop_desc_file) || length(reference_pop_desc_file) != 1 || !nzchar(reference_pop_desc_file)) {
    stop("'reference_pop_desc_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_pop_desc_file)) {
    stop(sprintf("Reference population description file does not exist: %s", reference_pop_desc_file), call. = FALSE)
  }

  if (!is.numeric(block_size) || length(block_size) != 1 || is.na(block_size) || block_size <= 0) {
    stop("'block_size' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(maf_cutoff) || length(maf_cutoff) != 1 || is.na(maf_cutoff) || maf_cutoff < 0 || maf_cutoff >= 0.5) {
    stop("'maf_cutoff' must be in [0, 0.5).", call. = FALSE)
  }
  if (!is.numeric(missing_cutoff) || length(missing_cutoff) != 1 || is.na(missing_cutoff) || missing_cutoff < 0 || missing_cutoff > 1) {
    stop("'missing_cutoff' must be in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(shrinkage) || length(shrinkage) != 1 || is.na(shrinkage) || shrinkage < 0 || shrinkage >= 1) {
    stop("'shrinkage' must be in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(min_eigenvalue) || length(min_eigenvalue) != 1 || is.na(min_eigenvalue) || min_eigenvalue <= 0) {
    stop("'min_eigenvalue' must be positive.", call. = FALSE)
  }
  if (!is.null(panel_name) && (!is.character(panel_name) || length(panel_name) != 1 || !nzchar(panel_name))) {
    stop("'panel_name' must be NULL or a non-empty string.", call. = FALSE)
  }
  if (!is.null(max_blocks) && (!is.numeric(max_blocks) || length(max_blocks) != 1 || is.na(max_blocks) || max_blocks < 1)) {
    stop("'max_blocks' must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.logical(write_snpinfo_sidecar) || length(write_snpinfo_sidecar) != 1 || is.na(write_snpinfo_sidecar)) {
    stop("'write_snpinfo_sidecar' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
}

.validate_prscs_panels_args <- function(reference_index_file,
                                        reference_data_file,
                                        reference_pop_desc_file,
                                        superpopulation_labels,
                                        output_dir,
                                        chromosomes,
                                        block_size,
                                        maf_cutoff,
                                        missing_cutoff,
                                        shrinkage,
                                        min_eigenvalue,
                                        max_blocks_per_chr,
                                        verbose) {
  if (!is.character(reference_index_file) || length(reference_index_file) != 1 || !nzchar(reference_index_file)) {
    stop("'reference_index_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_index_file)) {
    stop(sprintf("Reference index file does not exist: %s", reference_index_file), call. = FALSE)
  }
  if (!is.character(reference_data_file) || length(reference_data_file) != 1 || !nzchar(reference_data_file)) {
    stop("'reference_data_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_data_file)) {
    stop(sprintf("Reference genotype file does not exist: %s", reference_data_file), call. = FALSE)
  }
  if (!is.character(reference_pop_desc_file) || length(reference_pop_desc_file) != 1 || !nzchar(reference_pop_desc_file)) {
    stop("'reference_pop_desc_file' must be a non-empty file path.", call. = FALSE)
  }
  if (!file.exists(reference_pop_desc_file)) {
    stop(sprintf("Reference population description file does not exist: %s", reference_pop_desc_file), call. = FALSE)
  }
  if (!is.character(superpopulation_labels) || length(superpopulation_labels) < 2 ||
      anyNA(superpopulation_labels) || any(!nzchar(superpopulation_labels))) {
    stop("'superpopulation_labels' must be a character vector of 2 or more non-empty labels.", call. = FALSE)
  }
  if (anyDuplicated(toupper(superpopulation_labels)) != 0) {
    stop("'superpopulation_labels' must not contain duplicate labels.", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1 || !nzchar(output_dir)) {
    stop("'output_dir' must be a non-empty directory path.", call. = FALSE)
  }
  if (!is.numeric(chromosomes) || length(chromosomes) == 0 || anyNA(chromosomes)) {
    stop("'chromosomes' must be a non-empty integer vector.", call. = FALSE)
  }
  if (!is.numeric(block_size) || length(block_size) != 1 || is.na(block_size) || block_size <= 0) {
    stop("'block_size' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(maf_cutoff) || length(maf_cutoff) != 1 || is.na(maf_cutoff) || maf_cutoff < 0 || maf_cutoff >= 0.5) {
    stop("'maf_cutoff' must be in [0, 0.5).", call. = FALSE)
  }
  if (!is.numeric(missing_cutoff) || length(missing_cutoff) != 1 || is.na(missing_cutoff) || missing_cutoff < 0 || missing_cutoff > 1) {
    stop("'missing_cutoff' must be in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(shrinkage) || length(shrinkage) != 1 || is.na(shrinkage) || shrinkage < 0 || shrinkage >= 1) {
    stop("'shrinkage' must be in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(min_eigenvalue) || length(min_eigenvalue) != 1 || is.na(min_eigenvalue) || min_eigenvalue <= 0) {
    stop("'min_eigenvalue' must be positive.", call. = FALSE)
  }
  if (!is.null(max_blocks_per_chr) && (!is.numeric(max_blocks_per_chr) || length(max_blocks_per_chr) != 1 ||
      is.na(max_blocks_per_chr) || max_blocks_per_chr < 1)) {
    stop("'max_blocks_per_chr' must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
}

.h5_write_attr <- function(h5obj, name, value) {
  h5obj$create_attr(name, value, dtype = hdf5r::guess_dtype(value))
}

.read_gauss_index <- function(reference_index_file) {
  con <- gzfile(reference_index_file, open = "rt")
  on.exit(close(con))
  idx <- utils::read.table(
    con,
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

.build_contiguous_blocks <- function(bp, block_size) {
  starts <- seq(from = min(bp), to = max(bp), by = block_size)
  data.frame(
    start_bp = as.integer(starts),
    end_bp = as.integer(starts + block_size - 1),
    stringsAsFactors = FALSE
  )
}

.block_counts <- function(bp, blocks) {
  n <- nrow(blocks)
  counts <- integer(n)
  for (i in seq_len(n)) {
    counts[i] <- sum(bp >= blocks$start_bp[i] & bp <= blocks$end_bp[i])
  }
  data.frame(block = seq_len(n), blocks, n_snps_raw = counts, stringsAsFactors = FALSE)
}

.empty_snpinfo <- function() {
  data.frame(chr = integer(0), rsid = character(0), bp = integer(0),
    a1 = character(0), a2 = character(0), maf = numeric(0), stringsAsFactors = FALSE)
}

# Rcpp::stop() messages from computePopLDWindow/computeSuperPopLDWindow that mean
# "this window has too few SNPs to compute LD" rather than a real error -- these
# are turned into an empty stub block so block indices stay aligned across every
# panel built against the same reference_index_file/block_size (see
# write_prscs_ld_chr's roxygen for why alignment matters). Any other error
# message is rethrown as-is.
.insufficient_snps_pattern <- "No SNPs found in requested region|Fewer than 2 SNPs remain"

.compute_and_write_ld_block <- function(blk_group, chr, start_bp, end_bp, ld_window_fn, shrinkage, min_eigenvalue) {
  res <- tryCatch(
    ld_window_fn(chr, start_bp, end_bp),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    if (grepl(.insufficient_snps_pattern, conditionMessage(res))) {
      .h5_write_attr(blk_group, "start_bp", as.integer(start_bp))
      .h5_write_attr(blk_group, "end_bp", as.integer(end_bp))
      .h5_write_attr(blk_group, "m_snps", as.integer(0))
      blk_group[["ldblk"]] <- matrix(numeric(0), nrow = 0, ncol = 0)
      blk_group[["snplist"]] <- character(0)
      blk_group[["bp"]] <- integer(0)
      blk_group[["a1"]] <- character(0)
      blk_group[["a2"]] <- character(0)
      return(.empty_snpinfo())
    }
    stop(res)
  }

  ld <- .shrink_and_make_pd(res$cormat, shrinkage = shrinkage, min_eigenvalue = min_eigenvalue)
  snplist <- res$snplist

  .h5_write_attr(blk_group, "start_bp", as.integer(start_bp))
  .h5_write_attr(blk_group, "end_bp", as.integer(end_bp))
  .h5_write_attr(blk_group, "m_snps", as.integer(nrow(snplist)))
  blk_group[["ldblk"]] <- ld
  blk_group[["snplist"]] <- as.character(snplist$rsid)
  blk_group[["bp"]] <- as.integer(snplist$bp)
  blk_group[["a1"]] <- as.character(snplist$a1)
  blk_group[["a2"]] <- as.character(snplist$a2)

  data.frame(
    chr = snplist$chr,
    rsid = snplist$rsid,
    bp = snplist$bp,
    a1 = snplist$a1,
    a2 = snplist$a2,
    maf = pmin(snplist$af1pop, 1 - snplist$af1pop),
    stringsAsFactors = FALSE
  )
}

# Writes one chromosome's blocks for every population in `labels` together,
# reconciling each block's SNP set across populations so panels stay aligned
# (see generate_prscs_ld_panels' roxygen for why). Mirrors write_prscs_ld_chr's
# structure (HDF5 attrs, on.exit cleanup, empty stub blocks) but drives all
# populations' HDF5 files from a single shared block loop.
.write_prscs_ld_chr_aligned <- function(chr, labels, panel_dirs, panel_names,
                                        reference_index_file, reference_data_file,
                                        reference_pop_desc_file, block_size,
                                        maf_cutoff, missing_cutoff, shrinkage,
                                        min_eigenvalue, max_blocks, idx, verbose) {
  idx_chr <- idx[idx$chr == chr, , drop = FALSE]
  if (nrow(idx_chr) == 0) {
    stop(sprintf("No SNPs found for chromosome %s in index file: %s", chr, reference_index_file), call. = FALSE)
  }
  idx_chr <- idx_chr[order(idx_chr$bp), , drop = FALSE]

  blocks <- .build_contiguous_blocks(idx_chr$bp, block_size = block_size)
  if (!is.null(max_blocks)) {
    blocks <- blocks[seq_len(min(nrow(blocks), max_blocks)), , drop = FALSE]
  }

  h5_by_pop <- stats::setNames(vector("list", length(labels)), labels)
  on.exit({
    for (h5 in h5_by_pop) if (!is.null(h5)) h5$close_all()
  }, add = TRUE)

  for (label in labels) {
    output_h5 <- file.path(panel_dirs[[label]], sprintf("ldblk_1kg_chr%d.hdf5", chr))
    if (file.exists(output_h5)) {
      file.remove(output_h5)
    }
    h5 <- hdf5r::H5File$new(output_h5, mode = "w")
    h5_by_pop[[label]] <- h5

    .h5_write_attr(h5, "panel_name", panel_names[[label]])
    .h5_write_attr(h5, "chr", as.integer(chr))
    .h5_write_attr(h5, "block_bp", as.integer(block_size))
    .h5_write_attr(h5, "maf_cutoff", as.numeric(maf_cutoff))
    .h5_write_attr(h5, "missing_cutoff", as.numeric(missing_cutoff))
    .h5_write_attr(h5, "shrinkage_lambda", as.numeric(shrinkage))
    .h5_write_attr(h5, "eig_floor", as.numeric(min_eigenvalue))
    .h5_write_attr(h5, "superpopulation_label", label)
    .h5_write_attr(h5, "n_blocks", as.integer(nrow(blocks)))
  }

  chr_snpinfo_by_pop <- stats::setNames(
    lapply(labels, function(label) vector("list", nrow(blocks))),
    labels
  )

  for (i in seq_len(nrow(blocks))) {
    start_bp <- blocks$start_bp[i]
    end_bp <- blocks$end_bp[i]

    if (verbose) {
      message(sprintf("Processing chr%s block %s/%s [%s, %s] across %s",
        chr, i, nrow(blocks), start_bp, end_bp, paste(labels, collapse = ", ")
      ))
    }

    blk_groups <- stats::setNames(
      lapply(labels, function(label) h5_by_pop[[label]]$create_group(sprintf("blk_%d", i))),
      labels
    )

    per_pop <- .compute_and_write_aligned_ld_block(
      blk_groups = blk_groups,
      chr = chr,
      start_bp = start_bp,
      end_bp = end_bp,
      labels = labels,
      reference_index_file = reference_index_file,
      reference_data_file = reference_data_file,
      reference_pop_desc_file = reference_pop_desc_file,
      maf_cutoff = maf_cutoff,
      missing_cutoff = missing_cutoff,
      shrinkage = shrinkage,
      min_eigenvalue = min_eigenvalue
    )

    for (label in labels) {
      chr_snpinfo_by_pop[[label]][[i]] <- per_pop[[label]]
    }
  }

  stats::setNames(
    lapply(labels, function(label) {
      out <- do.call(rbind, chr_snpinfo_by_pop[[label]])
      if (is.null(out)) out <- .empty_snpinfo()
      out
    }),
    labels
  )
}

# Computes each population's raw SNP set/LD for one block independently (same
# contract as .compute_and_write_ld_block), then keeps only the SNPs common to
# every population -- in identical bp order -- before writing each
# population's HDF5 block. A population's "insufficient SNPs" error
# contributes an empty rsid set to the intersection, so a block that's empty
# in any one population comes out empty in every population's panel.
.compute_and_write_aligned_ld_block <- function(blk_groups, chr, start_bp, end_bp, labels,
                                                reference_index_file, reference_data_file,
                                                reference_pop_desc_file, maf_cutoff,
                                                missing_cutoff, shrinkage, min_eigenvalue) {
  raw <- stats::setNames(vector("list", length(labels)), labels)
  for (label in labels) {
    raw[[label]] <- tryCatch(
      computeSuperPopLDWindow(
        chr = chr,
        start_bp = start_bp,
        end_bp = end_bp,
        superpopulation_label = label,
        reference_index_file = reference_index_file,
        reference_data_file = reference_data_file,
        reference_pop_desc_file = reference_pop_desc_file,
        maf_cutoff = maf_cutoff,
        missing_cutoff = missing_cutoff
      ),
      error = function(e) e
    )
  }

  for (label in labels) {
    r <- raw[[label]]
    if (inherits(r, "error") && !grepl(.insufficient_snps_pattern, conditionMessage(r))) {
      stop(r)
    }
  }

  # Identify "the same physical SNP" across populations by (bp, a1, a2) --
  # the same tuple the reference index itself uses to key a SNP (see
  # ReadReferenceIndex's MapKey in src/gauss.cpp) -- rather than by rsid
  # alone. Real reference panels commonly reuse a placeholder rsid (e.g. ".")
  # for un-annotated variants, so rsid alone can't disambiguate two distinct
  # SNPs that happen to share it within the same block; a1/a2 are consistent
  # across populations here since every call reads them from the same
  # reference_index_file, independent of which population was selected.
  snp_key <- function(snplist) paste(snplist$bp, snplist$a1, snplist$a2, sep = ":")

  key_sets <- lapply(labels, function(label) {
    r <- raw[[label]]
    if (inherits(r, "error")) character(0) else snp_key(r$snplist)
  })
  common_keys <- Reduce(intersect, key_sets)

  result <- stats::setNames(vector("list", length(labels)), labels)

  if (length(common_keys) < 2) {
    for (label in labels) {
      blk_group <- blk_groups[[label]]
      .h5_write_attr(blk_group, "start_bp", as.integer(start_bp))
      .h5_write_attr(blk_group, "end_bp", as.integer(end_bp))
      .h5_write_attr(blk_group, "m_snps", as.integer(0))
      blk_group[["ldblk"]] <- matrix(numeric(0), nrow = 0, ncol = 0)
      blk_group[["snplist"]] <- character(0)
      blk_group[["bp"]] <- integer(0)
      blk_group[["a1"]] <- character(0)
      blk_group[["a2"]] <- character(0)
      result[[label]] <- .empty_snpinfo()
    }
    return(result)
  }

  # Every element of common_keys is present in every label's snplist (that's
  # what Reduce(intersect, ...) guarantees), so bp lookup/match below never
  # produces NA regardless of which label we read bp positions from.
  bp_lookup <- stats::setNames(raw[[labels[1]]]$snplist$bp, snp_key(raw[[labels[1]]]$snplist))
  common_keys <- common_keys[order(bp_lookup[common_keys])]

  for (label in labels) {
    blk_group <- blk_groups[[label]]
    snplist <- raw[[label]]$snplist
    pos <- match(common_keys, snp_key(snplist))
    snplist_sub <- snplist[pos, , drop = FALSE]
    ld_sub <- raw[[label]]$cormat[pos, pos, drop = FALSE]
    ld <- .shrink_and_make_pd(ld_sub, shrinkage = shrinkage, min_eigenvalue = min_eigenvalue)

    .h5_write_attr(blk_group, "start_bp", as.integer(start_bp))
    .h5_write_attr(blk_group, "end_bp", as.integer(end_bp))
    .h5_write_attr(blk_group, "m_snps", as.integer(nrow(snplist_sub)))
    blk_group[["ldblk"]] <- ld
    blk_group[["snplist"]] <- as.character(snplist_sub$rsid)
    blk_group[["bp"]] <- as.integer(snplist_sub$bp)
    blk_group[["a1"]] <- as.character(snplist_sub$a1)
    blk_group[["a2"]] <- as.character(snplist_sub$a2)

    result[[label]] <- data.frame(
      chr = snplist_sub$chr,
      rsid = snplist_sub$rsid,
      bp = snplist_sub$bp,
      a1 = snplist_sub$a1,
      a2 = snplist_sub$a2,
      maf = pmin(snplist_sub$af1pop, 1 - snplist_sub$af1pop),
      stringsAsFactors = FALSE
    )
  }

  result
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

.read_pop_desc <- function(reference_pop_desc_file) {
  pop_desc <- utils::read.table(
    reference_pop_desc_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_cols <- c("Population_Abbreviation", "Super_Population")
  if (!all(required_cols %in% names(pop_desc))) {
    stop("reference_pop_desc_file must include Population_Abbreviation and Super_Population columns.", call. = FALSE)
  }
  pop_desc
}
