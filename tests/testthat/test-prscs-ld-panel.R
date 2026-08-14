# Test doubles for computeSuperPopLDWindow -----------------------------------
#
# The real engine needs a BGZF-compressed reference_data_file with true
# bgzf_seek() virtual offsets baked into reference_index_file's fpos column --
# a plain-gzip fixture can't satisfy that (see src/bgzf.c's bgzf_seek()). So
# these tests exercise the orchestration layer (write_prscs_ld_chr,
# finalize_prscs_snpinfo, generate_prscs_ld_panel) against a mocked
# computeSuperPopLDWindow that mimics the real engine's contract: a data frame
# of SNPs in [start_bp, end_bp], and the same two "insufficient SNPs" error
# messages the real engine throws. Real-data verification against the actual
# 33KG reference happens separately (see plan's RedHawk verification step).

.mock_snp_pool <- function() {
  data.frame(
    rsid = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6"),
    chr = rep(1L, 6),
    bp = c(10L, 20L, 150L, 160L, 170L, 500L),
    a1 = c("A", "C", "G", "T", "A", "C"),
    a2 = c("G", "T", "A", "C", "G", "T"),
    af1pop = c(0.1, 0.2, 0.3, 0.4, 0.1, 0.2),
    stringsAsFactors = FALSE
  )
}

.mock_compute_superpop_ld_window <- function(chr, start_bp, end_bp, superpopulation_label,
                                             reference_index_file, reference_data_file,
                                             reference_pop_desc_file, maf_cutoff = NULL,
                                             missing_cutoff = NULL) {
  if (identical(superpopulation_label, "ZZZ")) {
    stop("Invalid superpopulation label: no reference populations selected.")
  }
  pool <- .mock_snp_pool()
  sub <- pool[pool$chr == chr & pool$bp >= start_bp & pool$bp <= end_bp, , drop = FALSE]
  if (nrow(sub) == 0) {
    stop("No SNPs found in requested region for provided reference index.")
  }
  if (nrow(sub) < 2) {
    stop("Fewer than 2 SNPs remain after missingness/MAF filtering.")
  }
  m <- nrow(sub)
  cormat <- matrix(0.3, m, m)
  diag(cormat) <- 1
  list(snplist = sub, cormat = cormat)
}

# Fixture helpers --------------------------------------------------------

.write_toy_index <- function(path, snp_pool = .mock_snp_pool()) {
  idx <- data.frame(
    rsid = snp_pool$rsid, chr = snp_pool$chr, bp = snp_pool$bp,
    a1 = snp_pool$a1, a2 = snp_pool$a2, af1ref = snp_pool$af1pop, fpos = 0
  )
  con <- gzfile(path, "wt")
  write.table(idx, file = con, row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)
}

.write_dummy_file <- function(path, lines = "dummy") {
  writeLines(lines, path)
}

.write_toy_pop_desc <- function(path) {
  writeLines(c(
    "Population_Abbreviation\tNumber_of_Subjects\tSuper_Population",
    "P1\t10\tEUR",
    "P2\t10\tEUR"
  ), path)
}

# write_prscs_ld_chr ------------------------------------------------------

test_that("write_prscs_ld_chr writes blk_1..blk_N at HDF5 root with ldblk/snplist", {
  skip_if_not_installed("hdf5r")

  td <- tempfile("gausstest"); dir.create(td)
  idx_path <- file.path(td, "toy_index.gz")
  geno_path <- file.path(td, "toy_geno.gz")
  pop_desc_path <- file.path(td, "toy_pop_desc.txt")
  panel_dir <- file.path(td, "ldblk_1kg_eur33kg")
  dir.create(panel_dir)

  .write_toy_index(idx_path)
  .write_dummy_file(geno_path)
  .write_toy_pop_desc(pop_desc_path)

  testthat::local_mocked_bindings(
    computeSuperPopLDWindow = .mock_compute_superpop_ld_window
  )

  # bp range 10..500, block_size=100 -> 5 contiguous blocks:
  # [10,109]=2 SNPs, [110,209]=3 SNPs, [210,309]=0, [310,409]=0, [410,509]=1
  result <- write_prscs_ld_chr(
    chr = 1,
    superpopulation_label = "EUR",
    panel_dir = panel_dir,
    reference_index_file = idx_path,
    reference_data_file = geno_path,
    reference_pop_desc_file = pop_desc_path,
    block_size = 100,
    verbose = FALSE
  )

  out_path <- file.path(panel_dir, "ldblk_1kg_chr1.hdf5")
  expect_true(file.exists(out_path))

  h5 <- hdf5r::H5File$new(out_path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  # Root must contain ONLY blk_1..blk_5, no wrapping "blocks" group.
  root_names <- names(h5)
  expect_setequal(root_names, sprintf("blk_%d", 1:5))

  expect_equal(h5$attr_open("n_blocks")$read(), 5L)
  expect_equal(h5$attr_open("chr")$read(), 1L)
  expect_equal(h5$attr_open("superpopulation_label")$read(), "EUR")

  # blk_1: 2 SNPs, real content
  b1 <- h5[["blk_1"]]
  expect_true(all(c("ldblk", "snplist") %in% names(b1)))
  ld1 <- b1[["ldblk"]]$read()
  expect_equal(dim(ld1), c(2L, 2L))
  expect_equal(diag(ld1), c(1, 1), tolerance = 1e-8)
  expect_equal(ld1, t(ld1), tolerance = 1e-8)
  snplist1 <- b1[["snplist"]]$read()
  expect_setequal(snplist1, c("rs1", "rs2"))

  # blk_2: 3 SNPs, real content
  b2 <- h5[["blk_2"]]
  ld2 <- b2[["ldblk"]]$read()
  expect_equal(dim(ld2), c(3L, 3L))
  snplist2 <- b2[["snplist"]]$read()
  expect_setequal(snplist2, c("rs3", "rs4", "rs5"))

  # blk_3, blk_4: 0 SNPs -> empty stub blocks
  for (blk in c("blk_3", "blk_4")) {
    bg <- h5[[blk]]
    expect_equal(bg$attr_open("m_snps")$read(), 0L)
    expect_equal(dim(bg[["ldblk"]]$read()), c(0L, 0L))
    expect_length(bg[["snplist"]]$read(), 0)
  }

  # blk_5: 1 SNP -> "fewer than 2" stub block
  b5 <- h5[["blk_5"]]
  expect_equal(b5$attr_open("m_snps")$read(), 0L)
  expect_equal(dim(b5[["ldblk"]]$read()), c(0L, 0L))

  # Returned snpinfo covers exactly the 5 SNPs found in real (non-stub) blocks
  expect_equal(nrow(result), 5L)
  expect_setequal(result$rsid, c("rs1", "rs2", "rs3", "rs4", "rs5"))
  expect_true(all(c("chr", "rsid", "bp", "a1", "a2", "maf") %in% names(result)))

  # snpinfo sidecar was written
  sidecar <- file.path(panel_dir, ".snpinfo_part_chr1.tsv")
  expect_true(file.exists(sidecar))
})

test_that("write_prscs_ld_chr rethrows errors that are not 'insufficient SNPs'", {
  skip_if_not_installed("hdf5r")

  td <- tempfile("gausstest"); dir.create(td)
  idx_path <- file.path(td, "toy_index.gz")
  geno_path <- file.path(td, "toy_geno.gz")
  pop_desc_path <- file.path(td, "toy_pop_desc.txt")
  panel_dir <- file.path(td, "ldblk_1kg_zzz33kg")
  dir.create(panel_dir)

  .write_toy_index(idx_path)
  .write_dummy_file(geno_path)
  .write_toy_pop_desc(pop_desc_path)

  testthat::local_mocked_bindings(
    computeSuperPopLDWindow = .mock_compute_superpop_ld_window
  )

  expect_error(
    write_prscs_ld_chr(
      chr = 1,
      superpopulation_label = "ZZZ",
      panel_dir = panel_dir,
      reference_index_file = idx_path,
      reference_data_file = geno_path,
      reference_pop_desc_file = pop_desc_path,
      block_size = 100,
      verbose = FALSE
    ),
    "Invalid superpopulation label"
  )
})

# .build_contiguous_blocks -------------------------------------------------

test_that("contiguous blocks are non-overlapping and every SNP falls in exactly one", {
  bp <- c(10L, 20L, 150L, 160L, 170L, 500L)
  blocks <- .build_contiguous_blocks(bp, block_size = 100)

  # Contiguous, non-overlapping by construction: end_bp[i] + 1 == start_bp[i+1]
  expect_equal(blocks$start_bp[-1], blocks$end_bp[-nrow(blocks)] + 1L)

  membership <- vapply(bp, function(x) {
    sum(x >= blocks$start_bp & x <= blocks$end_bp)
  }, integer(1))
  expect_true(all(membership == 1L))
})

# finalize_prscs_snpinfo ----------------------------------------------------

test_that("finalize_prscs_snpinfo combines an in-memory snpinfo_list", {
  td <- tempfile("gausstest"); dir.create(td)

  part1 <- data.frame(chr = 1L, rsid = "rs1", bp = 10L, a1 = "A", a2 = "G", maf = 0.1, stringsAsFactors = FALSE)
  part2 <- data.frame(chr = 2L, rsid = "rs2", bp = 20L, a1 = "C", a2 = "T", maf = 0.2, stringsAsFactors = FALSE)

  out_path <- finalize_prscs_snpinfo(
    panel_dir = td,
    chromosomes = c(1, 2),
    snpinfo_list = list(part1, part2),
    verbose = FALSE
  )

  expect_equal(out_path, file.path(td, "snpinfo_1kg_hm3"))
  expect_true(file.exists(out_path))

  written <- read.table(out_path, header = TRUE, stringsAsFactors = FALSE)
  expect_equal(names(written), c("CHR", "SNP", "BP", "A1", "A2", "MAF"))
  expect_equal(written$SNP, c("rs1", "rs2"))
})

test_that("finalize_prscs_snpinfo reads and cleans up on-disk sidecars", {
  td <- tempfile("gausstest"); dir.create(td)

  part1 <- data.frame(chr = 1L, rsid = "rs1", bp = 10L, a1 = "A", a2 = "G", maf = 0.1, stringsAsFactors = FALSE)
  write.table(part1, file = file.path(td, ".snpinfo_part_chr1.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  out_path <- finalize_prscs_snpinfo(
    panel_dir = td,
    chromosomes = 1,
    cleanup_sidecars = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(out_path))
  expect_false(file.exists(file.path(td, ".snpinfo_part_chr1.tsv")))

  written <- read.table(out_path, header = TRUE, stringsAsFactors = FALSE)
  expect_equal(written$SNP, "rs1")
})

test_that("finalize_prscs_snpinfo errors clearly on missing sidecar files", {
  td <- tempfile("gausstest"); dir.create(td)
  expect_error(
    finalize_prscs_snpinfo(panel_dir = td, chromosomes = 1, verbose = FALSE),
    "Missing snpinfo sidecar"
  )
})

# generate_prscs_ld_panel ----------------------------------------------------

test_that("generate_prscs_ld_panel builds a '1kg'-named directory with panel + snpinfo", {
  skip_if_not_installed("hdf5r")

  td <- tempfile("gausstest"); dir.create(td)
  idx_path <- file.path(td, "toy_index.gz")
  geno_path <- file.path(td, "toy_geno.gz")
  pop_desc_path <- file.path(td, "toy_pop_desc.txt")
  output_dir <- file.path(td, "prscs_ref")
  dir.create(output_dir)

  .write_toy_index(idx_path)
  .write_dummy_file(geno_path)
  .write_toy_pop_desc(pop_desc_path)

  testthat::local_mocked_bindings(
    computeSuperPopLDWindow = .mock_compute_superpop_ld_window
  )

  panel_dir <- generate_prscs_ld_panel(
    reference_index_file = idx_path,
    reference_data_file = geno_path,
    reference_pop_desc_file = pop_desc_path,
    superpopulation_label = "EUR",
    output_dir = output_dir,
    chromosomes = 1,
    block_size = 100,
    verbose = FALSE
  )

  expect_true(grepl("1kg", basename(panel_dir)))
  expect_equal(basename(panel_dir), "ldblk_1kg_eur33kg")
  expect_true(file.exists(file.path(panel_dir, "ldblk_1kg_chr1.hdf5")))
  expect_true(file.exists(file.path(panel_dir, "snpinfo_1kg_hm3")))

  snpinfo <- read.table(file.path(panel_dir, "snpinfo_1kg_hm3"), header = TRUE, stringsAsFactors = FALSE)
  expect_equal(names(snpinfo), c("CHR", "SNP", "BP", "A1", "A2", "MAF"))
  expect_setequal(snpinfo$SNP, c("rs1", "rs2", "rs3", "rs4", "rs5"))

  # No leftover sidecar files after finalize
  expect_length(Sys.glob(file.path(panel_dir, ".snpinfo_part_chr*.tsv")), 0)
})

test_that("dry_run returns block boundaries and raw SNP counts without touching genotypes", {
  td <- tempfile("gausstest"); dir.create(td)
  idx_path <- file.path(td, "toy_index.gz")
  .write_toy_index(idx_path)

  plan <- generate_prscs_ld_panel(
    reference_index_file = idx_path,
    reference_data_file = "does_not_exist.gz",
    reference_pop_desc_file = "does_not_exist.txt",
    superpopulation_label = "EUR",
    output_dir = file.path(td, "prscs_ref"),
    chromosomes = 1,
    block_size = 100,
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_s3_class(plan, "data.frame")
  expect_true(all(c("chr", "block", "start_bp", "end_bp", "n_snps_raw") %in% names(plan)))
  expect_equal(nrow(plan), 5L)
  expect_equal(plan$n_snps_raw, c(2L, 3L, 0L, 0L, 1L))
  # dry_run must not create the output directory
  expect_false(dir.exists(file.path(td, "prscs_ref")))
})

# Argument validation --------------------------------------------------------

test_that("generate_prscs_ld_panel validation reports a missing reference_index_file", {
  expect_error(
    generate_prscs_ld_panel(
      reference_index_file = "missing.gz",
      reference_data_file = "missing2.gz",
      reference_pop_desc_file = "missing3.txt",
      superpopulation_label = "EUR",
      output_dir = tempdir()
    ),
    "Reference index file does not exist"
  )
})

test_that("write_prscs_ld_chr validation requires panel_dir to already exist", {
  td <- tempfile("gausstest"); dir.create(td)
  idx_path <- file.path(td, "toy_index.gz")
  geno_path <- file.path(td, "toy_geno.gz")
  pop_desc_path <- file.path(td, "toy_pop_desc.txt")
  .write_toy_index(idx_path)
  .write_dummy_file(geno_path)
  .write_toy_pop_desc(pop_desc_path)

  expect_error(
    write_prscs_ld_chr(
      chr = 1,
      superpopulation_label = "EUR",
      panel_dir = file.path(td, "does_not_exist"),
      reference_index_file = idx_path,
      reference_data_file = geno_path,
      reference_pop_desc_file = pop_desc_path,
      verbose = FALSE
    ),
    "'panel_dir' does not exist"
  )
})
