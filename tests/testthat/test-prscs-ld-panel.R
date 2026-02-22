test_that("generate_prscs_ld_panel writes HDF5 blocks with metadata", {
  skip_if_not_installed("hdf5r")

  td <- tempdir()
  idx_path <- file.path(td, "toy_index.gz")
  geno_path <- file.path(td, "toy_geno.gz")
  out_path <- file.path(td, "toy_chr1.h5")

  idx <- data.frame(
    rsid = c("rs1", "rs2", "rs3"),
    chr = c(1, 1, 1),
    bp = c(100, 150, 220),
    a1 = c("A", "C", "G"),
    a2 = c("G", "T", "A"),
    af1ref = c(0.1, 0.2, 0.3),
    fpos = c(0, 0, 0)
  )

  con_idx <- gzfile(idx_path, "wt")
  write.table(idx, file = con_idx, row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con_idx)

  geno_lines <- c(
    "012 111 0.1 0.2",
    "120 101 0.2 0.3",
    "201 011 0.3 0.4"
  )
  con_geno <- gzfile(geno_path, "wt")
  writeLines(geno_lines, con_geno)
  close(con_geno)

  generate_prscs_ld_panel(
    reference_geno_file = geno_path,
    reference_index_file = idx_path,
    chromosome = 1,
    output_h5 = out_path,
    window_size = 100,
    step_size = 100,
    n_reference_populations = 2,
    shrinkage = 0.05,
    maf_min = 0,
    missing_max = 1,
    panel_name = "toy_panel",
    verbose = FALSE
  )

  expect_true(file.exists(out_path))

  h5 <- hdf5r::H5File$new(out_path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  expect_equal(h5$attr_open("panel_name")$read(), "toy_panel")
  expect_equal(h5$attr_open("chr")$read(), 1L)

  block_names <- names(h5[["blocks"]])
  expect_true(length(block_names) >= 1)

  b1 <- h5[["blocks"]][[block_names[[1]]]]
  expect_true("rsid" %in% names(b1))
  expect_true("a1" %in% names(b1))
  expect_true("a2" %in% names(b1))

  ld <- b1[["ld"]][]
  expect_equal(nrow(ld), ncol(ld))
  expect_true(all(is.finite(ld)))
  if (length(ld) > 0) {
    expect_equal(ld, t(ld), tolerance = 1e-8)
    expect_equal(diag(ld), rep(1, nrow(ld)), tolerance = 1e-8)
  }
})

test_that("dry_run returns window SNP counts", {
  td <- tempdir()
  idx_path <- file.path(td, "toy_index_dry.gz")
  geno_path <- file.path(td, "toy_geno_dry.gz")

  idx <- data.frame(
    rsid = c("rs1", "rs2"), chr = c(1, 1), bp = c(100, 200),
    a1 = c("A", "C"), a2 = c("G", "T"), af1ref = c(0.1, 0.2), fpos = c(0, 0)
  )
  con_idx <- gzfile(idx_path, "wt")
  write.table(idx, file = con_idx, row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con_idx)

  con_geno <- gzfile(geno_path, "wt")
  writeLines(c("012 111", "120 101"), con_geno)
  close(con_geno)

  plan <- generate_prscs_ld_panel(
    reference_geno_file = geno_path,
    reference_index_file = idx_path,
    chromosome = 1,
    output_h5 = tempfile(fileext = ".h5"),
    window_size = 100,
    step_size = 100,
    n_reference_populations = 2,
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_s3_class(plan, "data.frame")
  expect_true(all(c("block", "start_bp", "end_bp", "n_snps") %in% names(plan)))
})

test_that("argument validation returns informative errors", {
  expect_error(
    generate_prscs_ld_panel(
      reference_geno_file = "missing.gz",
      reference_index_file = "missing2.gz",
      chromosome = 1,
      output_h5 = tempfile(fileext = ".h5")
    ),
    "Reference genotype file does not exist"
  )
})
