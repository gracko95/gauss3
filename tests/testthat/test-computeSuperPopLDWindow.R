test_that("computeSuperPopLDWindow validates superpopulation label", {
  td <- tempdir()
  pop_desc <- file.path(td, "pop_desc_sup.txt")

  writeLines(c(
    "Population_Abbreviation\tNumber_of_Subjects\tSuper_Population",
    "P1\t10\tEUR",
    "P2\t10\tAFR"
  ), pop_desc)

  expect_error(
    computeSuperPopLDWindow(
      chr = 1,
      start_bp = 1,
      end_bp = 100,
      superpopulation_label = "XXX",
      reference_index_file = "missing_index.gz",
      reference_data_file = "missing_geno.gz",
      reference_pop_desc_file = pop_desc
    ),
    "invalid population name|Invalid superpopulation"
  )
})
