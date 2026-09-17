test_that("computePopLDWindow validates population index before LD computation", {
  td <- tempdir()
  pop_desc <- file.path(td, "pop_desc.txt")

  writeLines(c(
    "Population_Abbreviation\tNumber_of_Subjects\tSuper_Population",
    "POP1\t10\tSUP1",
    "POP2\t10\tSUP1"
  ), pop_desc)

  expect_error(
    computePopLDWindow(
      chr = 1,
      start_bp = 1,
      end_bp = 100,
      population_index = 3,
      reference_index_file = "missing_index.gz",
      reference_data_file = "missing_geno.gz",
      reference_pop_desc_file = pop_desc
    ),
    "Invalid population_index"
  )
})
