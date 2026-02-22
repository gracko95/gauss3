# Example: generate PRS-CS LD blocks for a single chromosome using GAUSS 33KG files.
#
# Usage:
# Rscript inst/scripts/generate_prscs_ld_chr_example.R /path/to/33kg_geno.gz /path/to/33kg_index.gz out_chr22.h5

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_prscs_ld_chr_example.R <33kg_geno.gz> <33kg_index.gz> <output.h5>")
}

geno_file <- args[[1]]
index_file <- args[[2]]
out_file <- args[[3]]

message("Running generate_prscs_ld_panel for chr22 with 3 windows...")

gauss::generate_prscs_ld_panel(
  reference_geno_file = geno_file,
  reference_index_file = index_file,
  chromosome = 22,
  output_h5 = out_file,
  window_size = 2e6,
  step_size = 1e6,
  n_reference_populations = 29,
  max_windows = 3,
  shrinkage = 0.01,
  verbose = TRUE
)

message("Done: ", out_file)
