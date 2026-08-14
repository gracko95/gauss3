# Example: generate a PRS-CS LD reference panel for one superpopulation using
# GAUSS 33KG files (chr22 only, capped block count for a fast example run).
#
# Usage:
# Rscript inst/scripts/generate_prscs_ld_chr_example.R \
#   /path/to/33kg_index.gz /path/to/33kg_geno.gz /path/to/33kg_pop_desc.txt \
#   /path/to/output_dir EUR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript generate_prscs_ld_chr_example.R <33kg_index.gz> <33kg_geno.gz> <33kg_pop_desc.txt> <output_dir> <superpopulation_label>")
}

index_file <- args[[1]]
geno_file <- args[[2]]
pop_desc_file <- args[[3]]
output_dir <- args[[4]]
superpop <- args[[5]]

message(sprintf("Running generate_prscs_ld_panel for chr22, superpopulation %s, 3 blocks...", superpop))

panel_dir <- gauss::generate_prscs_ld_panel(
  reference_index_file = index_file,
  reference_data_file = geno_file,
  reference_pop_desc_file = pop_desc_file,
  superpopulation_label = superpop,
  output_dir = output_dir,
  chromosomes = 22,
  block_size = 2e6,
  max_blocks_per_chr = 3,
  shrinkage = 0.01,
  verbose = TRUE
)

message("Done: ", panel_dir)
