# Example: compute one superpopulation-specific LD block from GAUSS 33KG references.
#
# Usage:
# Rscript inst/scripts/compute_superpop_ld_window_example.R \
#   /path/to/33kg_index.gz /path/to/33kg_geno.gz /path/to/ref_pop_desc.txt

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compute_superpop_ld_window_example.R <33kg_index.gz> <33kg_geno.gz> <ref_pop_desc.txt>")
}

idx <- args[[1]]
geno <- args[[2]]
pop_desc <- args[[3]]

res <- gauss::computeSuperPopLDWindow(
  chr = 22,
  start_bp = 17000000,
  end_bp = 17100000,
  superpopulation_label = "EUR",
  reference_index_file = idx,
  reference_data_file = geno,
  reference_pop_desc_file = pop_desc,
  maf_cutoff = 0.01,
  missing_cutoff = 0
)

print(utils::head(res$snplist))
print(dim(res$cormat))
