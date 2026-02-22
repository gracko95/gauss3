\name{generate_prscs_ld_panel}
\alias{generate_prscs_ld_panel}
\title{Generate PRS-CS compatible LD panels from 33KG GAUSS reference files}
\usage{
generate_prscs_ld_panel(
  reference_geno_file,
  reference_index_file,
  chromosome,
  output_h5,
  window_size = 3e+06,
  step_size = 1e+06,
  n_reference_populations = 29,
  population_indices = NULL,
  shrinkage = 0,
  min_eigenvalue = 1e-08,
  maf_min = 0,
  missing_max = 1,
  max_windows = NULL,
  dry_run = FALSE,
  panel_name = "GAUSS_33KG",
  dosage_allele = c("A1", "A2"),
  verbose = TRUE
)
}
\arguments{
\item{reference_geno_file}{Path to \code{33kg_geno.gz}.}
\item{reference_index_file}{Path to \code{33kg_index.gz}.}
\item{chromosome}{Integer chromosome to process.}
\item{output_h5}{Output HDF5 file path.}
\item{window_size}{Window size in base pairs.}
\item{step_size}{Step size in base pairs.}
\item{n_reference_populations}{Number of population genotype strings per line in reference geno file.}
\item{population_indices}{Indices of populations to include.}
\item{shrinkage}{Shrinkage parameter in \code{[0, 1)}.}
\item{min_eigenvalue}{Minimum eigenvalue used when enforcing positive definiteness.}
\item{maf_min}{Minimum minor allele frequency threshold.}
\item{missing_max}{Maximum missing genotype fraction allowed per SNP.}
\item{max_windows}{Optional cap on number of processed windows.}
\item{dry_run}{If \code{TRUE}, return only planned windows and SNP counts (no HDF5 write).}
\item{panel_name}{Panel identifier stored in HDF5 file attributes.}
\item{dosage_allele}{Allele represented by genotype dosage values.}
\item{verbose}{Logical; print progress messages.}
}
\value{If \code{dry_run = FALSE}, invisibly returns \code{output_h5}. If \code{dry_run = TRUE}, returns a data frame with per-window SNP counts.}
\description{
Builds per-chromosome LD blocks directly from GAUSS 33KG reference files without converting to PLINK/BGEN/VCF.
}
\examples{
\dontrun{
generate_prscs_ld_panel(
  reference_geno_file = "33kg_geno.gz",
  reference_index_file = "33kg_index.gz",
  chromosome = 22,
  output_h5 = "chr22_prscs_ld.h5",
  window_size = 2e6,
  step_size = 1e6,
  n_reference_populations = 29,
  max_windows = 3,
  shrinkage = 0.01
)
}
}
