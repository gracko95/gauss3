\name{computeSuperPopLDWindow}
\alias{computeSuperPopLDWindow}
\title{Compute superpopulation-specific LD for one reference window}
\usage{
computeSuperPopLDWindow(
  chr,
  start_bp,
  end_bp,
  superpopulation_label,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  maf_cutoff = NULL,
  missing_cutoff = NULL
)
}
\arguments{
\item{chr}{Chromosome number.}
\item{start_bp}{Start base-pair position of the region.}
\item{end_bp}{End base-pair position of the region.}
\item{superpopulation_label}{Superpopulation label (AFR/AMR/ASN/EUR/SAS).}
\item{reference_index_file}{Reference panel index file.}
\item{reference_data_file}{Reference panel genotype data file.}
\item{reference_pop_desc_file}{Reference panel population description file.}
\item{maf_cutoff}{Optional MAF threshold (default 0.01).}
\item{missing_cutoff}{Optional missingness threshold in [0,1] (default 1.0).}
}
\value{A list with \code{snplist} and \code{cormat}.}
\description{Computes a signed LD correlation matrix for one genomic region and one selected superpopulation directly from GAUSS reference files.}
\examples{
\dontrun{
res <- computeSuperPopLDWindow(
  chr = 22,
  start_bp = 17000000,
  end_bp = 17100000,
  superpopulation_label = "EUR",
  reference_index_file = "33kg_index.gz",
  reference_data_file = "33kg_geno.gz",
  reference_pop_desc_file = "33kg_pop_desc.txt"
)
str(res$snplist)
dim(res$cormat)
}
}
