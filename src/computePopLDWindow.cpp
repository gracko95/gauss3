#include <Rcpp.h>

using namespace Rcpp;

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "gauss.h"
#include "snp.h"
#include "util.h"

namespace {

struct PopLDWindowCleanup {
  std::map<MapKey, Snp*, LessThanMapKey>& snp_map;
  std::vector<Snp*>& snp_vec;
  bool genotype_loaded;

  PopLDWindowCleanup(std::map<MapKey, Snp*, LessThanMapKey>& map_ref,
                     std::vector<Snp*>& vec_ref)
      : snp_map(map_ref), snp_vec(vec_ref), genotype_loaded(false) {}

  ~PopLDWindowCleanup() {
    if (genotype_loaded) {
      FreeGenotype(snp_vec);
    }
    for (std::map<MapKey, Snp*, LessThanMapKey>::iterator it = snp_map.begin(); it != snp_map.end();) {
      (it->second)->ClearSnp();
      delete it->second;
      snp_map.erase(it++);
    }
  }
};

}  // namespace

//' Compute population-specific LD for one reference window
//'
//' Computes a signed LD correlation matrix for one genomic region and one selected
//' reference population directly from GAUSS 33KG reference files, without GWAS
//' summary statistics input.
//'
//' @param chr Chromosome number.
//' @param start_bp Start base-pair position of the region.
//' @param end_bp End base-pair position of the region.
//' @param population_index 1-based index of the selected population in
//'   `reference_pop_desc_file`.
//' @param reference_index_file Reference panel index file.
//' @param reference_data_file Reference panel genotype data file.
//' @param reference_pop_desc_file Reference panel population description file.
//' @param maf_cutoff Optional MAF threshold. Defaults to 0.01.
//' @param missing_cutoff Optional missingness threshold in [0,1]. Defaults to 1.0.
//' @return A list with `snplist` (data frame: rsid, chr, bp, a1, a2, af1pop)
//'   and `cormat` (signed correlation matrix).
// [[Rcpp::export]]
Rcpp::List computePopLDWindow(
    int chr,
    long long int start_bp,
    long long int end_bp,
    int population_index,
    std::string reference_index_file,
    std::string reference_data_file,
    std::string reference_pop_desc_file,
    Rcpp::Nullable<double> maf_cutoff = R_NilValue,
    Rcpp::Nullable<double> missing_cutoff = R_NilValue) {

  if (start_bp > end_bp) {
    Rcpp::stop("Invalid region: 'start_bp' must be <= 'end_bp'.");
  }

  const double maf_min = maf_cutoff.isNotNull() ? Rcpp::as<double>(maf_cutoff) : 0.01;
  const double missing_max = missing_cutoff.isNotNull() ? Rcpp::as<double>(missing_cutoff) : 1.0;

  if (!std::isfinite(maf_min) || maf_min < 0.0 || maf_min >= 0.5) {
    Rcpp::stop("Invalid 'maf_cutoff': must be finite and in [0, 0.5)." );
  }
  if (!std::isfinite(missing_max) || missing_max < 0.0 || missing_max > 1.0) {
    Rcpp::stop("Invalid 'missing_cutoff': must be finite and in [0, 1].");
  }

  Arguments args;
  args.chr = chr;
  args.start_bp = start_bp;
  args.end_bp = end_bp;
  args.wing_size = 0;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;

  read_ref_desc(args);

  if (population_index < 1 || population_index > args.num_pops) {
    Rcpp::stop("Invalid population_index: must be between 1 and %d (inclusive).", args.num_pops);
  }

  const std::string selected_pop = args.ref_pop_vec[population_index - 1];
  args.pop_wgt_map[selected_pop] = 1.0;
  init_pop_flag_wgt_vec(args);

  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  std::vector<Snp*> snp_vec;
  PopLDWindowCleanup cleanup(snp_map, snp_vec);

  ReadReferenceIndex(snp_map, args);

  for (std::map<MapKey, Snp*, LessThanMapKey>::iterator it = snp_map.begin(); it != snp_map.end(); ++it) {
    snp_vec.push_back(it->second);
  }

  if (snp_vec.empty()) {
    Rcpp::stop("No SNPs found in requested region for provided reference index.");
  }

  ReadGenotype(snp_vec, args);
  cleanup.genotype_loaded = true;

  std::vector<Snp*> qc_snps;
  std::vector<double> qc_af1;
  qc_snps.reserve(snp_vec.size());
  qc_af1.reserve(snp_vec.size());

  for (size_t i = 0; i < snp_vec.size(); ++i) {
    Snp* snp = snp_vec[i];
    std::vector<std::string>& geno_vec = snp->GetGenotypeVec();

    if (geno_vec.size() != 1) {
      Rcpp::stop("Single-population mode expected genotype_vec size 1, but got " +
                 std::to_string(geno_vec.size()) + " for SNP " + snp->GetRsid() + ".");
    }

    const std::string& geno = geno_vec[0];
    if (geno.empty()) {
      continue;
    }

    int total_n = static_cast<int>(geno.size());
    int nonmissing_n = 0;
    int missing_n = 0;
    double dosage_sum = 0.0;

    for (int j = 0; j < total_n; ++j) {
      char c = geno[j];
      if (c == '0' || c == '1' || c == '2') {
        dosage_sum += static_cast<double>(c - '0');
        nonmissing_n++;
      } else {
        missing_n++;
      }
    }

    if (nonmissing_n == 0) {
      continue;
    }

    const double missing_frac = static_cast<double>(missing_n) / static_cast<double>(total_n);
    const double af1pop = dosage_sum / (2.0 * static_cast<double>(nonmissing_n));
    const double maf = std::min(af1pop, 1.0 - af1pop);

    if (missing_frac <= missing_max && maf >= maf_min) {
      qc_snps.push_back(snp);
      qc_af1.push_back(af1pop);
    }
  }

  if (qc_snps.size() < 2) {
    Rcpp::stop("Fewer than 2 SNPs remain after missingness/MAF filtering.");
  }

  std::vector<Snp*> kept_snps;
  std::vector<double> kept_af1;
  std::vector<double> snp_sd;
  kept_snps.reserve(qc_snps.size());
  kept_af1.reserve(qc_snps.size());
  snp_sd.reserve(qc_snps.size());

  for (size_t i = 0; i < qc_snps.size(); ++i) {
    Snp* snp = qc_snps[i];
    std::vector<std::string>& geno_vec = snp->GetGenotypeVec();
    const std::string& geno = geno_vec[0];

    // CalWgtCov does not support non-(0/1/2) characters; enforce this strictly.
    bool has_missing = false;
    for (size_t j = 0; j < geno.size(); ++j) {
      char c = geno[j];
      if (!(c == '0' || c == '1' || c == '2')) {
        has_missing = true;
        break;
      }
    }
    if (has_missing) {
      continue;
    }

    const double v = CalWgtCov(geno_vec, geno_vec, args.pop_wgt_vec);
    if (!std::isfinite(v) || v <= 0.0) {
      continue;
    }

    kept_snps.push_back(snp);
    kept_af1.push_back(qc_af1[i]);
    snp_sd.push_back(std::sqrt(v));
  }

  if (kept_snps.size() < 2) {
    Rcpp::stop("Fewer than 2 SNPs remain after variance/non-missing genotype checks.");
  }

  const int m = static_cast<int>(kept_snps.size());
  Rcpp::NumericMatrix cor_mat(m, m);

  for (int i = 0; i < m; ++i) {
    cor_mat(i, i) = 1.0;
    for (int j = i + 1; j < m; ++j) {
      const double cov = CalWgtCov(kept_snps[i]->GetGenotypeVec(), kept_snps[j]->GetGenotypeVec(), args.pop_wgt_vec);
      const double cor = cov / (snp_sd[i] * snp_sd[j]);
      if (!std::isfinite(cor)) {
        Rcpp::stop("Non-finite correlation encountered for SNP pair (" +
                   kept_snps[i]->GetRsid() + ", " + kept_snps[j]->GetRsid() + ").");
      }
      cor_mat(i, j) = cor;
      cor_mat(j, i) = cor;
    }
  }

  Rcpp::StringVector rsid_vec;
  Rcpp::IntegerVector chr_vec;
  Rcpp::IntegerVector bp_vec;
  Rcpp::StringVector a1_vec;
  Rcpp::StringVector a2_vec;
  Rcpp::NumericVector af1pop_vec;


  for (int i = 0; i < m; ++i) {
    Snp* snp = kept_snps[i];
    rsid_vec.push_back(snp->GetRsid());
    chr_vec.push_back(snp->GetChr());
    bp_vec.push_back(snp->GetBp());
    a1_vec.push_back(snp->GetA1());
    a2_vec.push_back(snp->GetA2());
    af1pop_vec.push_back(kept_af1[i]);
  }

  Rcpp::DataFrame snplist = Rcpp::DataFrame::create(
    Rcpp::Named("rsid") = rsid_vec,
    Rcpp::Named("chr") = chr_vec,
    Rcpp::Named("bp") = bp_vec,
    Rcpp::Named("a1") = a1_vec,
    Rcpp::Named("a2") = a2_vec,
    Rcpp::Named("af1pop") = af1pop_vec
  );

  return Rcpp::List::create(
    Rcpp::Named("snplist") = snplist,
    Rcpp::Named("cormat") = cor_mat
  );
}
