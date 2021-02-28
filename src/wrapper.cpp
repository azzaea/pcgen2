#include "gemma.h"
#include <gsl/gsl_sf_bessel.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <Rcpp.h>
#define STRICT_R_HEADERS

using namespace Rcpp;
using namespace std;

// source using Rcpp::sourceCpp

//' R interface to 'Genome-wide Efficient Mixed Model Association' (GEMMA)
//'
//' GEMMA is an efficient association testing software of Linear Mixed Models
//' (LMM)s and related GWAS models. This interface is for the Multivariate
//' linear mixed model (mvLMM) functionality which also corrects for population
//' structure and sample (non)exchangability
//' @param genoinputs Character vector of input files. It should be either the
//'   prefix name of PLINK binary ped file, or the 3 BIMBAM files' names (the
//'   menan genotypes, the phenotyes and SNP annotaiton files)
//' @param kfile Character variable for the relatedness matrix file name (can be
//'   in gzip compressed format)
//' @param colnums Numeric vector for specifying column number from the
//'   phenotypes file used for association testing
//' @param miss,maf,r2,hwe Floating point variables for filtering SNPs at
//'   missingness cutt-off (default 0.05), minor allele frequency (default
//'   0.01), r2 threshold (default 0.9999) and HWE test p value threshold
//'   (default 0; no test)
//' @param notsnp Boolean variable to use real values as covariates (minor
//'   allele frequency cutoff is not used)
//' @param lmmMode Integer specifying the frequenist test to use: 1 for the Wald
//'   test, 2 for the likelihood ratio test, 3 for score test; and 4 for
//'   performing all the three tests
//' @param gxe Character variable for the file containing a column of
//'   environmental variables. If provided, GEMMA fits a LMM controlling for
//'   both the SNP main effect and the environmental main effect while testing
//'   for interaction effect.
//' @param outprefix Character variable for output file prefix
//' @param outdir Character variable for output directory path
//' @param license Boolean variable for printing GEMMA's license information
//'
// [[Rcpp::export]]
bool gemmaMVLMM(CharacterVector genoinputs,
                std::string kfile,
                NumericVector colnums,
                double miss = 0.05,
                double maf = 0.01,
                double r2 = 0.9999,
                double hwe = 0,
                bool notsnp = false,
                int lmmMode = 1,
                std::string gxe = "",
                std::string outprefix = "out",
                std::string outdir = "output",
                bool license = false){
  // @param predit Boolean variable to impute missing phenotypes before
  //   association testing (if a small proportion is missing)

  // int main(int argc, char *argv[])

  GEMMA cGemma;
  PARAM cPar;

  gsl_set_error_handler (&gemma_gsl_error_handler);

  if (license) {
    cGemma.PrintHeader();
    cGemma.PrintLicense();
    return EXIT_SUCCESS;
  }

  //cGemma.Assign(argc, argv, cPar);

  if (genoinputs.size() == 1){
    cPar.file_bfile = genoinputs[0];
    // Rcout << "-bfile: " << cPar.file_bfile << "\n";
  }
  if (genoinputs.size() >= 2) {
    cPar.file_geno = genoinputs[0];
    cPar.file_pheno = genoinputs[1];
  }
  if (genoinputs.size() == 3)
    cPar.file_anno = genoinputs[2];
  if (genoinputs.size() > 3){
    cGemma.PrintHeader();
    cGemma.PrintHelp(0);
    return EXIT_SUCCESS;
  }
  // Rcout << "-g: " << cPar.file_geno << "\t -p: " << cPar.file_pheno << "\t -a: " <<
  //cPar.file_anno << "\n";

  // set pheno column (list/range)
  (cPar.p_column).clear();
  for (int i = 0; i < colnums.size(); i++) {
    (cPar.p_column).push_back(colnums[i]);
    // Rcout << "-n:" << cPar.p_column[i] << "\n";
  }

  cPar.file_gxe = gxe;
  // Rcout << "-gxe: " << cPar.file_gxe << "\n";
  cPar.file_kin = kfile;
  // Rcout << "-k: " << cPar.file_kin << "\n";

  cPar.path_out = outdir;
  // Rcout << "-outdir:" << cPar.path_out << "\n";
  cPar.file_out = outprefix;
  // Rcout << "-o: " << cPar.file_out << "\n";

  // SNP QC Options:
  cPar.miss_level = miss;
  if (cPar.maf_level != -1) {
    cPar.maf_level = maf;
  }
  cPar.hwe_level = hwe;
  cPar.r2_level = r2;
  if (notsnp)
    cPar.maf_level = -1;

  // lmm options:
  if (cPar.a_mode != 0) {
    cPar.error = true;
    cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
    "-predict -calccor options is allowed."
    << endl;
    //break;
  }
  cPar.a_mode = lmmMode;
  // Rcout << "-lmm:" << cPar.a_mode << "\n";



  ifstream check_dir((cPar.path_out).c_str());
  if (!check_dir) {
#ifdef WINDOWS
    mkdir((cPar.path_out).c_str());
#else
    mkdir((cPar.path_out).c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif
  }

  if (!is_quiet_mode())
    cGemma.PrintHeader();

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  if (is_quiet_mode()) {
    stringstream ss;
    cout.rdbuf(ss.rdbuf());
  }

  cPar.CheckParam();

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  cGemma.BatchRun(cPar);

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }
/*
  cGemma.WriteLog(argc, argv, cPar);

  return EXIT_SUCCESS;
*/
  return true;
}

// for testing and development, this R code will be automatically
// run after the compilation.

/*** R

gemmaMVLMM(genoinputs = c("/home/p287664/github_repos/GEMMA/example/mouse_hs1940.geno.txt.gz",
                          "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.pheno.txt",
                          "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.anno.txt"),
           kfile = "/home/p287664/github_repos/GEMMA/output/mouse_hs1940.cXX.txt",
           colnums = c(1, 6), outprefix = "mouse_hs1940_CD8MCH_lmm", outdir = "rgemma")


*/
