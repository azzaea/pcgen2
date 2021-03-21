//#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include "lapack.h"

#include "gemma.h"
#include "gemma_io.h"
#include "lmm.h"
#include "mathfunc.h"
#include "mvlmm.h"


#include <sys/stat.h>
#include <sys/types.h>


#define STRICT_R_HEADERS

using namespace Rcpp;
using namespace std;

// source using Rcpp::sourceCpp
/* PKG_CPPFLAGS = -isystem/home/p287664/software/OpenBLAS/include/ -Icontrib/catch-1.9.7 -Isrc
 * PKG_LIBS = -lopenblas -lgsl -lz -pthread -O3
 * PKG_CXXFLAGS = -DOPENBLAS -DHAVE_INLINE -std=c++11 -fPIC
 Sys.setenv("PKG_CXXFLAGS"="-DOPENBLAS -DHAVE_INLINE -std=c++11 -fPIC")
 Sys.setenv("PKG_CPPFLAGS" = "-isystem/home/p287664/software/OpenBLAS/include/ -Icontrib/catch-1.9.7 -Isrc")
 Sys.setenv("PKG_LIBS"="-lopenblas -lgsl -lz -pthread -O3")
 Rcpp::sourceCpp("src/internalwrappers.cpp", verbose = T, rebuild = T)
 */



//' R interface to GEMMA's mvlmm.cpp
//'
//' GEMMA, Genome-wide Efficient Mixed Model Association, is an efficient
//' association testing software of Linear Mixed Models (LMM)s and related GWAS
//' models. This is an internal interface for the Multivariate linear mixed
//' model (mvLMM) functionality which also corrects for population structure and
//' sample (non)exchangability
//' @param genoinputs Character vector of input files. It should be either the
//'   prefix name of PLINK binary ped file, or the 3 BIMBAM files' names (the
//'   menan genotypes, the phenotyes and SNP annotaiton files)
//' @param kfile Character variable for the relatedness matrix file name (can be
//'   in gzip compressed format)
//' @param k_mode Integer variable for the type of the kinship/relatedness
//'   matrix type
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
Rcpp::List MvLMM(CharacterVector genoinputs,
                 std::string kfile,
                 NumericVector colnums,
                 NumericMatrix Gmat,
                 int k_mode = 1,
                 double miss = 0.05,
                 double maf = 0.01,
                 double r2 = 0.9999,
                 double hwe = 0,
                 bool notsnp = false,
                 int lmmMode = 1){

  // @param predit Boolean variable to impute missing phenotypes before
  //   association testing (if a small proportion is missing)

  GEMMA cGemma;
  PARAM cPar;

  //gsl_set_error_handler (&gemma_gsl_error_handler);
  gsl_set_error_handler_off();

  // Reading input parameters:

  if (genoinputs.size() == 1){
    cPar.file_bfile = genoinputs[0];
  }
  if (genoinputs.size() >= 2) {
    cPar.file_geno = genoinputs[0];
    cPar.file_pheno = genoinputs[1];
  }
  if (genoinputs.size() == 3)
    cPar.file_anno = genoinputs[2];
  if (genoinputs.size() > 3)
    Rcpp::stop("Incorrect inputs: Please user either PLINK or BIMBAM standard files");

  // set pheno column (list/range)
  (cPar.p_column).clear();
  for (int i = 0; i < colnums.size(); i++) {
    (cPar.p_column).push_back(colnums[i]);
  }

  cPar.file_kin = kfile;
  cPar.file_gxe = ""; cPar.path_out = "rgemma"; cPar.file_out = "mouse_hs1940";

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
  cPar.a_mode = lmmMode;

  cPar.CheckParam();

  //cGemma.BatchRun(cPar);

  if (is_check_mode()) enable_segfpe(); // fast NaN checking by default

  // Read Files.
  Rcout << "Reading Files ... " << endl;
  cPar.ReadFiles();
  cPar.CheckData();

  // LMM or mvLMM or Eigen-Decomposition
  RcppGSL::Matrix Y(cPar.ni_test, cPar.n_ph); // Advantage of RcppGSL types
  RcppGSL::Matrix W(Y->size1, cPar.n_cvt);    // is ease of exchange with R
  RcppGSL::Matrix G(Y->size1, Y->size1);      // compared with gsl_matrix
  gsl_matrix *U = gsl_matrix_safe_alloc(Y->size1, Y->size1);
  gsl_matrix *UtW = gsl_matrix_calloc(Y->size1, W->size2);
  gsl_matrix *UtY = gsl_matrix_calloc(Y->size1, Y->size2);
  gsl_vector *eval = gsl_vector_calloc(Y->size1);
  assert_issue(is_issue(26), UtY->data[0] == 0.0);

  // Remove missing individuals from W, Y, and G:

  // set covariates matrix W and phenotype matrix Y
  // an intercept should be included in W,
  cPar.CopyCvtPhen(W, Y, 0);
  //gsl_matrix *G = gsl_matrix_alloc(Y->size1, Y->size1);

  int row=0; int col=0;
  for (int i=0; i<Gmat.nrow(); i++) {// row index
    if (cPar.indicator_idv[i] == 0) continue;
    for (int j=0; j<Gmat.ncol(); ++j){ // col index
      if (cPar.indicator_idv[j] == 0) continue; // both row and column should not be put in G
      G(row, col) =  Gmat(i,j);
      col++;
    }
    row++;
    col = 0;
  }

  CenterMatrix(G); // center matrix G
  validate_K(G);

  // eigen-decomposition and calculate trace_G - main track
  // Azza: Eignvalues should be passed directly instead. This does
  // not work in this version of GEMMA anyway!
  cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 0);

  // calculate UtW and Uty
  CalcUtX(U, W, UtW);
  CalcUtX(U, Y, UtY);
  assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

  // Fit LMM or mvLMM (w. LOCO)
  Rcout << "fit mvLMM (multiple phenotypes)" << endl;
  MVLMM cMvlmm;
  cMvlmm.CopyFromParam(cPar); // set parameters
  cMvlmm.AnalyzeBimbam2(U, eval, UtW, UtY);
  cMvlmm.WriteFiles();
  cMvlmm.CopyToParam(cPar);

  // release all matrices and vectors
  gsl_matrix_safe_free(U);
  gsl_matrix_safe_free(UtW);
  gsl_matrix_safe_free(UtY);
  gsl_vector_safe_free(eval);

  Rcout << "Time spent on optimization iterations: " <<
    cPar.time_opt << " min " << endl;

  return Rcpp::List::create(Rcpp::Named("l_mle_null") = cPar.l_mle_null,
                            Rcpp::Named("l_remle_null") = cPar.l_remle_null,
                            Rcpp::Named("logl_mle_H0") = cPar.logl_mle_H0,
                            Rcpp::Named("logl_remle_H0") = cPar.logl_remle_H0,
                            Rcpp::Named("Vg_remle_null") = cPar.Vg_remle_null,
                            Rcpp::Named("Ve_remle_null") = cPar.Ve_remle_null,
                            Rcpp::Named("Vg_mle_null") = cPar.Vg_mle_null,
                            Rcpp::Named("Ve_mle_null") = cPar.Ve_mle_null,
                            Rcpp::Named("pve_null") = cPar.pve_null,
                            Rcpp::Named("pve_total") = cPar.pve_total,
                            // Rcpp::Named("pve_se_null") = cPar.pve_se_null,
                            // Rcpp::Named("se_pve_total") = cPar.se_pve_total
                            Rcpp::Named("v_beta") = cPar.v_beta, // REML estimator for beta.
                            Rcpp::Named("beta_remle_null") = cPar.beta_remle_null,
                            Rcpp::Named("beta_mle_null") = cPar.beta_mle_null
                              // Does not work: Rcpp::Named("mvlmm") = cMvlmm.sumStat.data()
                              // Rcpp::Named("se_beta_remle_null") = cPar.se_beta_remle_null
                              // Rcpp::Named("VVg_remle_null") = cPar.VVg_remle_null,
                              // Rcpp::Named("VVe_remle_null") = cPar.VVe_remle_null,
                              // Rcpp::Named("VVe_mle_null") = cPar.VVe_mle_null,
                              // Rcpp::Named("VVg_mle_null") = cPar.VVg_mle_null
  );
}

// for testing and development, this R code will be automatically
// run after the compilation.

/*** R
tictoc::tic()
file.remove("rgemma/mouse_hs1940.assoc.txt")
Gmat <- readRDS(file = "tmpscr/mouse_hs1940.cXX.rds")
mv <- MvLMM(genoinputs = c("/home/p287664/github_repos/GEMMA/example/mouse_hs1940.1snp.txt2",
                           "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.pheno.txt",
                           "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.anno.txt"),
            kfile = "/home/p287664/github_repos/GEMMA/output/mouse_hs1940.cXX.txt",
            Gmat = Gmat,
            colnums = c(1, 6))
tictoc::toc()

*/


//' R interface to 'Genome-wide Efficient Mixed Model Association' (GEMMA)
//'
//' GEMMA is an efficient association testing software of Linear Mixed Models
//' (LMM)s and related GWAS models. This interface is for the Relatedness matrix
//' calculation functionality in centered or standardized matrix form
//' @param gk Integer variable for calculating the relatedness matrix. A value
//'   of 1 (the default) generates a centered matrix. A value of 2 generates a
//'   standardized matrix
//' @inheritParams gemmaMVLMM
//'
//[[Rcpp::export]]
Rcpp::List CalcKin(CharacterVector genoinputs,
                   int gk,
                   NumericVector colnums,
                   double miss = 0.05,
                   double maf = 0.01,
                   double r2 = 0.9999,
                   double hwe = 0,
                   bool notsnp = false,
                   std::string outprefix = "result",
                   std::string outdir = "output",
                   bool license = false){

  GEMMA cGemma;
  PARAM cPar;

  //gsl_set_error_handler (&gemma_gsl_error_handler);
  gsl_set_error_handler_off();

  if (license) {
    cGemma.PrintHeader();
    cGemma.PrintLicense();
    Rcpp::stop("License printed");
  }

  /* Reading input parameters: cGemma.Assign(argc, argv, cPar);
   */

  if (genoinputs.size() == 1){
    cPar.file_bfile = genoinputs[0];
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
    Rcpp::stop("Incorrect inputs: Please user either PLINK or BIMBAM standard files");
  }

  // set pheno column (list/range)
  (cPar.p_column).clear();
  for (int i = 0; i < colnums.size(); i++) {
    (cPar.p_column).push_back(colnums[i]);
  }

  cPar.path_out = outdir;
  cPar.file_out = outprefix;

  // SNP QC Options:
  cPar.miss_level = miss;
  if (cPar.maf_level != -1)
    cPar.maf_level = maf;
  cPar.hwe_level = hwe;
  cPar.r2_level = r2;
  if (notsnp)
    cPar.maf_level = -1;
  cPar.a_mode = 20 + gk;

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

  if (cPar.error == true)
    Rcpp::stop("Error! Can not compute the Relatedness matrix");

  if (is_quiet_mode()) {
    stringstream ss;
    cout.rdbuf(ss.rdbuf());
  }

  cPar.CheckParam();
  if (cPar.error == true)
    Rcpp::stop("Error in checking parameters");

  /* cGemma.BatchRun(cPar); */

  if (is_check_mode()) enable_segfpe(); // fast NaN checking by default

  // Read Files.
  Rcout << "Reading Files ... " << endl;
  cPar.ReadFiles();
  if (cPar.error == true)
    Rcpp::stop("error! failed to read files. ");

  cPar.CheckData();
  if (cPar.error == true)
    Rcpp::stop("error! failed to check data. ");

  Rcout << "Calculating Relatedness Matrix ... " << endl;

  //gsl_matrix *G = gsl_matrix_safe_alloc(cPar.ni_total, cPar.ni_total);
  RcppGSL::matrix<double> G(cPar.ni_total, cPar.ni_total);
  enforce_msg(G, "allocate G"); // just to be sure

  clock_t time_start = clock();

  cPar.CalcKin(G);

  cPar.time_G = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

  if (cPar.error == true)
    Rcpp::stop("error! failed to calculate relatedness matrix. ");

  // Now we have the Kinship matrix test it
  validate_K(G);

  if (cPar.a_mode == M_KIN) {
    cPar.WriteMatrix(G, "cXX");
  } else { // M_KIN2
    cPar.WriteMatrix(G, "sXX");
  }

  Rcout << "Relatedness matrix calculations done in: " <<
    cPar.time_G << " min " << endl;

  return Rcpp::List::create(Rcpp::Named("kinship") = G);
}



/*** R
tictoc::tic()
k = gemmaGK(genoinputs = c("/home/p287664/github_repos/GEMMA/example/mouse_hs1940.geno.txt.gz",
                           "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.pheno.txt",
                           "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.anno.txt"),
            gk = 1, colnums = c(1), outprefix = "mouse_hs1940_k", outdir = "rgemma")
tictoc::toc()
saveRDS(k$kinship, file = "tmpscr/mouse_hs1940.cXX.rds")
str(k)

*/


// [[Rcpp::export]]
int sum_gsl_vector_int(const RcppGSL::vector<int> & vec){
  int res = std::accumulate(vec.begin(), vec.end(), 0);
  return res;
}

// [[Rcpp::export]]
int sum_gsl_matrix_int( RcppGSL::matrix<int> & mat){
  int res = 0;
  for (int i = 0; i < mat.nrow(); i++ )
    for (int j = 0; j < mat.ncol(); j++)
      res = res + mat(i,j);
  return res;
}


// [[Rcpp::export]]
void predictPhenos(){
  // The bulk is in the gemma.BatchRun function, corresponding to cPar.a_mode=43
  // This will need the MVLMM::CalcMvLmmVgVeBeta (for multi-traits) or LMM::CalcLambda
  return;
  }
