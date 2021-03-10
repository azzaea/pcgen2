#include <Rcpp.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include "lapack.h"

#include "gemma.h"
#include "gemma_io.h"
#include "lmm.h"
#include "mathfunc.h"
#include "mvlmm.h"

// #include <fstream>
// #include <iostream>
// #include <sstream>
// #include <vector>

#include <sys/stat.h>
#include <sys/types.h>


#define STRICT_R_HEADERS

using namespace Rcpp;
using namespace std;

// source using Rcpp::sourceCpp
/* PKG_CPPFLAGS = -isystem/home/p287664/software/OpenBLAS/include/ -Icontrib/catch-1.9.7 -Isrc
 * PKG_LIBS = -lopenblas -lgsl -lz -pthread -O3
 * PKG_CXXFLAGS = -DOPENBLAS -DHAVE_INLINE -std=c++11 -fPIC
 * Sys.setenv("PKG_CXXFLAGS"="-DOPENBLAS -DHAVE_INLINE -std=c++11 -fPIC")
 * Sys.setenv("PKG_CPPFLAGS" = "-isystem/home/p287664/software/OpenBLAS/include/ -Icontrib/catch-1.9.7 -Isrc")
 * Sys.setenv("PKG_LIBS"="-lopenblas -lgsl -lz -pthread -O3")
 * Rcpp::sourceCpp("src/wrappers.cpp", verbose = T, rebuild = T)
 */



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
bool gemmaMVLMM(CharacterVector genoinputs,
                std::string kfile,
                NumericVector colnums,
                int k_mode = 1,
                double miss = 0.05,
                double maf = 0.01,
                double r2 = 0.9999,
                double hwe = 0,
                bool notsnp = false,
                int lmmMode = 1,
                std::string gxe = "",
                std::string outprefix = "result",
                std::string outdir = "output",
                bool license = false){

  // @param predit Boolean variable to impute missing phenotypes before
  //   association testing (if a small proportion is missing)

  // int main(int argc, char *argv[])

  GEMMA cGemma;
  PARAM cPar;
  clock_t time_begin, time_start;
  time_begin = clock();

  gsl_set_error_handler (&gemma_gsl_error_handler);

  if (license) {
    cGemma.PrintHeader();
    cGemma.PrintLicense();
    return EXIT_SUCCESS;
  }

  // Reading input parameters:
   // cGemma.Assign(argc, argv, cPar);
   //

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
    return EXIT_SUCCESS;
  }

  // set pheno column (list/range)
  (cPar.p_column).clear();
  for (int i = 0; i < colnums.size(); i++) {
    (cPar.p_column).push_back(colnums[i]);
  }

  cPar.file_kin = kfile;
  cPar.file_gxe = gxe;
  cPar.path_out = outdir;
  cPar.file_out = outprefix;

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
    Rcout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
    "-predict -calccor options is allowed."
    << endl;
  }
  cPar.a_mode = lmmMode;

  // Rcout << "-bfile: " << cPar.file_bfile << "\n";
  // Rcout << "-g: " << cPar.file_geno << "\n";
  // Rcout << "-p: " << cPar.file_pheno << "\n";
  // Rcout << "-a: " << cPar.file_anno << "\n";
  // Rcout << "-n: " << cPar.p_column[i] << "\n";
  // Rcout << "-gxe: " << cPar.file_gxe << "\n";
  // Rcout << "-k: " << cPar.file_kin << "\n";
  // Rcout << "-outdir:" << cPar.path_out << "\n";
  // Rcout << "-o: " << cPar.file_out << "\n";
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
    Rcout << "Error in checking parameters \n" ;
    return EXIT_FAILURE;
  }

  // The real work:
   //
  //cGemma.BatchRun(cPar);

  if (is_check_mode()) enable_segfpe(); // fast NaN checking by default

  // Read Files.
  Rcout << "Reading Files ... " << endl;
  cPar.ReadFiles();
  if (cPar.error == true) {
    Rcout << "error! fail to read files. " << endl;
    return EXIT_FAILURE;
  }
  cPar.CheckData();
  if (cPar.error == true) {
    Rcout << "error! fail to check data. " << endl;
    return EXIT_FAILURE;
  }






  // LMM or mvLMM or Eigen-Decomposition
  if (cPar.a_mode == M_LMM1 || cPar.a_mode == M_LMM2 || cPar.a_mode == M_LMM3 ||
      cPar.a_mode == M_LMM4 || cPar.a_mode == M_LMM5 || cPar.a_mode == M_LMM9 ||
      cPar.a_mode == M_EIGEN  ) { // Fit LMM or mvLMM or eigen
    gsl_matrix *Y = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_ph);
    enforce_msg(Y, "allocate Y"); // just to be sure
    gsl_matrix *W = gsl_matrix_safe_alloc(Y->size1, cPar.n_cvt);
    gsl_matrix *B = gsl_matrix_safe_alloc(Y->size2, W->size2); // B is a d by c
    // matrix
    gsl_matrix *se_B = gsl_matrix_safe_alloc(Y->size2, W->size2);
    gsl_matrix *G = gsl_matrix_safe_alloc(Y->size1, Y->size1);
    gsl_matrix *U = gsl_matrix_safe_alloc(Y->size1, Y->size1);
    gsl_matrix *UtW = gsl_matrix_calloc(Y->size1, W->size2);
    gsl_matrix *UtY = gsl_matrix_calloc(Y->size1, Y->size2);
    gsl_vector *eval = gsl_vector_calloc(Y->size1);
    gsl_vector *env = gsl_vector_safe_alloc(Y->size1);
    gsl_vector *weight = gsl_vector_safe_alloc(Y->size1);
    debug_msg("Started on LMM");
    assert_issue(is_issue(26), UtY->data[0] == 0.0);

    // set covariates matrix W and phenotype matrix Y
    // an intercept should be included in W,
    cPar.CopyCvtPhen(W, Y, 0);
    if (!cPar.file_gxe.empty()) {
      cPar.CopyGxe(env);
    }

    // read relatedness matrix G
    if (!(cPar.file_kin).empty()) {
      ReadFile_kin(cPar.file_kin, cPar.indicator_idv, cPar.mapID2num,
                   cPar.k_mode, cPar.error, G);
      debug_msg("Read K/GRM file");
      if (cPar.error == true) {
        Rcout << "error! fail to read kinship/relatedness file. " << endl;
        return EXIT_FAILURE;
      }

      // center matrix G
      CenterMatrix(G);
      validate_K(G);

      // is residual weights are provided, then
      if (!cPar.file_weight.empty()) {
        cPar.CopyWeight(weight);
        double d, wi, wj;
        for (size_t i = 0; i < G->size1; i++) {
          wi = gsl_vector_get(weight, i);
          for (size_t j = i; j < G->size2; j++) {
            wj = gsl_vector_get(weight, j);
            d = gsl_matrix_get(G, i, j);
            if (wi <= 0 || wj <= 0) {
              d = 0;
            } else {
              d /= safe_sqrt(wi * wj);
            }
            gsl_matrix_set(G, i, j, d);
            if (j != i) {
              gsl_matrix_set(G, j, i, d);
            }
          }
        }
      }

      // eigen-decomposition and calculate trace_G - main track
      Rcout << "Start Eigen-Decomposition..." << endl;
      time_start = clock();

      if (cPar.a_mode == M_EIGEN) {
        cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 1);
      } else {
        cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 0);
      }
      // write(eval,"eval");

      if (!cPar.file_weight.empty()) {
        double wi;
        for (size_t i = 0; i < U->size1; i++) {
          wi = gsl_vector_get(weight, i);
          if (wi <= 0) {
            wi = 0;
          } else {
            wi = safe_sqrt(wi);
          }
          gsl_vector_view Urow = gsl_matrix_row(U, i);
          gsl_vector_scale(&Urow.vector, wi);
        }
      }

      cPar.time_eigen =
        (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    } else {
      ReadFile_eigenU(cPar.file_ku, cPar.error, U);
      if (cPar.error == true) {
        Rcout << "error! fail to read the U file. " << endl;
        return EXIT_FAILURE;
      }

      ReadFile_eigenD(cPar.file_kd, cPar.error, eval);
      if (cPar.error == true) {
        Rcout << "error! fail to read the D file. " << endl;
        return EXIT_FAILURE;
      }

      cPar.trace_G = 0.0;
      for (size_t i = 0; i < eval->size; i++) {
        if (gsl_vector_get(eval, i) < 1e-10) {
          gsl_vector_set(eval, i, 0);
        }
        cPar.trace_G += gsl_vector_get(eval, i);
      }
      cPar.trace_G /= (double)eval->size;
    }
    // write(eval,"eval2");

    if (cPar.a_mode == M_EIGEN) {
      cPar.WriteMatrix(U, "eigenU");
      cPar.WriteVector(eval, "eigenD");
    } else if (!cPar.file_gene.empty()) { // Run with gene file
      // calculate UtW and Uty
      CalcUtX(U, W, UtW);
      CalcUtX(U, Y, UtY);

      assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

      LMM cLmm;
      cLmm.CopyFromParam(cPar);

      gsl_vector_view Y_col = gsl_matrix_column(Y, 0);
      gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

      cLmm.AnalyzeGene(U, eval, UtW, &UtY_col.vector, W,
                       &Y_col.vector); // y is the predictor, not the phenotype

      cLmm.WriteFiles();
      cLmm.CopyToParam(cPar);
    } else {
      debug_msg("Main LMM track");
      // calculate UtW and Uty
      CalcUtX(U, W, UtW);
      CalcUtX(U, Y, UtY);
      assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

      // calculate REMLE/MLE estimate and pve for univariate model
      if (cPar.n_ph == 1) { // one phenotype
        gsl_vector_view beta = gsl_matrix_row(B, 0);
        gsl_vector_view se_beta = gsl_matrix_row(se_B, 0);
        gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

        assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

        CalcLambda('L', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max,
                   cPar.n_region, cPar.l_mle_null, cPar.logl_mle_H0);
        //assert(!isnan(UtY->data[0]));
        if (ISNAN(UtY->data[0]))
          Rcpp::stop("Nan detected: UtY->data[0]");

        CalcLmmVgVeBeta(eval, UtW, &UtY_col.vector, cPar.l_mle_null,
                        cPar.vg_mle_null, cPar.ve_mle_null, &beta.vector,
                        &se_beta.vector);

        //assert(!isnan(UtY->data[0]));
        if (ISNAN(UtY->data[0]))
          Rcpp::stop("Nan detected: UtY->data[0]");

        cPar.beta_mle_null.clear();
        cPar.se_beta_mle_null.clear();
        //assert(!isnan(B->data[0]));
        if (ISNAN(B->data[0]))
          Rcpp::stop("Nan detected: B->data[0]");

        //assert(!isnan(se_B->data[0]));
        if (ISNAN(se_B->data[0]))
          Rcpp::stop("Nan detected: se_B->data[0]");

        for (size_t i = 0; i < B->size2; i++) {
          cPar.beta_mle_null.push_back(gsl_matrix_get(B, 0, i));
          cPar.se_beta_mle_null.push_back(gsl_matrix_get(se_B, 0, i));
        }
        //assert(!isnan(UtY->data[0]));
        if (ISNAN(UtY->data[0]))
          Rcpp::stop("Nan detected: UtY->data[0]");
        // assert(!isnan(cPar.beta_mle_null.front()));
        if (ISNAN(cPar.beta_mle_null.front()))
          Rcpp::stop("Nan detected: cPar.beta_mle_null.front()");
        //assert(!isnan(cPar.se_beta_mle_null.front()));
        if (ISNAN(cPar.se_beta_mle_null.front()))
          Rcpp::stop("Nan detected: cPar.se_beta_mle_null.front()");

        // the following functions do not modify eval
        CalcLambda('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max,
                   cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
        CalcLmmVgVeBeta(eval, UtW, &UtY_col.vector, cPar.l_remle_null,
                        cPar.vg_remle_null, cPar.ve_remle_null, &beta.vector,
                        &se_beta.vector);

        cPar.beta_remle_null.clear();
        cPar.se_beta_remle_null.clear();
        //assert(!isnan(B->data[0]));
        if (ISNAN(B->data[0]))
          Rcpp::stop("Nan detected: B->data[0]");
        //assert(!isnan(se_B->data[0]));
        if (ISNAN(se_B->data[0]))
          Rcpp::stop("Nan detected: se_B->data[0]");

        for (size_t i = 0; i < B->size2; i++) {
          cPar.beta_remle_null.push_back(gsl_matrix_get(B, 0, i));
          cPar.se_beta_remle_null.push_back(gsl_matrix_get(se_B, 0, i));
        }

        CalcPve(eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.trace_G,
                cPar.pve_null, cPar.pve_se_null);
        debug_msg("main print summary");
        cPar.PrintSummary();

        // calculate and output residuals
        if (cPar.a_mode == M_LMM5) {
          gsl_vector *Utu_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *Ute_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *u_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *e_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *y_hat = gsl_vector_safe_alloc(Y->size1);

          // obtain Utu and Ute
          gsl_vector_safe_memcpy(y_hat, &UtY_col.vector);
          gsl_blas_dgemv(CblasNoTrans, -1.0, UtW, &beta.vector, 1.0, y_hat);

          double d, u, e;
          for (size_t i = 0; i < eval->size; i++) {
            d = gsl_vector_get(eval, i);
            u = cPar.l_remle_null * d / (cPar.l_remle_null * d + 1.0) *
              gsl_vector_get(y_hat, i);
            e = 1.0 / (cPar.l_remle_null * d + 1.0) * gsl_vector_get(y_hat, i);
            gsl_vector_set(Utu_hat, i, u);
            gsl_vector_set(Ute_hat, i, e);
          }

          // obtain u and e
          gsl_blas_dgemv(CblasNoTrans, 1.0, U, Utu_hat, 0.0, u_hat);
          gsl_blas_dgemv(CblasNoTrans, 1.0, U, Ute_hat, 0.0, e_hat);

          // output residuals
          cPar.WriteVector(u_hat, "residU");
          cPar.WriteVector(e_hat, "residE");

          gsl_vector_safe_free(u_hat);
          gsl_vector_safe_free(e_hat);
          gsl_vector_safe_free(y_hat);
        } // output residuals
      }

      // Fit LMM or mvLMM (w. LOCO)
      if (cPar.a_mode == M_LMM1 || cPar.a_mode == M_LMM2 || cPar.a_mode == M_LMM3 ||
          cPar.a_mode == M_LMM4 || cPar.a_mode == M_LMM9) {
        if (cPar.n_ph == 1) {
          debug_msg("fit LMM (one phenotype)");
          LMM cLmm;
          cLmm.CopyFromParam(cPar); // set parameters

          // if (is_check_mode()) disable_segfpe(); // disable fast NaN checking for now

          gsl_vector_view Y_col = gsl_matrix_column(Y, 0);
          gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

          if (!cPar.file_bfile.empty()) {
            // PLINK analysis
            if (cPar.file_gxe.empty()) {
              cLmm.AnalyzePlink(U, eval, UtW, &UtY_col.vector, W,
                                &Y_col.vector, cPar.setGWASnps);
            }
            else {
              cLmm.AnalyzePlinkGXE(U, eval, UtW, &UtY_col.vector, W,
                                   &Y_col.vector, env);
            }
          }
          else {
            // BIMBAM analysis

            if (cPar.file_gxe.empty()) {
              cLmm.AnalyzeBimbam(U, eval, UtW, &UtY_col.vector, W,
                                 &Y_col.vector, cPar.setGWASnps);
            } else {
              cLmm.AnalyzeBimbamGXE(U, eval, UtW, &UtY_col.vector, W,
                                    &Y_col.vector, env);
            }
          }
          cLmm.WriteFiles();
          cLmm.CopyToParam(cPar);
        } else {
          debug_msg("fit mvLMM (multiple phenotypes)");
          MVLMM cMvlmm;
          cMvlmm.CopyFromParam(cPar); // set parameters

          // if (is_check_mode()) disable_segfpe(); // disable fast NaN checking

          // write(eval,"eval3");

          if (!cPar.file_bfile.empty()) {
            if (cPar.file_gxe.empty()) {
              cMvlmm.AnalyzePlink(U, eval, UtW, UtY);
            } else {
              cMvlmm.AnalyzePlinkGXE(U, eval, UtW, UtY, env);
            }
          } else {
            if (cPar.file_gxe.empty()) {
              cMvlmm.AnalyzeBimbam(U, eval, UtW, UtY);
            } else {
              cMvlmm.AnalyzeBimbamGXE(U, eval, UtW, UtY, env);
            }
          }

          cMvlmm.WriteFiles();
          cMvlmm.CopyToParam(cPar);
        }
      }
    }

    // release all matrices and vectors
    gsl_matrix_safe_free(Y);
    gsl_matrix_safe_free(W);
    gsl_matrix_warn_free(B); // sometimes unused
    gsl_matrix_warn_free(se_B);
    gsl_matrix_warn_free(G);
    gsl_matrix_safe_free(U);
    gsl_matrix_safe_free(UtW);
    gsl_matrix_safe_free(UtY);
    gsl_vector_safe_free(eval);
    gsl_vector_free(env); // sometimes unused


  }


/////////
  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  Rcout << "Exit: " << EXIT_SUCCESS ;

  return EXIT_SUCCESS;
}

// for testing and development, this R code will be automatically
// run after the compilation.

/*** R
tictoc::tic()
gemmaMVLMM(genoinputs = c("/home/p287664/github_repos/GEMMA/example/mouse_hs1940.geno.txt.gz",
                          "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.pheno.txt",
                          "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.anno.txt"),
           kfile = "/home/p287664/github_repos/GEMMA/output/mouse_hs1940.cXX.txt",
           colnums = c(1, 6), outprefix = "mouse_hs1940_CD8MCH_lmm", outdir = "rgemma")
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
bool gemmaGK(CharacterVector genoinputs,
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
  clock_t time_start;

  gsl_set_error_handler (&gemma_gsl_error_handler);

  if (license) {
    cGemma.PrintHeader();
    cGemma.PrintLicense();
    return EXIT_SUCCESS;
  }

  /* Reading input parameters:
   * cGemma.Assign(argc, argv, cPar);
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
    return EXIT_SUCCESS;
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
  if (cPar.maf_level != -1) {
    cPar.maf_level = maf;
  }
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

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  if (is_quiet_mode()) {
    stringstream ss;
    cout.rdbuf(ss.rdbuf());
  }

  cPar.CheckParam();
  if (cPar.error == true) {
    Rcout << "Error in checking parameters \n" ;
    return EXIT_FAILURE;
  }

  /* The real work:
  cGemma.BatchRun(cPar); */

  if (is_check_mode()) enable_segfpe(); // fast NaN checking by default

  // Read Files.
  Rcout << "Reading Files ... " << endl;
  cPar.ReadFiles();
  if (cPar.error == true) {
    Rcout << "error! fail to read files. " << endl;
    return EXIT_FAILURE;
  }
  cPar.CheckData();
  if (cPar.error == true) {
    Rcout << "error! fail to check data. " << endl;
    return EXIT_FAILURE;
  }

  Rcout << "Calculating Relatedness Matrix ... " << endl;

  gsl_matrix *G = gsl_matrix_safe_alloc(cPar.ni_total, cPar.ni_total);
  enforce_msg(G, "allocate G"); // just to be sure

  time_start = clock();

  cPar.CalcKin(G);

  cPar.time_G = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
  if (cPar.error == true) {
    Rcout << "error! fail to calculate relatedness matrix. " << endl;
    return EXIT_FAILURE;
  }

  // Now we have the Kinship matrix test it
  validate_K(G);

  if (cPar.a_mode == M_KIN) {
    cPar.WriteMatrix(G, "cXX");
  } else { // M_KIN2
    cPar.WriteMatrix(G, "sXX");
  }

  gsl_matrix_safe_free(G);

  Rcout << "Relatedness matrix calculations done. \n"  ;

  Rcout << "##      time on calculating relatedness matrix = "
          << cPar.time_G << " min " << endl;

  return true;

}



/*** R
tictoc::tic()
gemmaGK(genoinputs = c("/home/p287664/github_repos/GEMMA/example/mouse_hs1940.geno.txt.gz",
                       "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.pheno.txt",
                       "/home/p287664/github_repos/GEMMA/example/mouse_hs1940.anno.txt"),
        gk = 1, colnums = c(1), outprefix = "mouse_hs1940_k", outdir = "rgemma")
tictoc::toc()

*/
