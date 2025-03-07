// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// MvLMM
Rcpp::List MvLMM(CharacterVector genoinputs, std::string kfile, NumericVector colnums, NumericMatrix Gmat, int k_mode, double miss, double maf, double r2, double hwe, bool notsnp, int lmmMode);
RcppExport SEXP _pcgen2_MvLMM(SEXP genoinputsSEXP, SEXP kfileSEXP, SEXP colnumsSEXP, SEXP GmatSEXP, SEXP k_modeSEXP, SEXP missSEXP, SEXP mafSEXP, SEXP r2SEXP, SEXP hweSEXP, SEXP notsnpSEXP, SEXP lmmModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type genoinputs(genoinputsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kfile(kfileSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colnums(colnumsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< int >::type k_mode(k_modeSEXP);
    Rcpp::traits::input_parameter< double >::type miss(missSEXP);
    Rcpp::traits::input_parameter< double >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< double >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< double >::type hwe(hweSEXP);
    Rcpp::traits::input_parameter< bool >::type notsnp(notsnpSEXP);
    Rcpp::traits::input_parameter< int >::type lmmMode(lmmModeSEXP);
    rcpp_result_gen = Rcpp::wrap(MvLMM(genoinputs, kfile, colnums, Gmat, k_mode, miss, maf, r2, hwe, notsnp, lmmMode));
    return rcpp_result_gen;
END_RCPP
}
// CalcKin
Rcpp::List CalcKin(CharacterVector genoinputs, int gk, NumericVector colnums, double miss, double maf, double r2, double hwe, bool notsnp, std::string outprefix, std::string outdir, bool license);
RcppExport SEXP _pcgen2_CalcKin(SEXP genoinputsSEXP, SEXP gkSEXP, SEXP colnumsSEXP, SEXP missSEXP, SEXP mafSEXP, SEXP r2SEXP, SEXP hweSEXP, SEXP notsnpSEXP, SEXP outprefixSEXP, SEXP outdirSEXP, SEXP licenseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type genoinputs(genoinputsSEXP);
    Rcpp::traits::input_parameter< int >::type gk(gkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colnums(colnumsSEXP);
    Rcpp::traits::input_parameter< double >::type miss(missSEXP);
    Rcpp::traits::input_parameter< double >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< double >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< double >::type hwe(hweSEXP);
    Rcpp::traits::input_parameter< bool >::type notsnp(notsnpSEXP);
    Rcpp::traits::input_parameter< std::string >::type outprefix(outprefixSEXP);
    Rcpp::traits::input_parameter< std::string >::type outdir(outdirSEXP);
    Rcpp::traits::input_parameter< bool >::type license(licenseSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcKin(genoinputs, gk, colnums, miss, maf, r2, hwe, notsnp, outprefix, outdir, license));
    return rcpp_result_gen;
END_RCPP
}
// sum_gsl_vector_int
int sum_gsl_vector_int(const RcppGSL::vector<int>& vec);
RcppExport SEXP _pcgen2_sum_gsl_vector_int(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RcppGSL::vector<int>& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_gsl_vector_int(vec));
    return rcpp_result_gen;
END_RCPP
}
// sum_gsl_matrix_int
int sum_gsl_matrix_int(RcppGSL::matrix<int>& mat);
RcppExport SEXP _pcgen2_sum_gsl_matrix_int(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RcppGSL::matrix<int>& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_gsl_matrix_int(mat));
    return rcpp_result_gen;
END_RCPP
}
// predictPhenos
void predictPhenos();
RcppExport SEXP _pcgen2_predictPhenos() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    predictPhenos();
    return R_NilValue;
END_RCPP
}
// gemmaMVLMM
Rcpp::List gemmaMVLMM(CharacterVector genoinputs, std::string kfile, NumericVector colnums, int k_mode, double miss, double maf, double r2, double hwe, bool notsnp, int lmmMode, std::string gxe, std::string outprefix, std::string outdir, bool license);
RcppExport SEXP _pcgen2_gemmaMVLMM(SEXP genoinputsSEXP, SEXP kfileSEXP, SEXP colnumsSEXP, SEXP k_modeSEXP, SEXP missSEXP, SEXP mafSEXP, SEXP r2SEXP, SEXP hweSEXP, SEXP notsnpSEXP, SEXP lmmModeSEXP, SEXP gxeSEXP, SEXP outprefixSEXP, SEXP outdirSEXP, SEXP licenseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type genoinputs(genoinputsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kfile(kfileSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colnums(colnumsSEXP);
    Rcpp::traits::input_parameter< int >::type k_mode(k_modeSEXP);
    Rcpp::traits::input_parameter< double >::type miss(missSEXP);
    Rcpp::traits::input_parameter< double >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< double >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< double >::type hwe(hweSEXP);
    Rcpp::traits::input_parameter< bool >::type notsnp(notsnpSEXP);
    Rcpp::traits::input_parameter< int >::type lmmMode(lmmModeSEXP);
    Rcpp::traits::input_parameter< std::string >::type gxe(gxeSEXP);
    Rcpp::traits::input_parameter< std::string >::type outprefix(outprefixSEXP);
    Rcpp::traits::input_parameter< std::string >::type outdir(outdirSEXP);
    Rcpp::traits::input_parameter< bool >::type license(licenseSEXP);
    rcpp_result_gen = Rcpp::wrap(gemmaMVLMM(genoinputs, kfile, colnums, k_mode, miss, maf, r2, hwe, notsnp, lmmMode, gxe, outprefix, outdir, license));
    return rcpp_result_gen;
END_RCPP
}
// gemmaGK
Rcpp::List gemmaGK(CharacterVector genoinputs, int gk, NumericVector colnums, double miss, double maf, double r2, double hwe, bool notsnp, std::string outprefix, std::string outdir, bool license);
RcppExport SEXP _pcgen2_gemmaGK(SEXP genoinputsSEXP, SEXP gkSEXP, SEXP colnumsSEXP, SEXP missSEXP, SEXP mafSEXP, SEXP r2SEXP, SEXP hweSEXP, SEXP notsnpSEXP, SEXP outprefixSEXP, SEXP outdirSEXP, SEXP licenseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type genoinputs(genoinputsSEXP);
    Rcpp::traits::input_parameter< int >::type gk(gkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colnums(colnumsSEXP);
    Rcpp::traits::input_parameter< double >::type miss(missSEXP);
    Rcpp::traits::input_parameter< double >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< double >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< double >::type hwe(hweSEXP);
    Rcpp::traits::input_parameter< bool >::type notsnp(notsnpSEXP);
    Rcpp::traits::input_parameter< std::string >::type outprefix(outprefixSEXP);
    Rcpp::traits::input_parameter< std::string >::type outdir(outdirSEXP);
    Rcpp::traits::input_parameter< bool >::type license(licenseSEXP);
    rcpp_result_gen = Rcpp::wrap(gemmaGK(genoinputs, gk, colnums, miss, maf, r2, hwe, notsnp, outprefix, outdir, license));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pcgen2_MvLMM", (DL_FUNC) &_pcgen2_MvLMM, 11},
    {"_pcgen2_CalcKin", (DL_FUNC) &_pcgen2_CalcKin, 11},
    {"_pcgen2_sum_gsl_vector_int", (DL_FUNC) &_pcgen2_sum_gsl_vector_int, 1},
    {"_pcgen2_sum_gsl_matrix_int", (DL_FUNC) &_pcgen2_sum_gsl_matrix_int, 1},
    {"_pcgen2_predictPhenos", (DL_FUNC) &_pcgen2_predictPhenos, 0},
    {"_pcgen2_gemmaMVLMM", (DL_FUNC) &_pcgen2_gemmaMVLMM, 14},
    {"_pcgen2_gemmaGK", (DL_FUNC) &_pcgen2_gemmaGK, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_pcgen2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
