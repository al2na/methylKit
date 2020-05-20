// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// methCall
void methCall(std::string read1, std::string type, bool nolap, int minqual, int mincov, bool phred64, std::string CpGfile, std::string CHHfile, std::string CHGfile, size_t verbosity);
RcppExport SEXP _methylKit_methCall(SEXP read1SEXP, SEXP typeSEXP, SEXP nolapSEXP, SEXP minqualSEXP, SEXP mincovSEXP, SEXP phred64SEXP, SEXP CpGfileSEXP, SEXP CHHfileSEXP, SEXP CHGfileSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type read1(read1SEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type nolap(nolapSEXP);
    Rcpp::traits::input_parameter< int >::type minqual(minqualSEXP);
    Rcpp::traits::input_parameter< int >::type mincov(mincovSEXP);
    Rcpp::traits::input_parameter< bool >::type phred64(phred64SEXP);
    Rcpp::traits::input_parameter< std::string >::type CpGfile(CpGfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type CHHfile(CHHfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type CHGfile(CHGfileSEXP);
    Rcpp::traits::input_parameter< size_t >::type verbosity(verbositySEXP);
    methCall(read1, type, nolap, minqual, mincov, phred64, CpGfile, CHHfile, CHGfile, verbosity);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_methylKit_methCall", (DL_FUNC) &_methylKit_methCall, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_methylKit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
