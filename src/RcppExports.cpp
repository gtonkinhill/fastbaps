// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bhier
List bhier(List data, List partitions, NumericVector d_k);
RcppExport SEXP _fastbaps_bhier(SEXP dataSEXP, SEXP partitionsSEXP, SEXP d_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type partitions(partitionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d_k(d_kSEXP);
    rcpp_result_gen = Rcpp::wrap(bhier(data, partitions, d_k));
    return rcpp_result_gen;
END_RCPP
}
// bhier_parallel
List bhier_parallel(List data, List partitions, NumericVector d_k, int n_cores);
RcppExport SEXP _fastbaps_bhier_parallel(SEXP dataSEXP, SEXP partitionsSEXP, SEXP d_kSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type partitions(partitionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d_k(d_kSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(bhier_parallel(data, partitions, d_k, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// calc_ddk
List calc_ddk(List data, arma::imat merges);
RcppExport SEXP _fastbaps_calc_ddk(SEXP dataSEXP, SEXP mergesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type merges(mergesSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_ddk(data, merges));
    return rcpp_result_gen;
END_RCPP
}
// import_fasta_to_vector_each_nt
List import_fasta_to_vector_each_nt(std::string file);
RcppExport SEXP _fastbaps_import_fasta_to_vector_each_nt(SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    rcpp_result_gen = Rcpp::wrap(import_fasta_to_vector_each_nt(file));
    return rcpp_result_gen;
END_RCPP
}
// part_llks
List part_llks(List data, List partitions);
RcppExport SEXP _fastbaps_part_llks(SEXP dataSEXP, SEXP partitionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type partitions(partitionsSEXP);
    rcpp_result_gen = Rcpp::wrap(part_llks(data, partitions));
    return rcpp_result_gen;
END_RCPP
}
// summarise_clusters
IntegerVector summarise_clusters(NumericMatrix merge, NumericVector rk, double threshold, int n_isolates);
RcppExport SEXP _fastbaps_summarise_clusters(SEXP mergeSEXP, SEXP rkSEXP, SEXP thresholdSEXP, SEXP n_isolatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rk(rkSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type n_isolates(n_isolatesSEXP);
    rcpp_result_gen = Rcpp::wrap(summarise_clusters(merge, rk, threshold, n_isolates));
    return rcpp_result_gen;
END_RCPP
}
// tree_llk
List tree_llk(List data, arma::imat merges);
RcppExport SEXP _fastbaps_tree_llk(SEXP dataSEXP, SEXP mergesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type merges(mergesSEXP);
    rcpp_result_gen = Rcpp::wrap(tree_llk(data, merges));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastbaps_bhier", (DL_FUNC) &_fastbaps_bhier, 3},
    {"_fastbaps_bhier_parallel", (DL_FUNC) &_fastbaps_bhier_parallel, 4},
    {"_fastbaps_calc_ddk", (DL_FUNC) &_fastbaps_calc_ddk, 2},
    {"_fastbaps_import_fasta_to_vector_each_nt", (DL_FUNC) &_fastbaps_import_fasta_to_vector_each_nt, 1},
    {"_fastbaps_part_llks", (DL_FUNC) &_fastbaps_part_llks, 2},
    {"_fastbaps_summarise_clusters", (DL_FUNC) &_fastbaps_summarise_clusters, 4},
    {"_fastbaps_tree_llk", (DL_FUNC) &_fastbaps_tree_llk, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastbaps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
