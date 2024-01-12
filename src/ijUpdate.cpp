#include <Rcpp.h>
using namespace Rcpp;

// Define the C++ function
// [[Rcpp::export]]
NumericVector ijUpdate(DataFrame ijDat, NumericVector newPi, DataFrame iDatNew) {
  Environment skipTrack("package:skipTrack") ;
  Function postCij = skipTrack["postCij"] ;

  // Extract relevant columns from DataFrame
  NumericVector ys = ijDat["ys"];
  IntegerVector Individual = ijDat["Individual"];
  NumericVector cs = ijDat["cs"];
  NumericVector mus = iDatNew["mus"];
  NumericVector taus = iDatNew["taus"];

  // Loop through each row and update cs
  for (int k = 0; k < ys.size(); ++k) {
    Rcpp::LogicalVector subsetCondition = (Individual == Individual[k]);
    SEXP result = postCij(ys[k], newPi, mus[Individual == subsetCondition], taus[Individual == subsetCondition]);

    // Convert the Vector to NumericVector
    NumericVector resultVec(result);

    // Assign individual values to cs
    cs[k] = resultVec[0];  // Adjust this if the result is a vector
  }

  return cs;
}
