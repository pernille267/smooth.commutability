#include <Rcpp.h>
using namespace Rcpp;


//' Reconstruct a 4-banded matrix based on flattened data
//'
//' @title Reconstruct a 4-banded matrix based on flattened data
//' @name reconstruct_4_band_matrix
//'
//' @param x A \code{numeric} vector signifying the flattened version of the 4-banded matrix.
//'
//' @description Reconstructs a 4-banded matrix based on flattened data.
//'
//' @details This reconstruction assumes symmetrical 4-banded matrix. End-users should not bother with this function.
//'
//' @return A 4-banded \code{matrix}.
//'
//' @examples \dontrun{
//'   print(1)
//' }
// [[Rcpp::export]]
NumericMatrix reconstruct_4_band_matrix(NumericVector x) {

  int n = x.size() / 4;
  NumericMatrix output(n, n);
  int count = 0;

  // Diagonal
  for(int i = 0; i < n; ++i){
    output(i, i) = x[count];
    ++count;
  }

  // Band
  for(int j = 1; j < 4; ++j){
    for(int i = j; i < n; ++i){
      output(i, i - j) = x[count];
      output(i - j, i) = output(i, i - j);
      ++count;
    }
    count += j;
  }
  return output;
}

//' Reconstruct matrices based on smooth.spline output
//'
//' @title Reconstruct matrices based on smooth.spline output
//' @name get_matrices
//'
//' @param auxM A \code{List} containing the \code{auxM} object from the smooth.spline fit.
//'
//' @description Reconstructs the penalty and B-spline matrices based on flattened data.
//'
//' @details This reconstruction assumes symmetrical 4-banded matrices. End-users should not bother with this function.
//'
//' @return A list of two 4-banded \code{matrix} objects and one \code{numeric} vector.
//'
//' @examples \dontrun{
//'   print(1)
//' }
// [[Rcpp::export]]
List get_matrices(List auxM){
  NumericVector Omega_flattened = auxM["Sigma"];
  NumericVector BTB_flattened = auxM["XWX"];
  NumericVector By_flattened = auxM["XWy"];
  NumericMatrix Omega = reconstruct_4_band_matrix(Omega_flattened);
  NumericMatrix BTB = reconstruct_4_band_matrix(BTB_flattened);
  List out = List::create(Named("BTy") = By_flattened, Named("BTB") = BTB, Named("Omega") = Omega);
  return out;
}

