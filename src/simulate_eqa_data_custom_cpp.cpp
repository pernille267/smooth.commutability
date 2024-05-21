#include <Rcpp.h>
using namespace Rcpp;

//' Simulate EQA data according to custom non-linear functions
//'
//' @title Simulate EQA data according to custom non-linear functions
//' @name simulate_eqa_data_custom_cpp
//'
//' @param parameters A \code{list} containing the simulation parameters.
//' @param type An \code{integer} value between 1 and 3. Which custom function should be simulated from.
//' @param AR A non-missing \code{logical} value. Set to \code{TRUE} if all replicated measurements should be returned. Otherwise the mean of replicates are returned
//'
//' @description Simulates EQA data according to particular non-linear functions
//'
//' @details This function has no other use than in simulation studies.
//'
//' @return A \code{list} containing the simulated EQA data for the specified \code{type}.
//'
//' @examples \dontrun{
//'
//'   # Generate some fictional data
//'   simulate_eqa_data_custom_cpp(list(n = 25, R = 3, cil = 2, ciu = 10,
//'    cvx = 0.02, cvy = 0.03, cve = 0), type = 2, AR = FALSE)
//'
//' }



// [[Rcpp::export]]
List simulate_eqa_data_custom_cpp(List parameters, int type = 1, bool AR = false) {
  int n = parameters["n"];
  double cil = parameters["cil"];
  double ciu = parameters["ciu"];
  double cvx = parameters["cvx"];
  double cvy = parameters["cvy"];
  double cve = parameters["cve"];
  int R = parameters["R"];

  NumericVector x_L = runif(n, cil, ciu);
  double mu = (cil + ciu) / 2.0;
  double sdx = cvx * mu / sqrt(R);
  double sdy = cvy * mu / sqrt(R);
  double sde = cve * mu / sqrt(R);
  NumericVector f_kL(n);

  if (type == 1) {
    f_kL = x_L + 0.9 * sin(0.4 * pow(x_L, 1.06));
  } else if (type == 2) {
    f_kL = x_L + 0.05 * exp(0.16 * pow(x_L, 1.35));
  } else if (type == 3) {
    f_kL = x_L - exp(-0.5 * pow((x_L - 1.5)/2, 2));
  } else {
    stop("ERROR");
  }

  if (!AR) {
    NumericVector y = f_kL + rnorm(n, 0, sdy) + rnorm(n, 0, sde);
    NumericVector x = x_L + rnorm(n, 0, sdx);
    return List::create(Named("SampleID") = seq_len(n),
                        Named("MP_A") = y,
                        Named("MP_B") = x);
  } else {
    sdx *= sqrt(R);
    sdy *= sqrt(R);
    sde *= sqrt(R);
    NumericVector y(n * R);
    NumericVector x(n * R);
    IntegerVector r(n * R);
    IntegerVector s(n * R);

    for (int i = 0; i < n; ++i) {
      int start_index = R * i;
      for (int j = 0; j < R; ++j) {
        y[start_index + j] = rnorm(1, f_kL[i], sdy)[0] + rnorm(1, 0, sde)[0];
        x[start_index + j] = rnorm(1, x_L[i], sdx)[0];
        r[start_index + j] = j + 1;
        s[start_index + j] = i + 1;
      }
    }
    return List::create(Named("SampleID") = s,
                        Named("ReplicateID") = r,
                        Named("MP_A") = y,
                        Named("MP_B") = x);
  }
}
