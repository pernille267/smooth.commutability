#include <Rcpp.h>
#include <cmath>

//' Estimate \eqn{\mathrm{df}^2(c)}
//'
//' @title Estimate \eqn{\mathrm{df}^2(c)}
//' @name df_2
//'
//' @param df_grid A \code{numeric} vector containing the candidate degrees of freedom for smoothing splines. Missing values are allowed.
//' @param cv A \code{numeric} vector contaning cross-validation values for various smoothing splines fit associated with the degrees of freedom found in \code{df_grid}. Missing values are allowed.
//' @param df_0 A \code{double} between 2 and the maximum value of \code{df_grid} representing the traditional optimal degrees of freedom. This is typically the minimized of the cross-validation curve.
//' @param D A positive \code{double} representing the largest allowable angle of the slopes of the cross-validation curve.
//' @param tol A non-negative \code{double} representing the maximum allowable difference between the angle of each slope of the cross-validation curve and \code{D}.
//' @param silence A \code{logical} value. If set to \code{TRUE}, informative messages are printed to the console regarding the iterative process of choosing \eqn{\mathrm{df}^2(c)}.
//'
//' @description Estimate a more conservative degrees of freedom, \eqn{\mathrm{df}^2(c)}, through an iterative process: starting from the degress of freedom immediately to the left of \eqn{\mathrm{df}^0(c)} on the degrees of freedom grid, \code{df_grid}, we examine its associated angle. If this angle does not exceed the threshold, \code{D}, we proceed to the next closest degrees of freedom and evaluate its angle. This process is repeated, moving leftward along the grid, until we identify the first angle that surpasses the threshold. The value of \eqn{\mathrm{df}^2(c)} is then assigned as the degrees of freedom corresponding to the position immediately preceding the one where the angle first exceeded \code{D}.
//'
//' @details Caution is required. With a limited length of \code{df_grid}, it is possible that no angle fall below \code{D}. Even if the length of \code{df_grid}, it is not assured that there exists a cross-validation curve slope smaller than \code{D}. Therefore, ensuring that that \code{D} is not set too low, will increase the likelihood of obtaining a meaningful value of \eqn{\mathrm{df}^2(c)}. If the iterative process fail to find a suitable value of \eqn{\mathrm{df}^2(c)}, the input value of \code{df_0} will be returned.
//'
//' @return A \code{double} representing the resulting value of \eqn{\mathrm{df}^2(c)}.
//'
//' @examples \dontrun{
//'
//'   # Generate some fictional data
//'   test_df_grid <- seq(from = 2, to = 12, length.out = 100)
//'   test_cv <- 0.5 + 0.2 * (df_grid - 5)**2
//'   test_df_0 <- 5
//'
//'   # Using D = 2.5
//'   df_2(test_df_grid, test_cv, test_df_0, D = 2.5, tol = 0)
//'   df_2(test_df_grid, test_cv, test_df_0, D = 2.5, tol = 1/4)
//'
//'   # Using D = 5.0
//'   df_2(test_df_grid, test_cv, test_df_0, D = 5, tol = 0)
//'   df_2(test_df_grid, test_cv, test_df_0, D = 5, tol = 1/2)
//'
//' }

using namespace Rcpp;

// [[Rcpp::export]]
double df_2(NumericVector df_grid, NumericVector cv, double df_0, double D, double tol = 0, bool silence = true){

  // Check if lengths of df_grid and cv are equal
  if(df_grid.size() != cv.size()){
    Rcpp::stop("'df_grid' and 'cv' must have the same length");
  }

  // Raw input df_grid and cv
  NumericVector raw_df_grid = clone(df_grid);
  NumericVector raw_cv = clone(cv);

  // Lengths of df_grid and cv before removing NA values
  int N = raw_df_grid.size();

  // After removing NA-values, store them in these two
  NumericVector clean_df_grid;
  NumericVector clean_cv;

  // Remove NA values
  for(int i = 0; i < raw_df_grid.size(); ++i){
    bool na_df_grid = ISNAN(raw_df_grid[i]);
    bool na_cv = ISNAN(raw_cv[i]);
    if(na_df_grid || na_cv){
      N--;
      continue;
    }
    clean_df_grid.push_back(raw_df_grid[i]);
    clean_cv.push_back(raw_cv[i]);
  }

  // Check validity of df_0
  if(df_0 < 2){
    Rcpp::stop("'df_0' [%d] must be larger than or equal to 2.",
               df_0);
  }
  else if(df_0 > max(clean_df_grid)){
    Rcpp::stop("'df_0' [%d] must be less than or equal to maximum of 'df_grid' [%d].",
               df_0,
               max(clean_df_grid));
  }

  // Check validity of D
  if(D <= 0){
    Rcpp::stop("'D' [%d] must be positive.",
               D);
  }
  else if(D > 30){
    Rcpp::stop("'D' [%d] must be less than or equal to 30 degrees.",
               df_0);
  }

  // Check the validity of tol
  if(tol < 0){
    Rcpp::stop("The tolerance, 'tol' [%d], must be non-negative.",
               tol);
  }

  // Create slopes D_j
  NumericVector D_j(N - 1);
  for(int j = 0; j < N - 1; ++j){
    double slope = (clean_cv[j + 1] - clean_cv[j]) / (clean_df_grid[j + 1] - clean_df_grid[j]);
    D_j[j] = std::atan(slope) * (180.0 / M_PI);
    if(D_j[j] < 0){
      D_j[j] = D_j[j] * (-1);
    }
  }

  // Iterative process
  for(int j = N - 2; j >= 0; --j){
    if(clean_df_grid[j + 1] < df_0 && D_j[j] <= D){
      double D_j_error = D_j[j] - D;
      bool last_satisfying = false;
      if(j >= 1 && D_j[j - 1] > D){
        if(!silence){
          Rcout << "D_j[" << j << " - 1] = " << D_j[j - 1] << " is larger than D "<< "\n";
        }
        last_satisfying = true;
      }
      if(D_j_error < 0){
        D_j_error = D_j_error * (-1);
      }
      if(!silence){
        Rcout << "Current abs(D_j - D) is " << D_j_error << "\n";
        Rcout << "Current df(c) is " << clean_df_grid[j + 1] << "\n";
      }
      if(last_satisfying){
        if(!silence){
          Rcout << "UTILIZING LAST SATISFYING " << j << "\n";
        }
        return clean_df_grid[j + 1];
      }
      else if(D_j_error <= tol){
        return clean_df_grid[j + 1];
      }
    }
    else{
      //Rcout << "df_grid[j] = " << df_grid[j + 1] << "is larger than df_0 = " << df_0 << "\n";
    }
  }
  return df_0;
}
