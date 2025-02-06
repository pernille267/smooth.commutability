#include <Rcpp.h>
using namespace Rcpp;

// Inline helper function for better performance
inline double b1_deming2(const double msxx, const double msyy, const double msxy, const double lambda) {
  const double sub_expression_1 = msyy - lambda * msxx;
  const double sub_expression_2 = std::sqrt(sub_expression_1 * sub_expression_1 + 4.0 * lambda * msxy * msxy);
  return (sub_expression_1 + sub_expression_2) / (2.0 * msxy);
}

// Optimized sum of squares calculation
inline double calculate_sxy2(const NumericVector& x, const NumericVector& y, const double mx, const double my, const int n) {
  double sxy = 0.0;
  for(int i = 0; i < n; ++i) {
    sxy += (x[i] - mx) * (y[i] - my);
  }
  return sxy;
}

//' Checks if \eqn{y_{n+j}} for \eqn{j = 1, \ldots, m} are inside Deming prediction intervals.
//'
//' @title Check Whether New Observations Are Inside Deming Prediction Intervals
//' @name inside_deming2
//'
//' @param data A \code{list} or \code{data.table}. Must contain \code{MP_A} (response variable values)
//'        and \code{MP_B} (predictor variable values).
//' @param new_data A \code{list} or \code{data.table}. Must contain \code{MP_A} (new response variable values)
//'        and \code{MP_B} (new predictor variable values).
//' @param imprecision_estimates A \code{list} or \code{data.table}. Must contain \code{lambda}
//'        and \code{Var_B}. See details for more information.
//' @param R A \code{integer} value. The number of replicated measurements. If the number of
//'        replicated measurements is sample-dependent, use some summary statistic.
//'        Defaults to \code{3}.
//' @param R_ratio A \code{double} value. If a different number of replicated measurements are used
//'        in data than in new_data, specify the ratio here. Defaults to \code{1}.
//' @param level A \code{double} value. The desired nominal confidence level for
//'        the estimated Deming prediction intervals. Must be between \code{0} and \code{1}.
//'        Sensible values include values larger than or equal to \code{0.50}. Defaults to \code{0.99}.
//'
//'
//' @details
//' The Deming regression \eqn{\mathrm{level} \cdot 100\%} prediction intervals are calculated by
//'
//' \eqn{PI_D(y_{n+j} | x_{n+j}) \approx b_0 + b_1 \cdot x_{n+j} \pm t_{n - 2, 1 - \alpha/2} \sqrt{\hat{\sigma}_{b_1}^2\Big[(\hat{x}_{n+j}^L - \hat{\mu})^2 + \tilde{\sigma}_h^2 \cdot R_{\mathrm{rat}}\Big] + R_{\mathrm{rat}} \cdot \Big(\frac{n+j}{n}\Big) \cdot (\tilde{\sigma}_v^2 + b_1^2 \cdot \tilde{\sigma}_h^2)}}
//'
//' This function calculates \eqn{z_j = \mathbb{I}\big[y_{n+j} \in PI_D(y_{n+j} | x_{n+j})\big]}
//' for each \eqn{j = 1, \ldots, m}. Note that
//' \itemize{
//'   \item \eqn{b_1} is the Deming regression slope estimate of \eqn{\beta_1}.
//'   \item \eqn{\hat{\sigma}_{b_1}^2} is the estimated variance of the Deming regression slope estimator \eqn{b_1}.
//'   \item \eqn{\hat{x}_{n+j}^L} is an estimate of the latent \eqn{x_{n+j}^L}.
//'   \item \eqn{\hat{\mu}} is an estimate of the mean of \eqn{x_{i}^L}.
//'   \item \eqn{\tilde{\sigma}_h^2} is the method of moments estimator of the repeatability variance \eqn{\sigma_h^2}.
//'   \item \eqn{\tilde{\sigma}_v^2} is the method of moments estimator of the repeatability variance \eqn{\sigma_v^2}.
//'   \item \eqn{R_{\mathrm{rat}}} is the \code{R_ratio} value.
//' }
//'
//' Missing values in \code{MP_A} or \code{MP_B} in \code{data} are allowed and will be silently removed
//' before \eqn{z_j} values are calculated. Missing values in \code{MP_A} or \code{MP_B} in \code{new_data}
//' are allowed, but corresponding \eqn{z_j} values will not be calculated in such cases.
//'
//' Note that \code{imprecision_estimates} should contain the following variables:
//' \itemize{
//'   \item \code{lambda}: This an estimator of the ratio \eqn{\sigma_v^2 / \sigma_h^2}.
//'         How you estimate this ratio is up to you. However, a natural estimator is
//'         \eqn{\hat{\sigma}_v^2 / \hat{\sigma}_h^2}, where \eqn{\hat{\sigma}_v^2} and
//'         \eqn{\hat{\sigma}_h^2} are pooled variances based on replicated measurements.
//'   \item \code{Var_B}: This is an estimator of \eqn{\sigma_h^2}. How you estimate this
//'         ratio is up to you. However a natural estimator is \eqn{\hat{\sigma}_h^2},
//'         which is a pooled variance based on replicated measurements.
//' }
//' N.B., if either \code{lambda} or \code{Var_B} is missing or takes negative values
//' Ordinary least squares regression is used instead to calculate \eqn{z_j} values.
//'
//'
//' @return An \code{integer} vector of length equal to the number of rows in \code{new_data}.
//'
//' @examples \dontrun{
//'   print(1)
//' }
// [[Rcpp::export]]
IntegerVector inside_deming2(const List& data, const List& new_data, const List& imprecision_estimates, const int R = 3, const double R_ratio = 1.0, const double level = 0.99) {

  // Input validation
  if (!data.containsElementNamed("MP_A") || !data.containsElementNamed("MP_B")) {
    stop("MP_A and MP_B must be present in data");
  }
  if (!new_data.containsElementNamed("MP_A") || !new_data.containsElementNamed("MP_B")) {
    stop("MP_A and MP_B must be present in new_data");
  }
  if (!imprecision_estimates.containsElementNamed("lambda") ||
      !imprecision_estimates.containsElementNamed("Var_B")) {
      stop("lambda and Var_B must be present in imprecision_estimates");
  }

  // Extract raw data
  NumericVector x_raw = data["MP_B"];
  NumericVector y_raw = data["MP_A"];
  const int n_raw = x_raw.length();

  // Handle NA values in training data
  LogicalVector valid_pairs(n_raw);
  int valid_count = 0;
  for(int i = 0; i < n_raw; ++i) {
    valid_pairs[i] = !ISNAN(x_raw[i]) && !ISNAN(y_raw[i]);
    if(valid_pairs[i]) valid_count++;
  }

  if(valid_count < 3) stop("Need at least 3 valid observations after removing NA values");

  // Create clean vectors without NA values
  NumericVector x(valid_count);
  NumericVector y(valid_count);
  int j = 0;
  for(int i = 0; i < n_raw; ++i) {
    if(valid_pairs[i]) {
      x[j] = x_raw[i];
      y[j] = y_raw[i];
      j++;
    }
  }
  const int n = valid_count;  // New sample size after removing NA values

  // Extract new data
  const NumericVector new_x = new_data["MP_B"];
  const NumericVector new_y = new_data["MP_A"];
  const int m = new_x.length();

  // Extract and validate parameters
  const double lambda = as<double>(imprecision_estimates["lambda"]);
  const double Var_B = as<double>(imprecision_estimates["Var_B"]);

  // Determine if OLS should be used
  const bool ols = ISNAN(lambda) || lambda <= 0.0 || ISNAN(Var_B) || Var_B <= 0.0;

  // Pre-allocate output vectors
  IntegerVector inside(m, NA_INTEGER);  // Initialize with NA
  NumericVector ny(m), var_pred_error(m), lwr(m), upr(m);

  // Calculate basic statistics using clean data
  const double mx = mean(x);
  const double my = mean(y);
  const double msxx = var(x);
  const double msyy = var(y);
  const double sxy = calculate_sxy2(x, y, mx, my, n);
  const double msxy = sxy / (n - 1.0);
  const double sxx = msxx * (n - 1.0);

  const double t_quantile = R::qt((1.0 - level) / 2.0, n - 2, 0, 0);

  if (ols) {
    const double b1 = sxy / sxx;
    const double b0 = my - b1 * mx;

    double mse = 0.0;
    for(int i = 0; i < n; ++i) {
      const double resid = y[i] - (b0 + b1 * x[i]);
      mse += resid * resid;
    }
    mse = (mse / (n - 2.0)) * R_ratio;

    for(int j = 0; j < m; ++j) {
      // Skip calculations if either new_x[j] or new_y[j] is NA
      if(ISNAN(new_x[j]) || ISNAN(new_y[j])) {
        continue;  // inside[j] remains NA
      }

      ny[j] = b0 + b1 * new_x[j];
      const double dev = new_x[j] - mx;
      var_pred_error[j] = mse * (1.0 + (1.0 / n) + (dev * dev) / sxx);
      const double margin = t_quantile * std::sqrt(var_pred_error[j]);
      lwr[j] = ny[j] - margin;
      upr[j] = ny[j] + margin;
      inside[j] = (new_y[j] >= lwr[j] && new_y[j] <= upr[j]) ? 1 : 0;
    }
  }
  else {
    const double b1 = b1_deming2(msxx, msyy, msxy, lambda);
    const double b0 = my - b1 * mx;

    const double var_b1 = (b1 * b1 / (n * msxy * msxy)) * ((msxx * msyy) - msxy * msxy);
    const double sigma_h_squared = (msyy + lambda * msxx -
                                    std::sqrt(std::pow(msyy - lambda * msxx, 2) + 4.0 * lambda * msxy * msxy)) / (2.0 * lambda);
    const double sigma_v_squared = lambda * sigma_h_squared;

    const double mu_ratio = msxy / (b1 * msxx);
    double mu_hat = 0.0;

    for(int i = 0; i < n; ++i) {
      mu_hat += ((lambda / (lambda + b1 * b1)) * x[i] +
        (b1 / (lambda + b1 * b1)) * (y[i] - b0)) / n;
    }

    for(int j = 0; j < m; ++j) {
      // Skip calculations if either new_x[j] or new_y[j] is NA
      if(ISNAN(new_x[j]) || ISNAN(new_y[j])) {
        continue;  // inside[j] remains NA
      }

      ny[j] = b0 + b1 * new_x[j];
      const double nx_latent = mu_hat * (1.0 - mu_ratio) + ny[j] * mu_ratio;
      const double dev = nx_latent - mu_hat;

      var_pred_error[j] = var_b1 * dev * dev +
        var_b1 * sigma_h_squared * R_ratio +
        (1.0 + 1.0/n) * (b1 * b1 * sigma_h_squared + sigma_v_squared) * R_ratio;

      const double margin = t_quantile * std::sqrt(var_pred_error[j]);
      lwr[j] = ny[j] - margin;
      upr[j] = ny[j] + margin;
      inside[j] = (new_y[j] >= lwr[j] && new_y[j] <= upr[j]) ? 1 : 0;
    }
  }

  return inside;
}

//' Calculates the Deming regression estimate of \eqn{\sigma_h^2}.
//'
//' @title Calculates the Deming Regression Estimate of \eqn{\sigma_h^2}.
//' @name sigma_h_squared_deming2
//'
//' @param data A \code{list} or \code{data.table}. Must contain \code{MP_A} (response variable values)
//'        and \code{MP_B} (predictor variable values).
//' @param lambda A \code{double} value.
//'
//' @details
//' This function is typically not relevant for end-users.
//'
//'
//' @return An \code{double} that is the Deming regression estimate of \eqn{\sigma_h^2}.
//'
//' @examples \dontrun{
//'   print(1)
//' }
// [[Rcpp::export]]
double sigma_h_squared_deming2(const List& data, const double lambda = 1) {

  // Input validation
  if (!data.containsElementNamed("MP_A") || !data.containsElementNamed("MP_B")) {
    stop("MP_A and MP_B must be present in data");
  }

  // Extract raw data
  NumericVector x_raw = data["MP_B"];
  NumericVector y_raw = data["MP_A"];
  const int n_raw = x_raw.length();

  // Handle NA values in training data
  LogicalVector valid_pairs(n_raw);
  int valid_count = 0;
  for(int i = 0; i < n_raw; ++i) {
    valid_pairs[i] = !ISNAN(x_raw[i]) && !ISNAN(y_raw[i]);
    if(valid_pairs[i]) valid_count++;
  }

  if(valid_count < 3) stop("Need at least 3 valid observations after removing NA values");

  // Create clean vectors without NA values
  NumericVector x(valid_count);
  NumericVector y(valid_count);
  int j = 0;
  for(int i = 0; i < n_raw; ++i) {
    if(valid_pairs[i]) {
      x[j] = x_raw[i];
      y[j] = y_raw[i];
      j++;
    }
  }
  const int n = valid_count;  // New sample size after removing NA values

  // Calculate basic statistics using clean data
  const double mx = mean(x);
  const double my = mean(y);
  const double msxx = var(x);
  const double msyy = var(y);
  const double sxy = calculate_sxy2(x, y, mx, my, n);
  const double msxy = sxy / (n - 1.0);

  //const double b1 = b1_deming2(msxx, msyy, msxy, lambda);
  //const double b0 = my - b1 * mx;

  //const double var_b1 = (b1 * b1 / (n * msxy * msxy)) * ((msxx * msyy) - msxy * msxy);
  const double sigma_h_squared = (msyy + lambda * msxx -
                                  std::sqrt(std::pow(msyy - lambda * msxx, 2) + 4.0 * lambda * msxy * msxy)) / (2.0 * lambda);
  //const double sigma_v_squared = lambda * sigma_h_squared;

  return sigma_h_squared;
}


