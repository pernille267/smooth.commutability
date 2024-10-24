#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// (*) Helper function for estimating variance
double estimate_variance(const std::vector<double>& values) {
  double sum = 0.0, sq_sum = 0.0;
  const int n = values.size();

  for (double value : values) {
    sum += value;
    sq_sum += value * value;
  }

  double mean = sum / n;
  return (sq_sum - sum * mean) / (n - 1);
}

// (**) Helper function for obtaining unique elements maintaining their original order
CharacterVector unique_preserve_order(CharacterVector x) {
  std::unordered_set<std::string> seen;
  CharacterVector result;

  for(int i = 0; i < x.length(); i++) {
    std::string current = as<std::string>(x[i]);
    if(seen.find(current) == seen.end()) {
      seen.insert(current);
      result.push_back(x[i]);
    }
  }

  return result;
}

//' Calculate imprecision point estimates of measurements in a given IVD-MD comparison
//'
//' @title Calculate imprecision point estimates of measurements in a given IVD-MD comparison
//' @name global_precision_estimates2
//'
//' @param data \code{list} or \code{data.table}. Must contain the following variables:
//' \itemize{
//'   \item \code{SampleID}: Sample identifiers. Must be \code{character}.
//'   \item \code{ReplicateID}: Replicate measurement identifiers within samples. Must be \code{character}.
//'   \item \code{MP_A}: Measurement results for the IVD-MD used as response variable.
//'   \item \code{MP_B}: Measurement results for the IVD-MD used as predictor variable.
//' }
//'
//' @details
//' Calculates five relevant global imprecision estimates. The term global is
//' used because these imprecision estimates are based on the whole dataset.
//' These five imprecision estimates are calculated:
//' \itemize{
//'   \item{\code{Var_A: }}{Pooled variance of all sample-variances based on MP_A}
//'   \item{\code{Var_B: }}{Pooled variance of all sample-variances based on MP_B}
//'   \item{\code{CV_A: }}{CV estimate based on Var_A and the grand mean of all measurements from MP_A}
//'   \item{\code{CV_B: }}{CV estimate based on Var_B and the grand mean of all measurements from MP_B}
//'   \item{\code{lambda: }}{Ratio of pooled variances Var_A and Var_B}
//' }
//' CV values \code{CV_A} and \code{CV_B} can also be represented in percent.
//' To convert these to percentage values, one just multiply their raw results by 100%.
//'
//'
//' @return A \code{list} of length \code{5} with the following point imprecision estimates:
//'         \code{Var_A}, \code{Var_B}, \code{CV_A}, \code{CV_B} and \code{lambda}.
//'         See details for more information on these statistics.
//'
//' @examples \dontrun{
//'   library(data.table)
//'   data <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.02, cvy = 0.3), AR = TRUE)
//'   data$SampleID <- as.character(data$SampleID)
//'   data$ReplicateID <- as.character(data$ReplicateID)
//'   print(global_prcision_estimates2(data = data) |> as.data.table())
//' }

 // [[Rcpp::export]]
List global_precision_estimates2(List data) {
  // Extract data columns
  CharacterVector SampleID = data["SampleID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];

  // Get unique samples and initialize result vectors
  CharacterVector summary_SampleID = unique_preserve_order(SampleID);
  const int n = summary_SampleID.size();
  const int N = SampleID.size();
  NumericVector ith_var_MS_A(n, NA_REAL);
  NumericVector ith_var_MS_B(n, NA_REAL);
  const int replicate_number_requirement = 2;

  // Create hash map for faster sample lookup
  std::unordered_map<String, std::vector<int>> sample_indices;
  for (int i = 0; i < N; ++i) {
    sample_indices[String(SampleID[i])].push_back(i);
  }

  // Calculate variances for each sample
  for (int i = 0; i < n; ++i) {
    String current_sample = String(summary_SampleID[i]);
    const std::vector<int>& indices = sample_indices[current_sample];

    std::vector<double> valid_measurements_A;
    std::vector<double> valid_measurements_B;
    valid_measurements_A.reserve(indices.size());
    valid_measurements_B.reserve(indices.size());

    // Collect valid measurements
    for (int idx : indices) {
      if (!ISNAN(MS_A[idx])) {
        valid_measurements_A.push_back(MS_A[idx]);
      }
      if (!ISNAN(MS_B[idx])) {
        valid_measurements_B.push_back(MS_B[idx]);
      }
    }

    // Calculate variances if enough replicates
    if (valid_measurements_A.size() >= replicate_number_requirement) {
      ith_var_MS_A[i] = estimate_variance(valid_measurements_A);
    }
    if (valid_measurements_B.size() >= replicate_number_requirement) {
      ith_var_MS_B[i] = estimate_variance(valid_measurements_B);
    }
  }

  // Calculate global statistics
  double var_MS_A = 0, var_MS_B = 0;
  int effective_n_A = 0, effective_n_B = 0;

  for (int i = 0; i < n; ++i) {
    if (!ISNAN(ith_var_MS_A[i])) {
      var_MS_A += ith_var_MS_A[i];
      effective_n_A++;
    }
    if (!ISNAN(ith_var_MS_B[i])) {
      var_MS_B += ith_var_MS_B[i];
      effective_n_B++;
    }
  }

  var_MS_A = effective_n_A > 0 ? var_MS_A / effective_n_A : NA_REAL;
  var_MS_B = effective_n_B > 0 ? var_MS_B / effective_n_B : NA_REAL;

  // Calculate means
  double mean_MS_A = 0, mean_MS_B = 0;
  int valid_count_A = 0, valid_count_B = 0;

  for (int i = 0; i < N; ++i) {
    if (!ISNAN(MS_A[i])) {
      mean_MS_A += MS_A[i];
      valid_count_A++;
    }
    if (!ISNAN(MS_B[i])) {
      mean_MS_B += MS_B[i];
      valid_count_B++;
    }
  }

  mean_MS_A = valid_count_A > 0 ? mean_MS_A / valid_count_A : NA_REAL;
  mean_MS_B = valid_count_B > 0 ? mean_MS_B / valid_count_B : NA_REAL;

  // Calculate final statistics
  double cv_MS_A = (!ISNAN(var_MS_A) && !ISNAN(mean_MS_A) && mean_MS_A != 0) ?
  sqrt(var_MS_A) / mean_MS_A : NA_REAL;
  double cv_MS_B = (!ISNAN(var_MS_B) && !ISNAN(mean_MS_B) && mean_MS_B != 0) ?
  sqrt(var_MS_B) / mean_MS_B : NA_REAL;
  double lambda = (!ISNAN(var_MS_A) && !ISNAN(var_MS_B) && var_MS_B != 0) ?
  var_MS_A / var_MS_B : NA_REAL;

  return List::create(
    Named("Var_A") = var_MS_A,
    Named("Var_B") = var_MS_B,
    Named("CV_A") = cv_MS_A,
    Named("CV_B") = cv_MS_B,
    Named("lambda") = lambda
  );
}

//' Resample clustered EQA clinical sample data
//'
//' @title Resample clustered EQA clinical sample data
//' @name resample_samples2
//'
//' @param data A \code{list} or a \code{data.table}. Must contain \code{SampleID},
//'        \code{ReplicateID}, \code{MP_A} and \code{MP_B}. The ID variables \code{SampleID}
//'        and \code{ReplicateID} must be of character type for the function to operate correctly.
//'
//' @details
//' This function is a very efficient method to resample clinical sample data on sample-level.
//' It is convenient to combine this function with \code{fasteqa} functions such as
//' \itemize{
//'   \item \code{global_precision_estimates()} to estimate bootstrap imprecision confidence intervals
//'   \item \code{estimate_zeta()} to estimate bootstrap zeta confidence intervals
//' }
//' Alternatively it could be used to estimate \code{smooth.commutability} functions such as
//' \itemize{
//'   \item \code{estimate_zeta_ss()} to estimate bootstrap smoothing spline zeta confidence intervals
//'   \item \code{estimate_df()} to estimate bootstrap smoothing spline df confidence intervals
//' }
//'
//' @return A \code{list} containing the resampled EQA clinical sample data.
//'
//' @examples
//' library(data.table)
//' fictive_data <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.01, cvy = 0.01), AR = TRUE)
//' resampled_data <- resample_samples2(fictive_data)
//' setDT(resampled_data)
//' print(resampled_data)

// [[Rcpp::export]]
List resample_samples2(List data) {
  CharacterVector samples = data["SampleID"];
  CharacterVector replicates = data["ReplicateID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  CharacterVector unique_samples = unique_preserve_order(samples);
  int n = unique_samples.size();
  int N = samples.size();

  CharacterVector resampled_samples = sample(unique_samples, n, true);

  // Create a map for faster lookup
  std::unordered_map<String, std::vector<int>> sample_indices;
  for (int i = 0; i < N; ++i) {
    sample_indices[samples[i]].push_back(i);
  }

  // Calculate output size and prepare output vectors
  int output_size = 0;
  for (int i = 0; i < n; ++i) {
    output_size += sample_indices[resampled_samples[i]].size();
  }

  CharacterVector new_samples(output_size);
  CharacterVector new_replicates(output_size);
  NumericVector new_MS_A(output_size);
  NumericVector new_MS_B(output_size);

  // Fill output vectors
  int counter = 0;
  for (int i = 0; i < n; ++i) {
    String current_sample = resampled_samples[i];
    const std::vector<int>& indices = sample_indices[current_sample];
    for (int idx : indices) {
      new_samples[counter] = unique_samples[i];
      new_replicates[counter] = replicates[idx];
      new_MS_A[counter] = MS_A[idx];
      new_MS_B[counter] = MS_B[idx];
      ++counter;
    }
  }

  return List::create(
    Named("SampleID") = new_samples,
    Named("ReplicateID") = new_replicates,
    Named("MP_A") = new_MS_A,
    Named("MP_B") = new_MS_B
  );
}

//' Resample cluster statistics based on EQA clinical sample data
//'
//' @title Resample cluster statistics based on EQA clinical sample data
//' @name resample_fun_samples2
//'
//' @param data A \code{list} or a \code{data.table}. The mean-of-replicates
//'        clinical sample data. Must contain \code{SampleID}, \code{MP_A} and \code{MP_B}.
//'        The ID variable \code{SampleID} must be \code{character}.
//' @param weight_data A \code{list} or a \code{data.table}. The weight data based
//'        on the clinical sample data. Must contain \code{SampleID}, \code{MP_A} and \code{MP_B}.
//'        The ID variable \code{SampleID} must be \code{character}.
//'
//' @details
//' This function is a very efficient method to resample aggregated sample statistics
//' based on clinical sample data.
//'
//' It is convenient to combine this function with \code{smooth.commutability} functions such as
//' \itemize{
//'   \item \code{predict_smoothing_splines()} to estimate inside rates for a IVD-MD comparison
//'   \item \code{smoothing_spline()} to estimate bootstrap distribution of LOO-CV chosen effective degrees of freedom.
//' }
//' If you do not have weight data available that you seek to resample jointly,
//' just pass \code{data} as both first and second argument.
//' See example under documentation of \code{resample_fun_samples2_all}.
//'
//' @return
//' A \code{list} of length two containing the resampled cluster statistics based on
//' \code{data} and \code{weight_data}. The output \code{list} have the following structure:
//' \itemize{
//'   \item \code{resampled_cs_data}: Contains the resampled \code{data}
//'   \item \code{resampled_weight_data}: Contains the resampled \code{weight_data}
//' }
//'
//'
//' @examples
//' library(data.table)
//' fictive_data <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.01, cvy = 0.02))
//' fictive_weight_data <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.03, cvy = 0.02))
//' resampled_data <- resample_fun_samples2(fictive_data, fictive_weight_data)
//' setDT(resampled_data)
//' print(resampled_data)

// [[Rcpp::export]]
List resample_fun_samples2(List data, List weight_data) {
  CharacterVector samples = data["SampleID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  NumericVector W_A = weight_data["MP_A"];
  NumericVector W_B = weight_data["MP_B"];
  int n = samples.size();

  CharacterVector resampled_samples = sample(samples, n, true);

  // Create a map for faster lookup
  std::unordered_map<String, std::vector<int>> sample_indices;
  for (int i = 0; i < n; ++i) {
    sample_indices[samples[i]].push_back(i);
  }

  // Calculate output size and prepare output vectors
  int output_size = 0;
  for (int i = 0; i < n; ++i) {
    output_size += sample_indices[resampled_samples[i]].size();
  }

  // Set up ...
  CharacterVector new_samples(output_size);
  NumericVector new_MS_A(output_size);
  NumericVector new_MS_B(output_size);
  NumericVector new_W_A(output_size);
  NumericVector new_W_B(output_size);

  // Fill output vectors
  int counter = 0;
  for (int i = 0; i < n; ++i) {
    String current_sample = resampled_samples[i];
    const std::vector<int>& indices = sample_indices[current_sample];
    for (int idx : indices) {
      new_samples[counter] = samples[i];
      new_MS_A[counter] = MS_A[idx];
      new_MS_B[counter] = MS_B[idx];
      new_W_A[counter] = W_A[idx];
      new_W_B[counter] = W_B[idx];
      ++counter;
    }
  }

  // Output 1
  List resampled_cs_data = List::create(Named("SampleID") = new_samples,
                                        Named("MP_A") = new_MS_A,
                                        Named("MP_B") = new_MS_B);

  // Output 2
  List resampled_weight_data = List::create(Named("SampleID") = new_samples,
                                            Named("MP_A") = new_W_A,
                                            Named("MP_B") = new_W_B);

  // Merged outputs
  return List::create(
    Named("resampled_cs_data") = resampled_cs_data,
    Named("resampled_weight_data") = resampled_weight_data
  );
}

//' Resample cluster statistics based on EQA clinical sample data for each IVD-MD comparison
//'
//' @title Resample cluster statistics based on EQA clinical sample data for each IVD-MD comparison
//' @name resample_fun_samples2_all
//'
//' @param data A \code{list} or a \code{data.table}. The mean-of-replicates
//'        clinical sample data. Must contain \code{comparison} \code{SampleID},
//'        \code{MP_A} and \code{MP_B}. The ID variables \code{comparison}
//'        and \code{SampleID} must be \code{character}.
//' @param weight_data A \code{list} or a \code{data.table}. The weight data based
//'        on the clinical sample data. Must contain \code{comparison} \code{SampleID},
//'        \code{MP_A} and \code{MP_B}.The ID variables \code{comparison}
//'        and \code{SampleID} must be \code{character}.
//'
//' @details
//' This function is a very efficient method to resample aggregated sample statistics
//' based on clinical sample data grouped by IVD-MD \code{comparison}.
//'
//' It is convenient to combine this function with \code{smooth.commutability} functions such as
//' \itemize{
//'   \item \code{predict_smoothing_splines()} to estimate inside rates for each IVD-MD comparison
//'   \item \code{smoothing_spline()} to estimate bootstrap distribution of LOO-CV chosen effective degrees of freedom.
//' }
//' If you do not have weight data available that you seek to resample jointly,
//' just pass \code{data} as both first and second argument. See example.
//'
//' @return
//' A \code{list} of length two containing the resampled cluster statistics based on
//' \code{data} and \code{weight_data}. The output \code{list} have the following structure:
//' \itemize{
//'   \item \code{resampled_cs_data}: Contains the resampled \code{data}
//'   \item \code{resampled_weight_data}: Contains the resampled \code{weight_data}
//' }
//'
//'
//' @examples
//' library(data.table)
//' fictive_data1 <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.01, cvy = 0.02))
//' fictive_data2 <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.03, cvy = 0.02))
//' fictive_data1$comparison <- rep("A - B", length(fictive_data1$MP_B))
//' fictive_data2$comparison <- rep("A - C", length(fictive_data2$MP_B))
//' fictive_data <- rbindlist(list(fictive_data1, fictive_data2))
//' resampled_data <- resample_fun_samples2_all(fictive_data, fictive_data)
//' setDT(resampled_data)
//' print(resampled_data)

// [[Rcpp::export]]
List resample_fun_samples2_all(List data, List weight_data) {
 CharacterVector comparison = data["comparison"];
 CharacterVector samples = data["SampleID"];
 NumericVector MS_A = data["MP_A"];
 NumericVector MS_B = data["MP_B"];
 NumericVector W_A = weight_data["MP_A"];
 NumericVector W_B = weight_data["MP_B"];
 CharacterVector unique_comparison = unique_preserve_order(comparison);
 int n = unique_comparison.size();
 int N = comparison.size();

 // Create a map for faster lookup
 std::unordered_map<String, std::vector<int>> comparison_indices;
 for (int i = 0; i < N; ++i) {
   comparison_indices[comparison[i]].push_back(i);
 }

 CharacterVector new_comparison(N);
 CharacterVector new_samples(N);
 NumericVector new_MS_A(N);
 NumericVector new_MS_B(N);
 NumericVector new_W_A(N);
 NumericVector new_W_B(N);

 // Fill output vectors
 int counter = 0;
 for (int i = 0; i < n; ++i) {
   int sub_counter = 0;
   String current_comparison = comparison[i];
   const std::vector<int>& indices = comparison_indices[current_comparison];
   int output_size = indices.size();
   CharacterVector filter_samples(output_size);
   NumericVector filter_MS_A(output_size);
   NumericVector filter_MS_B(output_size);
   NumericVector filter_W_A(output_size);
   NumericVector filter_W_B(output_size);
   for (int idx : indices) {
     filter_samples[sub_counter] = samples[idx];
     filter_MS_A[sub_counter] = MS_A[idx];
     filter_MS_B[sub_counter] = MS_B[idx];
     filter_W_A[sub_counter] = W_A[idx];
     filter_W_B[sub_counter] = W_B[idx];
     ++sub_counter;
   }
   List data_in = List::create(Named("SampleID") = filter_samples,
                               Named("MP_A") = filter_MS_A,
                               Named("MP_B") = filter_MS_B);

   List weight_data_in = List::create(Named("SampleID") = filter_samples,
                                      Named("MP_A") = filter_W_A,
                                      Named("MP_B") = filter_W_B);

   List out_data_both = resample_fun_samples2(data_in, weight_data_in);
   List out_data = out_data_both["resampled_cs_data"];
   List out_weight_data = out_data_both["resampled_weight_data"];

   CharacterVector new_filter_samples = out_data["SampleID"];
   NumericVector new_filter_MS_A = out_data["MP_A"];
   NumericVector new_filter_MS_B = out_data["MP_B"];
   NumericVector new_filter_W_A = out_weight_data["MP_A"];
   NumericVector new_filter_W_B = out_weight_data["MP_B"];

   for (int j = 0; j < output_size; ++j) {
     new_comparison[counter] = unique_comparison[i];
     new_samples[counter] = new_filter_samples[j];
     new_MS_A[counter] = new_filter_MS_A[j];
     new_MS_B[counter] = new_filter_MS_B[j];
     new_W_A[counter] = new_filter_W_A[j];
     new_W_B[counter] = new_filter_W_B[j];
     ++counter;
   }


 }

 List resampled_cs_data = List::create(Named("comparison") = new_comparison,
                                       Named("SampleID") = new_samples,
                                       Named("MP_A") = new_MS_A,
                                       Named("MP_B") = new_MS_B);

 List resampled_weight_data = List::create(Named("comparison") = new_comparison,
                                           Named("SampleID") = new_samples,
                                           Named("MP_A") = new_W_A,
                                           Named("MP_B") = new_W_B);


 return List::create(
   Named("resampled_cs_data") = resampled_cs_data,
   Named("resampled_weight_data") = resampled_weight_data
 );
}

//' Resample imprecision estimates based on clustered EQA clinical sample data
//'
//' @title Resample imprecision estimates based on clustered EQA clinical sample data
//' @name resample_imprecision2
//'
//' @param data A \code{list} or a \code{data.table}. Must contain \code{SampleID},
//'        \code{ReplicateID}, \code{MP_A} and \code{MP_B}. The ID variables \code{SampleID}
//'        and \code{ReplicateID} must be of character type for the function to operate correctly.
//'
//' @details
//' This function is a very efficient method to resample repeatability in clinical sample data.
//'
//' @return A \code{list} containing the resampled imprecision.
//'
//' @examples
//' library(data.table)
//' fictive_data <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.01, cvy = 0.01), AR = TRUE)
//' impr <- replicate(n = 5, expr = resample_imprecision2(fictive_data), simplify = FALSE)
//' print(impr)


// [[Rcpp::export]]
List resample_imprecision2(List data){
  List resampled_data = resample_samples2(data);
  List output = global_precision_estimates2(resampled_data);
  return output;
}

