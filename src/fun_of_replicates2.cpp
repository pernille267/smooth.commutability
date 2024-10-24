#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// (*) Helper function to calculate the specified statistic
double calculate_stat(const std::vector<double>& values, const std::string& fun) {
  if (values.empty()) {
    return NA_REAL;
  }

  double result = 0.0;
  size_t n = values.size();

  if (fun == "mean") {
    result = std::accumulate(values.begin(), values.end(), 0.0) / n;
  } else if (fun == "var") {
    double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
    double sqSum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
    result = (sqSum - n * mean * mean) / (n - 1);
  } else if (fun == "median") {
    std::vector<double> sortedValues = values;
    std::sort(sortedValues.begin(), sortedValues.end());
    result = (n % 2 == 0) ? (sortedValues[n/2 - 1] + sortedValues[n/2]) / 2 : sortedValues[n/2];
  } else if (fun == "sd") {
    double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
    double sqSum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
    result = std::sqrt((sqSum - n * mean * mean) / (n - 1));
  } else if (fun == "cv") {
      double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
      double sqSum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
      double stDev = std::sqrt((sqSum - n * mean * mean) / (n - 1));
      result = stDev / mean;
  } else if (fun == "min") {
    result = *std::min_element(values.begin(), values.end());
  } else if (fun == "max") {
    result = *std::max_element(values.begin(), values.end());
  } else {
    Rcpp::stop("Invalid statistic type specified");
  }

  return result;
}

//' Calculate a Statistic for Each Sample
//'
//' @title Calculate a Statistic for Each Sample
//' @name fun_of_replicates2
//'
//' @param data A \code{data.table} or \code{list} object. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//' @param fun A \code{character} string. Which statistic is to be calculated for
//'        each SampleID. Possible choices include \code{mean}, \code{var} (variance),
//'        \code{sd} (standard deviaton), \code{cv} (coefficient of variation),
//'        \code{median}, \code{min} (minimum) and \code{max} (maximum).
//'
//' @description Calculates a chosen statistic over the replicates.
//'
//' @details This function handles NA-values automatically.
//'
//' @return A \code{list} with elements \code{SampleID}, \code{MP_A} and \code{MP_B}.
//'         \code{MP_A} and \code{MP_B} contains the sample-wise statistics.
//'
//' @examples \dontrun{
//'   print(1)
//' }
// [[Rcpp::export]]
List fun_of_replicates2(List data, std::string fun = "mean") {
  CharacterVector sampleID = data["SampleID"];
  NumericVector mpA = data["MP_A"];
  NumericVector mpB = data["MP_B"];

  std::unordered_map<std::string, std::vector<double>> groupedDataA;
  std::unordered_map<std::string, std::vector<double>> groupedDataB;
  std::vector<std::string> uniqueIDsOrdered;

  int n = sampleID.size();

  for (int i = 0; i < n; ++i) {
    std::string id = as<std::string>(sampleID[i]);
    if (groupedDataA.find(id) == groupedDataA.end()) {
      uniqueIDsOrdered.push_back(id);
    }
    if (!NumericVector::is_na(mpA[i])) {
      groupedDataA[id].push_back(mpA[i]);
    }
    if (!NumericVector::is_na(mpB[i])) {
      groupedDataB[id].push_back(mpB[i]);
    }
  }

  CharacterVector uniqueIDs(uniqueIDsOrdered.size());
  NumericVector statA(uniqueIDsOrdered.size());
  NumericVector statB(uniqueIDsOrdered.size());

  for (size_t i = 0; i < uniqueIDsOrdered.size(); ++i) {
    const std::string& id = uniqueIDsOrdered[i];
    uniqueIDs[i] = id;

    statA[i] = calculate_stat(groupedDataA[id], fun);
    statB[i] = calculate_stat(groupedDataB[id], fun);
  }

  return List::create(
    Named("SampleID") = uniqueIDs,
    Named("MP_A") = statA,
    Named("MP_B") = statB
  );
}
