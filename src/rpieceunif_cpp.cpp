// File: src/rpieceunif_cpp.cpp
#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

//' Generate random numbers from piecewise uniform distribution (C++ version)
 //'
 //' This function generates random numbers from a piecewise uniform distribution
 //' using C++ for maximum performance.
 //'
 //' @param n A positive integer specifying the number of observations to generate.
 //' @param time A numeric vector of time points defining the intervals.
 //' @param intensity A numeric vector of intensities (or NULL).
 //' @param proportion A numeric vector of proportions (or NULL).
 //' @return A numeric vector of random numbers from the piecewise uniform distribution.
 //' @export
 // [[Rcpp::export]]
 NumericVector rpieceunif_cpp(int n,
                              NumericVector time,
                              Nullable<NumericVector> intensity = R_NilValue,
                              Nullable<NumericVector> proportion = R_NilValue) {

   // Input validation
   bool has_intensity = !intensity.isNull();
   bool has_proportion = !proportion.isNull();

   if (!has_intensity && !has_proportion) {
     stop("Either intensity or proportion must be specified");
   }

   if (has_intensity && has_proportion) {
     stop("Cannot specify both intensity and proportion arguments");
   }

   // Convert to std::vector for faster access
   std::vector<double> time_vec = as<std::vector<double>>(time);
   int n_intervals = time_vec.size() - 1;

   // Calculate cumulative proportions
   std::vector<double> cumulative_props(n_intervals);

   if (has_proportion) {
     // Using proportions
     NumericVector prop = as<NumericVector>(proportion);

     if (prop.size() != n_intervals) {
       stop("Length of time must be one more than length of proportion");
     }

     // Check for non-negative proportions
     for (int i = 0; i < prop.size(); ++i) {
       if (prop[i] < 0) {
         stop("All proportions must be non-negative");
       }
     }

     // Check if proportions sum to 1
     double prop_sum = 0.0;
     for (int i = 0; i < prop.size(); ++i) {
       prop_sum += prop[i];
     }

     if (std::abs(prop_sum - 1.0) > 1e-10) {
       stop("Proportions must sum to 1");
     }

     // Calculate cumulative proportions
     double running_sum = 0.0;
     for (int i = 0; i < n_intervals; ++i) {
       running_sum += prop[i];
       cumulative_props[i] = running_sum;
     }

   } else {
     // Using intensities
     NumericVector intens = as<NumericVector>(intensity);

     if (intens.size() != n_intervals) {
       stop("Length of time must be one more than length of intensity");
     }

     // Check for positive intensities
     for (int i = 0; i < intens.size(); ++i) {
       if (intens[i] <= 0) {
         stop("All intensities must be positive");
       }
     }

     // Check for strictly increasing time vector
     for (int i = 1; i < time_vec.size(); ++i) {
       if (time_vec[i] <= time_vec[i-1]) {
         stop("Time vector must be strictly increasing");
       }
     }

     // Calculate expected counts and cumulative proportions
     std::vector<double> expected_counts(n_intervals);
     double total_expected = 0.0;

     for (int i = 0; i < n_intervals; ++i) {
       double interval_length = time_vec[i+1] - time_vec[i];
       expected_counts[i] = intens[i] * interval_length;
       total_expected += expected_counts[i];
     }

     // Calculate cumulative proportions
     double running_sum = 0.0;
     for (int i = 0; i < n_intervals; ++i) {
       running_sum += expected_counts[i] / total_expected;
       cumulative_props[i] = running_sum;
     }
   }

   // Generate random numbers
   NumericVector result(n);
   NumericVector u_interval = runif(n);  // For interval selection
   NumericVector u_position = runif(n);  // For position within interval

   // Process each random number
   for (int i = 0; i < n; ++i) {
     // Find interval using binary search-like approach
     int interval_idx = 0;
     for (int j = 0; j < n_intervals; ++j) {
       if (u_interval[i] <= cumulative_props[j]) {
         interval_idx = j;
         break;
       }
     }

     // Generate position within interval
     double start_time = time_vec[interval_idx];
     double end_time = time_vec[interval_idx + 1];
     result[i] = start_time + u_position[i] * (end_time - start_time);
   }

   return result;
 }
