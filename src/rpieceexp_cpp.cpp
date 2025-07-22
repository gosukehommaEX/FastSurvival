// File: src/rpieceexp_cpp.cpp
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

//' Generate random numbers from piecewise exponential distribution (C++ version)
 //'
 //' This function generates random numbers from a piecewise exponential distribution
 //' using C++ for maximum performance.
 //'
 //' @param n A positive integer specifying the number of observations to generate.
 //' @param time A numeric vector of time points defining the intervals.
 //'   The last element must be Inf.
 //' @param hazard A numeric vector of hazard rates for each interval.
 //' @return A numeric vector of random numbers from the piecewise exponential distribution.
 //' @export
 // [[Rcpp::export]]
 NumericVector rpieceexp_cpp(int n, NumericVector time, NumericVector hazard) {

   // Input validation
   if (time.size() != hazard.size() + 1) {
     stop("Length of time must be one more than length of hazard");
   }

   if (is_false(any(hazard > 0))) {
     stop("All hazard rates must be positive");
   }

   // Check if last element is Inf
   if (!R_IsNA(time[time.size() - 1]) && !R_PosInf == time[time.size() - 1]) {
     if (std::isfinite(time[time.size() - 1])) {
       stop("Last element of time must be Inf");
     }
   }

   // Convert to std::vector for faster access
   std::vector<double> time_vec = as<std::vector<double>>(time);
   std::vector<double> hazard_vec = as<std::vector<double>>(hazard);

   int n_intervals = hazard_vec.size();

   // Check for strictly increasing finite time points
   for (int i = 1; i < n_intervals; ++i) {
     if (time_vec[i] <= time_vec[i-1]) {
       stop("Finite time points must be strictly increasing");
     }
   }

   // Calculate cumulative hazard at time change points
   std::vector<double> cumulative_hazard(n_intervals + 1);
   cumulative_hazard[0] = 0.0;

   for (int i = 1; i <= n_intervals; ++i) {
     if (i < n_intervals) {
       // Finite interval
       double interval_length = time_vec[i] - time_vec[i-1];
       cumulative_hazard[i] = cumulative_hazard[i-1] + hazard_vec[i-1] * interval_length;
     } else {
       // Infinite interval - set to a large value for comparison
       cumulative_hazard[i] = R_PosInf;
     }
   }

   // Generate uniform random numbers
   NumericVector u = runif(n);

   // Calculate target hazards
   NumericVector result(n);

   for (int i = 0; i < n; ++i) {
     double target_hazard = -log(u[i]);

     // Find which interval the target falls into
     int interval_idx = 0;
     for (int j = 1; j <= n_intervals; ++j) {
       if (target_hazard <= cumulative_hazard[j]) {
         interval_idx = j - 1;
         break;
       }
     }

     // Ensure interval_idx is within bounds
     if (interval_idx >= n_intervals) {
       interval_idx = n_intervals - 1;
     }

     // Calculate result using inverse CDF method
     double start_time = time_vec[interval_idx];
     double current_hazard = hazard_vec[interval_idx];
     double prev_cumulative_hazard = cumulative_hazard[interval_idx];

     double remaining_hazard = target_hazard - prev_cumulative_hazard;
     result[i] = start_time + remaining_hazard / current_hazard;
   }

   return result;
 }
