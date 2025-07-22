// File: src/simData_cpp.cpp (修正版)
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Internal function: Fast piecewise uniform generation
std::vector<double> rpieceunif_internal(int n,
                                        const std::vector<double>& time,
                                        const std::vector<double>& rates,
                                        bool use_intensity) {

  int n_intervals = time.size() - 1;
  std::vector<double> cumulative_props(n_intervals);

  if (use_intensity) {
    // Calculate expected counts and cumulative proportions
    std::vector<double> expected_counts(n_intervals);
    double total_expected = 0.0;

    for (int i = 0; i < n_intervals; ++i) {
      double interval_length = time[i+1] - time[i];
      expected_counts[i] = rates[i] * interval_length;
      total_expected += expected_counts[i];
    }

    // Calculate cumulative proportions
    double running_sum = 0.0;
    for (int i = 0; i < n_intervals; ++i) {
      running_sum += expected_counts[i] / total_expected;
      cumulative_props[i] = running_sum;
    }
  } else {
    // Using proportions directly
    double running_sum = 0.0;
    for (int i = 0; i < n_intervals; ++i) {
      running_sum += rates[i];
      cumulative_props[i] = running_sum;
    }
  }

  // Generate random numbers
  std::vector<double> result(n);

  // Generate all random numbers at once
  NumericVector u_interval = runif(n);
  NumericVector u_position = runif(n);

  for (int i = 0; i < n; ++i) {
    // Find interval
    int interval_idx = 0;
    for (int j = 0; j < n_intervals; ++j) {
      if (u_interval[i] <= cumulative_props[j]) {
        interval_idx = j;
        break;
      }
    }

    // Generate position within interval
    double start_time = time[interval_idx];
    double end_time = time[interval_idx + 1];
    result[i] = start_time + u_position[i] * (end_time - start_time);
  }

  return result;
}

// Internal function: Fast piecewise exponential generation
std::vector<double> rpieceexp_internal(int n,
                                       const std::vector<double>& time,
                                       const std::vector<double>& hazard) {

  int n_intervals = hazard.size();

  // Calculate cumulative hazard at time change points
  std::vector<double> cumulative_hazard(n_intervals + 1);
  cumulative_hazard[0] = 0.0;

  for (int i = 1; i <= n_intervals; ++i) {
    if (i < n_intervals) {
      double interval_length = time[i] - time[i-1];
      cumulative_hazard[i] = cumulative_hazard[i-1] + hazard[i-1] * interval_length;
    } else {
      cumulative_hazard[i] = R_PosInf;
    }
  }

  // Generate uniform random numbers
  NumericVector u = runif(n);

  std::vector<double> result(n);

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

    if (interval_idx >= n_intervals) {
      interval_idx = n_intervals - 1;
    }

    // Calculate result
    double start_time = time[interval_idx];
    double current_hazard = hazard[interval_idx];
    double prev_cumulative_hazard = cumulative_hazard[interval_idx];

    double remaining_hazard = target_hazard - prev_cumulative_hazard;
    result[i] = start_time + remaining_hazard / current_hazard;
  }

  return result;
}

//' Fast simulation data generation using C++ (Complete Implementation)
 //'
 //' This function generates simulation datasets for clinical trials with time-to-event endpoints
 //' using C++ for maximum performance.
 //'
 //' @param nsim A positive integer specifying the number of simulation iterations.
 //' @param N A positive integer specifying the total sample size per simulation.
 //' @param a_time A numeric vector of time points defining the accrual intervals.
 //' @param e_time A numeric vector of time points defining the survival time intervals.
 //' @param e_hazard A numeric vector of hazard rates for survival time.
 //' @param d_time A numeric vector of time points defining the dropout time intervals.
 //' @param d_hazard A numeric vector of hazard rates for dropout time.
 //' @param intensity A numeric vector of accrual intensities (or NULL).
 //' @param proportion A numeric vector of accrual proportions (or NULL).
 //' @return A list containing vectors for each column of the result.
 //' @export
 // [[Rcpp::export]]
 List simData_cpp(int nsim, int N,
                  NumericVector a_time,
                  NumericVector e_time, NumericVector e_hazard,
                  NumericVector d_time, NumericVector d_hazard,
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
   std::vector<double> a_time_vec = as<std::vector<double>>(a_time);
   std::vector<double> e_time_vec = as<std::vector<double>>(e_time);
   std::vector<double> e_hazard_vec = as<std::vector<double>>(e_hazard);
   std::vector<double> d_time_vec = as<std::vector<double>>(d_time);
   std::vector<double> d_hazard_vec = as<std::vector<double>>(d_hazard);

   std::vector<double> accrual_rates;
   bool use_intensity = false;

   if (has_intensity) {
     accrual_rates = as<std::vector<double>>(intensity);
     use_intensity = true;
   } else {
     accrual_rates = as<std::vector<double>>(proportion);
     use_intensity = false;
   }

   // Pre-calculate total sample size
   int total_n = nsim * N;

   // Pre-allocate result vectors
   std::vector<int> sim_ids(total_n);
   std::vector<double> accrual_times(total_n);
   std::vector<double> surv_times(total_n);
   std::vector<double> dropout_times(total_n);
   std::vector<double> tte(total_n);
   std::vector<double> total_times(total_n);
   std::vector<int> dropout_flags(total_n);

   // Generate all data at once for maximum efficiency
   std::vector<double> all_accrual = rpieceunif_internal(total_n, a_time_vec, accrual_rates, use_intensity);
   std::vector<double> all_surv = rpieceexp_internal(total_n, e_time_vec, e_hazard_vec);
   std::vector<double> all_dropout = rpieceexp_internal(total_n, d_time_vec, d_hazard_vec);

   // Fill result vectors
   for (int i = 0; i < total_n; ++i) {
     int sim_id = (i / N) + 1;

     sim_ids[i] = sim_id;
     accrual_times[i] = all_accrual[i];
     surv_times[i] = all_surv[i];
     dropout_times[i] = all_dropout[i];

     // Calculate derived values
     tte[i] = std::min(surv_times[i], dropout_times[i]);
     total_times[i] = accrual_times[i] + tte[i];
     dropout_flags[i] = (dropout_times[i] < surv_times[i]) ? 1 : 0;
   }

   // Return as R list
   return List::create(
     Named("simID") = wrap(sim_ids),
     Named("accrual.time") = wrap(accrual_times),
     Named("surv.time") = wrap(surv_times),
     Named("dropout.time") = wrap(dropout_times),
     Named("tte") = wrap(tte),
     Named("total") = wrap(total_times),
     Named("dropout") = wrap(dropout_flags)
   );
 }
