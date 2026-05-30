// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' Core restricted mean survival time computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the restricted mean survival time
//' (RMST) and its Greenwood-type variance for a single group up to a horizon
//' \code{tau}. A single left-to-right scan over the sorted data accumulates
//' the Kaplan-Meier product and the area under the survival step function,
//' storing one summary per distinct event time (the cumulative area to the
//' left of that time and the Greenwood increment). A second short pass over
//' the saved summaries forms the variance as a sum of squared right-side
//' areas weighted by the Greenwood increments. Accepts the event vector as a
//' numeric (double) vector to avoid an \code{as.integer()} copy in R. Not
//' intended to be called directly by users; use \code{rmst_fast()} instead.
//'
//' @param t_sorted A numeric vector of follow-up times sorted in ascending
//'   order.
//' @param e_sorted A numeric vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{t_sorted}. Accepted as numeric to
//'   avoid an integer-conversion copy on the R side.
//' @param tau A single numeric value specifying the restriction horizon.
//'
//' @return A numeric vector of length 2: \code{c(rmst, var_rmst)}, where
//'   \code{rmst} is the restricted mean survival time and \code{var_rmst} is
//'   its variance, so that SE[RMST] = sqrt(var_rmst). Returns
//'   \code{c(tau, 0.0)} when no events occur up to \code{tau}.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector rmst_core(
    const NumericVector& t_sorted,
    const NumericVector& e_sorted,
    double tau
) {
  const int n = t_sorted.size();
  if (n == 0) return NumericVector::create(R_NaReal, R_NaReal);

  // Per-event-time summary stored as a struct for cache locality.
  // cum_area: area under S from 0 to this event time (entering survival
  // level). g: Greenwood increment d / (n_risk * (n_risk - d)).
  // Reserve n: an upper bound on distinct event times. No prepass needed.
  struct EvSummary {
    double cum_area;
    double g;
  };
  std::vector<EvSummary> ev;
  ev.reserve(n);

  double surv   = 1.0;        // KM survival over the interval entering this time
  double rmst   = 0.0;        // area under S accumulated so far
  double prev_t = 0.0;        // left endpoint of the current interval
  double n_risk = (double)n;

  // Single scan: integrate the step function and store per-event summaries
  int i = 0;
  while (i < n) {
    const double t = t_sorted[i];

    // Stop integrating beyond the horizon
    if (t > tau) break;

    // Consume the tied block at time t: count deaths and block size
    int d = 0, c_blk = 0;
    int j = i;
    while (j < n && t_sorted[j] == t) {
      ++c_blk;
      if (e_sorted[j] >= 0.5) ++d;
      ++j;
    }

    // Area over [prev_t, t) at the survival level entering this time
    rmst   += surv * (t - prev_t);
    prev_t  = t;

    if (d > 0) {
      const double dd = (double)d;
      if (n_risk > dd) {
        ev.push_back(EvSummary{rmst, dd / (n_risk * (n_risk - dd))});
        surv *= 1.0 - dd / n_risk;
      } else {
        // n_risk == d: survival drops to zero; Greenwood increment is zero
        ev.push_back(EvSummary{rmst, 0.0});
        surv = 0.0;
      }
    }

    n_risk -= (double)c_blk;
    i       = j;
  }

  // Tail area from the last processed time to tau
  if (prev_t < tau) {
    rmst += surv * (tau - prev_t);
  }

  // Second pass: variance = sum of (area to the right of t_i)^2 * g_i
  double var_rmst = 0.0;
  const std::size_t K = ev.size();
  for (std::size_t k = 0; k < K; ++k) {
    const double a_right = rmst - ev[k].cum_area;
    var_rmst += a_right * a_right * ev[k].g;
  }

  return NumericVector::create(rmst, var_rmst);
}
