// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Core weighted log-rank computation on pooled sorted vectors (C++ backend)
//'
//' @description
//' Internal C++ function that computes the weighted log-rank numerator U and
//' variance V directly from a pooled sorted dataset and a group indicator.
//' Four weight schemes are supported, selected by \code{scheme}: Fleming-
//' Harrington G(rho, gamma), the modestly-weighted log-rank test, Gehan-
//' Breslow, and Tarone-Ware. The weights are functions of the left-continuous
//' pooled Kaplan-Meier estimate (for Fleming-Harrington and the modestly-
//' weighted test) or of the at-risk count (for Gehan-Breslow and Tarone-Ware).
//' Tied event times are processed atomically. Not intended to be called
//' directly by users; use \code{survdiff_fast()} with the \code{weight}
//' argument instead.
//'
//' @details
//' The scheme codes are: 0 = Fleming-Harrington G(rho, gamma) with weight
//' \code{S_minus^rho * (1 - S_minus)^gamma}; 1 = modestly-weighted log-rank
//' with weight \code{min(1 / S_minus, max_weight)}, where \code{max_weight}
//' is the reciprocal of the smallest right-continuous pooled Kaplan-Meier
//' value at or after \code{t_star} (and is 1 when \code{t_star = 0});
//' 2 = Gehan-Breslow with weight \code{n_j}; 3 = Tarone-Ware with weight
//' \code{sqrt(n_j)}. Here \code{S_minus} is the left-continuous pooled
//' Kaplan-Meier estimate just prior to each event time, initialized at 1.
//' The modestly-weighted scheme requires the largest weight to be known before
//' accumulation, so its \code{max_weight} is obtained in a first pass over the
//' event times; the other schemes accumulate in a single pass.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param scheme An integer weight-scheme code: 0 = Fleming-Harrington,
//'   1 = modestly-weighted, 2 = Gehan-Breslow, 3 = Tarone-Ware.
//' @param rho A numeric Fleming-Harrington first parameter (used when
//'   \code{scheme = 0}).
//' @param gamma A numeric Fleming-Harrington second parameter (used when
//'   \code{scheme = 0}).
//' @param t_star A numeric timepoint for the modestly-weighted scheme (used
//'   when \code{scheme = 1}).
//'
//' @return A numeric vector of length 3: \code{c(O1, U, V)}, where O1 is the
//'   observed number of events in the treatment group, U is the weighted
//'   log-rank numerator sum w (O1 - E1), and V is the weighted variance.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector weighted_logrank_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    int scheme,
    double rho,
    double gamma,
    double t_star
) {
  const int n = time_sorted.size();

  // At-risk counts per group, initialized to the group totals
  int n1_init = 0, n0_init = 0;
  for (int k = 0; k < n; ++k) {
    if (j_sorted[k] == 1) ++n1_init; else ++n0_init;
  }

  // ---------------------------------------------------------------- //
  //  First pass for the modestly-weighted scheme: determine max_weight
  //  as the reciprocal of the smallest right-continuous pooled KM value
  //  at or after t_star. When t_star = 0, max_weight is 1.
  // ---------------------------------------------------------------- //
  double max_weight = 1.0;
  if (scheme == 1 && t_star > 0.0) {
    int nrisk = n1_init + n0_init;   // pooled at-risk count
    double s = 1.0;                  // right-continuous pooled KM
    double s_min_post = R_PosInf;
    int i = 0;
    while (i < n) {
      const double t = time_sorted[i];
      int d = 0, c = 0;
      int jj = i;
      while (jj < n && time_sorted[jj] == t) {
        ++c;
        if (event_sorted[jj] == 1) ++d;
        ++jj;
      }
      if (d > 0 && nrisk > 0) {
        s *= (1.0 - (double)d / (double)nrisk);
        if (t >= t_star && s < s_min_post) s_min_post = s;
      }
      nrisk -= c;
      i = jj;
    }
    if (R_finite(s_min_post) && s_min_post > 0.0) {
      max_weight = 1.0 / s_min_post;
    } else if (s > 0.0) {
      // No events at or after t_star: fall back to the final KM value
      max_weight = 1.0 / s;
    }
  }

  // ---------------------------------------------------------------- //
  //  Main accumulation pass
  // ---------------------------------------------------------------- //
  int n1 = n1_init, n0 = n0_init;
  double s_minus = 1.0;   // left-continuous pooled KM just prior to t
  double s       = 1.0;   // right-continuous pooled KM

  double O1 = 0.0, U = 0.0, V = 0.0;

  int i = 0;
  while (i < n) {
    const double t = time_sorted[i];

    // Consume tied block at time t, splitting deaths and censorings per group
    int d1 = 0, d0 = 0, c1 = 0, c0 = 0;
    int jj = i;
    while (jj < n && time_sorted[jj] == t) {
      if (j_sorted[jj] == 1) {
        ++c1;
        if (event_sorted[jj] == 1) ++d1;
      } else {
        ++c0;
        if (event_sorted[jj] == 1) ++d0;
      }
      ++jj;
    }

    const int d  = d1 + d0;
    const int nj = n1 + n0;

    if (d > 0 && nj > 1) {
      const double dn1 = (double)n1;
      const double dnj = (double)nj;
      const double dd  = (double)d;

      // Weight from the left-continuous KM (s_minus) or at-risk count
      double w;
      if (scheme == 0) {
        w = std::pow(s_minus, rho) * std::pow(1.0 - s_minus, gamma);
      } else if (scheme == 1) {
        w = 1.0 / s_minus;
        if (w > max_weight) w = max_weight;
      } else if (scheme == 2) {
        w = dnj;
      } else {
        w = std::sqrt(dnj);
      }

      const double e1  = dd * dn1 / dnj;
      const double v1  = dd * dn1 * (dnj - dn1) * (dnj - dd) /
                         (dnj * dnj * (dnj - 1.0));

      O1 += d1;
      U  += w * ((double)d1 - e1);
      V  += w * w * v1;
    }

    // Update KM after using s_minus at this time: s_minus -> s -> next s_minus
    if (d > 0 && nj > 0) {
      s *= (1.0 - (double)d / (double)nj);
    }
    s_minus = s;

    // Decrement at-risk counts by this block's size
    n1 -= c1;
    n0 -= c0;
    i   = jj;
  }

  return NumericVector::create(O1, U, V);
}
