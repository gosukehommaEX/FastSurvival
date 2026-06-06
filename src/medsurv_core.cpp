#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// Single-scan core for non-parametric median survival time, supporting two
// variance methods.
//
// For each group g in 0, ..., ngroup - 1 the function walks the globally
// time-sorted data once, recording every distinct event time together with its
// event count and at-risk count, and the running Kaplan-Meier product. It
// returns, per group:
//   column 0: median survival time, inf{t : S(t) <= 0.5} on the Kaplan-Meier
//             estimate (NA if undefined). This is the point estimate for both
//             variance methods.
//   column 1: Kaplan-Meier estimate S at the median
//   column 2: Greenwood sum up to and including the median,
//             sum_{s <= median} d(s) / (Y(s) (Y(s) - d(s)))   [method "km"]
//   column 3: Ramlau-Hansen Epanechnikov kernel hazard at the median, with
//             bandwidth bw[g]                                   [method "km"]
//   column 4: nph variance increment sum up to the median,
//             sum_{s <= median} sum_{j = 0}^{d(s) - 1} 1 / (Y(s) - j)^2
//                                                               [method "nph"]
//   column 5: nph local constant hazard at the median, total events divided by
//             total person-time over a window of ceil(sqrt(m)) * 2 event rows
//             on each side of the median                        [method "nph"]
//
// The R wrapper turns these into a variance: method "km" uses
// column 2 / column 3^2, method "nph" uses column 4 / column 5^2. The nph
// branch reproduces the computations in nph::nphparams (param_type "Q",
// haz_method "local").
//
// Inputs time, event and grp must be sorted by time in ascending order; grp
// must be coded as 0, ..., ngroup - 1.

// [[Rcpp::export]]
NumericMatrix medsurv_core(NumericVector time, IntegerVector event,
                           IntegerVector grp, int ngroup, NumericVector bw) {
  int n_obs = time.size();
  NumericMatrix out(ngroup, 6);

  for (int g = 0; g < ngroup; ++g) {
    int atrisk = 0;
    for (int i = 0; i < n_obs; ++i) {
      if (grp[i] == g) ++atrisk;
    }

    double surv = 1.0;
    double median = NA_REAL;
    double surv_med = NA_REAL;
    int pos_med = -1;
    bool found = false;

    std::vector<double> event_time;
    std::vector<int> event_count;
    std::vector<double> at_risk;
    std::vector<double> hazard_incr;

    int i = 0;
    while (i < n_obs) {
      double t0 = time[i];
      int n_event = 0;
      int n_record = 0;
      int j = i;
      while (j < n_obs && time[j] == t0) {
        if (grp[j] == g) {
          ++n_record;
          if (event[j] == 1) ++n_event;
        }
        ++j;
      }
      int at_risk_here = atrisk;
      if (n_event > 0 && at_risk_here > 0) {
        event_time.push_back(t0);
        event_count.push_back(n_event);
        at_risk.push_back((double) at_risk_here);
        hazard_incr.push_back((double) n_event / (double) at_risk_here);
        surv *= (1.0 - (double) n_event / (double) at_risk_here);
        if (!found && surv <= 0.5) {
          found = true;
          median = t0;
          surv_med = surv;
          pos_med = (int) event_time.size() - 1;
        }
      }
      atrisk -= n_record;
      i = j;
    }

    double greenwood = NA_REAL;
    double kernel_hazard = NA_REAL;
    double nph_var = NA_REAL;
    double local_hazard = NA_REAL;

    if (found) {
      // Greenwood sum up to the median (km method).
      double gw = 0.0;
      for (int k = 0; k <= pos_med; ++k) {
        double y = at_risk[k];
        int d = event_count[k];
        if (y - d > 0) gw += (double) d / (y * (y - (double) d));
      }
      greenwood = gw;

      // Epanechnikov kernel hazard at the median (km method).
      double b = bw[g];
      if (b > 0 && !event_time.empty()) {
        double acc = 0.0;
        for (size_t k = 0; k < event_time.size(); ++k) {
          double u = (event_time[k] - median) / b;
          if (std::fabs(u) < 1.0) {
            acc += 0.75 * (1.0 - u * u) / b * hazard_incr[k];
          }
        }
        kernel_hazard = acc;
      }

      // nph variance increment sum up to the median (nph method).
      double vs = 0.0;
      for (int k = 0; k <= pos_med; ++k) {
        double y = at_risk[k];
        int d = event_count[k];
        for (int jj = 0; jj < d; ++jj) {
          double denom = y - (double) jj;
          vs += 1.0 / (denom * denom);
        }
      }
      nph_var = vs;

      // nph local constant hazard at the median (nph method).
      int n_rows = (int) event_time.size();
      double total_event = 0.0;
      for (int k = 0; k < n_rows; ++k) total_event += event_count[k];
      int delta = (int) (std::ceil(std::sqrt(total_event)) * 2.0);
      int pos1 = pos_med + 1;
      int low = pos1 - delta;
      if (low < 1) low = 1;
      int up = pos1 + delta;
      if (up > n_rows) up = n_rows;
      double event_sum = 0.0;
      double persontime_sum = 0.0;
      double prev_t = 0.0;
      for (int k = 0; k < n_rows; ++k) {
        double persontime = (event_time[k] - prev_t) * at_risk[k];
        prev_t = event_time[k];
        if (k + 1 >= low && k + 1 <= up) {
          event_sum += event_count[k];
          persontime_sum += persontime;
        }
      }
      if (persontime_sum > 0) local_hazard = event_sum / persontime_sum;
    }

    out(g, 0) = median;
    out(g, 1) = surv_med;
    out(g, 2) = greenwood;
    out(g, 3) = kernel_hazard;
    out(g, 4) = nph_var;
    out(g, 5) = local_hazard;
  }

  return out;
}
