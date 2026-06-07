#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Single-scan core for the window mean survival time (WMST).
//
// For each group g in 0, ..., ngroup - 1 the function walks the time-sorted
// data once, builds the Kaplan-Meier estimate up to tau2, and returns:
//   column 0: WMST over the window [tau1, tau2], the area under the
//             Kaplan-Meier curve between tau1 and tau2
//   column 1: variance of the WMST estimate
//
// The variance uses the counting process (Greenwood) form
//   var(WMST) = sum_{e_k <= tau2} RW(e_k)^2 d(e_k) / (Y(e_k) (Y(e_k) - d(e_k))),
// where RW(e_k) is the remaining window area integral_{max(e_k, tau1)}^{tau2}
// S(u) du. For tau1 = 0 this reduces to the restricted mean survival time and
// its variance (matching survRM2).
//
// Inputs time, event and grp must be sorted by time in ascending order; grp is
// coded as 0, ..., ngroup - 1; event is 1 for an event and 0 for censoring.

// [[Rcpp::export]]
NumericMatrix wmst_core(NumericVector time, IntegerVector event, IntegerVector grp,
                        int ngroup, double tau1, double tau2) {
  int n = time.size();
  NumericMatrix out(ngroup, 2);

  for (int g = 0; g < ngroup; ++g) {
    int atrisk = 0;
    for (int i = 0; i < n; ++i) {
      if (grp[i] == g) ++atrisk;
    }

    double surv = 1.0;
    std::vector<double> e_time;
    std::vector<double> at_risk;
    std::vector<double> surv_after;
    std::vector<int> n_event;

    int i = 0;
    while (i < n) {
      double t0 = time[i];
      if (t0 > tau2) break;
      int d = 0, n_record = 0, j = i;
      while (j < n && time[j] == t0) {
        if (grp[j] == g) {
          ++n_record;
          if (event[j] == 1) ++d;
        }
        ++j;
      }
      int y = atrisk;
      if (d > 0 && y > 0) {
        surv *= (1.0 - (double) d / y);
        e_time.push_back(t0);
        at_risk.push_back((double) y);
        surv_after.push_back(surv);
        n_event.push_back(d);
      }
      atrisk -= n_record;
      i = j;
    }

    int m = (int) e_time.size();

    double wmst = 0.0;
    double prev_t = 0.0, prev_s = 1.0;
    for (int k = 0; k < m; ++k) {
      double lo = std::max(prev_t, tau1);
      double hi = std::min(e_time[k], tau2);
      if (hi > lo) wmst += prev_s * (hi - lo);
      prev_t = e_time[k];
      prev_s = surv_after[k];
    }
    {
      double lo = std::max(prev_t, tau1);
      double hi = tau2;
      if (hi > lo) wmst += prev_s * (hi - lo);
    }

    std::vector<double> tail(m, 0.0);
    if (m >= 1) {
      tail[m - 1] = surv_after[m - 1] * (tau2 - e_time[m - 1]);
      for (int k = m - 2; k >= 0; --k) {
        tail[k] = surv_after[k] * (e_time[k + 1] - e_time[k]) + tail[k + 1];
      }
    }

    double variance = 0.0;
    for (int k = 0; k < m; ++k) {
      double rw = (e_time[k] >= tau1) ? tail[k] : wmst;
      double y = at_risk[k];
      int d = n_event[k];
      if (y - d > 0) {
        variance += rw * rw * (double) d / (y * (y - (double) d));
      }
    }

    out(g, 0) = wmst;
    out(g, 1) = variance;
  }

  return out;
}
