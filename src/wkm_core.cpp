#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// Single-scan core for the weighted Kaplan-Meier (Pepe-Fleming) two-sample test.
//
// This reproduces the computation in nphsim::wkm.Stat. The integrand is a step
// function that is constant between consecutive observation times, so the
// midpoint Riemann sum over the sorted observation times is the exact value of
// the step-function integral.
//
// For each interval between consecutive sorted observation times the core
// evaluates, on the open interval, the per-group Kaplan-Meier survival
// estimates (events only), the per-group Kaplan-Meier censoring survival
// estimates, and the pooled Kaplan-Meier survival estimate. The Pepe-Fleming
// combined weight is
//   w(t) = n G1(t) G2(t) / (n1 G1(t) + n2 G2(t)),
// where Gj is the censoring survival estimate in group j. The weight type
// argument selects w (0), sqrt(w) (1) or a constant 1 (2).
//
// Returns a numeric vector with:
//   [0] num_raw  = sum_k w_k (S1_k - S2_k) width_k, the weighted integrated
//                  survival difference (treatment minus control)
//   [1] variance = Pepe-Fleming variance estimate of the score statistic
//   [2] n1       = number of treatment subjects (group 1)
//   [3] n2       = number of control subjects (group 0)
//
// Inputs time, event and grp must be sorted by time in ascending order; grp is
// coded as 1 for treatment and 0 for control; event is 1 for an event and 0
// for censoring.

// [[Rcpp::export]]
NumericVector wkm_core(NumericVector time, IntegerVector event, IntegerVector grp,
                       int weight_type) {
  int n = time.size();
  int n1 = 0, n2 = 0;
  for (int i = 0; i < n; ++i) {
    if (grp[i] == 1) ++n1; else ++n2;
  }

  std::vector<double> width(n), wt(n), s1(n), s2(n), sp(n);

  double surv1 = 1.0, surv2 = 1.0, cen1 = 1.0, cen2 = 1.0, survp = 1.0;
  int y1 = n1, y2 = n2, yp = n;
  int absorbed = 0;

  for (int k = 0; k < n; ++k) {
    double t_prev = (k >= 1) ? time[k - 1] : 0.0;
    double t_cur = time[k];

    if (k >= 1) {
      while (absorbed < n && time[absorbed] <= t_prev) {
        double tau = time[absorbed];
        int d1ev = 0, d1cen = 0, d2ev = 0, d2cen = 0, blk = 0;
        int j = absorbed;
        while (j < n && time[j] == tau) {
          if (grp[j] == 1) {
            if (event[j] == 1) ++d1ev; else ++d1cen;
          } else {
            if (event[j] == 1) ++d2ev; else ++d2cen;
          }
          ++blk;
          ++j;
        }
        int dpev = d1ev + d2ev;
        if (y1 > 0) {
          if (d1ev > 0) surv1 *= (1.0 - (double) d1ev / y1);
          if (d1cen > 0) cen1 *= (1.0 - (double) d1cen / y1);
        }
        if (y2 > 0) {
          if (d2ev > 0) surv2 *= (1.0 - (double) d2ev / y2);
          if (d2cen > 0) cen2 *= (1.0 - (double) d2cen / y2);
        }
        if (yp > 0 && dpev > 0) survp *= (1.0 - (double) dpev / yp);
        y1 -= (d1ev + d1cen);
        y2 -= (d2ev + d2cen);
        yp -= blk;
        absorbed = j;
      }
    }

    width[k] = t_cur - t_prev;
    s1[k] = surv1;
    s2[k] = surv2;
    sp[k] = survp;

    double denom = (double) n1 * cen1 + (double) n2 * cen2;
    double w;
    if (cen1 + cen2 == 0.0 || denom <= 0.0) {
      w = 0.0;
    } else {
      w = ((double) n * cen1 * cen2) / denom;
    }
    if (weight_type == 1) {
      w = std::sqrt(w);
    } else if (weight_type == 2) {
      w = 1.0;
    }
    wt[k] = w;
  }

  // Absorb any remaining observations to obtain the final pooled survival.
  while (absorbed < n) {
    double tau = time[absorbed];
    int dpev = 0, blk = 0, j = absorbed;
    while (j < n && time[j] == tau) {
      if (event[j] == 1) ++dpev;
      ++blk;
      ++j;
    }
    if (yp > 0 && dpev > 0) survp *= (1.0 - (double) dpev / yp);
    yp -= blk;
    absorbed = j;
  }
  double sp_final = survp;

  double num_raw = 0.0;
  for (int k = 0; k < n; ++k) {
    num_raw += wt[k] * (s1[k] - s2[k]) * width[k];
  }

  std::vector<double> a_seq(n);
  double cum = 0.0;
  for (int k = 0; k < n; ++k) {
    cum += width[k] * wt[k] * sp[k];
    a_seq[k] = cum;
  }
  double total = (n > 0) ? a_seq[n - 1] : 0.0;

  double variance = 0.0;
  for (int k = 0; k < n; ++k) {
    double a_rem = total - a_seq[k];
    double sp_next = (k < n - 1) ? sp[k + 1] : sp_final;
    double d_sm = sp_next - sp[k];
    if (wt[k] > 0.0 && sp[k] > 0.0) {
      double term = (a_rem * a_rem) / (sp[k] * sp[k] * wt[k]) * d_sm;
      variance -= term;
    }
  }

  NumericVector out(4);
  out[0] = num_raw;
  out[1] = variance;
  out[2] = n1;
  out[3] = n2;
  return out;
}
