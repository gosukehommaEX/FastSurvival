#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cfloat>
using namespace Rcpp;

// Reusable scratch for the median survival core. Defined at file scope so the
// pointer-based implementation can take it by reference from the caller and
// avoid reallocating the per-event-time buffers per cell in a fused simulation
// loop.
struct MedsurvScratch {
  std::vector<double> event_time;
  std::vector<int> event_count;
  std::vector<double> at_risk;
  std::vector<double> hazard_incr;
};

// Forward declaration of the pointer-based implementation (defined below).
void medsurv_core_impl(const double*, const int*, int, double,
                       MedsurvScratch&, double*);

// Single-scan core for non-parametric median survival time.
//
// For each group g in 0, ..., ngroup - 1 the function walks the globally
// time-sorted data once, recording every distinct event time together with its
// event count and at-risk count, and the running Kaplan-Meier product. It
// returns, per group:
//   column 0: median survival time on the Kaplan-Meier estimate, following
//             the survival::survfit convention: the first event time at which
//             S(t) <= 0.5 (within a tolerance of sqrt(machine epsilon)), and,
//             when S(t) equals 0.5 on a flat stretch, the midpoint between that
//             event time and the next event time (or the last observed time if
//             no later event exists). NA if undefined, including a stretch
//             at 0.5 that ends at an infinite observed time, in which case
//             all the columns are NA.
//   column 1: Kaplan-Meier estimate S at the median
//   column 2: Greenwood sum up to and including the median,
//             sum_{s <= median} d(s) / (Y(s) (Y(s) - d(s)))
//   column 3: Ramlau-Hansen Epanechnikov kernel hazard at the median, with
//             bandwidth bw[g]
//
// The R wrapper turns these into the variance column 2 / column 3^2.
//
// Inputs time, event and grp must be sorted by time in ascending order; grp
// must be coded as 0, ..., ngroup - 1.

// [[Rcpp::export]]
NumericMatrix medsurv_core(NumericVector time, IntegerVector event,
                           IntegerVector grp, int ngroup, NumericVector bw) {
  int n_obs = time.size();
  NumericMatrix out(ngroup, 4);

  MedsurvScratch sc;
  std::vector<double> tg;
  std::vector<int> eg;
  tg.reserve(n_obs);
  eg.reserve(n_obs);

  for (int g = 0; g < ngroup; ++g) {
    tg.clear();
    eg.clear();
    for (int i = 0; i < n_obs; ++i) {
      if (grp[i] == g) {
        tg.push_back(time[i]);
        eg.push_back(event[i]);
      }
    }
    double res[4];
    medsurv_core_impl(tg.data(), eg.data(), (int) tg.size(), bw[g], sc, res);
    for (int c = 0; c < 4; ++c) out(g, c) = res[c];
  }

  return out;
}

// Pointer-based implementation with external linkage. The input must be a
// single group, sorted by time in ascending order. Writes
//   median, surv_med, greenwood, kernel_hazard
// into out[0..3]. The per-event-time buffers are taken from the caller-supplied
// scratch sc and reused across cells. Algorithm identical to the per-group body
// of the exported wrapper above.
void medsurv_core_impl(const double* time, const int* event, int n, double bw,
                       MedsurvScratch& sc, double* out) {
  std::vector<double>& event_time = sc.event_time;
  std::vector<int>& event_count = sc.event_count;
  std::vector<double>& at_risk = sc.at_risk;
  std::vector<double>& hazard_incr = sc.hazard_incr;
  event_time.clear();
  event_count.clear();
  at_risk.clear();
  hazard_incr.clear();

  // Tolerance used by survival::survfit when locating quantiles.
  const double tol = std::sqrt(DBL_EPSILON);

  int atrisk = n;
  double surv = 1.0;
  double median = NA_REAL;
  double surv_med = NA_REAL;
  int pos_med = -1;
  bool found = false;

  int i = 0;
  while (i < n) {
    double t0 = time[i];
    int n_event = 0;
    int n_record = 0;
    int j = i;
    while (j < n && time[j] == t0) {
      ++n_record;
      if (event[j] == 1) ++n_event;
      ++j;
    }
    int at_risk_here = atrisk;
    if (n_event > 0 && at_risk_here > 0) {
      event_time.push_back(t0);
      event_count.push_back(n_event);
      at_risk.push_back((double) at_risk_here);
      hazard_incr.push_back((double) n_event / (double) at_risk_here);
      surv *= (1.0 - (double) n_event / (double) at_risk_here);
      if (!found && surv <= 0.5 + tol) {
        found = true;
        median = t0;
        surv_med = surv;
        pos_med = (int) event_time.size() - 1;
      }
    }
    atrisk -= n_record;
    i = j;
  }

  // Flat stretch exactly at 0.5: survfit reports the midpoint between the
  // event time where S reaches 0.5 and the next event time, or the last
  // observed time when the curve stays at 0.5 to the end.
  if (found && std::fabs(surv_med - 0.5) < tol) {
    const int n_ev = (int) event_time.size();
    const double t_next = (pos_med + 1 < n_ev) ? event_time[pos_med + 1]
                                               : time[n - 1];
    if (std::isfinite(t_next)) {
      median = 0.5 * (event_time[pos_med] + t_next);
    } else {
      // The stretch ends at an infinite observed time (a subject without a
      // finite event or dropout time), so the midpoint is not defined. The
      // median and its variance components are NA, as when the curve does
      // not reach 0.5.
      found = false;
      median = NA_REAL;
      surv_med = NA_REAL;
    }
  }

  double greenwood = NA_REAL;
  double kernel_hazard = NA_REAL;

  if (found) {
    // Greenwood sum up to the median.
    double gw = 0.0;
    for (int k = 0; k <= pos_med; ++k) {
      double y = at_risk[k];
      int d = event_count[k];
      if (y - d > 0) gw += (double) d / (y * (y - (double) d));
    }
    greenwood = gw;

    // Epanechnikov kernel hazard at the median (Ramlau-Hansen smoother).
    double b = bw;
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
  }

  out[0] = median;
  out[1] = surv_med;
  out[2] = greenwood;
  out[3] = kernel_hazard;
}
