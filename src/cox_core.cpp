// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

// Counts at one distinct event time within a stratum: the numbers of
// treatment and control subjects at risk just before the time, the number of
// treatment-group events, and the total number of events. Defined at file
// scope so that the fused analysis loop can pass a reusable buffer of it.
struct CoxEvent {
  double nT;
  double nC;
  double OT;
  double O;
};

// Forward declaration of the pointer-based implementation (defined below).
void cox_core_impl(const double*, const int*, const int*, const int*, int,
                   bool, int, std::vector<CoxEvent>&, double*);

//' Core two-group Cox model fit (C++ backend)
//'
//' @description
//' Internal C++ function that fits the Cox proportional hazards model with a
//' single treatment indicator by Newton-Raphson iteration on the risk-set
//' counts at the distinct event times, as described in
//' \code{\link{coxph_fast}}. Not intended to be called directly by users; use
//' \code{coxph_fast()} instead.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order (by stratum, then time, when \code{strata_sorted} is not
//'   empty).
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param strata_sorted An integer vector of stratum codes aligned with
//'   \code{time_sorted}, with the rows of each stratum contiguous, or an
//'   integer vector of length 0 for the unstratified model.
//' @param efron \code{TRUE} for the Efron handling of tied event times,
//'   \code{FALSE} for the Breslow handling.
//' @param iter_max The maximum number of iterations (20 in
//'   \code{coxph_fast()}).
//'
//' @return A numeric vector of length 6: the estimate of the log hazard
//'   ratio, its standard error, the log partial likelihood at 0 and at the
//'   estimate, the number of iterations, and a status code: 0 converged;
//'   1 not converged within \code{iter_max} iterations, in which case the
//'   estimate is the last iterate; 2 the estimate is infinite or cannot be
//'   computed, in which case the other elements except the number of
//'   iterations are \code{NA}.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector cox_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    const IntegerVector& strata_sorted,
    bool efron,
    int iter_max = 20
) {
  const int n = time_sorted.size();
  const int* s = strata_sorted.size() == 0 ? nullptr : strata_sorted.begin();
  std::vector<CoxEvent> ev;
  ev.reserve(n);
  double out[6];
  cox_core_impl(time_sorted.begin(), event_sorted.begin(), j_sorted.begin(),
                s, n, efron, iter_max, ev, out);
  return NumericVector(out, out + 6);
}

// One left-to-right scan over data sorted by time (by stratum, then time,
// when strata is not a null pointer; rows of a stratum contiguous). Subjects
// censored at an event time are at risk at that time. Times are tied only
// when they are exactly equal.
static void cox_event_table(const double* time, const int* event,
                            const int* trt, const int* strata, int n,
                            std::vector<CoxEvent>& ev) {
  ev.clear();
  int b = 0;
  while (b < n) {
    int e = n;
    if (strata != nullptr) {
      const int st = strata[b];
      e = b;
      while (e < n && strata[e] == st) ++e;
    }
    int nT = 0, nC = 0;
    for (int k = b; k < e; ++k) {
      if (trt[k] == 1) ++nT; else ++nC;
    }
    int i = b;
    while (i < e) {
      const double t = time[i];
      int dT = 0, dC = 0, cT = 0, cC = 0;
      int j = i;
      while (j < e && time[j] == t) {
        if (trt[j] == 1) {
          ++cT;
          if (event[j] == 1) ++dT;
        } else {
          ++cC;
          if (event[j] == 1) ++dC;
        }
        ++j;
      }
      if (dT + dC > 0) {
        ev.push_back(CoxEvent{(double) nT, (double) nC, (double) dT,
                              (double) (dT + dC)});
      }
      nT -= cT;
      nC -= cC;
      i = j;
    }
    b = e;
  }
}

// The score U(beta) = O_T - (expected treatment-group events) decreases in
// beta. As beta -> +Inf it tends to O_T - H, where H is the number of events
// at times with a treatment-group subject at risk; as beta -> -Inf it tends
// to O_T - L, where L is the number of events at times with no control
// subject at risk. The same limits hold under the Efron handling of ties.
// The maximum partial likelihood estimate is therefore finite if and only if
// L < O_T < H.
static bool cox_finite(const std::vector<CoxEvent>& ev) {
  double OT = 0.0, L = 0.0, H = 0.0;
  for (const CoxEvent& x : ev) {
    OT += x.OT;
    if (x.nC == 0.0) L += x.O;
    if (x.nT > 0.0) H += x.O;
  }
  return L < OT && OT < H;
}

// Sum of the logarithms of positive numbers, kept as a mantissa in [0.5, 1)
// and a binary exponent so that std::log() is called once instead of once
// per term.
struct CoxLogSum {
  double m = 1.0;
  long e = 0;
  void add(double x) {
    int ex;
    m = std::frexp(m * x, &ex);
    e += ex;
  }
  double value() const {
    return std::log(m) + (double) e * 0.69314718055994530942;
  }
};

// Log partial likelihood, score, and observed information at beta. Under the
// Breslow handling, the term of an event time with O events is
// O_T beta - O log(n_C + n_T r), r = exp(beta). Under the Efron handling, the
// r-th of the O tied events (r = 0, ..., O - 1) has a fraction r / O of the
// treatment and control events at that time removed from the risk set.
// Returns false when a value is not finite.
static bool cox_eval(const std::vector<CoxEvent>& ev, bool efron, double beta,
                     double& ll, double& U, double& I) {
  const double r = std::exp(beta);
  if (!std::isfinite(r)) return false;
  double u = 0.0, inf = 0.0, OT = 0.0;
  CoxLogSum ls;
  for (const CoxEvent& x : ev) {
    OT += x.OT;
    const int d = (int) x.O;
    if (!efron || d == 1) {
      const double a = x.nT * r;
      const double inv = 1.0 / (x.nC + a);
      const double p = a * inv;
      const double q = x.nC * inv;
      for (int k = 0; k < d; ++k) ls.add(x.nC + a);
      u -= x.O * p;
      inf += x.O * p * q;
    } else {
      const double OC = x.O - x.OT;
      for (int k = 0; k < d; ++k) {
        const double f = (double) k / x.O;
        const double nT = x.nT - f * x.OT;
        const double nC = x.nC - f * OC;
        const double a = nT * r;
        const double den = nC + a;
        const double inv = 1.0 / den;
        const double p = a * inv;
        const double q = nC * inv;
        ls.add(den);
        u -= p;
        inf += p * q;
      }
    }
  }
  ll = OT * beta - ls.value();
  U = u + OT;
  I = inf;
  return std::isfinite(ll) && std::isfinite(U) && std::isfinite(I);
}

// Newton-Raphson iteration from beta = 0. A step after which the log partial
// likelihood is lower, or not finite, is halved, and the convergence check is
// skipped right after a halved step. The points evaluated so far bracket the
// root (a: the largest point with a positive score, b: the smallest point
// with a negative score); from the first halved step on, a trial value
// outside (a, b) is replaced by the midpoint of the bracket. Iteration stops
// when the relative change in the log partial likelihood is at most 1e-9.
// The estimate is the last iterate and the information is evaluated there,
// also when the iteration stops at iter_max without converging.
//
// out: coef, se, loglik at 0, loglik at coef, iterations, status
// (0 converged, 1 not converged, 2 infinite or undefined estimate).
static void cox_fit(const std::vector<CoxEvent>& ev, bool efron, int iter_max,
                    double* out) {
  const double eps = 1e-9;
  const double inf = std::numeric_limits<double>::infinity();
  for (int k = 0; k < 6; ++k) out[k] = NA_REAL;
  out[4] = 0.0;
  out[5] = 2.0;
  if (ev.empty() || !cox_finite(ev)) return;

  double beta = 0.0, ll = 0.0, U = 0.0, I = 0.0;
  if (!cox_eval(ev, efron, beta, ll, U, I) || !(I > 0.0)) return;
  out[2] = ll;

  double a = -inf, b = inf;
  bool guarded = false;
  auto bracket = [&](double x, double u) {
    if (u > 0.0 && x > a) a = x;
    if (u < 0.0 && x < b) b = x;
  };
  auto clamp = [&](double x) {
    if (guarded && !(x > a && x < b) && std::isfinite(a) && std::isfinite(b)) {
      return 0.5 * (a + b);
    }
    return x;
  };

  bracket(beta, U);
  double newbeta = beta + U / I;
  double ll_new = NA_REAL, U_new = NA_REAL, I_new = NA_REAL;
  bool halving = false;
  int status = 1, iter = 0;
  for (iter = 1; iter <= iter_max; ++iter) {
    const bool ok = cox_eval(ev, efron, newbeta, ll_new, U_new, I_new);
    if (ok) bracket(newbeta, U_new);
    if (ok && !halving && std::fabs(1.0 - ll / ll_new) <= eps) {
      status = 0;
      break;
    }
    if (iter == iter_max) break;
    if (!ok || ll_new < ll) {
      halving = true;
      guarded = true;
      newbeta = clamp(0.5 * (newbeta + beta));
    } else {
      halving = false;
      ll = ll_new;
      beta = newbeta;
      double step = U_new / I_new;
      if (!std::isfinite(step)) {
        step = (U_new > 0.0 ? 1.0 : -1.0) * std::max(1.0, std::fabs(beta));
      }
      newbeta = clamp(beta + step);
    }
  }
  out[0] = newbeta;
  out[1] = 1.0 / std::sqrt(I_new);
  out[3] = ll_new;
  out[4] = (double) (iter > iter_max ? iter_max : iter);
  out[5] = (double) status;
}

// Pointer-based implementation with external linkage, used by cox_core() and
// by the fused analysis loop. Writes coef, se, loglik at 0, loglik at coef,
// iterations, and status into out[0..5]. The per-event-time counts are stored
// in the caller-supplied buffer ev, which is cleared on entry and reused
// across cells. When strata is a null pointer the whole input is one stratum.
void cox_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    const int* strata_sorted,
    int n,
    bool efron,
    int iter_max,
    std::vector<CoxEvent>& ev,
    double* out
) {
  cox_event_table(time_sorted, event_sorted, j_sorted, strata_sorted, n, ev);
  cox_fit(ev, efron, iter_max, out);
}
