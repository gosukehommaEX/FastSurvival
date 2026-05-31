// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ------------------------------------------------------------------ //
//  Internal: administrative censoring at a calendar cutoff + time sort
// ------------------------------------------------------------------ //
// For a single simulation, applies an interim/sequential look at calendar
// time 'cutoff' and returns the look data sorted by observed time, ready to
// be passed to the analysis cores (logrank_core, pihe_core, rmst_core,
// km_core) with presorted = TRUE.
//
// For each subject i with accrual time a_i, observed time-from-entry tte_i
// (= min(survival, dropout)), and full-follow-up event indicator event_i:
//
//   enrolled_i = (a_i <= cutoff)
//   t_look_i   = min(tte_i, cutoff - a_i)                 (enrolled only)
//   e_look_i   = event_i * (a_i + tte_i <= cutoff)        (enrolled only)
//
// A subject whose natural event/dropout occurs after the cutoff is
// administratively censored at cutoff (t_look = cutoff - a_i, e_look = 0).
// Subjects not yet enrolled (a_i > cutoff) are dropped from the look.
//
// The returned vectors contain only enrolled subjects, sorted in ascending
// order of t_look. Ties in t_look are left in arbitrary order, which is
// sufficient because the downstream cores process all rows sharing an event
// time as one atomic block. The element 'idx' gives, for each output
// position, the 1-based index of that subject within the input vectors,
// allowing callers to realign per-subject covariates (e.g. subgroup labels)
// to the sorted order.
//
// [[Rcpp::export]]
List analysis_cut_core(
    const NumericVector& accrual,  // a_i, one simulation
    const NumericVector& tte,      // observed time from entry, one simulation
    const IntegerVector& event,    // full-follow-up event indicator (1/0)
    const IntegerVector& j,        // treatment indicator (1 = trt, 0 = ctrl)
    double cutoff                  // look calendar time
) {
  const int n = accrual.size();

  // Collect enrolled subjects with their look time/event/group and original
  // (within-simulation) index, so callers can realign per-subject covariates
  // such as subgroup labels to the sorted output order.
  std::vector<int>    orig;
  std::vector<double> t_look;
  std::vector<int>    e_look;
  std::vector<int>    j_look;
  orig.reserve(n);
  t_look.reserve(n);
  e_look.reserve(n);
  j_look.reserve(n);

  int n_event = 0;
  for (int i = 0; i < n; ++i) {
    const double a = accrual[i];
    if (a > cutoff) continue;  // not yet enrolled at this look

    const double full   = tte[i];
    const double avail  = cutoff - a;
    const bool   before = (a + full <= cutoff);
    const double t      = before ? full : avail;
    const int    e      = before ? event[i] : 0;

    orig.push_back(i);
    t_look.push_back(t);
    e_look.push_back(e);
    j_look.push_back(j[i]);
    if (e == 1) ++n_event;
  }

  const int m = static_cast<int>(t_look.size());

  // Sort indices by t_look ascending
  std::vector<int> ord(m);
  for (int k = 0; k < m; ++k) ord[k] = k;
  std::sort(ord.begin(), ord.end(),
            [&t_look](int p, int q) { return t_look[p] < t_look[q]; });

  NumericVector time_out(m);
  IntegerVector event_out(m);
  IntegerVector j_out(m);
  IntegerVector idx_out(m);  // 1-based original within-simulation indices
  for (int k = 0; k < m; ++k) {
    const int p   = ord[k];
    time_out[k]   = t_look[p];
    event_out[k]  = e_look[p];
    j_out[k]      = j_look[p];
    idx_out[k]    = orig[p] + 1;
  }

  return List::create(
    Named("time")    = time_out,
    Named("event")   = event_out,
    Named("j")       = j_out,
    Named("idx")     = idx_out,
    Named("n")       = m,
    Named("n_event") = n_event
  );
}

// ------------------------------------------------------------------ //
//  Internal: n-th smallest event calendar time (event-driven look)
// ------------------------------------------------------------------ //
// Given the calendar times at which events occur in one simulation
// (cal_event = accrual + tte for subjects with event == 1), returns the
// calendar time at which the cumulative event count first reaches 'target'.
// This is the target-th order statistic of cal_event, used as the look
// cutoff for an information-based (event-count) interim analysis.
//
// Uses std::nth_element for average O(n) selection rather than a full sort.
// Returns R_PosInf when 'target' exceeds the number of events available, so
// the caller can detect that the target was never reached (the full data are
// then used and the look is flagged).
//
// [[Rcpp::export]]
double nth_event_time_core(
    const NumericVector& cal_event,  // calendar times of events, one simulation
    int target                       // target cumulative event count (1-based)
) {
  const int n = cal_event.size();
  if (target < 1 || target > n) {
    return R_PosInf;
  }

  std::vector<double> v(cal_event.begin(), cal_event.end());
  std::nth_element(v.begin(), v.begin() + (target - 1), v.end());
  return v[target - 1];
}
