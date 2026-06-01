#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Compute the Kaplan-Meier survival probability and Greenwood variance at a
// milestone timepoint for two groups in a single pass over the data.
//
// The two groups must be coded as 0 (control) and 1 (treatment). When
// presorted is TRUE the input is assumed to be sorted by time in ascending
// order, in which case each group's records remain in ascending order and no
// internal sort is performed. Censored observations sharing an event time are
// treated as being at risk for that event (the standard Kaplan-Meier tie
// convention used by survfit).
//
// [[Rcpp::export]]
List milestone_core(NumericVector time, IntegerVector status,
                    IntegerVector grp, double tau, bool presorted) {
  int n = time.size();

  // Build the time-ascending traversal order.
  std::vector<int> ord(n);
  for (int i = 0; i < n; i++) ord[i] = i;
  if (!presorted) {
    std::stable_sort(ord.begin(), ord.end(),
                     [&](int a, int b) { return time[a] < time[b]; });
  }

  // Total sample size per group, used to initialise the at-risk counts.
  int n0 = 0, n1 = 0;
  for (int i = 0; i < n; i++) {
    if (grp[i] == 0) n0++;
    else n1++;
  }

  int r0 = n0, r1 = n1;          // current number at risk per group
  double surv0 = 1.0, surv1 = 1.0;
  double vsum0 = 0.0, vsum1 = 0.0;  // Greenwood sum d / (Y (Y - d))
  int de0 = 0, de1 = 0;          // number of events with time <= tau per group

  int i = 0;
  while (i < n) {
    double t_cur = time[ord[i]];
    if (t_cur > tau) break;

    // Risk sets just before processing the current time.
    int yb0 = r0, yb1 = r1;

    // Tally events and total records at the current time for each group.
    int d0 = 0, d1 = 0, m0 = 0, m1 = 0;
    int j = i;
    while (j < n && time[ord[j]] == t_cur) {
      int g = grp[ord[j]];
      int s = status[ord[j]];
      if (g == 0) {
        m0++;
        if (s == 1) d0++;
      } else {
        m1++;
        if (s == 1) d1++;
      }
      j++;
    }

    // Kaplan-Meier and Greenwood updates for control.
    if (d0 > 0 && yb0 > 0) {
      surv0 *= 1.0 - (double) d0 / yb0;
      if (yb0 > d0) vsum0 += (double) d0 / ((double) yb0 * (yb0 - d0));
      else vsum0 = R_PosInf;
      de0 += d0;
    }
    // Kaplan-Meier and Greenwood updates for treatment.
    if (d1 > 0 && yb1 > 0) {
      surv1 *= 1.0 - (double) d1 / yb1;
      if (yb1 > d1) vsum1 += (double) d1 / ((double) yb1 * (yb1 - d1));
      else vsum1 = R_PosInf;
      de1 += d1;
    }

    // Remove the current time's records from the risk sets.
    r0 -= m0;
    r1 -= m1;
    i = j;
  }

  double v0 = surv0 * surv0 * vsum0;
  double v1 = surv1 * surv1 * vsum1;

  return List::create(
    _["surv0"] = surv0, _["surv1"] = surv1,
    _["var0"] = v0, _["var1"] = v1,
    _["n0"] = n0, _["n1"] = n1,
    _["n_event0"] = de0, _["n_event1"] = de1,
    _["n_risk0"] = r0, _["n_risk1"] = r1
  );
}
