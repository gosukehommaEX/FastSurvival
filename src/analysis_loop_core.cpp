// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Forward declarations of the existing statistic cores. These are ordinary
// C++ functions with external linkage defined in the other source files of
// this package; RcppExports.cpp already relies on the same cross-translation-
// unit linkage, so resolving them here is guaranteed. Signatures must match
// the definitions exactly (const-reference arguments).
NumericVector logrank_core(const NumericVector&, const IntegerVector&,
                           const IntegerVector&);
NumericVector weighted_logrank_core(const NumericVector&, const IntegerVector&,
                                    const IntegerVector&, int, double, double,
                                    double);
NumericVector stratified_logrank_core(const NumericVector&, const IntegerVector&,
                                       const IntegerVector&, const IntegerVector&);
NumericVector stratified_weighted_logrank_core(const NumericVector&,
                                               const IntegerVector&,
                                               const IntegerVector&,
                                               const IntegerVector&, int, double,
                                               double, double);
NumericVector pihe_core(const NumericVector&, const IntegerVector&,
                        const IntegerVector&);
NumericVector rmst_core(const NumericVector&, const NumericVector&, double);
NumericVector km_core(const NumericVector&, const NumericVector&, double);
List combo_logrank_core(const NumericVector&, const IntegerVector&,
                        const IntegerVector&, const NumericVector&,
                        const NumericVector&);
NumericVector ahsw_core(const NumericVector&, const IntegerVector&, double);

// Fused analysis loop. Walks simulations, looks, and populations in C++,
// applies the administrative cut inline, calls the requested statistic cores
// on each cell, and returns the raw core outputs. All transformation of those
// raw outputs into z-scores, confidence intervals, and p-values is done in the
// R wrapper, where it mirrors the existing wrapper formulas. The max-combo
// multivariate-normal p-value is also left to the R wrapper.
//
// The data must be sorted by simulation, with sim_ptr giving 0-based row
// offsets so that simulation si occupies rows [sim_ptr[si], sim_ptr[si + 1]).
// Output cells are ordered (sim, look, population) to match analysis_fast.
//
// [[Rcpp::export]]
List analysis_loop_core(
    const IntegerVector& sim_ptr,
    const NumericVector& accrual,
    const NumericVector& tte,
    const IntegerVector& event,
    const IntegerVector& j,
    int look_type,                  // 0 = event-driven, 1 = calendar time
    const NumericVector& look_values,
    const IntegerVector& pop_col,   // -1 = overall, else 0-based column in sub_mat
    const IntegerVector& pop_level,
    const IntegerMatrix& sub_mat,   // N x n_subcols (subgroup labels)
    const IntegerVector& strata,    // length N or length 0
    bool use_strata,
    bool do_logrank, bool do_coxph, bool do_rmst, bool do_km,
    bool do_maxcombo, bool do_ahsw,
    int weight_scheme,              // -1 plain, 0 FH, 1 mwlrt, 2 gehan, 3 tarone-ware
    double rho, double gamma, double t_star,
    const NumericVector& mc_rho, const NumericVector& mc_gamma,
    double tau, double t_eval
) {
  const int nsim    = sim_ptr.size() - 1;
  const int n_looks = look_values.size();
  const int n_pop   = pop_col.size();
  const int ncells  = nsim * n_looks * n_pop;
  const int nw      = mc_rho.size();

  NumericVector cutoff_out(ncells, NA_REAL);
  IntegerVector reached_out(ncells);
  IntegerVector nen_out(ncells);
  IntegerVector nev_out(ncells);

  NumericMatrix lr_mat, cox_mat, rmst_mat, km_mat, ahsw_mat, mc_Uout, mc_Vout;
  if (do_logrank) { lr_mat   = NumericMatrix(ncells, 2); std::fill(lr_mat.begin(),   lr_mat.end(),   NA_REAL); }
  if (do_coxph)   { cox_mat  = NumericMatrix(ncells, 4); std::fill(cox_mat.begin(),  cox_mat.end(),  NA_REAL); }
  if (do_rmst)    { rmst_mat = NumericMatrix(ncells, 4); std::fill(rmst_mat.begin(), rmst_mat.end(), NA_REAL); }
  if (do_km)      { km_mat   = NumericMatrix(ncells, 2); std::fill(km_mat.begin(),   km_mat.end(),   NA_REAL); }
  if (do_ahsw)    { ahsw_mat = NumericMatrix(ncells, 8); std::fill(ahsw_mat.begin(), ahsw_mat.end(), NA_REAL); }
  if (do_maxcombo) {
    mc_Uout = NumericMatrix(ncells, nw);      std::fill(mc_Uout.begin(), mc_Uout.end(), NA_REAL);
    mc_Vout = NumericMatrix(ncells, nw * nw); std::fill(mc_Vout.begin(), mc_Vout.end(), NA_REAL);
  }

  // Reusable per-simulation and per-cut buffers.
  std::vector<double> cal_ev;
  std::vector<int>    orig_cut;
  std::vector<double> t_cut;
  std::vector<int>    e_cut;
  std::vector<int>    j_cut;

  int pos = 0;
  for (int si = 0; si < nsim; ++si) {
    const int g0 = sim_ptr[si];
    const int g1 = sim_ptr[si + 1];

    // Calendar event times for event-driven cutoffs (whole simulation).
    if (look_type == 0) {
      cal_ev.clear();
      for (int g = g0; g < g1; ++g) {
        if (event[g] == 1) cal_ev.push_back(accrual[g] + tte[g]);
      }
    }

    for (int l = 0; l < n_looks; ++l) {
      double cut_cutoff;
      bool   reached;
      if (look_type == 0) {
        const int target = (int) look_values[l];
        const int n_ev_sim = (int) cal_ev.size();
        if (target < 1 || target > n_ev_sim) {
          cut_cutoff = R_PosInf;
          reached    = false;
        } else {
          std::vector<double> tmp(cal_ev);
          std::nth_element(tmp.begin(), tmp.begin() + (target - 1), tmp.end());
          cut_cutoff = tmp[target - 1];
          reached    = true;
        }
      } else {
        cut_cutoff = look_values[l];
        reached    = true;
      }

      // Administrative cut for this (sim, look): enrolled subjects, censored
      // at the cutoff, sorted by observed time. orig_cut holds the within-sim
      // 0-based index so subgroup and stratum labels can be realigned.
      orig_cut.clear(); t_cut.clear(); e_cut.clear(); j_cut.clear();
      for (int g = g0; g < g1; ++g) {
        const double a = accrual[g];
        if (a > cut_cutoff) continue;
        const double full   = tte[g];
        const bool   before = (a + full <= cut_cutoff);
        const double t      = before ? full : (cut_cutoff - a);
        const int    e      = before ? event[g] : 0;
        orig_cut.push_back(g - g0);
        t_cut.push_back(t);
        e_cut.push_back(e);
        j_cut.push_back(j[g]);
      }
      const int m = (int) t_cut.size();

      // Sort the cut by observed time.
      std::vector<int> ord(m);
      for (int k = 0; k < m; ++k) ord[k] = k;
      std::sort(ord.begin(), ord.end(),
                [&t_cut](int p, int q) { return t_cut[p] < t_cut[q]; });

      const double cutoff_report = reached ? cut_cutoff : NA_REAL;

      for (int pi = 0; pi < n_pop; ++pi) {
        // Select the rows of this population, preserving the time order.
        const int col = pop_col[pi];
        const int lev = pop_level[pi];

        std::vector<int> sel;
        sel.reserve(m);
        for (int k = 0; k < m; ++k) {
          const int kk = ord[k];
          if (col < 0) {
            sel.push_back(kk);
          } else if (sub_mat(g0 + orig_cut[kk], col) == lev) {
            sel.push_back(kk);
          }
        }
        const int sz = (int) sel.size();

        // Pooled population vectors (time-sorted).
        NumericVector t_sel(sz);
        IntegerVector ei_sel(sz);
        IntegerVector j_sel(sz);
        int n0 = 0, n1 = 0, n_ev = 0;
        for (int k = 0; k < sz; ++k) {
          const int kk = sel[k];
          t_sel[k]  = t_cut[kk];
          ei_sel[k] = e_cut[kk];
          j_sel[k]  = j_cut[kk];
          if (j_cut[kk] == 1) ++n1; else ++n0;
          if (e_cut[kk] == 1) ++n_ev;
        }
        const bool both = (n0 > 0 && n1 > 0);

        cutoff_out[pos]  = cutoff_report;
        reached_out[pos] = reached ? 1 : 0;
        nen_out[pos]     = sz;
        nev_out[pos]     = n_ev;

        // ---- log-rank (plain / weighted, optionally stratified) ----------
        if (do_logrank && n_ev > 0 && both) {
          double num = NA_REAL, var = NA_REAL;
          if (!use_strata) {
            if (weight_scheme < 0) {
              NumericVector r = logrank_core(t_sel, ei_sel, j_sel);
              num = r[0] - r[1];
              var = r[2];
            } else {
              NumericVector r = weighted_logrank_core(t_sel, ei_sel, j_sel,
                                                      weight_scheme, rho, gamma,
                                                      t_star);
              num = r[1];
              var = r[2];
            }
          } else {
            // Reorder the population subset by (stratum, time).
            IntegerVector st_sel(sz);
            for (int k = 0; k < sz; ++k) {
              st_sel[k] = strata[g0 + orig_cut[sel[k]]];
            }
            std::vector<int> so(sz);
            for (int k = 0; k < sz; ++k) so[k] = k;
            std::sort(so.begin(), so.end(), [&](int p, int q) {
              if (st_sel[p] != st_sel[q]) return st_sel[p] < st_sel[q];
              return t_sel[p] < t_sel[q];
            });
            NumericVector ts(sz); IntegerVector es(sz), js(sz), ss(sz);
            for (int k = 0; k < sz; ++k) {
              ts[k] = t_sel[so[k]]; es[k] = ei_sel[so[k]];
              js[k] = j_sel[so[k]]; ss[k] = st_sel[so[k]];
            }
            if (weight_scheme < 0) {
              NumericVector r = stratified_logrank_core(ts, es, js, ss);
              num = r[0] - r[1];
              var = r[2];
            } else {
              NumericVector r = stratified_weighted_logrank_core(
                ts, es, js, ss, weight_scheme, rho, gamma, t_star);
              num = r[1];
              var = r[2];
            }
          }
          lr_mat(pos, 0) = num;
          lr_mat(pos, 1) = var;
        }

        // ---- Cox (PiHE) ---------------------------------------------------
        if (do_coxph && n_ev > 0 && both) {
          NumericVector r = pihe_core(t_sel, ei_sel, j_sel);
          cox_mat(pos, 0) = r[0]; cox_mat(pos, 1) = r[1];
          cox_mat(pos, 2) = r[2]; cox_mat(pos, 3) = r[3];
        }

        // Group splits reused by rmst / km / ahsw.
        const bool need_split = (do_rmst && both) || do_km ||
                                (do_ahsw && n_ev > 0 && both);
        NumericVector t0d, t1d, e0d, e1d;
        IntegerVector e0i, e1i;
        if (need_split) {
          t0d = NumericVector(n0); e0d = NumericVector(n0); e0i = IntegerVector(n0);
          t1d = NumericVector(n1); e1d = NumericVector(n1); e1i = IntegerVector(n1);
          int a0 = 0, a1 = 0;
          for (int k = 0; k < sz; ++k) {
            if (j_sel[k] == 0) {
              t0d[a0] = t_sel[k]; e0d[a0] = (double) ei_sel[k]; e0i[a0] = ei_sel[k]; ++a0;
            } else {
              t1d[a1] = t_sel[k]; e1d[a1] = (double) ei_sel[k]; e1i[a1] = ei_sel[k]; ++a1;
            }
          }
        }

        // ---- RMST ---------------------------------------------------------
        if (do_rmst && both) {
          NumericVector r0 = rmst_core(t0d, e0d, tau);
          NumericVector r1 = rmst_core(t1d, e1d, tau);
          rmst_mat(pos, 0) = r0[0]; rmst_mat(pos, 1) = r0[1];
          rmst_mat(pos, 2) = r1[0]; rmst_mat(pos, 3) = r1[1];
        }

        // ---- Kaplan-Meier at t_eval (per available group) -----------------
        if (do_km) {
          if (n0 > 0) { NumericVector k0 = km_core(t0d, e0d, t_eval); km_mat(pos, 0) = k0[0]; }
          if (n1 > 0) { NumericVector k1 = km_core(t1d, e1d, t_eval); km_mat(pos, 1) = k1[0]; }
        }

        // ---- max-combo (component numerators + covariance) ----------------
        if (do_maxcombo && n_ev > 0 && both) {
          List r = combo_logrank_core(t_sel, ei_sel, j_sel, mc_rho, mc_gamma);
          NumericVector U = r["U"];
          NumericMatrix V = r["V"];
          for (int a = 0; a < nw; ++a) mc_Uout(pos, a) = U[a];
          for (int a = 0; a < nw; ++a)
            for (int b = 0; b < nw; ++b)
              mc_Vout(pos, a * nw + b) = V(a, b);
        }

        // ---- AHSW ---------------------------------------------------------
        if (do_ahsw && n_ev > 0 && both) {
          NumericVector a0c = ahsw_core(t0d, e0i, tau);
          NumericVector a1c = ahsw_core(t1d, e1i, tau);
          ahsw_mat(pos, 0) = a0c[2]; ahsw_mat(pos, 1) = a0c[3];
          ahsw_mat(pos, 2) = a0c[4]; ahsw_mat(pos, 3) = (double) n0;
          ahsw_mat(pos, 4) = a1c[2]; ahsw_mat(pos, 5) = a1c[3];
          ahsw_mat(pos, 6) = a1c[4]; ahsw_mat(pos, 7) = (double) n1;
        }

        ++pos;
      }
    }
  }

  return List::create(
    _["cutoff"]     = cutoff_out,
    _["reached"]    = reached_out,
    _["n_enrolled"] = nen_out,
    _["n_event"]    = nev_out,
    _["logrank"]    = do_logrank  ? (SEXP) lr_mat   : R_NilValue,
    _["coxph"]      = do_coxph    ? (SEXP) cox_mat  : R_NilValue,
    _["rmst"]       = do_rmst     ? (SEXP) rmst_mat : R_NilValue,
    _["km"]         = do_km       ? (SEXP) km_mat   : R_NilValue,
    _["ahsw"]       = do_ahsw     ? (SEXP) ahsw_mat : R_NilValue,
    _["mc_U"]       = do_maxcombo ? (SEXP) mc_Uout  : R_NilValue,
    _["mc_V"]       = do_maxcombo ? (SEXP) mc_Vout  : R_NilValue
  );
}
