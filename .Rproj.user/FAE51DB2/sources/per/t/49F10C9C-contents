#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "specfunc.h"

// --- 1. Von Mises Density (dvm) ---
SEXP C_dvm(SEXP x, SEXP mu, SEXP kappa, SEXP log_prob) {
  int n = length(x);
  double *dx = REAL(x), *dmu = REAL(mu), *dk = REAL(kappa);
  int give_log = asLogical(log_prob);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  for (int i = 0; i < n; i++) {
    double k = dk[i];
    double d_angle = (dx[i] - dmu[i]) / 2.0;
    double sin_val = sin(d_angle);
    double log_bessel = get_log_bessel_scaled(k);
    double log_dens = -2.0 * k * sin_val * sin_val - M_LN_2PI - log_bessel;

    if (give_log) dout[i] = log_dens;
    else dout[i] = exp(log_dens);
  }

  UNPROTECT(1);
  return out;
}

// --- 2. Link-Adjusted Density (dlavm) ---
SEXP C_dlavm(SEXP x, SEXP eta, SEXP kappa, SEXP log_prob) {
  int n = length(x);
  double *dx = REAL(x), *deta = REAL(eta), *dk = REAL(kappa);
  int give_log = asLogical(log_prob);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  for (int i = 0; i < n; i++) {
    double y = tan(dx[i] / 2.0);
    double u = y - deta[i];
    double z = 2.0 * atan(u);

    double log_J = log1p(y * y) - log1p(u * u);

    double k = dk[i];
    double log_bessel = get_log_bessel_scaled(k);

    double d_angle = z / 2.0;
    double sin_val = sin(d_angle);
    double log_vm = -2.0 * k * sin_val * sin_val - M_LN_2PI - log_bessel;

    double total_log = log_vm + log_J;

    if (give_log) dout[i] = total_log;
    else dout[i] = exp(total_log);
  }

  UNPROTECT(1);
  return out;
}

// --- Helper: Fast Linear Interpolation via Binary Search ---
double approx_linear(double x, const double *grid_x, const double *grid_y, int n_grid) {
  if (x <= grid_x[0]) return grid_y[0];
  if (x >= grid_x[n_grid - 1]) return grid_y[n_grid - 1];

  int l = 0, r = n_grid - 1;
  while (l <= r) {
    int m = l + (r - l) / 2;
    if (grid_x[m] == x) return grid_y[m];
    if (grid_x[m] < x) l = m + 1;
    else r = m - 1;
  }

  // now r < l, and x is between grid_x[r] and grid_x[l]
  double x0 = grid_x[r], x1 = grid_x[l];
  double y0 = grid_y[r], y1 = grid_y[l];
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}


// --- 3. Cumulative Distribution (pvm) ---
SEXP C_pvm(SEXP q, SEXP mu, SEXP strategy, SEXP log_prob, SEXP grid_x, SEXP grid_y) {
  int n = length(q);
  int n_grid = length(grid_x);
  double *dq = REAL(q), *dmu = REAL(mu);
  double *gx = REAL(grid_x), *gy = REAL(grid_y);
  int strat = asInteger(strategy);
  int give_log = asLogical(log_prob);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  for (int i = 0; i < n; i++) {
    double z;
    if (strat == 0) { // circular strategy
      z = dq[i] - dmu[i] + M_PI;
      z = z - 2.0 * M_PI * floor(z / (2.0 * M_PI));
      z = z - M_PI;
    } else { // linear strategy
      z = dq[i] - dmu[i];
    }

    double logit_val = approx_linear(z, gx, gy, n_grid);
    dout[i] = plogis(logit_val, 0.0, 1.0, 1, give_log);
  }
  UNPROTECT(1);
  return out;
}

// --- 4. Quantile Function (qvm) ---
SEXP C_qvm(SEXP p, SEXP mu, SEXP grid_x, SEXP grid_y) {
  int n = length(p);
  int n_grid = length(grid_x);
  double *dp = REAL(p), *dmu = REAL(mu);
  double *gx = REAL(grid_x), *gy = REAL(grid_y);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  double min_logit = gy[0];
  double max_logit = gy[n_grid - 1];

  for (int i = 0; i < n; i++) {
    double p_val = dp[i];
    if (ISNAN(p_val)) { dout[i] = NA_REAL; continue; }
    if (p_val <= 0.0) { dout[i] = dmu[i] - M_PI; continue; }
    if (p_val >= 1.0) { dout[i] = dmu[i] + M_PI; continue; }

    double logit_p = qlogis(p_val, 0.0, 1.0, 1, 0);
    if (logit_p < min_logit) logit_p = min_logit;
    if (logit_p > max_logit) logit_p = max_logit;

    double z = approx_linear(logit_p, gy, gx, n_grid);
    dout[i] = z + dmu[i];
  }
  UNPROTECT(1);
  return out;
}

// --- 5. Random Generation (rvm) ---
SEXP C_rvm(SEXP n_sxp, SEXP mu, SEXP grid_x, SEXP grid_y) {
  int n = asInteger(n_sxp);
  int n_grid = length(grid_x);
  double *dmu = REAL(mu);
  double *gx = REAL(grid_x), *gy = REAL(grid_y);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  double min_logit = gy[0];
  double max_logit = gy[n_grid - 1];

  GetRNGstate();
  for (int i = 0; i < n; i++) {
    double p_val = unif_rand();

    double logit_p = qlogis(p_val, 0.0, 1.0, 1, 0);
    if (logit_p < min_logit) logit_p = min_logit;
    if (logit_p > max_logit) logit_p = max_logit;

    double z = approx_linear(logit_p, gy, gx, n_grid);
    dout[i] = z + dmu[i];
  }
  PutRNGstate();
  UNPROTECT(1);
  return out;
}

// --- 6. PLAVM ---
SEXP C_plavm(SEXP q, SEXP eta, SEXP log_prob, SEXP grid_x, SEXP grid_y) {
  int n = length(q);
  int n_grid = length(grid_x);
  double *dq = REAL(q), *deta = REAL(eta);
  double *gx = REAL(grid_x), *gy = REAL(grid_y);
  int give_log = asLogical(log_prob);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  for (int i = 0; i < n; i++) {
    double u = tan(dq[i] / 2.0) - deta[i];
    double z = 2.0 * atan(u);

    double logit_val = approx_linear(z, gx, gy, n_grid);
    dout[i] = plogis(logit_val, 0.0, 1.0, 1, give_log);
  }
  UNPROTECT(1);
  return out;
}

// --- 7. QLAVM ---
SEXP C_qlavm(SEXP p, SEXP eta, SEXP grid_x, SEXP grid_y) {
  int n = length(p);
  int n_grid = length(grid_x);
  double *dp = REAL(p), *deta = REAL(eta);
  double *gx = REAL(grid_x), *gy = REAL(grid_y);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  double min_logit = gy[0];
  double max_logit = gy[n_grid - 1];

  for (int i = 0; i < n; i++) {
    double p_val = dp[i];
    if (ISNAN(p_val)) { dout[i] = NA_REAL; continue; }
    if (p_val <= 0.0) p_val = 1e-15;
    if (p_val >= 1.0) p_val = 1.0 - 1e-15;

    double logit_p = qlogis(p_val, 0.0, 1.0, 1, 0);
    if (logit_p < min_logit) logit_p = min_logit;
    if (logit_p > max_logit) logit_p = max_logit;

    double z = approx_linear(logit_p, gy, gx, n_grid);
    dout[i] = 2.0 * atan(tan(z / 2.0) + deta[i]);
  }
  UNPROTECT(1);
  return out;
}

// --- 8. RLAVM ---
SEXP C_rlavm(SEXP n_sxp, SEXP eta, SEXP grid_x, SEXP grid_y) {
  int n = asInteger(n_sxp);
  int n_grid = length(grid_x);
  double *deta = REAL(eta);
  double *gx = REAL(grid_x), *gy = REAL(grid_y);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  double min_logit = gy[0];
  double max_logit = gy[n_grid - 1];

  GetRNGstate();
  for (int i = 0; i < n; i++) {
    double p_val = unif_rand();

    double logit_p = qlogis(p_val, 0.0, 1.0, 1, 0);
    if (logit_p < min_logit) logit_p = min_logit;
    if (logit_p > max_logit) logit_p = max_logit;

    double z = approx_linear(logit_p, gy, gx, n_grid);
    dout[i] = 2.0 * atan(tan(z / 2.0) + deta[i]);
  }
  PutRNGstate();
  UNPROTECT(1);
  return out;
}
