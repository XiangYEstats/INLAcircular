#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "cgeneric.h"
#include "specfunc.h"

#define Malloc(n_, type_) (type_ *)malloc((n_) * sizeof(type_))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif
#ifndef M_LN_2PI
#define M_LN_2PI 1.83787706640934548356
#endif

// =========================================================================
// ANALYTICAL LAMBDA EVALUATORS (Thread-safe, no static caching)
// =========================================================================

static double get_lavm_pc_inf_lambda(double u, double alpha) {
  return -log(1.0 - alpha) / sqrt(1.0 - u);
}

static double get_lavm_pc_0_lambda(double u, double alpha) {
  if (u <= 1e-5) return -log(alpha) / sqrt(u);

  double k;
  if (u < 0.53) {
    k = 2.0 * u + u * u * u + 5.0 * pow(u, 5) / 6.0;
  } else if (u < 0.85) {
    k = -0.4 + 1.39 * u + 0.43 / (1.0 - u);
  } else {
    k = 1.0 / (2.0 * (1.0 - u));
  }

  for (int i = 0; i < 7; i++) {
    double i0 = cc_bessel_i(k, 0, 1);
    double i1 = cc_bessel_i(k, 1, 1);
    double a = i1 / i0;
    double diff = a - u;
    if (fabs(diff) < 1e-8) break;
    double deriv = 1.0 - a * a - a / k;
    k = k - diff / deriv;
  }

  double i0_exp = cc_bessel_i(k, 0, 1);
  double dis_sq = k * u - k - log(i0_exp);
  double T_val = sqrt(fmax(1e-16, dis_sq));
  return -log(alpha) / T_val;
}

// =========================================================================
// LIKELIHOODS
// =========================================================================
void lavm_inverse_tangent(double *res, double x, double eta, double th) {
  if (x > M_PI - 1e-12) x = M_PI - 1e-12;
  if (x < -M_PI + 1e-12) x = -M_PI + 1e-12;
  double y = tan(x / 2.0), u = y - eta;
  double term1 = (u*u > 0) ? -2.0 * exp(th + log(u*u) - log1p(u*u)) : 0.0;
  res[0] = term1 - M_LN_2PI - get_log_bessel_scaled(th) + (log1p(y*y) - log1p(u*u));
}

void lavm_scaled_logit(double *res, double x, double eta, double th) {
  double x_n = (x + M_PI) / (2.0 * M_PI);
  double u = c_qlogis(x_n) - eta, p = c_plogis(u);
  double sin2 = pow(sin((2.0 * M_PI * p - M_PI) / 2.0), 2.0);
  res[0] = -2.0 * exp(th + log(sin2)) - M_LN_2PI - get_log_bessel_scaled(th) + (log(p*(1.0-p)) - log(x_n*(1.0-x_n)));
}

void lavm_scaled_probit(double *res, double x, double eta, double th) {
  double x_n = (x + M_PI) / (2.0 * M_PI);
  double y_v = c_qnorm(x_n), u = y_v - eta;
  double sin2 = pow(sin((2.0 * M_PI * c_pnorm(u) - M_PI) / 2.0), 2.0);
  res[0] = -2.0 * exp(th + log(sin2)) - M_LN_2PI - get_log_bessel_scaled(th) + 0.5 * (y_v*y_v - u*u);
}

// =========================================================================
// PRIORS (Lambda is now passed locally, making it thread-safe)
// =========================================================================
double lavm_pc_prior_inf(double th, double lam) {
  double dist, log_j;
  if (th <= 11.512925) {
    double k = exp(th), i0 = cc_bessel_i(k,0,1), i1 = cc_bessel_i(k,1,1), i2 = cc_bessel_i(k,2,1);
    double r = i1/i0; dist = sqrt(fmax(1e-16, 1.0-r));
    log_j = -M_LN2 - log(dist) + log(fmax(1e-16, (i0+i2)/(2.0*i0) - r*r));
  } else {
    double ik = exp(-th), ik2 = ik*ik;
    double log_dist = -M_LN2/2.0 - th/2.0 + 0.125*ik + (7.0/64.0)*ik2;
    dist = exp(log_dist);
    log_j = -M_LN2 - log_dist + (-M_LN2 - 2.0*th + 0.5*ik + 0.625*ik2);
  }
  return log(lam) - lam * dist + log_j + th;
}

double lavm_pc_prior_0(double th, double lam) {
  double k = exp(th);

  if (k < 1e-4) {
    double dens = lam/2.0 - (lam*lam * k)/4.0 + (-9.0*lam/64.0 + lam*lam*lam/16.0) * k*k;
    if (dens <= 1e-16) dens = 1e-16;
    return log(dens) + th;
  }

  double dist, log_j_partial;
  if (th <= 11.512925) {
    double i0 = cc_bessel_i(k, 0, 1);
    double i1 = cc_bessel_i(k, 1, 1);
    double i2 = cc_bessel_i(k, 2, 1);
    double r = i1 / i0;
    double dis_sq = k * r - k - log(i0);
    dist = sqrt(fmax(1e-16, dis_sq));
    log_j_partial = log(k) + log(fmax(1e-16, (i0 + i2)/(2.0 * i0) - r*r));
  } else {
    double ik = exp(-th);
    double ik2 = ik * ik;
    double dis_sq = 0.5 * M_LN_2PI - 0.5 + 0.5 * th - 0.25 * ik - 0.1875 * ik2
    - (25.0/96.0) * ik2*ik - (65.0/128.0) * ik2*ik2 - (3219.0/2560.0) * ik2*ik2*ik;
    dist = sqrt(fmax(1e-16, dis_sq));
    log_j_partial = -M_LN2 - th + 0.5 * ik + 0.625 * ik2 + (59.0/48.0) * ik2*ik
    + (203.0/64.0) * ik2*ik2 + (12743.0/1280.0) * ik2*ik2*ik;
  }

  double log_jac = -M_LN2 - log(dist) + log_j_partial;
  return log(lam) - lam * dist + log_jac + th;
}

// =========================================================================
// CLOGLIKE INTERFACE
// =========================================================================
double *inla_cloglike_lavm(inla_cloglike_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data, int ny, double *y, int nx, double *x, double *result) {
  double *ret = NULL, init_th = 6.0, prior_u = 0.5, prior_alpha = 0.5;
  int link = 0, prior_code = 0, fixed_th = 0;

  if (data) {
    for (int i = 0; i < data->n_doubles; i++) {
      if (strcmp(data->doubles[i]->name, "lavm.link") == 0) link = (int)data->doubles[i]->doubles[0];
      else if (strcmp(data->doubles[i]->name, "lavm.prior") == 0) prior_code = (int)data->doubles[i]->doubles[0];
      else if (strcmp(data->doubles[i]->name, "lavm.u") == 0) prior_u = data->doubles[i]->doubles[0];
      else if (strcmp(data->doubles[i]->name, "lavm.alpha") == 0) prior_alpha = data->doubles[i]->doubles[0];
      else if (strcmp(data->doubles[i]->name, "lavm.initial.theta") == 0) init_th = data->doubles[i]->doubles[0];
      else if (strcmp(data->doubles[i]->name, "lavm.fixed.theta") == 0) fixed_th = (int)data->doubles[i]->doubles[0];
    }
  }

  double th = fixed_th ? init_th : 0.0;
  if (!fixed_th && theta) {
    th = theta[0];
  }
  if (th > 30.0) th = 30.0;
  if (th < -30.0) th = -30.0;

  switch (cmd) {
  case INLA_CLOGLIKE_INITIAL:
    if (fixed_th) {
      ret = Malloc(1, double);
      ret[0] = 0;
    } else {
      ret = Malloc(2, double);
      ret[0] = 1;
      ret[1] = init_th;
    }
    break;

  case INLA_CLOGLIKE_LOG_PRIOR: {
    // Determine lambda locally within the thread
    double current_lambda = (prior_code == 1) ? get_lavm_pc_0_lambda(prior_u, prior_alpha) : get_lavm_pc_inf_lambda(prior_u, prior_alpha);

    ret = Malloc(1, double);
    if (!theta) {
      ret[0] = 0.0;
    } else {
      ret[0] = (prior_code == 1) ? lavm_pc_prior_0(th, current_lambda) : lavm_pc_prior_inf(th, current_lambda);
    }
    break;
  }
  case INLA_CLOGLIKE_LOGLIKE:
    if ((!fixed_th && !theta) || y[0] != y[0]) { for(int i=0; i<nx; i++) result[i] = 0.0; break; }
    for (int i = 0; i < nx; i++) {
      double r[1];
      if (link == 0) lavm_inverse_tangent(r, y[0], x[i], th);
      else if (link == 1) lavm_scaled_logit(r, y[0], x[i], th);
      else lavm_scaled_probit(r, y[0], x[i], th);
      result[i] = r[0];
    }
    break;

  case INLA_CLOGLIKE_CDF:
    for (int i = 0; i < nx; i++) result[i] = NAN;
    break;
  case INLA_CLOGLIKE_QUIT:
    break;
  }
  return ret;
}
