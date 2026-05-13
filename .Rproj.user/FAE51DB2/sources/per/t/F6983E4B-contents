#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "specfunc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN_2PI
#define M_LN_2PI 1.83787706640934548356
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524401
#endif

// =========================================================================
// 1. LINK FUNCTIONS & INVERSE CDFS
// =========================================================================
double c_qlogis(double p) { return log(p / (1.0 - p)); }
double c_plogis(double x) { return 1.0 / (1.0 + exp(-x)); }
double c_pnorm(double x) { return 0.5 * (1.0 + erf(x * M_SQRT1_2)); }

double c_qnorm(double p) {
  double a1 = -39.69683028665376, a2 = 220.9460984245205, a3 = -275.9285104469687;
  double a4 = 138.3577518672690, a5 = -30.66479806614716, a6 = 2.506628277459239;
  double b1 = -54.47609879822406, b2 = 161.5858368580409, b3 = -155.6989798598866;
  double b4 = 66.80131188771972, b5 = -13.28068155288572;
  double c1 = -0.007784894002430293, c2 = -0.3223964580411365, c3 = -2.400758277161838;
  double c4 = -2.549732539343734, c5 = 4.374664141464968, c6 = 2.938163982698783;
  double d1 = 0.007784695709041462, d2 = 0.3224671290700398, d3 = 2.445134137142996, d4 = 3.754408661907416;
  double q, r;
  if (p < 0.02425) {
    q = sqrt(-2.0 * log(p));
    return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
      ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
  } else if (p <= 0.97575) {
    q = p - 0.5; r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
      (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
  } else {
    q = sqrt(-2.0 * log(1.0 - p));
    return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
      ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
  }
}

// =========================================================================
// 2. BESSEL FUNCTIONS (GSL CHEBYSHEV PORT)
// =========================================================================
static const double bi0_data[12] = { -.07660547252839144951, 1.92733795399380827000, .22826445869203013390, .01304891466707290428, .00043442709008164874, .00000942265768600193, .00000014340062895106, .00000000161384906966, .00000000001396650044, .00000000000009579451, .00000000000000053339, .00000000000000000245 };
static const double ai0_data[21] = { .07575994494023796, .00759138081082334, .00041531313389237, .00001070076463439, -.00000790117997921, -.00000078261435014, .00000027838499429, .00000000825247260, -.00000001204463945, .00000000155964859, .00000000022925563, -.00000000011916228, .00000000001757854, .00000000000112822, -.00000000000114684, .00000000000027155, -.00000000000002415, -.00000000000000608, .00000000000000314, -.00000000000000071, .00000000000000007 };
static const double ai02_data[22] = { .05449041101410882, .00336911647825569, .00006889758346918, .00000289137052082, .00000020489185893, .00000002266668991, .00000000339623203, .00000000049406022, .00000000001188914, -.00000000003149915, -.00000000001321580, -.00000000000179419, .00000000000071801, .00000000000038529, .00000000000001539, -.00000000000004151, -.00000000000000954, .00000000000000382, .00000000000000176, -.00000000000000034, -.00000000000000027, .00000000000000003 };
static const double bi1_data[11] = { -0.001971713261099859, 0.407348876675464810, 0.034838994299959456, 0.001545394556300123, 0.000041888521098377, 0.000000764902676483, 0.000000010042493924, 0.000000000099322077, 0.000000000000766380, 0.000000000000004741, 0.000000000000000024 };
static const double ai1_data[21] = { -0.02846744181881479, -0.01922953231443221, -0.00061151858579437, -0.00002069971253350, 0.00000858561914581, 0.00000104949824671, -0.00000029183389184, -0.00000001559378146, 0.00000001318012367, -0.00000000144842341, -0.00000000029085122, 0.00000000012663889, -0.00000000001664947, -0.00000000000166665, 0.00000000000124260, -0.00000000000027315, 0.00000000000002023, 0.00000000000000730, -0.00000000000000333, 0.00000000000000071, -0.00000000000000006 };
static const double ai12_data[22] = { 0.02857623501828014, -0.00976109749136147, -0.00011058893876263, -0.00000388256480887, -0.00000025122362377, -0.00000002631468847, -0.00000000383538039, -0.00000000055897433, -0.00000000001897495, 0.00000000003252602, 0.00000000001412580, 0.00000000000203564, -0.00000000000071985, -0.00000000000040836, -0.00000000000002101, 0.00000000000004273, 0.00000000000001041, -0.00000000000000382, -0.00000000000000186, 0.00000000000000033, 0.00000000000000028, -0.00000000000000003 };

static double cc_cheb_eval(const double *c, const int order, const double x) {
  double d = 0.0, dd = 0.0, y2 = 2.0 * x;
  for (int j = order; j >= 1; j--) {
    double temp = d;
    d = y2 * d - dd + c[j];
    dd = temp;
  }
  return x * d - dd + 0.5 * c[0];
}

static double cc_bessel_I0_scaled(double x) {
  double y = fabs(x);
  if (y < 2.0 * 1.490116e-08) return 1.0 - y;
  if (y <= 3.0) return exp(-y) * (2.75 + cc_cheb_eval(bi0_data, 11, y * y / 4.5 - 1.0));
  double sy = sqrt(y);
  if (y <= 8.0) return (0.375 + cc_cheb_eval(ai0_data, 20, (48.0 / y - 11.0) / 5.0)) / sy;
  return (0.375 + cc_cheb_eval(ai02_data, 21, 16.0 / y - 1.0)) / sy;
}

static double cc_bessel_I0_unscaled(double x) {
  double y = fabs(x);
  if (y < 2.0 * 1.490116e-08) return 1.0;
  if (y <= 3.0) return 2.75 + cc_cheb_eval(bi0_data, 11, y * y / 4.5 - 1.0);
  return exp(y) * cc_bessel_I0_scaled(x);
}

static double cc_bessel_I1_scaled(double x) {
  double y = fabs(x), s = (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
  if (y == 0.0) return 0.0;
  if (y < 2.828427 * 1.490116e-08) return 0.5 * y;
  if (y <= 3.0) return y * exp(-y) * (0.875 + cc_cheb_eval(bi1_data, 10, y * y / 4.5 - 1.0));
  double sy = sqrt(y);
  if (y <= 8.0) return s * (0.375 + cc_cheb_eval(ai1_data, 20, (48.0 / y - 11.0) / 5.0)) / sy;
  return s * (0.375 + cc_cheb_eval(ai12_data, 21, 16.0 / y - 1.0)) / sy;
}

static double cc_bessel_I1_unscaled(double x) {
  double y = fabs(x);
  if (y == 0.0) return 0.0;
  if (y < 2.828427 * 1.490116e-08) return 0.5 * y;
  if (y <= 3.0) return y * (0.875 + cc_cheb_eval(bi1_data, 10, y * y / 4.5 - 1.0));
  return exp(y) * cc_bessel_I1_scaled(y);
}

static double cc_bessel_I_CF1(double nu, double x) {
  double tk = 1.0, sum = 1.0, rhok = 0.0;
  for (int k = 1; k < 20000; k++) {
    double ak = 0.25 * (x / (nu + k)) * x / (nu + k + 1.0);
    rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
    tk *= rhok;
    sum += tk;
    if (fabs(tk / sum) < 2e-16) break;
  }
  return (x / (2.0 * (nu + 1.0))) * sum;
}

double cc_bessel_i(double x, double nu, int scaled) {
  double ax = fabs(x);
  int n = (int)(nu + 0.5);
  double sgn = (x < 0.0 && (n % 2 != 0)) ? -1.0 : 1.0;

  double res_scaled;
  if (n == 0) res_scaled = cc_bessel_I0_scaled(ax);
  else if (n == 1) res_scaled = cc_bessel_I1_scaled(ax);
  else if (ax == 0.0) return 0.0;
  else {
    double rat = cc_bessel_I_CF1((double)n, ax);
    double Ikp1 = rat * 1e-100, Ik = 1e-100, Ikm1;
    for (int k = n; k >= 1; k--) {
      Ikm1 = Ikp1 + (2.0 * k / ax) * Ik;
      Ikp1 = Ik;
      Ik = Ikm1;
    }
    res_scaled = cc_bessel_I0_scaled(ax) * (1e-100 / Ik);
  }

  double final_res = sgn * res_scaled;
  if (!scaled) final_res *= exp(ax);
  return final_res;
}

double get_log_bessel_scaled(double th) {
  if (th > 11.512925) {
    return -0.5 * (M_LN_2PI + th) + 0.125 * exp(-th);
  } else {
    return log(cc_bessel_i(exp(th), 0.0, 1));
  }
}

// =========================================================================
// 3. R C-API INTERFACE (.Call wrapper)
// =========================================================================
SEXP C_bessel_i(SEXP x, SEXP nu, SEXP expon_scaled) {
  int n = length(x);
  double *dx = REAL(x);
  double dnu = asReal(nu);
  int scaled = asLogical(expon_scaled);

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *dout = REAL(out);

  for (int i = 0; i < n; i++) {
    if (ISNAN(dx[i])) dout[i] = NA_REAL;
    else dout[i] = cc_bessel_i(dx[i], dnu, scaled);
  }

  UNPROTECT(1);
  return out;
}
