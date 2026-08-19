#ifndef INLACIRC_COMMON_H
#       define INLACIRC_COMMON_H
#       include <stdlib.h>

/*
 * Numerical functions shared by the R API and the INLA cloglike module.
 *
 * This header deliberately depends only on the ISO C mathematical library.
 * In particular, it must never include R.h, Rinternals.h, Rmath.h, or any
 * other R header: src-cloglike is compiled directly into the INLA code base.
 * All definitions have internal linkage and all names use the INLAcirc_
 * prefix to avoid collisions in a large native executable.
 */

#       include <float.h>
#       include <math.h>

#       define INLACIRC_LOG_2PI 1.837877066409345483560659472811235279
#       define INLACIRC_LOG_1E5 11.512925464970228420089957273421821038

static inline double INLAcirc_qlogis(double p)
{
	return log(p / (1.0 - p));
}

static inline double INLAcirc_plogis(double x)
{
	return 1.0 / (1.0 + exp(-x));
}

static inline double INLAcirc_pnorm(double x)
{
	return 0.5 * (1.0 + erf(x * M_SQRT1_2));
}

/* Acklam's rational approximation to the standard-normal quantile. */
static inline double INLAcirc_qnorm(double p)
{
	const double a1 = -39.69683028665376;
	const double a2 = 220.9460984245205;
	const double a3 = -275.9285104469687;
	const double a4 = 138.3577518672690;
	const double a5 = -30.66479806614716;
	const double a6 = 2.506628277459239;
	const double b1 = -54.47609879822406;
	const double b2 = 161.5858368580409;
	const double b3 = -155.6989798598866;
	const double b4 = 66.80131188771972;
	const double b5 = -13.28068155288572;
	const double c1 = -0.007784894002430293;
	const double c2 = -0.3223964580411365;
	const double c3 = -2.400758277161838;
	const double c4 = -2.549732539343734;
	const double c5 = 4.374664141464968;
	const double c6 = 2.938163982698783;
	const double d1 = 0.007784695709041462;
	const double d2 = 0.3224671290700398;
	const double d3 = 2.445134137142996;
	const double d4 = 3.754408661907416;
	double q;
	double r;

	if (p <= 0.0) {
		return -INFINITY;
	}
	if (p >= 1.0) {
		return INFINITY;
	}
	if (p < 0.02425) {
		q = sqrt(-2.0 * log(p));
		return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
	}
	if (p <= 0.97575) {
		q = p - 0.5;
		r = q * q;
		return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
	}

	q = sqrt(-2.0 * log(1.0 - p));
	return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
}

/* Chebyshev coefficients from the GSL implementation of modified Bessel I. */
static const double INLAcirc_bi0_data[12] = {
	-0.07660547252839144951, 1.92733795399380827000,
	0.22826445869203013390, 0.01304891466707290428,
	0.00043442709008164874, 0.00000942265768600193,
	0.00000014340062895106, 0.00000000161384906966,
	0.00000000001396650044, 0.00000000000009579451,
	0.00000000000000053339, 0.00000000000000000245
};

static const double INLAcirc_ai0_data[21] = {
	0.07575994494023796, 0.00759138081082334, 0.00041531313389237,
	0.00001070076463439, -0.00000790117997921, -0.00000078261435014,
	0.00000027838499429, 0.00000000825247260, -0.00000001204463945,
	0.00000000155964859, 0.00000000022925563, -0.00000000011916228,
	0.00000000001757854, 0.00000000000112822, -0.00000000000114684,
	0.00000000000027155, -0.00000000000002415, -0.00000000000000608,
	0.00000000000000314, -0.00000000000000071, 0.00000000000000007
};

static const double INLAcirc_ai02_data[22] = {
	0.05449041101410882, 0.00336911647825569, 0.00006889758346918,
	0.00000289137052082, 0.00000020489185893, 0.00000002266668991,
	0.00000000339623203, 0.00000000049406022, 0.00000000001188914,
	-0.00000000003149915, -0.00000000001321580, -0.00000000000179419,
	0.00000000000071801, 0.00000000000038529, 0.00000000000001539,
	-0.00000000000004151, -0.00000000000000954, 0.00000000000000382,
	0.00000000000000176, -0.00000000000000034, -0.00000000000000027,
	0.00000000000000003
};

static const double INLAcirc_bi1_data[11] = {
	-0.001971713261099859, 0.407348876675464810, 0.034838994299959456,
	0.001545394556300123, 0.000041888521098377, 0.000000764902676483,
	0.000000010042493924, 0.000000000099322077, 0.000000000000766380,
	0.000000000000004741, 0.000000000000000024
};

static const double INLAcirc_ai1_data[21] = {
	-0.02846744181881479, -0.01922953231443221, -0.00061151858579437,
	-0.00002069971253350, 0.00000858561914581, 0.00000104949824671,
	-0.00000029183389184, -0.00000001559378146, 0.00000001318012367,
	-0.00000000144842341, -0.00000000029085122, 0.00000000012663889,
	-0.00000000001664947, -0.00000000000166665, 0.00000000000124260,
	-0.00000000000027315, 0.00000000000002023, 0.00000000000000730,
	-0.00000000000000333, 0.00000000000000071, -0.00000000000000006
};

static const double INLAcirc_ai12_data[22] = {
	0.02857623501828014, -0.00976109749136147, -0.00011058893876263,
	-0.00000388256480887, -0.00000025122362377, -0.00000002631468847,
	-0.00000000383538039, -0.00000000055897433, -0.00000000001897495,
	0.00000000003252602, 0.00000000001412580, 0.00000000000203564,
	-0.00000000000071985, -0.00000000000040836, -0.00000000000002101,
	0.00000000000004273, 0.00000000000001041, -0.00000000000000382,
	-0.00000000000000186, 0.00000000000000033, 0.00000000000000028,
	-0.00000000000000003
};

static inline double INLAcirc_chebyshev_eval(const double *coefficients, int order, double x)
{
	double d = 0.0;
	double dd = 0.0;
	const double two_x = 2.0 * x;

	for (int j = order; j >= 1; --j) {
		const double previous_d = d;
		d = two_x * d - dd + coefficients[j];
		dd = previous_d;
	}
	return x * d - dd + 0.5 * coefficients[0];
}

static inline double INLAcirc_bessel_i0_scaled(double x)
{
	const double y = fabs(x);

	if (y < 2.0 * 1.490116e-08) {
		return 1.0 - y;
	}
	if (y <= 3.0) {
		return exp(-y) * (2.75 + INLAcirc_chebyshev_eval(INLAcirc_bi0_data, 11, y * y / 4.5 - 1.0));
	}
	if (y <= 8.0) {
		return (0.375 + INLAcirc_chebyshev_eval(INLAcirc_ai0_data, 20, (48.0 / y - 11.0) / 5.0)) / sqrt(y);
	}
	return (0.375 + INLAcirc_chebyshev_eval(INLAcirc_ai02_data, 21, 16.0 / y - 1.0)) / sqrt(y);
}

static inline double INLAcirc_bessel_i1_scaled(double x)
{
	const double y = fabs(x);
	const double sign = (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);

	if (y == 0.0) {
		return 0.0;
	}
	if (y < 2.828427 * 1.490116e-08) {
		return 0.5 * x * exp(-y);
	}
	if (y <= 3.0) {
		return x * exp(-y) * (0.875 + INLAcirc_chebyshev_eval(INLAcirc_bi1_data, 10, y * y / 4.5 - 1.0));
	}
	if (y <= 8.0) {
		return sign * (0.375 + INLAcirc_chebyshev_eval(INLAcirc_ai1_data, 20, (48.0 / y - 11.0) / 5.0)) / sqrt(y);
	}
	return sign * (0.375 + INLAcirc_chebyshev_eval(INLAcirc_ai12_data, 21, 16.0 / y - 1.0)) / sqrt(y);
}

static inline double INLAcirc_bessel_i_cf1(double nu, double x)
{
	double term = 1.0;
	double sum = 1.0;
	double rho = 0.0;

	for (int k = 1; k < 20000; ++k) {
		const double a = 0.25 * (x / (nu + k)) * x / (nu + k + 1.0);
		rho = -a * (1.0 + rho) / (1.0 + a * (1.0 + rho));
		term *= rho;
		sum += term;
		if (fabs(term / sum) < 2e-16) {
			break;
		}
	}
	return (x / (2.0 * (nu + 1.0))) * sum;
}

static inline double INLAcirc_bessel_i(double x, double nu, int scaled)
{
	const double absolute_x = fabs(x);
	const int order = (int) (nu + 0.5);
	const double sign = (x < 0.0 && order % 2 != 0) ? -1.0 : 1.0;
	double scaled_result;

	if (order == 0) {
		scaled_result = INLAcirc_bessel_i0_scaled(absolute_x);
	} else if (order == 1) {
		scaled_result = INLAcirc_bessel_i1_scaled(absolute_x);
	} else if (absolute_x == 0.0) {
		return 0.0;
	} else {
		const double ratio = INLAcirc_bessel_i_cf1((double) order, absolute_x);
		double i_k_plus_1 = ratio * 1e-100;
		double i_k = 1e-100;
		double i_k_minus_1 = 0.0;

		for (int k = order; k >= 1; --k) {
			i_k_minus_1 = i_k_plus_1 + (2.0 * k / absolute_x) * i_k;
			i_k_plus_1 = i_k;
			i_k = i_k_minus_1;
		}
		scaled_result = INLAcirc_bessel_i0_scaled(absolute_x) * (1e-100 / i_k);
	}

	scaled_result *= sign;
	if (!scaled) {
		scaled_result *= exp(absolute_x);
	}
	return scaled_result;
}

/*
 * PC prior with the circular uniform distribution (kappa = 0) as base.
 *
 *     d_0(kappa)^2 = kappa I1(kappa) / I0(kappa) - log I0(kappa).
 *
 * The small-kappa density polynomial and the large-kappa expansions are the
 * same approximations used in the original INLA reference implementation.
 */
static inline double INLAcirc_pc_vm0_small_density(double kappa, double lambda)
{
	const double lambda2 = lambda * lambda;
	const double lambda3 = lambda2 * lambda;
	const double lambda4 = lambda2 * lambda2;
	const double lambda5 = lambda4 * lambda;
	const double lambda6 = lambda3 * lambda3;
	const double coefficient2 = -9.0 * lambda / 64.0 + lambda3 / 16.0;
	const double coefficient3 = 3.0 * lambda2 / 32.0 - lambda4 / 96.0;
	const double coefficient4 = 1195.0 * lambda / 36864.0 + 3.0 * lambda3 / 512.0 + lambda5 / 768.0;
	const double coefficient5 = -79.0 * lambda2 / 6144.0 - lambda6 / 7680.0;

	return lambda / 2.0 + kappa *
	    (-lambda2 / 4.0 + kappa * (coefficient2 + kappa * (coefficient3 + kappa * (coefficient4 + kappa * coefficient5))));
}

static inline void INLAcirc_pc_vm0_geometry_from_log_kappa(double log_kappa, double *distance, double *log_abs_derivative)
{
	if (isnan(log_kappa)) {
		*distance = NAN;
		*log_abs_derivative = NAN;
		return;
	}
	if (log_kappa == -INFINITY) {
		*distance = 0.0;
		*log_abs_derivative = -M_LN2;
		return;
	}
	if (log_kappa == INFINITY) {
		*distance = INFINITY;
		*log_abs_derivative = -INFINITY;
		return;
	}

	if (log_kappa < -9.21034037197618273607) {
		const double kappa = exp(log_kappa);
		const double kappa_squared = kappa * kappa;
		const double derivative = 0.5 - (9.0 / 64.0) * kappa_squared;

		*distance = kappa * (0.5 - (3.0 / 64.0) * kappa_squared);
		*log_abs_derivative = log(derivative);
		return;
	}

	if (log_kappa <= INLACIRC_LOG_1E5) {
		const double kappa = exp(log_kappa);
		const double i0 = INLAcirc_bessel_i(kappa, 0.0, 1);
		const double i1 = INLAcirc_bessel_i(kappa, 1.0, 1);
		const double ratio = i1 / i0;
		const double ratio_derivative = fmax(0.0, 1.0 - ratio * ratio - ratio / kappa);
		const double distance_squared = fmax(0.0, kappa * ratio - kappa - log(i0));

		*distance = sqrt(distance_squared);
		if (*distance == 0.0 || ratio_derivative == 0.0) {
			*log_abs_derivative = -INFINITY;
		} else {
			*log_abs_derivative = log(kappa) + log(ratio_derivative) - M_LN2 - log(*distance);
		}
		return;
	}

	{
		const double inverse_kappa = exp(-log_kappa);
		const double inverse_kappa_squared = inverse_kappa * inverse_kappa;
		const double inverse_kappa_cubed = inverse_kappa_squared * inverse_kappa;
		const double inverse_kappa_fourth = inverse_kappa_squared * inverse_kappa_squared;
		const double inverse_kappa_fifth = inverse_kappa_fourth * inverse_kappa;
		const double distance_squared =
		    0.5 * INLACIRC_LOG_2PI - 0.5 + 0.5 * log_kappa -
		    0.25 * inverse_kappa -
		    (3.0 / 16.0) * inverse_kappa_squared -
		    (25.0 / 96.0) * inverse_kappa_cubed - (65.0 / 128.0) * inverse_kappa_fourth - (3219.0 / 2560.0) * inverse_kappa_fifth;
		const double log_derivative_partial =
		    -M_LN2 - log_kappa +
		    0.5 * inverse_kappa +
		    (5.0 / 8.0) * inverse_kappa_squared +
		    (59.0 / 48.0) * inverse_kappa_cubed + (203.0 / 64.0) * inverse_kappa_fourth + (12743.0 / 1280.0) * inverse_kappa_fifth;

		*distance = sqrt(fmax(0.0, distance_squared));
		*log_abs_derivative = -M_LN2 - log(*distance) + log_derivative_partial;
	}
}

static inline void INLAcirc_pc_vm0_geometry(double kappa, double *distance, double *log_abs_derivative)
{
	if (isnan(kappa) || kappa < 0.0) {
		*distance = NAN;
		*log_abs_derivative = NAN;
		return;
	}
	INLAcirc_pc_vm0_geometry_from_log_kappa((kappa == 0.0) ? -INFINITY : log(kappa), distance, log_abs_derivative);
}

static inline double INLAcirc_pc_vm0_log_density(double kappa, double lambda)
{
	double distance;
	double log_abs_derivative;

	if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0) {
		return NAN;
	}
	if (kappa < 0.0 || kappa == INFINITY) {
		return -INFINITY;
	}
	if (kappa < 1e-4) {
		return log(INLAcirc_pc_vm0_small_density(kappa, lambda));
	}

	INLAcirc_pc_vm0_geometry(kappa, &distance, &log_abs_derivative);
	return log(lambda) - lambda * distance + log_abs_derivative;
}

static inline double INLAcirc_pc_vm0_log_density_from_log_kappa(double log_kappa, double lambda)
{
	double distance;
	double log_abs_derivative;

	if (isnan(log_kappa) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0) {
		return NAN;
	}
	if (log_kappa == INFINITY) {
		return -INFINITY;
	}
	if (log_kappa < -9.21034037197618273607) {
		const double kappa = exp(log_kappa);
		return log(INLAcirc_pc_vm0_small_density(kappa, lambda)) + log_kappa;
	}

	INLAcirc_pc_vm0_geometry_from_log_kappa(log_kappa, &distance, &log_abs_derivative);
	return log(lambda) - lambda * distance + log_abs_derivative + log_kappa;
}

static inline double INLAcirc_pc_vm0_log_cdf(double kappa, double lambda)
{
	double distance;
	double unused_log_abs_derivative;

	if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0) {
		return NAN;
	}
	if (kappa < 0.0) {
		return -INFINITY;
	}
	if (kappa == INFINITY) {
		return 0.0;
	}

	INLAcirc_pc_vm0_geometry(kappa, &distance, &unused_log_abs_derivative);
	return log(-expm1(-lambda * distance));
}

static inline double INLAcirc_pc_vm0_quantile(double probability, double lambda)
{
	const double lower_log_kappa = log(DBL_MIN);
	const double upper_log_kappa = log(DBL_MAX);
	double target_log_distance;
	double lower = lower_log_kappa;
	double upper = upper_log_kappa;
	double distance;
	double unused_log_abs_derivative;

	if (isnan(probability) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0 || probability < 0.0 || probability > 1.0) {
		return NAN;
	}
	if (probability == 0.0) {
		return 0.0;
	}
	if (probability == 1.0) {
		return INFINITY;
	}

	target_log_distance = log(-log1p(-probability)) - log(lambda);

	INLAcirc_pc_vm0_geometry_from_log_kappa(upper, &distance, &unused_log_abs_derivative);
	if (log(distance) < target_log_distance) {
		return INFINITY;
	}

	for (int iteration = 0; iteration < 100; ++iteration) {
		const double middle = lower + 0.5 * (upper - lower);

		INLAcirc_pc_vm0_geometry_from_log_kappa(middle, &distance, &unused_log_abs_derivative);
		if (log(distance) < target_log_distance) {
			lower = middle;
		} else {
			upper = middle;
		}
	}

	return exp(lower + 0.5 * (upper - lower));
}

/*
 * PC prior with the point-mass limit (kappa -> infinity) as base.
 *
 *     d_inf(kappa) = sqrt(1 - I1(kappa) / I0(kappa)).
 *
 * The large-kappa expansions avoid cancellation in 1 - I1/I0.
 */
static inline void INLAcirc_pc_vminf_geometry_from_log_kappa(double log_kappa, double *distance, double *log_abs_derivative)
{
	if (isnan(log_kappa)) {
		*distance = NAN;
		*log_abs_derivative = NAN;
		return;
	}
	if (log_kappa == INFINITY) {
		*distance = 0.0;
		*log_abs_derivative = -INFINITY;
		return;
	}

	if (log_kappa > INLACIRC_LOG_1E5) {
		const double inverse_kappa = exp(-log_kappa);
		const double inverse_kappa_squared = inverse_kappa * inverse_kappa;
		const double inverse_kappa_cubed = inverse_kappa_squared * inverse_kappa;
		const double inverse_kappa_fourth = inverse_kappa_squared * inverse_kappa_squared;
		const double inverse_kappa_fifth = inverse_kappa_fourth * inverse_kappa;
		const double log_distance =
		    -0.5 * M_LN2 - 0.5 * log_kappa +
		    0.125 * inverse_kappa +
		    (7.0 / 64.0) * inverse_kappa_squared +
		    (1.0 / 6.0) * inverse_kappa_cubed + (715.0 / 2048.0) * inverse_kappa_fourth + (293.0 / 320.0) * inverse_kappa_fifth;
		const double log_derivative_partial =
		    -M_LN2 - 2.0 * log_kappa +
		    0.5 * inverse_kappa +
		    (5.0 / 8.0) * inverse_kappa_squared +
		    (59.0 / 48.0) * inverse_kappa_cubed + (203.0 / 64.0) * inverse_kappa_fourth + (12743.0 / 1280.0) * inverse_kappa_fifth;

		*distance = exp(log_distance);
		*log_abs_derivative = -M_LN2 - log_distance + log_derivative_partial;
		return;
	}

	{
		const double kappa = exp(log_kappa);
		const double i0 = INLAcirc_bessel_i(kappa, 0.0, 1);
		const double i1 = INLAcirc_bessel_i(kappa, 1.0, 1);
		const double ratio = i1 / i0;
		const double distance_squared = fmax(0.0, 1.0 - ratio);
		const double ratio_derivative = (kappa == 0.0)
		    ? 0.5 : fmax(0.0, 1.0 - ratio * ratio - ratio / kappa);

		*distance = sqrt(distance_squared);
		if (*distance == 0.0 || ratio_derivative == 0.0) {
			*log_abs_derivative = -INFINITY;
		} else {
			*log_abs_derivative = -M_LN2 - log(*distance) + log(ratio_derivative);
		}
	}
}

static inline void INLAcirc_pc_vminf_geometry(double kappa, double *distance, double *log_abs_derivative)
{
	if (isnan(kappa) || kappa < 0.0) {
		*distance = NAN;
		*log_abs_derivative = NAN;
		return;
	}
	INLAcirc_pc_vminf_geometry_from_log_kappa((kappa == 0.0) ? -INFINITY : log(kappa), distance, log_abs_derivative);
}

static inline double INLAcirc_pc_vminf_log_density(double kappa, double lambda)
{
	double distance;
	double log_abs_derivative;

	if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0) {
		return NAN;
	}
	if (kappa < 0.0 || kappa == INFINITY) {
		return -INFINITY;
	}

	INLAcirc_pc_vminf_geometry(kappa, &distance, &log_abs_derivative);
	return log(lambda) - lambda * distance + log_abs_derivative;
}

static inline double INLAcirc_pc_vminf_log_density_from_log_kappa(double log_kappa, double lambda)
{
	double distance;
	double log_abs_derivative;

	if (isnan(log_kappa) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0) {
		return NAN;
	}
	if (log_kappa == INFINITY) {
		return -INFINITY;
	}

	INLAcirc_pc_vminf_geometry_from_log_kappa(log_kappa, &distance, &log_abs_derivative);
	return log(lambda) - lambda * distance + log_abs_derivative + log_kappa;
}

/* Includes the boundary mass exp(-lambda) at kappa = 0. */
static inline double INLAcirc_pc_vminf_log_cdf(double kappa, double lambda)
{
	double distance;
	double unused_log_abs_derivative;

	if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0) {
		return NAN;
	}
	if (kappa < 0.0) {
		return -INFINITY;
	}
	if (kappa == INFINITY) {
		return 0.0;
	}

	INLAcirc_pc_vminf_geometry(kappa, &distance, &unused_log_abs_derivative);
	return -lambda * distance;
}

static inline double INLAcirc_pc_vminf_quantile(double probability, double lambda)
{
	const double lower_log_kappa = log(DBL_MIN);
	const double upper_log_kappa = log(DBL_MAX);
	double target_log_distance;
	double lower = lower_log_kappa;
	double upper = upper_log_kappa;
	double distance;
	double unused_log_abs_derivative;

	if (isnan(probability) || isnan(lambda) || !isfinite(lambda) || lambda <= 0.0 || probability < 0.0 || probability > 1.0) {
		return NAN;
	}
	if (probability <= exp(-lambda)) {
		return 0.0;
	}
	if (probability == 1.0) {
		return INFINITY;
	}

	target_log_distance = log(-log(probability)) - log(lambda);
	if (target_log_distance >= 0.0) {
		return 0.0;
	}

	INLAcirc_pc_vminf_geometry_from_log_kappa(upper, &distance, &unused_log_abs_derivative);
	if (log(distance) > target_log_distance) {
		return INFINITY;
	}

	for (int iteration = 0; iteration < 100; ++iteration) {
		const double middle = lower + 0.5 * (upper - lower);

		INLAcirc_pc_vminf_geometry_from_log_kappa(middle, &distance, &unused_log_abs_derivative);
		if (log(distance) > target_log_distance) {
			lower = middle;
		} else {
			upper = middle;
		}
	}

	return exp(lower + 0.5 * (upper - lower));
}

/* log{I0(kappa) exp(-kappa)} for kappa supplied on its natural scale. */
static inline double INLAcirc_log_bessel_i0_scaled(double kappa)
{
	return log(INLAcirc_bessel_i(kappa, 0.0, 1));
}

/*
 * log{I0(kappa) exp(-kappa)} for log(kappa).  The asymptotic branch avoids
 * evaluating exp(log_kappa) once the concentration is extremely large.
 */
static inline double INLAcirc_log_bessel_i0_scaled_from_log_kappa(double log_kappa)
{
	if (log_kappa > 11.512925) {
		return -0.5 * (INLACIRC_LOG_2PI + log_kappa) + 0.125 * exp(-log_kappa);
	}
	return INLAcirc_log_bessel_i0_scaled(exp(log_kappa));
}

#endif
