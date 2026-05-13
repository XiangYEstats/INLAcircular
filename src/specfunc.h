#ifndef SPECFUNC_H
#define SPECFUNC_H

/* Link and inverse CDF functions */
double c_qlogis(double p);
double c_plogis(double x);
double c_pnorm(double x);
double c_qnorm(double p);

/* Bessel functions */
double cc_bessel_i(double x, double nu, int scaled);
double get_log_bessel_scaled(double th);

#endif
