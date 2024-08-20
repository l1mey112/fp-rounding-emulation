#pragma once

double nextafter_1_reg_f(double x);
double nextafter_2_reg_f(double x);
double nextafter_3_reg_f(double x);
double nextafter_1_reg_e(double x);
double nextafter_2_reg_e(double x);
double nextafter_3_reg_e(double x);

double frexp_reg_e_nozero_noinf(double x, int *eptr);
double ldexp_reg_e_nozero_noinf(double x, int n);
