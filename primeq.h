#pragma once
#include <math.h>
inline double f64high24bit(
	const double x, const double scalar) {
	const double d = scalar * (1 << 24);
	const double shr24 = x / d;
	const double hi24 = floor(shr24) * d;
	return hi24;
}
inline double f64_mul_48(const double ah,
	const double al, const double bh,
	const double bl) {
	const double ll = al * bl;
	const double m1 = ah * bl;
	const double m2 = al * bh;
	const double ll_hi = f64high24bit(ll, 1);
	const double ll_lo = ll - ll_hi;
	const double mm = m1 + m2 + ll_hi;
	const double mm_lo
		= mm - f64high24bit(mm, (1 << 24));
	return mm_lo + ll_lo;
}
inline double f48_mul_mod(const double a,
	const double b, const double mod_hi,
	const double mod_lo) {
	const double prod = a * b;
	const double scalar
		= floor(prod / (mod_hi + mod_lo));
	const double a_hi = f64high24bit(a, 1);
	const double a_lo = a - a_hi;
	const double b_hi = f64high24bit(b, 1);
	const double b_lo = b - b_hi;
	const double scalar_hi = f64high24bit(scalar, 1);
	const double scalar_lo = scalar - scalar_hi;
	const double diff_hi
		= a_hi * b_hi - scalar_hi * mod_hi;
	const double diff_m1
		= a_hi * b_lo - scalar_hi * mod_lo;
	const double diff_m2
		= a_lo * b_hi - scalar_lo * mod_hi;
	const double diff_lo
		= a_lo * b_lo - scalar_lo * mod_lo;
	return diff_hi + diff_m1 + diff_m2 + diff_lo;
}
double remove_factor2(const double x) {
	for(double hx = x;;) {
		const double hhx = hx / 2;
		if(floor(hhx) != hhx) return hx;
		hx = hhx;
	}
}
void primeq_rabin_miller(const double n,
	double bases[], double answs[],
	const int ans_len) {
	const double d = remove_factor2(n - 1);
	const double s = (n - 1) / d;
	typedef long long i64;
	const double n_hi = f64high24bit(n, 1);
	const double n_lo = n - n_hi;
	for(int i = 0; i < ans_len; i++) { answs[i] = 1; }
	for(i64 exp = (i64)d; exp; exp /= 2) {
		if(exp & 1) {
			for(int i = 0; i < ans_len; i++) {
				answs[i] = f48_mul_mod(
					answs[i], bases[i], n_hi, n_lo);
			}
			for(int i = 0; i < ans_len; i++) {
				answs[i]
					= (answs[i] < 0 ? (answs[i] + n)
									: answs[i]);
			}
		}
		for(int i = 0; i < ans_len; i++) {
			bases[i] = f48_mul_mod(
				bases[i], bases[i], n_hi, n_lo);
		}
		for(int i = 0; i < ans_len; i++) {
			bases[i] = (bases[i] < 0 ? (bases[i] + n)
									 : bases[i]);
		}
	}
	for(int i = 0; i < ans_len; i++) {
		bases[i] = ((answs[i] == 1) ? 0 : answs[i]);
	}
	for(i64 k = (i64)s; k; k >>= 1) {
		for(int i = 0; i < ans_len; i++) {
			answs[i]
				= (bases[i] == (n - 1) ? 1 : answs[i]);
		}
		for(int i = 0; i < ans_len; i++) {
			bases[i] = f48_mul_mod(
				bases[i], bases[i], n_hi, n_lo);
		}
		for(int i = 0; i < ans_len; i++) {
			bases[i] = (bases[i] < 0 ? (bases[i] + n)
									 : bases[i]);
		}
	}
}

int primeq(const double n) {
	/* Ref: https://primes.utm.edu/prove/prove2_3.html */
	double Primes[] = {2, 3, 5, 7, 11, 13, 17, 19};
	double Ans[sizeof(Primes) / sizeof(double)];
	if(n <= 1 || n >= (1ll << 48)) return 0;
	for(int i = 0; i < sizeof(Primes) / sizeof(double);
		i++) {
		if(n == Primes[i]) return (1);
	}
	primeq_rabin_miller(
		n, Primes, Ans, sizeof(Ans) / sizeof(double));
	for(int i = 0; i < sizeof(Ans) / sizeof(double);
		i++) {
		if(Ans[i] != 1) return 0;
	}
	return 1;
}