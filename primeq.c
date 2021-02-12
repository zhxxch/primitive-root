#include "primeq.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int ac, char *av[]) {
	if(ac <= 1) { return 1; }
	const uint64_t N = strtoll(av[1], NULL, 10);
	if(ac <= 2) {
		const int prime_bit = primeq(N);
		printf("%i\n", prime_bit);
		return prime_bit;
	}
	const uint64_t N1 = strtoll(av[2], NULL, 10);
	if(ac <= 3) {
		for(uint64_t i = N; i < N1; i++) {
			if(primeq(i)) { printf("%lli\n", i); }
		}
		return 1;
	}
	const double a = strtod(av[1], NULL);
	const double b = strtod(av[2], NULL);
	const double N2 = strtod(av[3], NULL);
	if(ac <= 4) {
		const double a_hi = f64high24bit(a, 1);
		const double a_lo = a - a_hi;
		const double b_hi = f64high24bit(b, 1);
		const double b_lo = b - b_hi;
		const double mod_hi = f64high24bit(N2, 1);
		const double mod_lo = N2 - mod_hi;
		const double ans
			= f48_mul_mod(a, b, mod_hi, mod_lo);
		const double prod
			= f64_mul_48(a_hi, a_lo, b_hi, b_lo);
		printf("%f\n", ans);
		return 0;
	}
	return 1;
}