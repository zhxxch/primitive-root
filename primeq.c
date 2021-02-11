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
	for(uint64_t i = N; i < N1; i++) {
		if(primeq(i)){
			printf("%lli\n", i);
		}
	}
	return 1;
}