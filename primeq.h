#pragma once
void primeq_rabin_miller(const unsigned long long n,
	unsigned long long bases[],
	unsigned long long answs[], const int ans_len) {
	typedef unsigned long long u64;
	const u64 m = n - 1;
	const u64 s = m & (~m + 1);
	const u64 d = m / s;
	for(int i = 0; i < ans_len; i++) { answs[i] = 1; }
	for(u64 exp = d; exp; exp /= 2) {
		if(exp & 1) {
			for(int i = 0; i < ans_len; i++) {
				answs[i] = (answs[i] * bases[i]) % n;
			}
		}
		for(int i = 0; i < ans_len; i++) {
			bases[i] = (bases[i] * bases[i]) % n;
		}
	}
	for(int i = 0; i < ans_len; i++) {
		bases[i] = ((answs[i] == 1) ? 0 : answs[i]);
	}
	for(u64 k = s; k; k >>= 1) {
		for(int i = 0; i < ans_len; i++) {
			answs[i] = (bases[i] == m ? 1 : answs[i]);
		}
		for(int i = 0; i < ans_len; i++) {
			bases[i] = (bases[i] * bases[i]) % n;
		}
	}
}

int primeq(unsigned long long n) {
	/* Ref: https://primes.utm.edu/prove/prove2_3.html */
	typedef unsigned long long u64;
	u64 Primes[] = {2, 3, 5, 7, 11, 13, 17, 19};
	u64 Ans[sizeof(Primes) / sizeof(u64)];
	if(n <= 1) return 0;
	for(int i = 0; i < sizeof(Primes) / sizeof(u64);
		i++) {
		if(n == Primes[i]) return (1);
		if((n % Primes[i]) == 0) return (0);
	}
	primeq_rabin_miller(
		n, Primes, Ans, sizeof(Ans) / sizeof(u64));
	for(int i = 0; i < sizeof(Ans) / sizeof(u64);
		i++) {
		if(Ans[i] != 1) return 0;
	}
	return 1;
}