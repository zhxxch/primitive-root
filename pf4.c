#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include "speck.h"
const uint32_t pf4 = 0x10001u;
const uint32_t pf4_factors[] = {2};
const uint32_t pf4_orders[] = {0x10000 / 2};
inline uint32_t add_mod_pf4(
	const uint32_t a, const uint32_t b) {
	return (a + b) % 0x10001u;
}
inline uint32_t mul_mod_pf4(
	const uint32_t a, const uint32_t b) {
	return (uint32_t)(
		((uint64_t)a * (uint64_t)b) % 0x10001ull);
}
inline uint32_t add_mod_pf4_shr(
	const uint32_t a, const uint32_t b) {
	const uint32_t acc = a + b;
	const uint32_t srem
		= (acc & 0xffffu) - (acc >> 16);
	const uint32_t rem = srem
		+ ((srem & 0x80000000u) >> 15)
		+ ((srem & 0x80000000u) >> 31);
	return rem;
}
inline uint32_t mul_mod_pf4_shr(
	const uint32_t a, const uint32_t b) {
	const uint32_t pd_lo = a * b;
	const uint32_t pd_hi = a & b & 0x10000u;
	const uint32_t srem = (pd_lo & 0xffffu) - pd_hi
		- ((pd_lo & 0xffff0000u) >> 16);
	const uint32_t rem = srem
		+ ((srem & 0x80000000u) >> 15)
		+ ((srem & 0x80000000u) >> 31);
	return rem;
}
#if 0
void vec_modp3013(const int Num, int *restrict Iter) {
	for(int i = 0; i < Num; i++) {
		Iter[i] = mod_p3013(Iter[i]);
	}
}
#endif
uint64_t pow_Zp_ring(const uint64_t p,
	const uint64_t base, const uint64_t exponent) {
	uint64_t ans = 1;
	for(uint64_t exp = exponent, g = base; exp;) {
		if(exp & 1) { ans = (ans * g) % p; }
		exp = (exp >> 1);
		g = (g * g) % p;
	}
	return ans;
}
bool is_pf4_gen(const int32_t a) {
	for(int i = 0; i < 3; i++) {
		const uint64_t r
			= pow_Zp_ring(0x10001u, a, pf4_orders[i]);
		if(r <= 1) return false;
	}
	return true;
}
int search_pf4_gens(const int Num) {
	for(int i = 0; i < Num; i++) {
		if(is_pf4_gen(i)) {
			printf("Generator: %i\n", i);
		}
	}
	return 0;
}
int test_mod_pf4(const int Num, const int Nonce) {
	if(Num <= 0) return 0;
	const uint64_t a
		= speck64u96(Nonce, Nonce + 3, Nonce + 5, 99)
		% 0x10001u;
	const uint64_t b
		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
		% 0x10001u;
	const uint64_t add_ans = (a + b) % 0x10001ull;
	const uint64_t add_res = add_mod_pf4_shr(a, b);
	const uint64_t mul_ans = (a * b) % 0x10001ull;
	const uint64_t mul_res = mul_mod_pf4_shr(a, b);
	if(add_ans != add_res)
		printf("%s:\t0x%.8llx, 0x%.8llx (ans, res)\n",
			(add_ans == add_res ? "Add Eq"
								: "Add Neq"),
			add_ans, add_res);
	if(mul_ans != mul_res)
		printf("%s:\t0x%.8llx, 0x%.8llx (ans, res)\n",
			(mul_ans == mul_res ? "Mul Eq"
								: "Mul Neq"),
			mul_ans, mul_res);
	return test_mod_pf4(Num - 1, Nonce + 100001);
}
void vec_mod_pf4(const int Num,
	int32_t *restrict IterRes, int32_t *restrict IterA,
	int32_t *restrict IterB) {
	for(int i = 0; i < Num; i++) {
#if 1
		IterRes[i]
			= mul_mod_pf4_shr(IterA[i], IterB[i]);
#else
		IterRes[i] = (IterA[i] * IterB[i]) % 0x10001u;
#endif
	}
}
void fill_list_rand(const int Num, int *restrict List,
	const int Nonce) {
	for(int i = 0; i < Num; i++) {
		List[i] = (int)speck64u96(
			100, i, Nonce, i + 999999);
	}
}
int main(void) {
	const int VecLen = 1000000;
	int32_t *Vecs
		= malloc(VecLen * 3 * sizeof Vecs[0]);
	fill_list_rand(VecLen * 3, Vecs, 223344);
	for(int i = 0; i < 1000; i++) {
		vec_mod_pf4(VecLen, Vecs, Vecs + VecLen,
			Vecs + VecLen * 2);
	}
	return 0;
}