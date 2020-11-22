#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include "speck.h"
const int p3013_factors[] = {2, 3, 43691};
const int p3013_orders[] = {0x40002000 / 2,
	0x40002000 / 3, 0x40002000 / 43691};
inline int32_t mod_p3013(const uint64_t v) {
	const uint64_t f1 = v & (0x3ffull << 51);
	const uint64_t m1 = f1 + (f1 >> 17) + (f1 >> 30);
	const uint64_t vi1 = v - m1;
	const uint64_t v1
		= vi1 + (vi1 > v ? (0x40002001ull << 20) : 0);
	const uint64_t f2 = v1 & (0x3ffull << 41);
	const uint64_t m2 = f2 + (f2 >> 17) + (f2 >> 30);
	const uint64_t vi2 = v1 - m2;
	const uint64_t v2
		= vi2 + (vi2 > v1 ? (0x40002001ull << 10) : 0);
	const uint64_t f3 = v2 & (0x7ffull << 30);
	const uint64_t m3 = f3 + (f3 >> 17) + (f3 >> 30);
	const uint64_t vi3 = v2 - m3;
	const uint64_t v3
		= vi3 + (vi3 > v2 ? 0x40002001ull : 0);
	return (int32_t)(v3 & 0xffffffffull);
}
inline int32_t mod_p3013_2(const uint64_t v) {
	const uint64_t c2 = 0x1ffffffffull/0x10001ull+1;
	const uint64_t c = UINT64_MAX / 0x40002001ull + 1;
	const uint64_t m = (c * v * 0x40002001ull) >> 32;
	return (int32_t)m;
}
int test_modp3013(const int Num, const int Nonce) {
	if(Num <= 0) return 0;
	const uint64_t a
		= speck64u96(Nonce, Nonce + 3, Nonce + 5, 99)
		% 0x40002001;
	const uint64_t b
		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
		% 0x40002001;
	const int32_t ans = (a * b) % 0x40002001;
	const int32_t res = mod_p3013_2(a * b);
	if(ans != res)
		printf("%s:\t0x%.8x, 0x%.8x (ans, res)\n",
			(ans == res ? "Eq" : "Neq"), ans, res);
	return test_modp3013(Num - 1, Nonce + 100001);
}
void vec_modp3013(const int Num, int *restrict Iter) {
	for(int i = 0; i < Num; i++) {
		Iter[i] = mod_p3013(Iter[i]);
	}
}
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
bool is_p3013_gen(const int32_t a) {
	for(int i = 0; i < 3; i++) {
		const uint64_t r = pow_Zp_ring(
			0x40002001, a, p3013_orders[i]);
		if(r <= 1) return false;
	}
	return true;
}
int search_p3013_gens(const int Num) {
	for(int i = 0; i < Num; i++) {
		if(is_p3013_gen(i)) {
			printf("Generator: %i\n", i);
		}
	}
	return 0;
}
void fill_list_rand(
	const int Num, int *restrict List) {
	for(int i = 0; i < Num; i++) {
		List[i] = (int)speck64u96(100, 200, i, i + 2);
	}
}
int main(void) {
	const int Num = 10000;
	int *List = malloc(Num * sizeof List[0]);
	fill_list_rand(Num, List);
	vec_modp3013(Num, List);
	test_modp3013(10000, 100);
	return search_p3013_gens(10);
}