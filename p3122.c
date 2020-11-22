#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include "speck.h"
const int p3133 = 0x80400001ul;
const int p3122_factors[] = {2, 3, 19};
const int p3122_orders[] = {
	0x80400000 / 2, 0x80400000 / 3, 0x80400000 / 19};
inline uint32_t mul_mod_p3122(
	const uint32_t a, const uint32_t b) {
	const uint64_t r62 = (uint64_t)a * b;
	const uint64_t rf1 = r62 & (0x1ffull << 54);
	const uint64_t f1 = rf1 + (rf1 >> 9) + (rf1 >> 31);
	const uint64_t r54
		= r62 - f1 + (0x80400001ull << 22);
	const uint64_t rf2 = r54 & (0x1ffull << 46);
	const uint64_t f2 = rf2 + (rf2 >> 9) + (rf2 >> 31);
	const uint64_t r46
		= r54 - f2 + (0x80400001ull << 14);
	const uint64_t rf3 = r46 & (0x1ffull << 38);
	const uint64_t f3 = rf3 + (rf3 >> 9) + (rf3 >> 31);
	const uint64_t r38
		= r46 - f3 + (0x80400001ull << 6);
	const uint64_t rf4 = r38 & (0xffull << 31);
	const uint64_t f4 = rf4 + (rf4 >> 9) + (rf4 >> 31);
	const uint64_t sr = r38 - f4;
	const uint64_t r
		= sr + (0x80400001ull & (sr >> 32));
	return (uint32_t)r;
}
int test_modp3122(const int Num, const int Nonce) {
	if(Num <= 0) return 0;
	const uint64_t a
		= speck64u96(Nonce, Nonce + 3, Nonce + 5, 99)
		% 0x80400001ull;
	const uint64_t b
		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
		% 0x80400001ull;
	const int32_t ans = (a * b) % 0x80400001ull;
	const int32_t res = mul_mod_p3122(a, b);
	if(ans != res)
		printf("%s:\t0x%.8x, 0x%.8x (ans, res)\n",
			(ans == res ? "Eq" : "Neq"), ans, res);
	return test_modp3122(Num - 1, Nonce + 100001);
}
#if 0
void vec_modp3122(const int Num, int *restrict Iter) {
	for(int i = 0; i < Num; i++) {
		Iter[i] = mul_mod_p3122(Iter[i]);
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
bool is_p3122_gen(const int32_t a) {
	for(int i = 0; i < 3; i++) {
		const uint64_t r = pow_Zp_ring(
			0x40002001, a, p3122_orders[i]);
		if(r <= 1) return false;
	}
	return true;
}
int search_p3122_gens(const int Num) {
	for(int i = 0; i < Num; i++) {
		if(is_p3122_gen(i)) {
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

	test_modp3122(10, 1000);
	return 0;
}