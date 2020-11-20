#include <stdio.h>
#include <stdint.h>
#include "speck.h"
inline int32_t mod_p3013(const uint64_t v) {
	const uint64_t f1 = v & (0x3ffull << 51);
	const uint64_t m1 = f1 + (f1 >> 17) + (f1 >> 30);
	const int64_t vi1 = v - m1;
	const uint64_t v1
		= vi1 + (vi1 < 0 ? (0x40002001ull << 20) : 0);
	const uint64_t f2 = v1 & (0x3ffull << 41);
	const uint64_t m2 = f2 + (f2 >> 17) + (f2 >> 30);
	const int64_t vi2 = v1 - m2;
	const uint64_t v2
		= vi2 + (vi2 < 0 ? (0x40002001ull << 10) : 0);
	const uint64_t f3 = v2 & (0x7ffull << 30);
	const uint64_t m3 = f3 + (f3 >> 17) + (f3 >> 30);
	const int64_t vi4 = v2 - m3;
	const uint64_t v4
		= vi4 + (vi4 < 0 ? 0x40002001ull : 0);
	return (int32_t)(v4 & 0xffffffffull);
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
	const int32_t res = mod_p3013(a * b);
	printf("%s:\t0x%.8x, 0x%.8x (ans, res)\n",
		(ans == res ? "Eq" : "Neq"), ans, res);
	return test_modp3013(Num - 1, Nonce + 100001);
}
int main(void) { return test_modp3013(10, 2222); }