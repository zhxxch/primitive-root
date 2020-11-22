#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "speck.h"
int adc_u32(const int n, const unsigned int a[],
	const unsigned int b[], unsigned int ans[]) {
	int carry = 0;
	for(int i = 0; i < n; i++) {
		const unsigned int a_i = a[i], b_i = b[i],
						   res_i = a_i + b_i + carry;
		carry = (res_i < a_i) | (res_i < b_i);
		ans[i] = res_i;
	}
	return carry;
}
int main(void) {
	const int Num = 100;
	unsigned int *List = malloc(3 * Num * sizeof List[0]);
	for(int i = 0; i < Num; i++) {
		List[i] = (uint32_t)speck64u96(i, 22, 33, i + 99999);
		List[i + Num]
			= (uint32_t)speck64u96(i + 1, 98, 87, i + 3);
	}
	adc_u32(Num, List, List + Num, List + 2 * Num);
}