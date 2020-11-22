#include <stdio.h>
#include <stdlib.h>
int cmp_bit_rev(const int*restrict a, const int*restrict b) {
	const int d = *a ^ *b;
	const int d_lsb1 = d ^ (d - 1);
	return (d_lsb1 & *a) - (d_lsb1 & *b);
}
void fill_iota(const int Num, int Iter[]) {
	for(int i = 0; i < Num; i++) { Iter[i] = i; }
}
void print_arr(const int Num, const int Iter[]) {
	for(int i = 0; i < Num; i++) {
		printf("%i\n", Iter[i]);
	}
}
int main(void) {
	const int Num = 8;
	int *Buffer = malloc(Num * sizeof(Buffer[0]));
	fill_iota(Num, Buffer);
	qsort(Buffer, Num, sizeof Buffer[0], cmp_bit_rev);
	print_arr(Num, Buffer);
	return 0;
}