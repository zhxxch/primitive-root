#include <iostream>
#include <memory>
#include <vector>
#include "fft707.hpp"

inline int fft_rev_bit(const int src_idx) {
	int src_rev_idx = 0;
	for(int src_shf = 3, res_shf = 0; src_shf--;
		res_shf++)
		src_rev_idx += ((src_idx & ((int)1 << src_shf))
						   >> src_shf)
			<< res_shf;
	return (src_rev_idx);
}
template<typename Iter_t>
void fill_iota(const int Num, Iter_t Iter) {
	for(int i = 0; i < Num; i++) { Iter[i] = i; }
}
template<typename Iter_t>
void print_arr(const int Num, const Iter_t Iter) {
	for(int i = 0; i < Num; i++) {
		std::cout << Iter[i] << std::endl;
	}
}
int main(void) {
	const int Num = 16;
	std::vector<std::complex<double>> Buffer(Num);
	std::vector<std::complex<double>> Res(Num);
	fill_iota(Num/2, Buffer.begin());
	/*
	qsort(Buffer, Num, sizeof Buffer[0], cmp_bit_rev);
	*/
	print_arr(Num, Buffer.begin());
	
	printf("%i\n",fft707::ilog2<uint32_t>(1<<12));
	return 0;
}