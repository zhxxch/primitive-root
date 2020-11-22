#pragma once
#ifndef __cplusplus
#error
#endif
#include<complex>
#include<iterator>
template<int forward, int size, typename R_t>
inline void fft_c()
template<int forward, typename out_iter_t, typename in_iter_t>
inline void fft(out_iter_t Out, in_iter_t In, const int Size){
	auto reverse_bit = [Size](int Idx)->int{
		int RevIdx = 0;
		for(int SrcIdx = Size, DstShf = 0; SrcIdx >>= 1; DstShf++)
			RevIdx += ((Idx&SrcIdx) / SrcIdx) << DstShf;
		return(RevIdx); };
	for(int i = 0; i < Size; i++)
		Out[i] = In[reverse_bit(i)];
	using complex = typename std::iterator_traits<out_iter_t>::value_type;
	complex U1Gen = complex(-1, 0);
	for(int SubFTSize = 2; SubFTSize <= Size; SubFTSize <<= 1){
		for(int SubFTIt = 0; SubFTIt < Size; SubFTIt += SubFTSize){
			complex W = complex(1, 0);
			for(int i = 0; i < SubFTSize / 2; i++){
				const complex Even = Out[SubFTIt + i],
				Odd = W * Out[SubFTIt + i + SubFTSize / 2];
				Out[SubFTIt + i] = Even + Odd;
				Out[SubFTIt + i + SubFTSize / 2] = Even - Odd;
				W *= forward > 0 ? std::conj(U1Gen) : U1Gen;
			}
		}
		U1Gen = sqrt(U1Gen);
	}
}