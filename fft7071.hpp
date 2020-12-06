#pragma once
/***__ |  \   /   __|
****  /     /    |
****___|  _/ _\ \___|
	fft7071.hpp (2020)
*/
#ifndef __cplusplus
#error "__cplusplus not defined"
#endif
#include <complex>
#include <vector>
#include <iterator>
#include <stdint.h>
namespace fft7071 {
using std::advance;
using std::distance;
const double pi = 3.141592653589793238462643383;
/*
0x1.921fb54442d18p+1<pi<0x1.921fb54442d19p+1
*/
const double pi_0 = 0x1.921fb54442d18p+1;
const double pi_1 = 0x1.921fb54442d19p+1;
const double pi_0e = 1.2246467991473531772e-16;
const double pi_1e = 3.2162452993532729845e-16;
/* 0x1.6a09e667f3bcdp-1 */
const double root7071 = 0.7071067811865475244008443621;

const uint64_t SINPI_4_TAYLOR[13][3]
	= {{0xc4c6628cull << 32, 0x2168c234ull << 32,
		   0xc90fdaa2ull << 32},
		{0xabb8e5eeull << 32, 0x25be52beull << 32,
			0x14abbce6ull << 32},
		{0x923f3422ull << 32, 0x3bad570eull << 32,
			0xa335e3ull << 32},
		{0xb7ccbdc3ull << 32, 0x99cc57b0ull << 32,
			0x265a5ull << 32},
		{0xe06e41ebull << 32, 0xe0d21fb9ull << 32,
			0x541ull << 32},
		{0x21c7d88full << 32, 0x8c1d3f7aull << 32,
			0x7ull << 32},
		{0x406358cdull << 32, 0x7a3d0d3ull << 32},
		{0xe7c555d1ull << 32, 0x5beb6ull << 32},
		{0xd8655f26ull << 32, 0x355ull << 32},
		{0x8a404212ull << 32, 0x1ull << 32},
		{0x943b81ull << 32}, {0x2e43ull << 32},
		{0xcull << 32}};
const uint64_t COSPI_4_TAYLOR[12][3]
	= {{0x2b71366dull << 32, 0xf9177969ull << 32,
		   0x4ef4f326ull << 32},
		{0xd4cc0780ull << 32, 0x06d6b0ecull << 32,
			0x40f07c2ull << 32},
		{0xfc54fadbull << 32, 0x7e3cbff9ull << 32,
			0x155d3cull << 32},
		{0x75e8c9d1ull << 32, 0xa0d12375ull << 32,
			0x3c3eull << 32},
		{0x2a2ea69full << 32, 0xb47ca881ull << 32,
			0x69ull << 32},
		{0xd8f30a37ull << 32, 0x7e74e28dull << 32},
		{0xd12c4a3dull << 32, 0x6db893ull << 32},
		{0x8b0bcb60ull << 32, 0x4831ull << 32},
		{0x418b235full << 32, 0x25ull << 32},
		{0xf7b7184ull << 32}, {0x54ab9ull << 32},
		{0x184ull << 32}};
inline uint64_t adc_iter_le64(const uint64_t carry,
	const uint64_t src1, const uint64_t src2,
	uint64_t *dst) {
	const uint64_t res = src1 + src2 + carry;
	*dst = res;
	return (res < src1) | (res < src2);
}
inline uint64_t neg_le64(const uint64_t num,
	const uint64_t src[], uint64_t dst[]) {
	uint64_t carry = 1;
	for(uint64_t i = 0; i < num; i++) {
		carry = adc_iter_le64(
			carry, ~src[i], 0, dst + i);
	}
	return carry;
}
inline uint64_t adc_le64(const uint64_t num,
	const uint64_t src1[], const uint64_t src2[],
	uint64_t dst[]) {
	uint64_t carry = 0;
	for(uint64_t i = 0; i < num; i++) {
		carry = adc_iter_le64(
			carry, src1[i], src2[i], dst + i);
	}
	return carry;
}
inline void round96_3le64(
	const uint64_t src[3], uint64_t dst[3]) {
	dst[0] = src[1] << 32;
	dst[2] = src[2] & (uint64_t)(-(1ll << 32));
	dst[1] = src[2] << 32;
}
inline uint64_t mul_3le64(const uint64_t src1[3],
	const uint64_t src2[3], uint64_t dst[3]) {
	uint64_t s1[3], s2[3];
	for(int i = 0; i < 3; i++) {
		s1[i] = src1[i] >> 32;
		s2[i] = src2[i] >> 32;
	}
	uint64_t im[3][3] = {0};
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			im[i][j] = s1[i] * s2[j];
		}
	}
	uint64_t carry_i = 0;
	carry_i = adc_iter_le64(
		0, im[0][1], im[1][0], &im[1][0]);
	carry_i = adc_iter_le64(
		carry_i, im[1][2], im[2][1], &im[2][1]);
	adc_iter_le64(
		0, carry_i << 32, im[2][2], &im[2][2]);
	carry_i = adc_iter_le64(
		0, im[0][2], im[1][1], &im[1][1]);
	adc_iter_le64(carry_i, im[2][2], 0, &im[2][2]);
	carry_i = adc_iter_le64(
		0, im[1][1], im[2][0], &im[2][0]);
	adc_iter_le64(carry_i, im[2][2], 0, &im[2][2]);
	carry_i = adc_iter_le64(
		0, im[0][0] >> 32, im[1][0], &im[1][0]);
	carry_i = adc_iter_le64(0,
		(carry_i << 32) + (im[1][0] >> 32), im[2][0],
		&im[2][0]);
	carry_i = adc_iter_le64(0,
		(carry_i << 32) + (im[2][0] >> 32), im[2][1],
		&im[2][1]);
	carry_i = adc_iter_le64(0,
		(carry_i << 32) + (im[2][1] >> 32), im[2][2],
		&im[2][2]);
	dst[0] = im[2][1] << 32;
	round96_3le64(im[2], dst);
	return carry_i;
}
inline double cr_sinpi_hq(const double frac) {
	const uint64_t theta = (uint64_t)(0x1p64 * frac);
	const uint64_t x[3] = {0, 0, theta};
	uint64_t ans[3] = {0};
	uint64_t x2[3] = {0};
	mul_3le64(x, x, x2);
	mul_3le64(x2, SINPI_4_TAYLOR[12], ans);
	for(int i = 11; i > 0; i--) {
		neg_le64(3, ans, ans);
		adc_le64(3, SINPI_4_TAYLOR[i], ans, ans);
		mul_3le64(x2, ans, ans);
	}
	neg_le64(3, ans, ans);
	adc_le64(3, SINPI_4_TAYLOR[0], ans, ans);
	mul_3le64(x, ans, ans);
	const double im2 = 0x1p-128 * (double)ans[0];
	const double im1 = im2 + 0x1p-96 * (double)ans[1];
	const double im0 = im1 + 0x1p-64 * (double)ans[2];
	return im0;
}
inline double cr_cospi_hq(const double frac) {
	const uint64_t theta = (uint64_t)(0x1p64 * frac);
	const uint64_t x[3] = {0, 0, theta};
	uint64_t ans[3] = {0};
	uint64_t x2[3] = {0};
	mul_3le64(x, x, x2);
	mul_3le64(x2, COSPI_4_TAYLOR[11], ans);
	for(int i = 11; i--;) {
		neg_le64(3, ans, ans);
		adc_le64(3, COSPI_4_TAYLOR[i], ans, ans);
		mul_3le64(x2, ans, ans);
	}
	neg_le64(3, ans, ans);
	for(int i = 0; i < 3; i++) {
		ans[i] = ans[i] & (uint64_t)(-(1ll << 32));
	}
	const double im2 = 0x1p-128 * (double)ans[0];
	const double im1 = im2 + 0x1p-96 * (double)ans[1];
	const double im0 = im1 + 0x1p-64 * (double)ans[2];
	return im0;
}
inline std::complex<double> c_unit_root_exp(
	const double phase, const double max_freq,
	const double k) {
	using namespace std::complex_literals;
	const double theta = k / max_freq;
	if(theta == 0) return 1;
	if(theta == 3. / 8)
		return (-root7071 + phase * root7071 * 1.0i);
	if(theta == 2. / 8) return (phase * 1.0i);
	if(theta == 1. / 8)
		return (root7071 + phase * root7071 * 1.0i);
	const double theta2
		= theta > 2. / 8 ? 0.5 - theta : theta;
	const double theta4
		= theta2 > 1. / 8 ? 0.25 - theta2 : theta2;
	const double s = cr_sinpi_hq(theta4 * 8);
	const double c = cr_cospi_hq(theta4 * 8);
	if(theta > 3. / 8) return (-c + phase * s * 1.0i);
	if(theta > 2. / 8) return (-s + phase * c * 1.0i);
	if(theta > 1. / 8) return (s + phase * c * 1.0i);
	return (c + phase * s * 1.0i);
}
template<typename T> class complex_unit_root_iter {
  public:
	typename T::value_type MaxFreq;
	typename T::value_type freq;
	typename T::value_type Phase;
	using iterator_category
		= std::forward_iterator_tag;
	using value_type = typename T;
	using difference_type = int;
	using pointer = typename T;
	using reference = typename T;
	complex_unit_root_iter(
		const typename T::value_type max_freq,
		const typename T::value_type phase) :
		MaxFreq(max_freq),
		freq(0), Phase(phase){};
	complex_unit_root_iter() :
		MaxFreq(1), freq(0), Phase(1){};
	T operator*() const;
	complex_unit_root_iter<T> operator++(int) {
		complex_unit_root_iter<T> R = *this;
		freq += 1;
		return R;
	};
	complex_unit_root_iter<T> &operator++() {
		freq += 1;
		return *this;
	};

	bool operator==(
		complex_unit_root_iter<T> const &other) const {
		return MaxFreq == other.MaxFreq && freq
			= other.freq && Phase == other.Phase;
	};
	bool operator!=(
		complex_unit_root_iter<T> const &other) const {
		return !(*this == other);
	};
};
template<typename T>
T complex_unit_root_iter<T>::operator*() const {
	return exp(
		Phase * 2 * freq * pi * T(0, 1) / MaxFreq);
}
template<>
std::complex<double> complex_unit_root_iter<
	std::complex<double>>::operator*() const {
	return c_unit_root_exp(Phase, MaxFreq, freq);
}

template<typename T>
void advance(
	complex_unit_root_iter<T> &iter, const size_t n) {
	iter.freq += n;
}

template<typename I> class strided_iterator {

  public:
	I BaseIter;
	size_t LDA;
	using iterator_category
		= std::forward_iterator_tag;
	using value_type = typename I::value_type;
	using difference_type =
		typename I::difference_type;
	using pointer = typename I::pointer;
	using reference = typename I::reference;
	strided_iterator(
		const I iter, const size_t stride) :
		BaseIter(iter),
		LDA(stride){};
	strided_iterator() : BaseIter((void *)0), LDA(0){};
	value_type &operator*() const {
		return *BaseIter;
	};
	strided_iterator<I> operator++(int) {
		strided_iterator<I> R = *this;
		advance(BaseIter, LDA);
		return R;
	};
	strided_iterator<I> &operator++() {
		advance(BaseIter, LDA);
		return *this;
	};

	bool operator==(
		strided_iterator<I> const &other) const {
		return BaseIter == other.BaseIter
			&& LDA == other.LDA;
	};
	bool operator!=(
		strided_iterator<I> const &other) const {
		return !(*this == other);
	};
};
template<typename I>
void advance(
	strided_iterator<I> &iter, const size_t n) {
	advance(iter.BaseIter, n * iter.LDA);
}
template<typename I>
typename I::difference_type distance(
	strided_iterator<I> i, strided_iterator<I> s) {
	return (s.BaseIter - i.BaseIter) / i.LDA;
}
#if __cplusplus > 201703L
template<typename I_x, typename S, typename I_w>
concept fft_in_situ_iters
	= std::sentinel_for<S, I_x> &&
		std::indirectly_movable<I_x, I_x> &&
			std::incrementable<I_x> &&
				std::forward_iterator<I_w> &&requires(
					typename I_x::value_type x,
					typename I_w::value_type w) {
	x + x *w;
	x - x *w;
};
#endif
template<typename x_iter_t, typename x_sentinel_t,
	typename w_iter_t>
#if __cplusplus > 201703L
requires fft_in_situ_iters<x_iter_t, x_sentinel_t,
	w_iter_t>
#endif
	inline void fft_in_situ(
		x_iter_t X_0, x_sentinel_t X_N, w_iter_t W_0) {
	using arith_t = typename x_iter_t::value_type;
	const size_t Length = distance(X_0, X_N);
	size_t sub_ft_size = 1;
	size_t num_sub_ft = Length / sub_ft_size;
	size_t num_sub_ft_pair = num_sub_ft / 2;
	x_iter_t parit00_it_2 = X_0;
	x_iter_t parit01_it_2 = X_0;
	++parit01_it_2;
	x_iter_t parit10_it_2 = X_0;
	advance(parit10_it_2, num_sub_ft_pair);
	x_iter_t parit11_it_2 = parit10_it_2;
	++parit11_it_2;
	while(parit10_it_2 != X_N) {
		const arith_t parit00 = *parit00_it_2;
		const arith_t parit01 = *parit10_it_2;
		const arith_t parit10 = *parit01_it_2;
		const arith_t parit11 = *parit11_it_2;
		*parit00_it_2 = parit00 + parit01;
		*parit01_it_2 = parit00 - parit01;
		*parit10_it_2 = parit10 + parit11;
		*parit11_it_2 = parit10 - parit11;
		advance(parit00_it_2, 2);
		advance(parit01_it_2, 2);
		advance(parit10_it_2, 2);
		advance(parit11_it_2, 2);
	}
	for(sub_ft_size *= 2, num_sub_ft /= 2,
		num_sub_ft_pair /= 2;
		sub_ft_size < num_sub_ft_pair;
		sub_ft_size *= 2, num_sub_ft /= 2,
		num_sub_ft_pair /= 2) {
		for(x_iter_t couple_group_it = X_0;
			couple_group_it != X_N;
			advance(couple_group_it,
				2 * num_sub_ft_pair)) {
			x_iter_t sub_ft_it = couple_group_it;
			x_iter_t sub_ft_s = sub_ft_it;
			advance(sub_ft_s, 2 * sub_ft_size);
			for(; sub_ft_it != sub_ft_s;
				advance(sub_ft_it, 2 * sub_ft_size)) {
				x_iter_t parit00_it = sub_ft_it;
				x_iter_t parit01_it = parit00_it;
				advance(parit01_it, sub_ft_size);
				x_iter_t parit10_it = sub_ft_it;
				advance(parit10_it, num_sub_ft_pair);
				x_iter_t parit11_it = parit10_it;
				advance(parit11_it, sub_ft_size);
				for(w_iter_t nth_pow
					= W_0;
					parit01_it != sub_ft_s;
					advance(nth_pow, num_sub_ft_pair),
					++parit00_it, ++parit01_it,
					++parit10_it, ++parit11_it) {
					typename w_iter_t::value_type W
						= *nth_pow;
					const arith_t parit00
						= *parit00_it;
					const arith_t parit01
						= *parit10_it * W;
					const arith_t parit10
						= *parit01_it;
					const arith_t parit11
						= *parit11_it * W;
					*parit00_it = parit00 + parit01;
					*parit01_it = parit00 - parit01;
					*parit10_it = parit10 + parit11;
					*parit11_it = parit10 - parit11;
				}
			}
		}
	}
	for(; sub_ft_size < Length; sub_ft_size *= 2,
		num_sub_ft /= 2, num_sub_ft_pair /= 2) {
		for(x_iter_t sub_ft_it = X_0; sub_ft_it != X_N;
			advance(sub_ft_it, 2 * sub_ft_size)) {
			x_iter_t parit0_it = sub_ft_it,
					 parit1_it = sub_ft_it;
			advance(parit1_it, sub_ft_size);
			const x_iter_t parit0_s = parit1_it;
			for(w_iter_t nth_pow
				= W_0;
				parit0_it != parit0_s;
				++parit0_it, ++parit1_it,
				advance(nth_pow, num_sub_ft_pair)) {
				const arith_t parit0 = *parit0_it;
				const arith_t parit1
					= (*parit1_it) * (*nth_pow);
				*parit0_it = parit0 + parit1;
				*parit1_it = parit0 - parit1;
			}
		}
	}
}
} // namespace fft7071
