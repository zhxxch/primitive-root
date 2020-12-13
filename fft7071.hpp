#pragma once
/***__ |  \   /   __|
****  /     /    |
****___|  _/ _\ \___|
	fft7071.hpp (2020) copyright, github.com/zhxxch
*/
#ifndef __cplusplus
#error "__cplusplus not defined"
#endif
#include <complex>
#include <iterator>
#include <cmath>
#include <cstdint>
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
	const double phase, const double N,
	const double k) {
	using namespace std::complex_literals;
	const double theta_l = k / N;
	const double theta = theta_l - floor(theta_l);
	if(theta == 0) return 1;
	if(theta == 1. / 2) return -1;
	const double theta2
		= theta > 1. / 2 ? 1 - theta : theta;
	const double p = theta > 1. / 2 ? -phase : phase;
	if(theta2 == 3. / 8)
		return (-root7071 + p * root7071 * 1.0i);
	if(theta2 == 2. / 8) return (p * 1.0i);
	if(theta2 == 1. / 8)
		return (root7071 + p * root7071 * 1.0i);
	const double theta4
		= theta2 > 2. / 8 ? 0.5 - theta2 : theta2;
	const double theta8
		= theta4 > 1. / 8 ? 0.25 - theta4 : theta4;
	const double s = cr_sinpi_hq(theta8 * 8);
	const double c = cr_cospi_hq(theta8 * 8);
	if(theta2 > 3. / 8) return (-c + p * s * 1.0i);
	if(theta2 > 2. / 8) return (-s + p * c * 1.0i);
	if(theta2 > 1. / 8) return (s + p * c * 1.0i);
	return (c + p * s * 1.0i);
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
	const typename T::value_type theta
		= Phase * 2 * freq * pi / MaxFreq;
	return T(cos(theta), sin(theta));
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
	using value_type =
		typename std::iterator_traits<I>::value_type;
	using difference_type =
		typename std::iterator_traits<
			I>::difference_type;
	using pointer =
		typename std::iterator_traits<I>::pointer;
	using reference =
		typename std::iterator_traits<I>::reference;
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
typename std::iterator_traits<I>::difference_type
distance(
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
					typename std::iterator_traits<
						I_x>::value_type x,
					typename std::iterator_traits<
						I_w>::value_type w) {
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
	using Component = typename std::iterator_traits<
		x_iter_t>::value_type;
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
		const Component parit00 = *parit00_it_2;
		const Component parit01 = *parit10_it_2;
		const Component parit10 = *parit01_it_2;
		const Component parit11 = *parit11_it_2;
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
			x_sentinel_t sub_ft_s = sub_ft_it;
			advance(sub_ft_s, num_sub_ft_pair);
			for(; sub_ft_it != sub_ft_s;
				advance(sub_ft_it, 2 * sub_ft_size)) {
				x_iter_t parit00_it = sub_ft_it;
				x_iter_t parit01_it = parit00_it;
				advance(parit01_it, sub_ft_size);
				const x_sentinel_t it_s = parit01_it;
				x_iter_t parit10_it = sub_ft_it;
				advance(parit10_it, num_sub_ft_pair);
				x_iter_t parit11_it = parit10_it;
				advance(parit11_it, sub_ft_size);
				for(w_iter_t nth_pow
					= W_0;
					parit00_it != it_s;
					advance(nth_pow, num_sub_ft_pair),
					++parit00_it, ++parit01_it,
					++parit10_it, ++parit11_it) {
					typename std::iterator_traits<
						w_iter_t>::value_type W
						= *nth_pow;
					const Component parit00
						= *parit00_it;
					const Component parit01
						= *parit10_it * W;
					const Component parit10
						= *parit01_it;
					const Component parit11
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
			const x_sentinel_t parit0_s = parit1_it;
			for(w_iter_t nth_pow
				= W_0;
				parit0_it != parit0_s;
				++parit0_it, ++parit1_it,
				advance(nth_pow, num_sub_ft_pair)) {
				const Component parit0 = *parit0_it;
				const Component parit1
					= (*parit1_it) * (*nth_pow);
				*parit0_it = parit0 + parit1;
				*parit1_it = parit0 - parit1;
			}
		}
	}
}
template<typename interleave_w_it, typename sep_w_it>
#if __cplusplus > 201703L
requires std::output_iterator<sep_w_it,
	typename std::iterator_traits<
		interleave_w_it>::value_type::value_type> &&
	std::forward_iterator<interleave_w_it>
#endif
	void copy_ex_w(interleave_w_it W_0,
		const long long W_LEN, sep_w_it ReW0,
		sep_w_it ImW0) {
	*ReW0++ = 0;
	*ImW0++ = 0;
	for(long long L = W_LEN; L > 0; L /= 2) {
		interleave_w_it w_it = W_0;
		for(long long i = 0; i < W_LEN; i += L) {
			typename std::iterator_traits<
				interleave_w_it>::value_type W
				= *w_it;
			*ReW0++ = std::real(W);
			*ImW0++ = std::imag(W);
			advance(w_it, L);
		}
	}
}
template<typename x_iter_t, typename w_iter_t>
#if __cplusplus > 201703L
requires std::contiguous_iterator<x_iter_t> &&
	std::contiguous_iterator<w_iter_t>
#endif
	void fft_in_situ_sep(x_iter_t ReX0, x_iter_t ImX0,
		w_iter_t ReW0, w_iter_t ImW0,
		const long long N) {
	using Component = typename std::iterator_traits<
		x_iter_t>::value_type;
	long long sub_ft_size = 1;
	long long num_sub_ft = N / sub_ft_size;
	long long num_sub_ft_pair = num_sub_ft / 2;
	++ReW0, ++ImW0;
	for(long long sub_ft_pos = 0;
		sub_ft_pos < num_sub_ft_pair;
		sub_ft_pos += 2) {
		const Component parit00re = ReX0[sub_ft_pos];
		const Component parit00im = ImX0[sub_ft_pos];
		const Component parit01re
			= ReX0[sub_ft_pos + num_sub_ft_pair];
		const Component parit01im
			= ImX0[sub_ft_pos + num_sub_ft_pair];
		const Component parit10re
			= ReX0[sub_ft_pos + 1];
		const Component parit10im
			= ImX0[sub_ft_pos + 1];
		const Component parit11re
			= ReX0[sub_ft_pos + 1 + num_sub_ft_pair];
		const Component parit11im
			= ImX0[sub_ft_pos + 1 + num_sub_ft_pair];
		ReX0[sub_ft_pos] = parit00re + parit01re;
		ImX0[sub_ft_pos] = parit00im + parit01im;
		ReX0[sub_ft_pos + 1] = parit00re - parit01re;
		ImX0[sub_ft_pos + 1] = parit00im - parit01im;
		ReX0[sub_ft_pos + num_sub_ft_pair]
			= parit10re + parit11re;
		ImX0[sub_ft_pos + num_sub_ft_pair]
			= parit10im + parit11im;
		ReX0[sub_ft_pos + 1 + num_sub_ft_pair]
			= parit10re - parit11re;
		ImX0[sub_ft_pos + 1 + num_sub_ft_pair]
			= parit10im - parit11im;
	}
	for(ReW0 += sub_ft_size, ImW0 += sub_ft_size,
		sub_ft_size *= 2, num_sub_ft /= 2,
		num_sub_ft_pair /= 2;
		sub_ft_size < num_sub_ft_pair;
		ReW0 += sub_ft_size, ImW0 += sub_ft_size,
		sub_ft_size *= 2, num_sub_ft /= 2,
		num_sub_ft_pair /= 2) {
		for(long long perm_pos = 0; perm_pos < N;
			perm_pos += 2 * num_sub_ft_pair) {
			for(long long sub_ft_pos = perm_pos;
				sub_ft_pos
				< perm_pos + num_sub_ft_pair;
				sub_ft_pos += 2 * sub_ft_size) {
				const x_iter_t parit00_re_it
					= ReX0 + sub_ft_pos;
				const x_iter_t parit00_im_it
					= ImX0 + sub_ft_pos;
				const x_iter_t parit01_re_it
					= ReX0 + sub_ft_pos + sub_ft_size;
				const x_iter_t parit01_im_it
					= ImX0 + sub_ft_pos + sub_ft_size;
				const x_iter_t parit10_re_it = ReX0
					+ sub_ft_pos + num_sub_ft_pair;
				const x_iter_t parit10_im_it = ImX0
					+ sub_ft_pos + num_sub_ft_pair;
				const x_iter_t parit11_re_it = ReX0
					+ sub_ft_pos + sub_ft_size
					+ num_sub_ft_pair;
				const x_iter_t parit11_im_it = ImX0
					+ sub_ft_pos + sub_ft_size
					+ num_sub_ft_pair;
#pragma omp simd
				for(long long i = 0; i < sub_ft_size;
					i++) {
					const Component w_re = ReW0[i];
					const Component w_im = ImW0[i];
					const Component parit10q_re
						= parit10_re_it[i];
					const Component parit10q_im
						= parit10_im_it[i];
					const Component parit00_re
						= parit00_re_it[i];
					const Component parit00_im
						= parit00_im_it[i];
					const Component parit01_re
						= parit10_re_it[i] * w_re
						- parit10_im_it[i] * w_im;
					const Component parit01_im
						= parit10_re_it[i] * w_im
						+ parit10_im_it[i] * w_re;
					const Component parit10_re
						= parit01_re_it[i];
					const Component parit10_im
						= parit01_im_it[i];
					const Component parit11_re
						= parit11_re_it[i] * w_re
						- parit11_im_it[i] * w_im;
					const Component parit11_im
						= parit11_re_it[i] * w_im
						+ parit11_im_it[i] * w_re;
					parit00_re_it[i]
						= parit00_re + parit01_re;
					parit00_im_it[i]
						= parit00_im + parit01_im;
					parit01_re_it[i]
						= parit00_re - parit01_re;
					parit01_im_it[i]
						= parit00_im - parit01_im;
					parit10_re_it[i]
						= parit10_re + parit11_re;
					parit10_im_it[i]
						= parit10_im + parit11_im;
					parit11_re_it[i]
						= parit10_re - parit11_re;
					parit11_im_it[i]
						= parit10_im - parit11_im;
				}
			}
		}
	}
	for(; sub_ft_size < N; ReW0 += sub_ft_size,
		ImW0 += sub_ft_size, sub_ft_size *= 2,
		num_sub_ft /= 2, num_sub_ft_pair /= 2) {
		for(long long sub_ft_pos = 0; sub_ft_pos < N;
			sub_ft_pos += 2 * sub_ft_size) {
			const x_iter_t parit0_re_it
				= ReX0 + sub_ft_pos;
			const x_iter_t parit0_im_it
				= ImX0 + sub_ft_pos;
			const x_iter_t parit1_re_it
				= ReX0 + sub_ft_pos + sub_ft_size;
			const x_iter_t parit1_im_it
				= ImX0 + sub_ft_pos + sub_ft_size;
#pragma omp simd
			for(long long i = 0; i < sub_ft_size;
				i++) {
				const Component w_re = ReW0[i];
				const Component w_im = ImW0[i];
				const Component parit0_re
					= parit0_re_it[i];
				const Component parit0_im
					= parit0_im_it[i];
				const Component parit1_re
					= parit1_re_it[i] * w_re
					- parit1_im_it[i] * w_im;
				const Component parit1_im
					= parit1_re_it[i] * w_im
					+ parit1_im_it[i] * w_re;
				parit0_re_it[i]
					= parit0_re + parit1_re;
				parit0_im_it[i]
					= parit0_im + parit1_im;
				parit1_re_it[i]
					= parit0_re - parit1_re;
				parit1_im_it[i]
					= parit0_im - parit1_im;
			}
		}
	}
}
template<typename x_iter_t, typename w_iter_t>
void real_conv(x_iter_t A_0, x_iter_t B_0,
	w_iter_t ReWpos, w_iter_t ImWpos, w_iter_t ReWneg,
	w_iter_t ImWneg, const long long N) {
	fft_in_situ_sep(A_0, B_0, ReWneg, ImWneg, N);
	using real_t = typename std::iterator_traits<
		x_iter_t>::value_type;
	A_0[0] = A_0[0] * B_0[0];
	B_0[0] = 0;
	A_0[N / 2] = A_0[N / 2] * B_0[N / 2];
	B_0[N / 2] = 0;
	for(long long i = 1; i < N / 2; i++) {
		const real_t a = A_0[i];
		const real_t b = B_0[i];
		const real_t c = A_0[N - i];
		const real_t d = B_0[N - i];
		A_0[i] = 0.5 * (a * b + c * d);
		B_0[i]
			= 0.25 * (b * b + c * c - d * d - a * a);
		A_0[N - i] = 0.5 * (a * b + c * d);
		B_0[N - i]
			= 0.25 * (d * d + a * a - b * b - c * c);
	}
	fft_in_situ_sep(A_0, B_0, ReWpos, ImWpos, N);
	const real_t s = 1. / (real_t)N;
#pragma omp simd
	for(long long i = 0; i < N; i++) { A_0[i] *= s; }
}
template<typename w_iter_t, typename chirp_iter_t>
void chirp_z_modulator(w_iter_t W_2czN_neg,
	const long long CZ_N, chirp_iter_t M_0) {
	for(long long i = 0; i < CZ_N; i++) {
		const long long k = (i * i) % (2 * CZ_N);
		w_iter_t w = W_2czN_neg;
		advance(w, k);
		*M_0++ = *w;
	}
}
template<typename w_iter_t, typename chirp_iter_t>
void chirp_z_filter(w_iter_t W_2czN_pos,
	const long long CZ_N, const long long FFT_N,
	chirp_iter_t F_0) {
	for(long long i = 0; i < CZ_N; i++) {
		const long long k = (i * i) % (2 * CZ_N);
		w_iter_t w = W_2czN_pos;
		advance(w, k);
		F_0[i] = *w;
	}
	for(long long i = 1; i < CZ_N; i++) {
		const long long k = (i * i) % (2 * CZ_N);
		w_iter_t w = W_2czN_pos;
		advance(w, k);
		F_0[FFT_N - i] = *w;
	}
}
template<typename x_iter_t, typename cz_m_it_t,
	typename cz_f_it_t, typename fft_w_it_t>
void chirp_z(x_iter_t X_0, const long long CZ_N,
	const long long FFT_N, cz_m_it_t cz_modulator,
	cz_f_it_t ft_cz_filter, fft_w_it_t Wneg,
	fft_w_it_t Wpos) {
	for(long long i = 0; i < CZ_N; i++) {
		X_0[i] = X_0[i] * cz_modulator[i];
	}
	fft_in_situ(X_0, X_0 + FFT_N, Wneg);
	for(long long i = 0; i < FFT_N; i++) {
		X_0[i] = X_0[i] * ft_cz_filter[i];
	}
	fft_in_situ(X_0, X_0 + FFT_N, Wpos);
	for(long long i = 0; i < CZ_N; i++) {
		X_0[i] = X_0[i] * cz_modulator[i];
	}
}
} // namespace fft7071
