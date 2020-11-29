#pragma once
/***__ |  \   /   __|
****  /     /    |
****___|  _/ _\ \___|
	fft7071.hpp
	zhxxch (2020) All rights reserved
*/
#ifndef __cplusplus
#error "__cplusplus not defined"
#endif
#include <complex>
#include <vector>
#include <stdint.h>
namespace fft7071 {
const double pi = 3.141592653589793238462643383;
/*
0x1.921fb54442d18p+1<pi<0x1.921fb54442d19p+1
*/
const double pi_0 = 0x1.921fb54442d18p+1;
const double pi_1 = 0x1.921fb54442d19p+1;
const double pi_0e = 1.2246467991473531772e-16;
const double pi_1e = 3.2162452993532729845e-16;
const double root7071 = 0.7071067811865475244008443621;

const uint64_t SINPI_4_TAYLOR[13][3] = {
	{
		0xc4c6628cull << 32,
		0x2168c234ull << 32,
		0xc90fdaa2ull << 32,
	},
	{
		0xabb8e5eeull << 32,
		0x25be52beull << 32,
		0x14abbce6ull << 32,
	},
	{
		0x923f3422ull << 32,
		0x3bad570eull << 32,
		0xa335e3ull << 32,
	},
	{
		0xb7ccbdc3ull << 32,
		0x99cc57b0ull << 32,
		0x265a5ull << 32,
	},
	{
		0xe06e41ebull << 32,
		0xe0d21fb9ull << 32,
		0x541ull << 32,
	},
	{
		0x21c7d88full << 32,
		0x8c1d3f7aull << 32,
		0x7ull << 32,
	},
	{
		0x406358cdull << 32,
		0x7a3d0d3ull << 32,
	},
	{
		0xe7c555d1ull << 32,
		0x5beb6ull << 32,
	},
	{
		0xd8655f26ull << 32,
		0x355ull << 32,
	},
	{
		0x8a404212ull << 32,
		0x1ull << 32,
	},
	{
		0x943b81ull << 32,
	},
	{
		0x2e43ull << 32,
	},
	{
		0xcull << 32,
	},
};
const uint64_t COSPI_4_TAYLOR[12][3]
	= {{
		   0x2b71366dull << 32,
		   0xf9177969ull << 32,
		   0x4ef4f326ull << 32,
	   },
		{
			0xd4cc0780ull << 32,
			0x06d6b0ecull << 32,
			0x40f07c2ull << 32,
		},
		{
			0xfc54fadbull << 32,
			0x7e3cbff9ull << 32,
			0x155d3cull << 32,
		},
		{
			0x75e8c9d1ull << 32,
			0xa0d12375ull << 32,
			0x3c3eull << 32,
		},
		{
			0x2a2ea69full << 32,
			0xb47ca881ull << 32,
			0x69ull << 32,
		},
		{
			0xd8f30a37ull << 32,
			0x7e74e28dull << 32,
		},
		{
			0xd12c4a3dull << 32,
			0x6db893ull << 32,
		},
		{
			0x8b0bcb60ull << 32,
			0x4831ull << 32,
		},
		{
			0x418b235full << 32,
			0x25ull << 32,
		},
		{
			0xf7b7184ull << 32,
		},
		{
			0x54ab9ull << 32,
		},
		{
			0x184ull << 32,
		}};
inline uint64_t adc_iter_le64(const uint64_t carry,
	const uint64_t src1, const uint64_t src2,
	uint64_t *dst) {
	const uint64_t res = src1 + src2 + carry;
	*dst = res;
	return (res < src1) | (res < src2);
}
inline uint64_t neg_iter_le64(const uint64_t carry,
	const uint64_t src, uint64_t *dst) {
	return adc_iter_le64(carry, (~src), 0, dst);
}
inline uint64_t neg_le64(const uint64_t num,
	const uint64_t src[], uint64_t dst[]) {
	uint64_t carry = 1;
	for(uint64_t i = 0; i < num; i++) {
		carry = neg_iter_le64(carry, src[i], dst + i);
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
	round96_3le64(im[2], dst);
	return carry_i;
}
inline double sinpihq(const double frac) {
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
inline double cospihq(const double frac) {
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
auto exp_c_unit_root(
	const double phase, const double max_freq) {
	return [phase, max_freq](
			   const std::complex<double> k)
			   -> std::complex<double> {
		using namespace std::complex_literals;
		const double theta = k.real() / max_freq;
		if(theta == 0) return 1;
		if(theta == 3. / 8)
			return (
				-root7071 + phase * root7071 * 1.0i);
		if(theta == 2. / 8) return (phase * 1.0i);
		if(theta == 1. / 8)
			return (
				root7071 + phase * root7071 * 1.0i);
		const double theta2
			= theta > 2. / 8 ? 0.5 - theta : theta;
		const double theta4
			= theta2 > 1. / 8 ? 0.25 - theta2 : theta2;
		const double s = sinpihq(theta4 * 8);
		const double c = cospihq(theta4 * 8);
		if(theta > 3. / 8)
			return (-c + phase * s * 1.0i);
		if(theta > 2. / 8)
			return (-s + phase * c * 1.0i);
		if(theta > 1. / 8)
			return (s + phase * c * 1.0i);
		return (c + phase * s * 1.0i);
	};
}
template<typename ft_vec_it_t, typename vec_iter_t>
inline void fft_in_situ(ft_vec_it_t VecIter,
	vec_iter_t IterAt1, vec_iter_t IterAtNeg1) {
	using ft_image_t = ft_vec_it_t::value_type;
	const size_t VecLen
		= 2 * std::distance(IterAt1, IterAtNeg1);
	size_t sub_ft_size = 1;
	size_t num_sub_ft = VecLen / sub_ft_size;
	size_t num_sub_ft_pair = num_sub_ft / 2;
	for(size_t sub_ft_it = 0;
		sub_ft_it < num_sub_ft_pair; sub_ft_it += 2) {
		const ft_image_t parit00 = VecIter[sub_ft_it];
		const ft_image_t parit01
			= VecIter[num_sub_ft_pair + sub_ft_it];
		const ft_image_t parit10
			= VecIter[sub_ft_it + 1];
		const ft_image_t parit11
			= VecIter[num_sub_ft_pair + sub_ft_it + 1];
		VecIter[sub_ft_it] = parit00 + parit01;
		VecIter[sub_ft_it + 1] = parit00 - parit01;
		VecIter[num_sub_ft_pair + sub_ft_it]
			= parit10 + parit11;
		VecIter[num_sub_ft_pair + sub_ft_it + 1]
			= parit10 - parit11;
	}
	for(sub_ft_size *= 2, num_sub_ft /= 2,
		num_sub_ft_pair /= 2;
		sub_ft_size < num_sub_ft_pair;
		sub_ft_size *= 2, num_sub_ft /= 2,
		num_sub_ft_pair /= 2) {
		for(size_t perm_it = 0; perm_it < VecLen;
			perm_it += 2 * num_sub_ft_pair) {
			for(size_t sub_ft_it = perm_it;
				sub_ft_it < perm_it + num_sub_ft_pair;
				sub_ft_it += 2 * sub_ft_size) {
				for(size_t i = sub_ft_it, nth_exp = 0;
					i < sub_ft_it + sub_ft_size;
					i++, nth_exp += num_sub_ft_pair) {
					const ft_image_t parit00
						= VecIter[i];
					const ft_image_t parit01
						= VecIter[num_sub_ft_pair + i]
						* IterAt1[nth_exp];
					const ft_image_t parit10
						= VecIter[i + sub_ft_size];
					const ft_image_t parit11
						= VecIter[num_sub_ft_pair + i
							  + sub_ft_size]
						* IterAt1[nth_exp];
					VecIter[i] = parit00 + parit01;
					VecIter[i + sub_ft_size]
						= parit00 - parit01;
					VecIter[num_sub_ft_pair + i]
						= parit10 + parit11;
					VecIter[num_sub_ft_pair + i
						+ sub_ft_size]
						= parit10 - parit11;
				}
			}
		}
	}
	for(; sub_ft_size < VecLen; sub_ft_size *= 2,
		num_sub_ft /= 2, num_sub_ft_pair /= 2) {
		for(size_t sub_ft_it = 0; sub_ft_it < VecLen;
			sub_ft_it += 2 * sub_ft_size) {
			for(size_t i = sub_ft_it, nth_exp = 0;
				i < sub_ft_it + sub_ft_size;
				i++, nth_exp += num_sub_ft_pair) {
				const ft_image_t parit1
					= IterAt1[nth_exp]
					* VecIter[i + sub_ft_size];
				const ft_image_t parit0 = VecIter[i];
				VecIter[i] = parit0 + parit1;
				VecIter[i + sub_ft_size]
					= parit0 - parit1;
			}
		}
	}
}
template<typename input_iter_t, typename A_t,
	typename vec_iter_t,
	typename vec_t
	= std::vector<vec_iter_t::value_type>>
inline vec_t fft(input_iter_t Series, A_t Scale,
	vec_iter_t IterAt1, vec_iter_t IterAtNeg1) {
	const size_t VecLen
		= 2 * std::distance(IterAt1, IterAtNeg1);
	vec_t FtVec(VecLen);
	std::transform(Series, Series + VecLen,
		FtVec.begin(),
		[Scale](auto x) { return Scale * x; });
	fft_in_situ(FtVec.begin(), IterAt1, IterAtNeg1);
	return FtVec;
}
} // namespace fft7071