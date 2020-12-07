/*
sse_f64_complex.hpp (2020) copyright, github.com/zhxxch
*/
#pragma once
#include <array>
#include <complex>
#include <immintrin.h>
namespace fft7071 {
class sse_f64_complex {
	std::array<double, 2> v;

  public:
	using value_type = double;
	sse_f64_complex() {
		v[0] = 0;
		v[1] = 0;
	};
	sse_f64_complex(const double re, const double im) {
		v[0] = re;
		v[1] = im;
	};
	sse_f64_complex(const std::complex<double> &z) {
		v[0] = z.real();
		v[1] = z.imag();
	};
	sse_f64_complex conj() const {
		return sse_f64_complex(v[0], -v[1]);
	};
	operator std::complex<double>() const {
		return std::complex<double>(v[0], v[1]);
	};
	template<typename A_t>
	sse_f64_complex &operator=(const A_t &other) {
		v[0] = other;
		return *this;
	}
	template<>
	sse_f64_complex &operator=<sse_f64_complex>(
		const sse_f64_complex &other) {
		v = other.v;
		return *this;
	}
	template<>
	sse_f64_complex &operator=<std::complex<double>>(
		const std::complex<double> &other) {
		v[0] = other.real();
		v[1] = other.imag();
		return *this;
	}
	sse_f64_complex operator+(
		const sse_f64_complex &other) const {
		sse_f64_complex ans;
		__m128d lhs = _mm_loadu_pd(v.data());
		__m128d rhs = _mm_loadu_pd(other.v.data());
		_mm_storeu_pd(
			ans.v.data(), _mm_add_pd(lhs, rhs));
		return ans;
	};
	sse_f64_complex operator-(
		const sse_f64_complex &other) const {
		sse_f64_complex ans;
		const __m128d lhs = _mm_loadu_pd(v.data());
		const __m128d rhs
			= _mm_loadu_pd(other.v.data());
		_mm_storeu_pd(
			ans.v.data(), _mm_sub_pd(lhs, rhs));
		return ans;
	};
	sse_f64_complex operator*(
		const sse_f64_complex &other) const {
		sse_f64_complex ans;
		const __m128d lhs = _mm_loadu_pd(v.data());
		const __m128d rhs
			= _mm_loadu_pd(other.v.data());
		const __m128d rhs_rere
			= _mm_unpacklo_pd(rhs, rhs);
		const __m128d rhs_imim
			= _mm_unpackhi_pd(rhs, rhs);
		const __m128d lhs_refl
			= _mm_shuffle_pd(lhs, lhs, 1);
		_mm_storeu_pd(ans.v.data(),
			_mm_fmaddsub_pd(lhs, rhs_rere,
				_mm_mul_pd(lhs_refl, rhs_imim)));
		return ans;
	};
};
} // namespace fft7071