#pragma once
#include<stdint.h>
typedef uint32_t GF_M31_t;
const GF_M31_t TotientM31 = ((1u << 31) - 2);
const GF_M31_t GF_M31_Factor[] = {
	2,3,7,11,31,151,331,
	/* 2 * 3^2 * 7 * 11 * 31 * 151 * 331
	= 2147483646 = (2^31 - 2) = M31 - 1
	(8 prime factors, 7 distinct)
	M31^2-1 = (2^31-1)^1-1 = 2^31*(M31-1)
	*/
};
const GF_M31_t GF_M31_FactorInv[] = {
#define M31 ((1ull << 31) - 1)
	(M31 - 1) / 2,
	(M31 - 1) / 3,
	(M31 - 1) / 7,
	(M31 - 1) / 11,
	(M31 - 1) / 31,
	(M31 - 1) / 151,
	(M31 - 1) / 331,
#undef M31
};
const uint64_t GF_M31_2_FactorInv[] = {
#define M31 ((1ull << 31) - 1)
#define N2M1(x) (x*x-1)
	N2M1(M31) / 2,
	N2M1(M31) / 3,
	N2M1(M31) / 7,
	N2M1(M31) / 11,
	N2M1(M31) / 31,
	N2M1(M31) / 151,
	N2M1(M31) / 331,
#undef N2M1
#undef M31
};

typedef union{
	GF_M31_t vec[2];
	struct{ GF_M31_t a0; GF_M31_t a1; }coeff;
	struct{ GF_M31_t Re; GF_M31_t Im; }Z;
}GF_M31_2_t;
const GF_M31_2_t GF_M31_Primitive_Root1
= {.Z.Re = 8,.Z.Im = 3};
const GF_M31_2_t GF_M31_Primitive_Root2
= {.Z.Re = 7,.Z.Im = 2};
const GF_M31_2_t GF_M31_GaloisTf_W[] = {
	{1,0},{0x7ffffffeu,0},
{0,1},{ 0x8000u,0x8000u },
{0x233668e2u,0x45abdd8au},
{0x39aea997u,0x49fb5248u},
{0x1b389fb1u,0x5d739c92u},
{0x705f4834u,0x1858e449u},
{0x405c70e7u,0x0bc8aa05u}
};
GF_M31_t add_GF_M31(
	const GF_M31_t a, const GF_M31_t b){
	const GF_M31_t M31 = (1u << 31) - 1;
	const GF_M31_t res = a + b;
	const GF_M31_t res_m = res - M31;
	return(res_m < M31 ? res_m : res);
}
GF_M31_t sub_GF_M31(
	const GF_M31_t a, const GF_M31_t b){
	const GF_M31_t M31 = (1u << 31) - 1;
	return add_GF_M31(a, (M31 - b));
}
GF_M31_t mul_GF_M31(
	const GF_M31_t a, const GF_M31_t b){
	const GF_M31_t M31 = (1u << 31) - 1;
	const uint64_t x = (uint64_t)a * (uint64_t)b;
	const GF_M31_t lo = (GF_M31_t)(x & M31);
	const GF_M31_t hi = (GF_M31_t)(x >> 31);
	return add_GF_M31(lo, hi);
}
GF_M31_t pow_GF_M31(
	const GF_M31_t base, const GF_M31_t exponent){
	GF_M31_t ans = 1;
	for(GF_M31_t exp = exponent, g = base; exp;){
		if(exp & 1){
			ans = mul_GF_M31(ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_M31(g, g);
	}
	return ans;
}
GF_M31_2_t mod_GF_M31_2_quad(const GF_M31_t a2,
	const GF_M31_2_t a, const GF_M31_2_t q_mod){
	GF_M31_2_t ans;
	ans.coeff.a0 = sub_GF_M31(a.coeff.a0,
		mul_GF_M31(q_mod.coeff.a0, a2));
	ans.coeff.a1 = sub_GF_M31(a.coeff.a1,
		mul_GF_M31(q_mod.coeff.a1, a2));
	return ans;
}

GF_M31_2_t mul_GF_M31_2_quad(const GF_M31_2_t q_mod,
	const GF_M31_2_t a, const GF_M31_2_t b){
	const GF_M31_t a2 = mul_GF_M31(
		a.coeff.a1, b.coeff.a1);
	GF_M31_2_t imm = {0};
	imm.coeff.a1 = add_GF_M31(
		mul_GF_M31(a.coeff.a1, b.coeff.a0),
		mul_GF_M31(a.coeff.a0, b.coeff.a1)
	);
	imm.coeff.a0 = mul_GF_M31(
		a.coeff.a0, b.coeff.a0);
	return mod_GF_M31_2_quad(a2, imm, q_mod);
}

GF_M31_2_t pow_GF_M31_2_quad(const GF_M31_2_t q_mod,
	const GF_M31_2_t base, const uint64_t exponent){
	GF_M31_2_t ans = {.coeff.a0 = 1,.coeff.a1 = 0};
	GF_M31_2_t g = base;
	for(uint64_t exp = exponent; exp;){
		if(exp & 1){
			ans = mul_GF_M31_2_quad(q_mod, ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_M31_2_quad(q_mod, g, g);
	}
	return ans;
}

int is_proot_GF_M31(const GF_M31_t x){
	for(int i = 0; i < 7; i++){
		const GF_M31_t ans = pow_GF_M31(
			x, GF_M31_FactorInv[i]);
		if(ans == 1)return 0;
	}
	return 1;
}
GF_M31_t JacobiS_GF_M31(const GF_M31_t x){
	return pow_GF_M31(x, ((1u << 31) - 1) / 2);
}
int is_realroot_GF_M31_2_quad(const GF_M31_2_t q){
	GF_M31_t a1 = q.coeff.a1;
	GF_M31_t Delta = sub_GF_M31(mul_GF_M31(a1, a1),
		mul_GF_M31(4, q.coeff.a0));
	GF_M31_t jacobi = JacobiS_GF_M31(Delta);
	return(jacobi == 1);
}
int is_ppoly_GF_M31_2(const GF_M31_2_t ppoly){
	const GF_M31_2_t x = {0, 1};
	if(is_realroot_GF_M31_2_quad(ppoly)){
		return 0;
	}
	const GF_M31_2_t id_elem
		= pow_GF_M31_2_quad(ppoly, x,
		((1ull << 31) - 1)*((1ull << 31) - 1) - 1);
	if(id_elem.coeff.a0 != 1 || id_elem.coeff.a1 != 0){
		return 0;
	}
	for(int i = 0; i < 7; i++){
		const GF_M31_2_t ans = pow_GF_M31_2_quad(
			ppoly, x, GF_M31_2_FactorInv[i]);
		if(ans.coeff.a0 == 1 && ans.coeff.a1 == 0){
			return 0;
		}
	}
	return 1;
}
GF_M31_2_t iter_GF_M31_2_ppoly_test(
	const GF_M31_2_t min){
	for(GF_M31_2_t x = min; x.coeff.a0 < 10;
		x.coeff.a0++){
		if(!is_proot_GF_M31(x.coeff.a0)){
			continue;
		}
		for(x.coeff.a1 = min.coeff.a1;
			x.coeff.a1 < 10; x.coeff.a1++){
			if(is_ppoly_GF_M31_2(x)){
				return x;
			}
		}
	}
	return (GF_M31_2_t){ 0, 0 };
}
GF_M31_2_t mul_GF_M31_2_Z(
	const GF_M31_2_t a, const GF_M31_2_t b){
	GF_M31_2_t ans = {0};
	const GF_M31_t
		ar_ai = add_GF_M31(a.Z.Re, a.Z.Im),
		br_bi = add_GF_M31(b.Z.Re, b.Z.Im),
		ar_br = mul_GF_M31(a.Z.Re, b.Z.Re),
		ai_bi = mul_GF_M31(a.Z.Im, b.Z.Im),
		aa_bb = mul_GF_M31(ar_ai, br_bi),
		a_b = add_GF_M31(ar_br, ai_bi);
	ans.Z.Re = sub_GF_M31(ar_br, ai_bi);
	ans.Z.Im = sub_GF_M31(aa_bb, a_b);
	return ans;
}
GF_M31_2_t add_GF_M31_2_Z(
	const GF_M31_2_t a, const GF_M31_2_t b){
	GF_M31_2_t ans = {0};
	ans.Z.Re = add_GF_M31(a.Z.Re, b.Z.Re);
	ans.Z.Im = add_GF_M31(a.Z.Im, b.Z.Im);
	return ans;
}
GF_M31_2_t sub_GF_M31_2_Z(
	const GF_M31_2_t a, const GF_M31_2_t b){
	GF_M31_2_t ans = {0};
	ans.Z.Re = sub_GF_M31(a.Z.Re, b.Z.Re);
	ans.Z.Im = sub_GF_M31(a.Z.Im, b.Z.Im);
	return ans;
}
GF_M31_2_t pow_GF_M31_2_Z(
	const GF_M31_2_t base, const uint64_t exponent){
	GF_M31_2_t ans = {.Z.Re = 1,.Z.Im = 0};
	GF_M31_2_t g = base;
	for(uint64_t exp = exponent; exp;){
		if(exp & 1){
			ans = mul_GF_M31_2_Z(ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_M31_2_Z(g, g);
	}
	return ans;
}