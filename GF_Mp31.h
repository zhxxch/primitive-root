#pragma once
#include<stdint.h>
typedef uint32_t GF_Mp31_t;
const GF_Mp31_t Mp31 = ((1u << 31) - 1);
const GF_Mp31_t TotientMp31 = ((1u << 31) - 2);
const GF_Mp31_t GF_Mp31_Factor[] = {
	2,3,7,11,31,151,331,
	/* 2 * 3^2 * 7 * 11 * 31 * 151 * 331
	= 2147483646 = (2^31 - 2) = Mp31 - 1
	(8 prime factors, 7 distinct)
	Mp31^2-1 = (2^31-1)^1-1 = 2^31*(Mp31-1)
	*/
};
const GF_Mp31_t GF_Mp31_FactorInv[] = {
#define MP31 ((1ull << 31) - 1)
	(MP31 - 1) / 2,
	(MP31 - 1) / 3,
	(MP31 - 1) / 7,
	(MP31 - 1) / 11,
	(MP31 - 1) / 31,
	(MP31 - 1) / 151,
	(MP31 - 1) / 331,
#undef MP31
};
const uint64_t GF_Mp31_2_FactorInv[] = {
#define MP31 ((1ull << 31) - 1)
#define N2M1(x) (x*x-1)
	N2M1(MP31) / 2,
	N2M1(MP31) / 3,
	N2M1(MP31) / 7,
	N2M1(MP31) / 11,
	N2M1(MP31) / 31,
	N2M1(MP31) / 151,
	N2M1(MP31) / 331,
#undef N2M1
#undef MP31
};

typedef union{
	GF_Mp31_t vec[2];
	struct{ GF_Mp31_t a0; GF_Mp31_t a1; }coeff;
	struct{ GF_Mp31_t Re; GF_Mp31_t Im; }Z;
}GF_Mp31_2_t;

GF_Mp31_t add_GF_Mp31(
	const GF_Mp31_t a, const GF_Mp31_t b){
	const GF_Mp31_t res = a + b;
	const GF_Mp31_t m = (res&((1u << 31) - 1))
		+ (res >> 31);
	return((m^((1u<<31)-1))?m:0);
}
GF_Mp31_t sub_GF_Mp31(
	const GF_Mp31_t a, const GF_Mp31_t b){
	const GF_Mp31_t res = a + ((1u << 31) - 1 - b);
	const GF_Mp31_t m = (res&((1u << 31) - 1))
		+ (res >> 31);
	return((m^((1u<<31)-1))?m:0);
}
GF_Mp31_t mul_GF_Mp31(
	const GF_Mp31_t a, const GF_Mp31_t b){
	const uint64_t x = (uint64_t)a * (uint64_t)b;
	const uint64_t x1
		= ((((1ull << 31) - 1) & x) << 31) + x;
	const uint64_t x2
		= (((1ull << 62) & x1) >> 31) + x1;
	const GF_Mp31_t x3 = (GF_Mp31_t)
		(((1ull << 31) - 1) & (x2 >> 31));
	return((x3^((1u<<31)-1))?x3:0);
}
GF_Mp31_t pow_GF_Mp31(
	const GF_Mp31_t base, const GF_Mp31_t exponent){
	GF_Mp31_t ans = 1;
	for(GF_Mp31_t exp = exponent, g = base; exp;){
		if(exp & 1){
			ans = mul_GF_Mp31(ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_Mp31(g, g);
	}
	return ans;
}
GF_Mp31_2_t mod_GF_Mp31_2_quad(const GF_Mp31_t a2,
	const GF_Mp31_2_t a, const GF_Mp31_2_t q_mod){
	GF_Mp31_2_t ans;
	ans.coeff.a0 = sub_GF_Mp31(a.coeff.a0,
		mul_GF_Mp31(q_mod.coeff.a0, a2));
	ans.coeff.a1 = sub_GF_Mp31(a.coeff.a1,
		mul_GF_Mp31(q_mod.coeff.a1, a2));
	return ans;
}

GF_Mp31_2_t mul_GF_Mp31_2_quad(const GF_Mp31_2_t q_mod,
	const GF_Mp31_2_t a, const GF_Mp31_2_t b){
	const GF_Mp31_t a2 = mul_GF_Mp31(
		a.coeff.a1, b.coeff.a1);
	GF_Mp31_2_t imm = {0};
	imm.coeff.a1 = add_GF_Mp31(
		mul_GF_Mp31(a.coeff.a1, b.coeff.a0),
		mul_GF_Mp31(a.coeff.a0, b.coeff.a1)
	);
	imm.coeff.a0 = mul_GF_Mp31(
		a.coeff.a0, b.coeff.a0);
	return mod_GF_Mp31_2_quad(a2, imm, q_mod);
}

GF_Mp31_2_t pow_GF_Mp31_2_quad(const GF_Mp31_2_t q_mod,
	const GF_Mp31_2_t base, const uint64_t exponent){
	GF_Mp31_2_t ans = {.coeff.a0 = 1,.coeff.a1 = 0};
	GF_Mp31_2_t g = base;
	for(uint64_t exp = exponent; exp;){
		if(exp & 1){
			ans = mul_GF_Mp31_2_quad(q_mod, ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_Mp31_2_quad(q_mod, g, g);
	}
	return ans;
}

int is_proot_GF_Mp31(const GF_Mp31_t x){
	for(int i = 0; i < 7; i++){
		const GF_Mp31_t ans = pow_GF_Mp31(
			x, GF_Mp31_FactorInv[i]);
		if(ans == 1)return 0;
	}
	return 1;
}
GF_Mp31_t JacobiS_GF_Mp31(const GF_Mp31_t x){
	return pow_GF_Mp31(x, ((1u << 31) - 1) / 2);
}
int is_realroot_GF_Mp31_2_quad(const GF_Mp31_2_t q){
	GF_Mp31_t a1 = q.coeff.a1;
	GF_Mp31_t Delta = sub_GF_Mp31(mul_GF_Mp31(a1, a1),
		mul_GF_Mp31(4, q.coeff.a0));
	GF_Mp31_t jacobi = JacobiS_GF_Mp31(Delta);
	return(jacobi == 1);
}
int is_ppoly_GF_Mp31_2(const GF_Mp31_2_t ppoly){
	const GF_Mp31_2_t x = {0, 1};
	if(is_realroot_GF_Mp31_2_quad(ppoly)){
		return 0;
	}
	const GF_Mp31_2_t id_elem
		= pow_GF_Mp31_2_quad(ppoly, x,
		((1ull << 31) - 1)*((1ull << 31) - 1) - 1);
	if(id_elem.coeff.a0 != 1 || id_elem.coeff.a1 != 0){
		return 0;
	}
	for(int i = 0; i < 7; i++){
		const GF_Mp31_2_t ans = pow_GF_Mp31_2_quad(
			ppoly, x, GF_Mp31_2_FactorInv[i]);
		if(ans.coeff.a0 == 1 && ans.coeff.a1 == 0){
			return 0;
		}
	}
	return 1;
}
GF_Mp31_2_t iter_GF_Mp31_2_ppoly_test(
	const GF_Mp31_2_t min){
	for(GF_Mp31_2_t x = min; x.coeff.a0 < 10;
		x.coeff.a0++){
		if(!is_proot_GF_Mp31(x.coeff.a0)){
			continue;
		}
		for(x.coeff.a1 = min.coeff.a1;
			x.coeff.a1 < 10; x.coeff.a1++){
			if(is_ppoly_GF_Mp31_2(x)){
				return x;
			}
		}
	}
	return (GF_Mp31_2_t){ 0, 0 };
}
GF_Mp31_2_t mul_GF_Mp31_2_Z(
	const GF_Mp31_2_t a, const GF_Mp31_2_t b){
	GF_Mp31_2_t ans = {0};
	ans.Z.Re = sub_GF_Mp31(
		mul_GF_Mp31(a.Z.Re, b.Z.Re),
		mul_GF_Mp31(a.Z.Im, b.Z.Im));
	ans.Z.Im = add_GF_Mp31(
		mul_GF_Mp31(a.Z.Re, b.Z.Im),
		mul_GF_Mp31(a.Z.Im, b.Z.Re));
	return ans;
}
GF_Mp31_2_t add_GF_Mp31_2_Z(
	const GF_Mp31_2_t a, const GF_Mp31_2_t b){
	GF_Mp31_2_t ans = {0};
	ans.Z.Re = add_GF_Mp31(a.Z.Re, b.Z.Re);
	ans.Z.Im = add_GF_Mp31(a.Z.Im, b.Z.Im);
	return ans;
}
GF_Mp31_2_t sub_GF_Mp31_2_Z(
	const GF_Mp31_2_t a, const GF_Mp31_2_t b){
	GF_Mp31_2_t ans = {0};
	ans.Z.Re = sub_GF_Mp31(a.Z.Re, b.Z.Re);
	ans.Z.Im = sub_GF_Mp31(a.Z.Im, b.Z.Im);
	return ans;
}
GF_Mp31_2_t pow_GF_Mp31_2_Z(
	const GF_Mp31_2_t base, const uint64_t exponent){
	GF_Mp31_2_t ans = {.Z.Re = 1,.Z.Im = 0};
	GF_Mp31_2_t g = base;
	for(uint64_t exp = exponent; exp;){
		if(exp & 1){
			ans = mul_GF_Mp31_2_Z(ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_Mp31_2_Z(g, g);
	}
	return ans;
}