#pragma once
#include<stdint.h>
typedef uint64_t Z64_t;
typedef const Z64_t Z64_ct;
typedef union{
	Z64_t vec[2];
	struct{ Z64_t a0; Z64_t a1; }coeff;
}GFp2_t;
typedef union{
	uint32_t vec[2];
	struct{ uint32_t lo; uint32_t hi; }coeff;
}GF_Mp61_t;
Z64_ct Fp30 = 3 * (1ull << 30) + 1;

Z64_ct FactorFp30[2] = {
	2,3,
};
Z64_t test_exps_arr1[7];
Z64_t test_exps_arr2[7];

inline Z64_t add_Zp(Z64_ct p, Z64_ct a, Z64_ct b){
	return((a + b) % p);
}

inline Z64_t sub_Zp(Z64_ct p, Z64_ct a, Z64_ct b){
	return((a + (p - (b % p))) % p);
}

inline Z64_t mul_Zp(Z64_ct p, Z64_ct a, Z64_ct b){
	return((a * b) % p);
}

Z64_t pow_Zp_ring(Z64_ct p,
	Z64_ct base, Z64_ct exponent){
	Z64_t ans = 1;
	for(Z64_t exp = exponent, g = base; exp;){
		if(exp & 1){
			ans = (ans * g) % p;
		}
		exp = (exp >> 1);
		g = (g * g) % p;
	}
	return ans;
}

inline Z64_t JacobiS(
	Z64_ct a, Z64_ct p){
	return pow_Zp_ring(p, a, (p - 1) / 2);
}

inline Z64_t quadratic_eval(
	Z64_ct p, Z64_ct x, const GFp2_t q){
	return((((x*x) % p) + ((q.coeff.a1*x) % p) + q.coeff.a0) % p);
}

GFp2_t quadratic_mod(Z64_ct p, Z64_ct a2,
	const GFp2_t a, const GFp2_t q_mod){
	GFp2_t ans;
	ans.coeff.a0 = sub_Zp(p,
		a.coeff.a0, q_mod.coeff.a0*a2);
	ans.coeff.a1 = sub_Zp(p,
		a.coeff.a1, q_mod.coeff.a1*a2);
	return ans;
}

inline Z64_t mul_quadratic_a2(
	Z64_ct p, Z64_ct a1_a, Z64_ct a1_b){
	return mul_Zp(p, a1_a, a1_b);
}
inline Z64_t mul_quadratic_a1(Z64_ct p,
	Z64_ct a1_a, Z64_ct a0_a,
	Z64_ct a1_b, Z64_ct a0_b){
	return add_Zp(p,
		mul_Zp(p, a1_a, a0_b), mul_Zp(p, a0_a, a1_b));
}
inline Z64_t mul_quadratic_a0(
	Z64_ct p, Z64_ct a0_a, Z64_ct a0_b){
	return mul_Zp(p, a0_a, a0_b);
}

GFp2_t quadratic_mul(Z64_ct p, const GFp2_t q_mod,
	const GFp2_t a, const GFp2_t b){
	Z64_ct a2 = mul_quadratic_a2(p,
		a.coeff.a1, b.coeff.a1);
	GFp2_t imm;
	imm.coeff.a1 = mul_quadratic_a1(p,
		a.coeff.a1, a.coeff.a0, b.coeff.a1, b.coeff.a0);
	imm.coeff.a0 = mul_quadratic_a0(p,
		a.coeff.a0, b.coeff.a0);
	const GFp2_t ans = quadratic_mod(p, a2, imm, q_mod);
	return ans;
}

GFp2_t quadratic_pow(Z64_ct p, const GFp2_t q_mod,
	const GFp2_t base, Z64_ct exponent){
	GFp2_t ans = {.coeff.a0 = 1,.coeff.a1 = 0};
	GFp2_t g = base;
	for(Z64_t exp = exponent; exp;){
		if(exp & 1){
			ans = quadratic_mul(p, q_mod, ans, g);
		}
		exp = (exp >> 1);
		g = quadratic_mul(p, q_mod, g, g);
	}
	return ans;
}

Z64_ct* init_test_exps(
	const int num_factor, Z64_ct totient,
	Z64_ct factors[], Z64_t out_arr[]){
	for(int i = 0; i < num_factor; i++){
		out_arr[i] = totient / factors[i];
	}
	return out_arr;
}
Z64_t is_proot(Z64_ct x, Z64_ct order,
	Z64_ct num_factors, Z64_ct exps_arr[]){
	for(int i = 0; i < num_factors; i++){
		Z64_ct ans = pow_Zp_ring(order, x, exps_arr[i]);
		if(ans == 1)return 0;
	}
	return 1;
}
Z64_t iter_proot_test(Z64_ct min, Z64_ct order,
	Z64_ct num_factors, Z64_ct exps_arr[]){
	for(Z64_t x = min; x < order; x++){
		if(is_proot(x, order, num_factors, exps_arr)){
			return x;
		}
	}
	return 0;
}
int is_realroot_quadratic(Z64_ct p, const GFp2_t q){
	Z64_ct a1 = q.coeff.a1;
	Z64_ct Delta = sub_Zp(p, mul_Zp(p, a1, a1),
		mul_Zp(p, 4, q.coeff.a0));
	Z64_ct jacobi = JacobiS(Delta, p);
	return(jacobi == 1);
}
Z64_t is_ppolynomial(const GFp2_t ppoly, Z64_ct p,
	Z64_ct num_factors, Z64_ct exps_arr[]){
	const GFp2_t x = {0, 1};
	if(is_realroot_quadratic(p, ppoly)){
		return 0;
	}
	const GFp2_t id_elem = quadratic_pow(
		p, ppoly, x, p*p - 1);
	if(id_elem.coeff.a0 != 1 || id_elem.coeff.a1 != 0){
		return 0;
	}
	for(int i = 0; i < num_factors; i++){
		const GFp2_t ans = quadratic_pow(p, ppoly, x, exps_arr[i]);
		if(ans.coeff.a0 == 1 && ans.coeff.a1 == 0){
			return 0;
		}
	}
	return 1;
}
GFp2_t iter_ppoly_test(
	const GFp2_t min, Z64_ct p,
	Z64_ct num_factors_GFp, Z64_ct num_factors_GFp2,
	Z64_ct exps_GFp[], Z64_ct exps_GFp2[]){
	for(GFp2_t x = min; x.coeff.a1 < p; x.coeff.a1++){
		for(x.coeff.a0 = min.coeff.a0;
			x.coeff.a0 < p; x.coeff.a0++){
			if(is_proot(x.coeff.a0, p,
				num_factors_GFp, exps_GFp)){
				if(is_ppolynomial(
					x, p, num_factors_GFp2, exps_GFp2)){
					return x;
				}
			}
		}
	}
	return (GFp2_t){ 0, 0 };
}