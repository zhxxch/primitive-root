#pragma once
#include<stdint.h>

typedef union{
	uint32_t vec[4];
	struct{ uint32_t w0, w1, w2, w3; }w;
}GF_M107_t;
const GF_M107_t GF_M107_MUL_Id = {.w={1,0,0,0}};
const GF_M107_t GF_M107_ADD_Id1 = {0};
const GF_M107_t GF_M107_ADD_Id2 = {.w={0x0fffffffu,0x0fffffffu,0x0fffffffu,((1u<<23)-1)}};
GF_M107_t GF_M107_MUL_Id(void){
	/* TODO */
	return{0};
}
const uint64_t GF_M107_Factors[] = {
	/*
	2*3*107*6361*69431*20394401*28059810762433
	= (2^107-1)-1
	2^108*3*107*6361*69431*20394401*28059810762433
	= (2^107-1)^2-1
	*/
	2,3,107,6361,69431,20394401ull,28059810762433ull
};
const GF_M107_t GF_M107_FactorInv[] = {
	{0}/* TODO */
};
GF_M107_t add_GF_M107(
	const GF_M107_t a, const GF_M107_t b){
	GF_M107_t vec_ans;
	int carry = 0;
	for(int i = 0; i < 3; i++){
		vec_ans.vec[i] = a.vec[i] + b.vec[i] + carry;
		carry = vec_ans.vec[i] >> 28;
		vec_ans.vec[i] &= 0x0fffffffu;
	}
	vec_ans.vec[0] += vec_ans.vec[4] >> 23;
	vec_ans.vec[4] &= ((1u << 23) - 1);
	return vec_ans;
}
GF_M107_t add_inv_GF_M107(const GF_M107_t a){
	GF_M107_t vec_ans;
	for(int i = 0; i < 4; i++){
		vec_ans.vec[i] = 0x0fffffffu - vec_ans.vec[i];
	}
	vec_ans.vec[4] &= ((1u << 23) - 1);
	return vec_ans;
}
GF_M107_t sub_GF_M107(
	const GF_M107_t a, const GF_M107_t b){
	return add_GF_M107(a, add_inv_GF_M107(b));
}
GF_M107_t mul_GF_M107(
	const GF_M107_t a, const GF_M107_t b){
	uint64_t res[4][4], semisum[4] = {0};
	for(int ai = 0; ai < 4; ai++){
		for(int bi = 0; bi < 4; bi++){
			const uint64_t num_res
				= (uint64_t)a.vec[ai]
				* (uint64_t)b.vec[bi];
			res[ai][bi] = res;
		}
	}
	semisum[0] = (res[0][0] & 0x0fffffffu)
		+ (res[4][4] >> 51);
	semisum[3] = (res[4][4] >> 23)
		+ (res[2][3] >> 51) + (res[3][2] >> 51);
	semisum[1] = (res[0][1] >> 28)
		+ (res[1][0] & 0x0ffffffu);
	semisum[2] = (res[0][1] >> 28) + (res[1][0] >> 28)
		+ ((res[2][3] >> 23) & 0x0fffffffu)
		+ ((res[2][1] >> 23) & 0x0fffffffu)
		+ (res[3][3] & ((1u << 23) - 1));
	for(int ai = 0; ai < 3; ai++){
		const uint64_t r = res[ai][2 - ai];
		const uint64_t r2 = r >> 28;
		const uint64_t l = res[ai + 1][3 - ai];
		const uint64_t l2 = l >> 23;
		semisum[2] += (r & 0x0fffffffu);
		semisum[3] += r2 & ((1u << 23) - 1));
		semisum[0] += r2 >> 23;
		semisum[0] += (l&((1u << 23) - 1)) << 5;
		semisum[1] += l2 & 0x0fffffffu;
		semisum[2] += l2 >> 28;
	}
	for(int ai = 0; ai < 4; ai++){
		const uint64_t r = res[ai][3 - ai];
		const uint64_t r2 = r >> 23;
		semisum[3] += (r & ((1u << 23) - 1));
		semisum[0] += r2 & 0x0fffffffu;
		semisum[1] += r2 >> 28;
	}
	uint64_t carry[4];
	for(int i = 0; i < 4; i++){
		carry[i] = semisum[i] >> 28;
		semisum[i] &= 0x0fffffffull;
	}
	for(int i = 0; i < 3; i++){
		semisum[i] += carry[i];
	}
	semisum[3] = ((1u << 23) - 1)
		&(semisum[3] + semisum[3] >> 23);
	GF_M107_t ans;
	for(int i = 0; i < 4; i++){
		ans.vec[i] = (GF_M107_t)semisum[i];
	}
}
GF_M107_t pow_GF_M107(
	const GF_M107_t base, const GF_M107_t exponent){
	GF_M107_t ans = {1,0,0,0}, g = base;
	for(GF_M107_t exp = exponent, g = base; exp;){
		for(int i = 0; i < 4; i++){
			for(int j = 1; j < (1u << 28); j = j << 1;){
				if(exp.vec[i] & j){
					ans = mul_GF_M107(ans, g);
				}
				g = mul_GF_M107(g, g);
			}
		}
	}
	return ans;
}
typedef union{
	uint32_t vecZ28[8];
	GF_M107_t vec[2];
	struct{ GF_M107_t a0; GF_M107_t a1; }coeff;
	struct{ GF_M107_t Re; GF_M107_t Im; }Z;
}GF_M107_2_t;
const GF_M107_2_t GF_M107_2_FactorInv[] = {
	{0}/* TODO */
};
GF_M107_2_t mod_GF_M107_2_quad(const GF_M107_t a2,
	const GF_M107_2_t a, const GF_M107_2_t q_mod){
	GF_M107_2_t ans;
	ans.coeff.a0 = sub_GF_M107(a.coeff.a0,
		mul_GF_M107(q_mod.coeff.a0, a2));
	ans.coeff.a1 = sub_GF_M107(a.coeff.a1,
		mul_GF_M107(q_mod.coeff.a1, a2));
	return ans;
}
GF_M107_2_t mul_GF_M107_2_quad(const GF_M107_2_t q_mod,
	const GF_M107_2_t a, const GF_M107_2_t b){
	const GF_M107_t a2 = mul_GF_M107(
		a.coeff.a1, b.coeff.a1);
	GF_M107_2_t imm = {0};
	imm.coeff.a1 = add_GF_M107(
		mul_GF_M107(a.coeff.a1, b.coeff.a0),
		mul_GF_M107(a.coeff.a0, b.coeff.a1));
	imm.coeff.a0 = mul_GF_M107(
		a.coeff.a0, b.coeff.a0);
	return mod_GF_M107_2_quad(a2, imm, q_mod);
}
GF_M107_2_t pow_GF_M107_2_quad(
	const GF_M107_2_t base, const GF_M107_2_t exponent){
	GF_M107_2_t ans = {0}, g = base;
	ans.coeff.a0.w.w0 = 1;
	for(GF_M107_2_t exp = exponent, g = base; exp;){
		for(int i = 0; i < 8; i++){
			for(int j = 1; j < (1u << 28); j = j << 1;){
				if(exp.vec[i] & j){
					ans = mul_GF_M107_2_quad(ans, g);
				}
				g = mul_GF_M107_2_quad(g, g);
			}
		}
	}
	return ans;
}
int is_mul_identity_GF_M107(const GF_M107_t x){
	return(x.w.w0 == 1
		&& x.w.w1 == 0 && x.w.w2 == 0 & x.w.w3 == 0);
}
int is_add_identity_GF_M107(const GF_M107_t x){
	int all0 = 1, all1 = 1;
	for(int i = 0; i < 4; i++){
		if(x.vec[i] != 0)all0 = 0;
	}
	int all1 = 1;
	for(int i = 0; i < 3; i++){
		if(x.vec[i] != 0x0fffffffu)all1 = 0;
	}
	if(x.w.w3 != ((1u << 23) - 1))all1 = 0;
	
	return(all0 || all1);
}
int is_mul_identity_GF_M107_2(const GF_M107_2_t x){
	return(is_mul_identity_GF_M107(x.Z.Re)
		&& is_add_identity_GF_M107(x.Z.Im));
}
int is_proot_GF_M107(const GF_M107_t x){
	for(int i = 0; i < 7; i++){
		const GF_M107_t ans = pow_GF_M107(
			x, GF_M107_FactorInv[i]);
		if(is_mul_identity_GF_M107(ans))return 0;
	}
	return 1;
}
GF_M107_t JacobiS_GF_M107(const GF_M107_t x){
	//TODO
	const GF_M107_t exp = {0/* (2^107-2)/2 */};
	return pow_GF_M107(x);
}
int is_realroot_GF_M107_2_quad(const GF_M107_2_t q){
	GF_M107_t a1 = q.coeff.a1;
	GF_M107_t Delta = sub_GF_M107(mul_GF_M107(a1, a1),
		mul_GF_M107(4, q.coeff.a0));
	GF_M107_t jacobi = JacobiS_GF_M107(Delta);
	return(is_mul_identity_GF_M107(jacobi);
}
int is_ppoly_GF_M107_2(const GF_M107_2_t ppoly){
	const GF_M107_2_t x = {0/* TODO */};
	if(is_realroot_GF_M107_2_quad(ppoly)){
		return 0;
	}
	const GF_M107_2_t id_elem
		= pow_GF_M107_2_quad(ppoly, x,
		0/* TODO: (2^107-1)-1 */);
	if(is_mul_identity_GF_M107_2(id_elem)){
		return 0;
	}
	for(int i = 0; i < 7; i++){
		const GF_M107_2_t ans = pow_GF_M107_2_quad(
			ppoly, x, GF_M107_2_FactorInv[i]);
		if(is_mul_identity_GF_M107_2(ans)){
			return 0;
		}
	}
	return 1;
}
GF_M107_2_t mul_GF_M107_2_Z(
	const GF_M107_2_t a, const GF_M107_2_t b){
	GF_M107_2_t ans = {0};
	ans.Z.Re = sub_GF_M107(
		mul_GF_M107(a.Z.Re, b.Z.Re),
		mul_GF_M107(a.Z.Im, b.Z.Im));
	ans.Z.Im = add_GF_M107(
		mul_GF_M107(a.Z.Re, b.Z.Im),
		mul_GF_M107(a.Z.Im, b.Z.Re));
	return ans;
}
GF_M107_2_t add_GF_M107_2_Z(
	const GF_M107_2_t a, const GF_M107_2_t b){
	GF_M107_2_t ans = {0};
	ans.Z.Re = add_GF_M107(a.Z.Re, b.Z.Re);
	ans.Z.Im = add_GF_M107(a.Z.Im, b.Z.Im);
	return ans;
}
GF_M107_2_t sub_GF_M107_2_Z(
	const GF_M107_2_t a, const GF_M107_2_t b){
	GF_M107_2_t ans = {0};
	ans.Z.Re = sub_GF_M107(a.Z.Re, b.Z.Re);
	ans.Z.Im = sub_GF_M107(a.Z.Im, b.Z.Im);
	return ans;
}
GF_M107_2_t pow_GF_M107_2_Z(
	const GF_M107_2_t base, const GF_M107_2_t exponent){
	/* TODO */

	/*
	GF_M107_2_t ans = {.Z.Re = 1,.Z.Im = 0};
	GF_M107_2_t g = base;
	for(uint64_t exp = exponent; exp;){
		if(exp & 1){
			ans = mul_GF_M107_2_Z(ans, g);
		}
		exp = (exp >> 1);
		g = mul_GF_M107_2_Z(g, g);
	}
	return ans;
	*/
}