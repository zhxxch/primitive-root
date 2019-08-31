#include<stdio.h>
#include<stdint.h>
#include<math.h>
#include"GF_M31.h"
#include"GF_p.h"
void print_quadratic(const GFp2_t q,
	Z64_ct p, const char delimiter[]){
	printf("x^2+%llix+%lli\t/ " "x^2-%llix-%lli%s",
		q.coeff.a1, q.coeff.a0,
		p - q.coeff.a1, p - q.coeff.a0,
		delimiter);
}
void print_GFp2_elem(
	const GFp2_t e, const char delimiter[]){
	printf("%llix+%lli%s",
		e.coeff.a1, e.coeff.a0, delimiter);
}
void print_GFp2_Z(const uint32_t Re, const uint32_t Im,
	const uint32_t p, const char delimiter[]){
	const uint32_t nRe = p - Re, nIm = p - Im;
	const uint32_t print_Re = Re < nRe ? Re : nRe;
	const uint32_t print_Im = Im < nIm ? Im : nIm;
	const char *signRe = Re < nRe ? "" : "-";
	const char *signIm = Im < nIm ? "+" : "-";
	printf("(%s0x%x%s0x%xi)%s",
		"", Re, "+", Im, delimiter);
}
int main(int argc, char* argv[]){
	const GF_M31_2_t id_elem = {1,0};
	const GF_M31_2_t poly_x = {0,1};
	const GF_M31_2_t poly_0 = {1,1};
	const GF_M31_t M31 = (1u << 31) - 1;
	GFp2_t ppoly_ans;
	GF_M31_2_t ppoly_ans_GFM31
		= iter_GF_M31_2_ppoly_test(poly_0);
	ppoly_ans.coeff.a0 = ppoly_ans_GFM31.coeff.a0;
	ppoly_ans.coeff.a1 = ppoly_ans_GFM31.coeff.a1;
	print_quadratic(ppoly_ans, M31, "\n");
	const GF_M31_2_t proot2 = {.Z.Re = 8,.Z.Im = 3};
	const uint64_t NTT_N_M31
		= (uint64_t)M31*(uint64_t)M31 - 1;
	for(int i = 0; i < 9; i++){
		printf("w_M31^(N/%i) = ", (1 << i));
		const GF_M31_2_t rs_complex = pow_GF_M31_2_Z(
			proot2, NTT_N_M31 / (1ull << i));
		print_GFp2_Z(rs_complex.Z.Re,
			rs_complex.Z.Im, M31, "\n");
	}
	for(int k = 0; k <= 256; k++){
		printf("w_M31^(%i*N/256) = ", k);
		const GF_M31_2_t rs_complex = pow_GF_M31_2_Z(
			proot2, k*(NTT_N_M31 / 256));
		print_GFp2_Z(rs_complex.Z.Re,
			rs_complex.Z.Im, M31, "\n");
	}
	for(int i = 3; i <= 27; i = i * i){
		printf("w_M31^(N/%i) = ", i);
		const GF_M31_2_t rs_complex = pow_GF_M31_2_Z(
			proot2, NTT_N_M31 / i);
		print_GFp2_Z(rs_complex.Z.Re,
			rs_complex.Z.Im, M31, "\n");
	}
	for(int k = 0; k <= 3; k++){
		printf("w_M31^(%i*N/3) = ", k);
		const GF_M31_2_t rs_complex = pow_GF_M31_2_Z(
			proot2, k*(NTT_N_M31 / 3));
		print_GFp2_Z(rs_complex.Z.Re,
			rs_complex.Z.Im, M31, "\n");
	}
#if 0
	Z64_ct F3_30 = 3 * (1u << 30) + 1;
	Z64_ct F3_30_Factors[] = {2,3};
	Z64_ct F3_30_2_Factors[] = {2,3,79,20387503};
	Z64_t F3_30_test_F[2], F3_30_2_test_F[4];
	init_test_exps(2, F3_30 - 1,
		F3_30_Factors, F3_30_test_F);
	init_test_exps(4, F3_30*F3_30 - 1,
		F3_30_2_Factors, F3_30_2_test_F);
	GFp2_t ppoly_F3_30_min = {1,1};
	for(int i = 2, progress_ctr = 2;
		i < (F3_30 - 1) / 2; i += 2){
		ppoly_F3_30_min.coeff.a1 = i;
		for(GFp2_t ppoly_F3_30_ans = ppoly_F3_30_min;
			ppoly_F3_30_ans.coeff.a0 < (F3_30 - 1) / 2; ){
			ppoly_F3_30_ans
				= iter_ppoly_test(ppoly_F3_30_ans, F3_30,
					2, 4, F3_30_test_F, F3_30_2_test_F);
			Z64_ct a = (ppoly_F3_30_ans.coeff.a1 / 2);
			Z64_ct b_sq = ppoly_F3_30_ans.coeff.a0 - a * a;
			Z64_ct b_sqrt = (Z64_t)floor(sqrt((double)b_sq));
			if(b_sq < ppoly_F3_30_ans.coeff.a0
				&& b_sq == b_sqrt * b_sqrt){
				print_quadratic(ppoly_F3_30_ans, F3_30, " / ");
				printf("(%lli+%llii)\n", a, b_sqrt);
			}
			ppoly_F3_30_ans.coeff.a0
				= ppoly_F3_30_ans.coeff.a0 + 1;
		}
		progress_ctr += 2;
		if(progress_ctr == (F3_30 - 1) / (2 * 1024)){
			progress_ctr = 0;
			printf("%lli/1024\n", i / ((F3_30 - 1) / (2 * 1024)));
		}
	}
#endif
	return(0);
}