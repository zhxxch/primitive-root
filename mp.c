#include<stdio.h>
#include<stdint.h>
#include<math.h>
#include"GF_Mp31.h"
#include"GF_p.h"

void print_quadratic(const GFp2_t q,
	Z64_ct p, const char delimiter[]){
	printf("x^2+%llix+%lli\t/\t" "x^2-%llix-%lli%s",
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
	printf("(%s%i%s%ii)%s",
		signRe, print_Re, signIm, print_Im, delimiter);
}
int main(int argc, char* argv[]){
	const GF_Mp31_2_t id_elem = {1,0};
	const GF_Mp31_2_t poly_x = {0,1};
	const GF_Mp31_2_t poly_0 = {1,1};
	GFp2_t ppoly_ans;
	GF_Mp31_2_t ppoly_ans_GFMp31
		= iter_GF_Mp31_2_ppoly_test(poly_0);
	ppoly_ans.coeff.a0 = ppoly_ans_GFMp31.coeff.a0;
	ppoly_ans.coeff.a1 = ppoly_ans_GFMp31.coeff.a1;
	print_quadratic(ppoly_ans, Mp31, "\n");
	const GF_Mp31_2_t proot2 = {.Z.Re = 8,.Z.Im = 3};
	const uint64_t NTT_N
		= (uint64_t)Mp31*(uint64_t)Mp31-1;
	for(int i = 0; i < 10; i++){
		printf("x^(N/%i) = ", (1<<i));
		const GF_Mp31_2_t rs_complex = pow_GF_Mp31_2_Z(
			proot2, NTT_N / (1ull << i));
			print_GFp2_Z(rs_complex.Z.Re,
				rs_complex.Z.Im, Mp31, "\n");
	}
	return(0);
}