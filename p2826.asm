; Listing generated by Microsoft (R) Optimizing Compiler Version 19.28.29333.0 

include listing.inc

INCLUDELIB OLDNAMES

PUBLIC	??_C@_02EAMLCBHB@Eq@				; `string'
PUBLIC	??_C@_03NNEIEOJB@Neq@				; `string'
PUBLIC	??_C@_0BP@LGDPFEOC@?$CFs?3?70x?$CF?48x?0?50x?$CF?48x?5?$CIans?0?5res?$CJ?6@ ; `string'
PUBLIC	??_C@_0P@DMHBINJL@Generator?3?5?$CFi?6@		; `string'
PUBLIC	p2826_orders
PUBLIC	p2826_factors
PUBLIC	p2826
EXTRN	__imp___stdio_common_vfprintf:PROC
EXTRN	__imp___acrt_iob_func:PROC
EXTRN	__security_check_cookie:PROC
p2826_orders DQ	000000008a000000H
	DQ	000000005c000000H
	DQ	000000000c000000H
p2826_factors DQ 0000000000000002H
	DQ	0000000000000003H
	DQ	0000000000000017H
p2826	DQ	0000000114000001H
CONST	ENDS
;	COMDAT ??_C@_0P@DMHBINJL@Generator?3?5?$CFi?6@
CONST	SEGMENT
??_C@_0P@DMHBINJL@Generator?3?5?$CFi?6@ DB 'Generator: %i', 0aH, 00H ; `string'
CONST	ENDS
;	COMDAT ??_C@_0BP@LGDPFEOC@?$CFs?3?70x?$CF?48x?0?50x?$CF?48x?5?$CIans?0?5res?$CJ?6@
CONST	SEGMENT
??_C@_0BP@LGDPFEOC@?$CFs?3?70x?$CF?48x?0?50x?$CF?48x?5?$CIans?0?5res?$CJ?6@ DB '%'
	DB	's:', 09H, '0x%.8x, 0x%.8x (ans, res)', 0aH, 00H ; `string'
CONST	ENDS
;	COMDAT ??_C@_03NNEIEOJB@Neq@
CONST	SEGMENT
??_C@_03NNEIEOJB@Neq@ DB 'Neq', 00H			; `string'
CONST	ENDS
;	COMDAT ??_C@_02EAMLCBHB@Eq@
CONST	SEGMENT
??_C@_02EAMLCBHB@Eq@ DB 'Eq', 00H			; `string'
CONST	ENDS
PUBLIC	main
PUBLIC	fill_list_rand
PUBLIC	search_p2826_gens
PUBLIC	is_p2826_gen
PUBLIC	pow_Zp_ring
PUBLIC	test_modp2826
PUBLIC	mod_p2826
PUBLIC	speck64u96
PUBLIC	printf
PUBLIC	_vfprintf_l
PUBLIC	__local_stdio_printf_options
COMM	?_OptionsStorage@?1??__local_stdio_printf_options@@9@9:QWORD							; `__local_stdio_printf_options'::`2'::_OptionsStorage
_DATA	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$main DD	imagerel $LN36
	DD	imagerel $LN36+205
	DD	imagerel $unwind$main
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$fill_list_rand DD imagerel $LN21
	DD	imagerel $LN21+157
	DD	imagerel $unwind$fill_list_rand
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$search_p2826_gens DD imagerel $LN34
	DD	imagerel $LN34+205
	DD	imagerel $unwind$search_p2826_gens
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$is_p2826_gen DD imagerel $LN23
	DD	imagerel $LN23+146
	DD	imagerel $unwind$is_p2826_gen
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$test_modp2826 DD imagerel $LN28
	DD	imagerel $LN28+18
	DD	imagerel $unwind$test_modp2826
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$0$test_modp2826 DD imagerel $LN28+18
	DD	imagerel $LN28+333
	DD	imagerel $chain$0$test_modp2826
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$1$test_modp2826 DD imagerel $LN28+333
	DD	imagerel $LN28+362
	DD	imagerel $chain$1$test_modp2826
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$printf DD imagerel $LN6
	DD	imagerel $LN6+85
	DD	imagerel $unwind$printf
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$_vfprintf_l DD imagerel $LN4
	DD	imagerel $LN4+68
	DD	imagerel $unwind$_vfprintf_l
pdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$_vfprintf_l DD 060f01H
	DD	09640fH
	DD	08340fH
	DD	0700b520fH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$printf DD 041b01H
	DD	07017521bH
	DD	030156016H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$chain$1$test_modp2826 DD 021H
	DD	imagerel $LN28
	DD	imagerel $LN28+18
	DD	imagerel $unwind$test_modp2826
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$chain$0$test_modp2826 DD 020521H
	DD	043405H
	DD	imagerel $LN28
	DD	imagerel $LN28+18
	DD	imagerel $unwind$test_modp2826
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$test_modp2826 DD 010401H
	DD	04204H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$is_p2826_gen DD 040a01H
	DD	02740aH
	DD	013405H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$search_p2826_gens DD 0a1a01H
	DD	09741aH
	DD	08641aH
	DD	07541aH
	DD	06341aH
	DD	0e016321aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$fill_list_rand DD 041201H
	DD	027412H
	DD	01340dH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$main DD	0a1a01H
	DD	09741aH
	DD	08641aH
	DD	07541aH
	DD	06341aH
	DD	0e016321aH
xdata	ENDS
; Function compile flags: /Ogtpy
; File C:\Program Files (x86)\Windows Kits\10\include\10.0.18362.0\ucrt\corecrt_stdio_config.h
;	COMDAT __local_stdio_printf_options
_TEXT	SEGMENT
__local_stdio_printf_options PROC			; COMDAT

; 87   :         static unsigned __int64 _OptionsStorage;
; 88   :         return &_OptionsStorage;

	lea	rax, OFFSET FLAT:?_OptionsStorage@?1??__local_stdio_printf_options@@9@9 ; `__local_stdio_printf_options'::`2'::_OptionsStorage

; 89   :     }

	ret	0
__local_stdio_printf_options ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Program Files (x86)\Windows Kits\10\include\10.0.18362.0\ucrt\stdio.h
;	COMDAT _vfprintf_l
_TEXT	SEGMENT
_Stream$ = 64
_Format$ = 72
_Locale$dead$ = 80
_ArgList$ = 88
_vfprintf_l PROC					; COMDAT

; 642  :     {

$LN4:
	mov	QWORD PTR [rsp+8], rbx
	mov	QWORD PTR [rsp+16], rsi
	push	rdi
	sub	rsp, 48					; 00000030H
	mov	rbx, r9
	mov	rdi, rdx
	mov	rsi, rcx

; 643  :         return __stdio_common_vfprintf(_CRT_INTERNAL_LOCAL_PRINTF_OPTIONS, _Stream, _Format, _Locale, _ArgList);

	call	__local_stdio_printf_options
	xor	r9d, r9d
	mov	QWORD PTR [rsp+32], rbx
	mov	r8, rdi
	mov	rdx, rsi
	mov	rcx, QWORD PTR [rax]
	call	QWORD PTR __imp___stdio_common_vfprintf

; 644  :     }

	mov	rbx, QWORD PTR [rsp+64]
	mov	rsi, QWORD PTR [rsp+72]
	add	rsp, 48					; 00000030H
	pop	rdi
	ret	0
_vfprintf_l ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Program Files (x86)\Windows Kits\10\include\10.0.18362.0\ucrt\stdio.h
;	COMDAT printf
_TEXT	SEGMENT
_Format$ = 80
printf	PROC						; COMDAT

; 954  :     {

$LN6:
	mov	QWORD PTR [rsp+8], rcx
	mov	QWORD PTR [rsp+16], rdx
	mov	QWORD PTR [rsp+24], r8
	mov	QWORD PTR [rsp+32], r9
	push	rbx
	push	rsi
	push	rdi
	sub	rsp, 48					; 00000030H
	mov	rdi, rcx

; 955  :         int _Result;
; 956  :         va_list _ArgList;
; 957  :         __crt_va_start(_ArgList, _Format);

	lea	rsi, QWORD PTR _Format$[rsp+8]

; 958  :         _Result = _vfprintf_l(stdout, _Format, NULL, _ArgList);

	mov	ecx, 1
	call	QWORD PTR __imp___acrt_iob_func
	mov	rbx, rax

; 643  :         return __stdio_common_vfprintf(_CRT_INTERNAL_LOCAL_PRINTF_OPTIONS, _Stream, _Format, _Locale, _ArgList);

	call	__local_stdio_printf_options
	xor	r9d, r9d
	mov	QWORD PTR [rsp+32], rsi
	mov	r8, rdi
	mov	rdx, rbx
	mov	rcx, QWORD PTR [rax]
	call	QWORD PTR __imp___stdio_common_vfprintf

; 959  :         __crt_va_end(_ArgList);
; 960  :         return _Result;
; 961  :     }

	add	rsp, 48					; 00000030H
	pop	rdi
	pop	rsi
	pop	rbx
	ret	0
printf	ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\speck.h
;	COMDAT speck64u96
_TEXT	SEGMENT
Pt$ = 8
KeyLo$ = 16
KeyHiLo$ = 24
KeyHiHi$ = 32
speck64u96 PROC						; COMDAT

; 11   : 	extern unsigned int _rotl(unsigned int, int);
; 12   : 	extern unsigned int _rotr(unsigned int, int);
; 13   : #define _SPECKXXR(x, y, k) \
; 14   : 	(y = _rotr(x ^ y, 3),  \
; 15   : 		x = _rotl(         \
; 16   : 			(unsigned long long)((x ^ (k)) - y), 8))
; 17   : #define _SPECKXX(Hi, Lo, k)          \
; 18   : 	(Hi = (_rotr(Hi, 8) + Lo) ^ (k), \
; 19   : 		Lo = _rotl(Lo, 3) ^ Hi)
; 20   : 	unsigned int CtHi = Pt >> 32,

	mov	rax, rcx
	shr	rax, 32					; 00000020H

; 21   : 				 CtLo = Pt & 0xffffffffu;
; 22   : 	for(int i = 0; i < 26;)

	xor	r10d, r10d
	npad	6
$LL2@speck64u96:

; 23   : 		(_SPECKXX(CtHi, CtLo, KeyLo),
; 24   : 			_SPECKXX(KeyHiLo, KeyLo, i++),
; 25   : 			_SPECKXX(CtHi, CtLo, KeyLo),

	ror	eax, 8
	add	eax, ecx
	ror	r8d, 8
	xor	eax, edx
	rol	ecx, 3
	add	r8d, edx
	ror	r9d, 8
	xor	ecx, eax
	rol	edx, 3
	xor	r8d, r10d
	ror	eax, 8
	xor	edx, r8d
	add	eax, ecx
	add	r9d, edx
	rol	ecx, 3
	xor	eax, edx
	inc	r10d
	xor	r9d, r10d
	rol	edx, 3
	xor	edx, r9d
	xor	ecx, eax
	inc	r10d
	cmp	r10d, 26
	jl	SHORT $LL2@speck64u96

; 26   : 			_SPECKXX(KeyHiHi, KeyLo, i++));
; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	shl	rax, 32					; 00000020H
	or	rax, rcx

; 28   : #undef _SPECKXXR
; 29   : #undef _SPECKXX
; 30   : }

	ret	0
speck64u96 ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT mod_p2826
_TEXT	SEGMENT
lo$ = 8
mod_p2826 PROC						; COMDAT

; 11   : 	return lo % 0x114000001ull;

	mov	rax, rcx
	xor	edx, edx
	mov	rcx, 4630511617				; 0000000114000001H
	div	rcx
	mov	rax, rdx

; 12   : }

	ret	0
mod_p2826 ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT test_modp2826
_TEXT	SEGMENT
Num$ = 48
Nonce$ = 56
test_modp2826 PROC					; COMDAT

; 13   : int test_modp2826(const int Num, const int Nonce) {

$LN28:
	sub	rsp, 40					; 00000028H
	movsxd	r8, edx

; 14   : 
; 15   : 	if(Num <= 0) return 0;

	test	ecx, ecx
	jle	$LN26@test_modp2

; 16   : 	const uint64_t a

	mov	rax, r8
	mov	QWORD PTR [rsp+32], rbx
; File C:\Users\user\Desktop\primitive-root\speck.h

; 20   : 	unsigned int CtHi = Pt >> 32,

	shr	rax, 32					; 00000020H
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 16   : 	const uint64_t a

	lea	r10d, DWORD PTR [r8+5]
; File C:\Users\user\Desktop\primitive-root\speck.h

; 22   : 	for(int i = 0; i < 26;)

	xor	r11d, r11d
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 16   : 	const uint64_t a

	lea	ecx, DWORD PTR [r8+3]
; File C:\Users\user\Desktop\primitive-root\speck.h

; 22   : 	for(int i = 0; i < 26;)

	mov	ebx, r11d
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 16   : 	const uint64_t a

	mov	r9d, 99					; 00000063H
; File C:\Users\user\Desktop\primitive-root\speck.h

; 21   : 				 CtLo = Pt & 0xffffffffu;

	mov	edx, r8d
$LL6@test_modp2:

; 23   : 		(_SPECKXX(CtHi, CtLo, KeyLo),
; 24   : 			_SPECKXX(KeyHiLo, KeyLo, i++),
; 25   : 			_SPECKXX(CtHi, CtLo, KeyLo),

	ror	eax, 8
	add	eax, edx
	ror	r10d, 8
	xor	eax, ecx
	rol	edx, 3
	add	r10d, ecx
	ror	r9d, 8
	xor	edx, eax
	rol	ecx, 3
	xor	r10d, ebx
	ror	eax, 8
	xor	ecx, r10d
	add	eax, edx
	add	r9d, ecx
	rol	edx, 3
	xor	eax, ecx
	inc	ebx
	xor	r9d, ebx
	rol	ecx, 3
	xor	ecx, r9d
	xor	edx, eax
	inc	ebx
	cmp	ebx, 26
	jl	SHORT $LL6@test_modp2

; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	mov	r10d, eax
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 19   : 	const uint64_t b

	lea	ecx, DWORD PTR [r8+16]
; File C:\Users\user\Desktop\primitive-root\speck.h

; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	mov	eax, edx
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 16   : 	const uint64_t a

	mov	rbx, 9223301660057661447		; 7fffbffe2001f007H
; File C:\Users\user\Desktop\primitive-root\speck.h

; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	shl	r10, 32					; 00000020H
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 19   : 	const uint64_t b

	mov	r9d, 111				; 0000006fH
; File C:\Users\user\Desktop\primitive-root\speck.h

; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	or	r10, rax
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 16   : 	const uint64_t a

	mov	rax, rbx
	mul	r10
	shr	rdx, 29
	imul	rax, rdx, 1073750017			; 40002001H

; 19   : 	const uint64_t b

	lea	edx, DWORD PTR [r8+2]
	sub	r10, rax
	movsxd	rax, edx
; File C:\Users\user\Desktop\primitive-root\speck.h

; 20   : 	unsigned int CtHi = Pt >> 32,

	shr	rax, 32					; 00000020H
$LL11@test_modp2:

; 23   : 		(_SPECKXX(CtHi, CtLo, KeyLo),
; 24   : 			_SPECKXX(KeyHiLo, KeyLo, i++),
; 25   : 			_SPECKXX(CtHi, CtLo, KeyLo),

	ror	eax, 8
	add	eax, edx
	ror	r9d, 8
	xor	eax, ecx
	rol	edx, 3
	add	r9d, ecx
	ror	r8d, 8
	xor	edx, eax
	rol	ecx, 3
	xor	r9d, r11d
	ror	eax, 8
	xor	ecx, r9d
	add	eax, edx
	add	r8d, ecx
	rol	edx, 3
	xor	eax, ecx
	inc	r11d
	xor	r8d, r11d
	rol	ecx, 3
	xor	ecx, r8d
	xor	edx, eax
	inc	r11d
	cmp	r11d, 26
	jl	SHORT $LL11@test_modp2

; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	mov	ecx, edx
	mov	r9d, eax
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 19   : 	const uint64_t b

	mov	rax, rbx
; File C:\Users\user\Desktop\primitive-root\speck.h

; 27   : 	return ((((unsigned long long)CtHi) << 32) | CtLo);

	shl	r9, 32					; 00000020H
	or	r9, rcx
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 19   : 	const uint64_t b

	mul	r9

; 20   : 		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
; 21   : 		% 0x40002001;
; 22   : 	const int32_t ans = (a * b) % 0x40002001;

	mov	rax, rbx

; 11   : 	return lo % 0x114000001ull;

	mov	rbx, QWORD PTR [rsp+32]

; 19   : 	const uint64_t b

	shr	rdx, 29
	imul	rcx, rdx, 1073750017			; 40002001H
	sub	r9, rcx

; 20   : 		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
; 21   : 		% 0x40002001;
; 22   : 	const int32_t ans = (a * b) % 0x40002001;

	imul	r9, r10
	mul	r9
	mov	r8, r9

; 11   : 	return lo % 0x114000001ull;

	mov	rax, r9

; 20   : 		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
; 21   : 		% 0x40002001;
; 22   : 	const int32_t ans = (a * b) % 0x40002001;

	shr	rdx, 29
	imul	rcx, rdx, 1073750017			; 40002001H

; 11   : 	return lo % 0x114000001ull;

	xor	edx, edx

; 20   : 		= speck64u96(Nonce + 2, Nonce + 16, 111, Nonce)
; 21   : 		% 0x40002001;
; 22   : 	const int32_t ans = (a * b) % 0x40002001;

	sub	r8, rcx

; 11   : 	return lo % 0x114000001ull;

	mov	rcx, 4630511617				; 0000000114000001H
	div	rcx

; 23   : 	const int32_t res = mod_p2826(a * b);
; 24   : 	if(ans != res)

	cmp	r8d, edx
	je	SHORT $LN26@test_modp2

; 25   : 		printf("%s:\t0x%.8x, 0x%.8x (ans, res)\n",

	mov	r9d, edx
	lea	rcx, OFFSET FLAT:??_C@_0BP@LGDPFEOC@?$CFs?3?70x?$CF?48x?0?50x?$CF?48x?5?$CIans?0?5res?$CJ?6@
	lea	rdx, OFFSET FLAT:??_C@_03NNEIEOJB@Neq@
	call	printf
$LN26@test_modp2:

; 26   : 			(ans == res ? "Eq" : "Neq"), ans, res);
; 27   : 	return 0;
; 28   : }

	xor	eax, eax
	add	rsp, 40					; 00000028H
	ret	0
test_modp2826 ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT pow_Zp_ring
_TEXT	SEGMENT
p$dead$ = 8
base$ = 16
exponent$ = 24
pow_Zp_ring PROC					; COMDAT

; 37   : 	const uint64_t base, const uint64_t exponent) {

	mov	r9, rdx

; 38   : 	uint64_t ans = 1;

	mov	ecx, 1

; 39   : 	for(uint64_t exp = exponent, g = base; exp;) {

	test	r8, r8
	je	SHORT $LN11@pow_Zp_rin
	mov	r10, 4630511617				; 0000000114000001H
$LL2@pow_Zp_rin:

; 40   : 		if(exp & 1) { ans = (ans * g) % p; }

	test	r8b, 1
	je	SHORT $LN5@pow_Zp_rin
	mov	rax, r9
	xor	edx, edx
	imul	rax, rcx
	div	r10
	mov	rcx, rdx
$LN5@pow_Zp_rin:

; 41   : 		exp = (exp >> 1);
; 42   : 		g = (g * g) % p;

	imul	r9, r9
	xor	edx, edx
	mov	rax, r9
	div	r10
	shr	r8, 1
	mov	r9, rdx
	jne	SHORT $LL2@pow_Zp_rin
$LN11@pow_Zp_rin:

; 43   : 	}
; 44   : 	return ans;
; 45   : }

	mov	rax, rcx
	ret	0
pow_Zp_ring ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT is_p2826_gen
_TEXT	SEGMENT
a$ = 8
is_p2826_gen PROC					; COMDAT

; 46   : bool is_p2826_gen(const int32_t a) {

$LN23:
	mov	QWORD PTR [rsp+8], rbx
	mov	QWORD PTR [rsp+16], rdi
	movsxd	r11, ecx
	lea	r10, OFFSET FLAT:p2826_orders
	mov	rbx, 4630511617				; 0000000114000001H
	lea	rdi, OFFSET FLAT:p2826_orders+24
	npad	11
$LL4@is_p2826_g:

; 39   : 	for(uint64_t exp = exponent, g = base; exp;) {

	mov	rcx, QWORD PTR [r10]
	mov	r9d, 1
	mov	r8, r11
	test	rcx, rcx
	je	SHORT $LN15@is_p2826_g
$LL8@is_p2826_g:

; 40   : 		if(exp & 1) { ans = (ans * g) % p; }

	test	cl, 1
	je	SHORT $LN11@is_p2826_g
	mov	rax, r8
	xor	edx, edx
	imul	rax, r9
	div	rbx
	mov	r9, rdx
$LN11@is_p2826_g:

; 41   : 		exp = (exp >> 1);
; 42   : 		g = (g * g) % p;

	imul	r8, r8
	xor	edx, edx
	mov	rax, r8
	div	rbx
	shr	rcx, 1
	mov	r8, rdx
	jne	SHORT $LL8@is_p2826_g

; 48   : 		const uint64_t r = pow_Zp_ring(
; 49   : 			0x114000001ull, a, p2826_orders[i]);
; 50   : 		if(r <= 1) return false;

	cmp	r9, 1
	jbe	SHORT $LN15@is_p2826_g

; 47   : 	for(int i = 0; i < 3; i++) {

	add	r10, 8
	cmp	r10, rdi
	jl	SHORT $LL4@is_p2826_g

; 51   : 	}
; 52   : 	return true;

	mov	al, 1

; 53   : }

	mov	rbx, QWORD PTR [rsp+8]
	mov	rdi, QWORD PTR [rsp+16]
	ret	0
$LN15@is_p2826_g:
	mov	rbx, QWORD PTR [rsp+8]
	xor	al, al
	mov	rdi, QWORD PTR [rsp+16]
	ret	0
is_p2826_gen ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT search_p2826_gens
_TEXT	SEGMENT
Num$dead$ = 48
search_p2826_gens PROC					; COMDAT

; 54   : int search_p2826_gens(const int Num) {

$LN34:
	mov	QWORD PTR [rsp+8], rbx
	mov	QWORD PTR [rsp+16], rbp
	mov	QWORD PTR [rsp+24], rsi
	mov	QWORD PTR [rsp+32], rdi
	push	r14
	sub	rsp, 32					; 00000020H

; 55   : 	for(int i = 0; i < Num; i++) {

	xor	edi, edi
	lea	r14, OFFSET FLAT:p2826_orders
	mov	ebx, edi
	lea	rbp, OFFSET FLAT:p2826_orders+24
	mov	rsi, 4630511617				; 0000000114000001H
	npad	10
$LL4@search_p28:

; 47   : 	for(int i = 0; i < 3; i++) {

	mov	r10, r14
	npad	13
$LL10@search_p28:

; 39   : 	for(uint64_t exp = exponent, g = base; exp;) {

	mov	rcx, QWORD PTR [r10]
	mov	r9d, 1
	mov	r8, rbx
	test	rcx, rcx
	je	SHORT $LN2@search_p28
$LL14@search_p28:

; 40   : 		if(exp & 1) { ans = (ans * g) % p; }

	test	cl, 1
	je	SHORT $LN17@search_p28
	mov	rax, r8
	xor	edx, edx
	imul	rax, r9
	div	rsi
	mov	r9, rdx
$LN17@search_p28:

; 41   : 		exp = (exp >> 1);
; 42   : 		g = (g * g) % p;

	imul	r8, r8
	xor	edx, edx
	mov	rax, r8
	div	rsi
	shr	rcx, 1
	mov	r8, rdx
	jne	SHORT $LL14@search_p28

; 50   : 		if(r <= 1) return false;

	cmp	r9, 1
	jbe	SHORT $LN2@search_p28

; 43   : 	}
; 44   : 	return ans;
; 45   : }
; 46   : bool is_p2826_gen(const int32_t a) {
; 47   : 	for(int i = 0; i < 3; i++) {

	add	r10, 8
	cmp	r10, rbp
	jl	SHORT $LL10@search_p28

; 56   : 		if(is_p2826_gen(i)) {
; 57   : 			printf("Generator: %i\n", i);

	mov	edx, edi
	lea	rcx, OFFSET FLAT:??_C@_0P@DMHBINJL@Generator?3?5?$CFi?6@
	call	printf
$LN2@search_p28:

; 55   : 	for(int i = 0; i < Num; i++) {

	inc	edi
	inc	rbx
	cmp	edi, 10
	jl	SHORT $LL4@search_p28

; 58   : 		}
; 59   : 	}
; 60   : 	return 0;
; 61   : }

	mov	rbx, QWORD PTR [rsp+48]
	xor	eax, eax
	mov	rbp, QWORD PTR [rsp+56]
	mov	rsi, QWORD PTR [rsp+64]
	mov	rdi, QWORD PTR [rsp+72]
	add	rsp, 32					; 00000020H
	pop	r14
	ret	0
search_p2826_gens ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
; File C:\Users\user\Desktop\primitive-root\speck.h
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT fill_list_rand
_TEXT	SEGMENT
Num$ = 8
List$ = 16
fill_list_rand PROC					; COMDAT

; 63   : 	const int Num, int *restrict List) {

$LN21:

; 64   : 	for(int i = 0; i < Num; i++) {

	test	ecx, ecx
	jle	$LN19@fill_list_
	mov	QWORD PTR [rsp+8], rbx
	mov	QWORD PTR [rsp+16], rdi

; 63   : 	const int Num, int *restrict List) {

	mov	rbx, rdx
	mov	edi, ecx

; 64   : 	for(int i = 0; i < Num; i++) {

	xor	r10d, r10d
	npad	6
$LL4@fill_list_:
; File C:\Users\user\Desktop\primitive-root\speck.h

; 20   : 	unsigned int CtHi = Pt >> 32,

	xor	eax, eax
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 65   : 		List[i] = (int)speck64u96(100, 200, i, i + 2);

	lea	r8d, DWORD PTR [r10+2]
; File C:\Users\user\Desktop\primitive-root\speck.h

; 22   : 	for(int i = 0; i < 26;)

	xor	r11d, r11d
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 65   : 		List[i] = (int)speck64u96(100, 200, i, i + 2);

	mov	r9d, r10d
	mov	edx, 200				; 000000c8H
; File C:\Users\user\Desktop\primitive-root\speck.h

; 21   : 				 CtLo = Pt & 0xffffffffu;

	mov	ecx, 100				; 00000064H
	npad	10
$LL7@fill_list_:

; 23   : 		(_SPECKXX(CtHi, CtLo, KeyLo),
; 24   : 			_SPECKXX(KeyHiLo, KeyLo, i++),
; 25   : 			_SPECKXX(CtHi, CtLo, KeyLo),

	ror	eax, 8
	add	eax, ecx
	ror	r9d, 8
	xor	eax, edx
	rol	ecx, 3
	add	r9d, edx
	ror	r8d, 8
	xor	ecx, eax
	rol	edx, 3
	xor	r9d, r11d
	ror	eax, 8
	xor	edx, r9d
	add	eax, ecx
	add	r8d, edx
	rol	ecx, 3
	xor	eax, edx
	inc	r11d
	xor	r8d, r11d
	rol	edx, 3
	xor	edx, r8d
	xor	ecx, eax
	inc	r11d
	cmp	r11d, 26
	jl	SHORT $LL7@fill_list_
; File C:\Users\user\Desktop\primitive-root\p2826.c

; 64   : 	for(int i = 0; i < Num; i++) {

	mov	DWORD PTR [rbx], ecx
	inc	r10d
	add	rbx, 4
	cmp	r10d, edi
	jl	SHORT $LL4@fill_list_

; 66   : 	}
; 67   : }

	mov	rbx, QWORD PTR [rsp+8]
	mov	rdi, QWORD PTR [rsp+16]
$LN19@fill_list_:
	ret	0
fill_list_rand ENDP
_TEXT	ENDS
; Function compile flags: /Ogtpy
; File C:\Users\user\Desktop\primitive-root\p2826.c
;	COMDAT main
_TEXT	SEGMENT
main	PROC						; COMDAT

; 68   : int main(void) { return search_p2826_gens(10); }

$LN36:
	mov	QWORD PTR [rsp+8], rbx
	mov	QWORD PTR [rsp+16], rbp
	mov	QWORD PTR [rsp+24], rsi
	mov	QWORD PTR [rsp+32], rdi
	push	r14
	sub	rsp, 32					; 00000020H

; 55   : 	for(int i = 0; i < Num; i++) {

	xor	edi, edi
	lea	r14, OFFSET FLAT:p2826_orders
	mov	ebx, edi
	lea	rbp, OFFSET FLAT:p2826_orders+24
	mov	rsi, 4630511617				; 0000000114000001H
	npad	10
$LL6@main:

; 47   : 	for(int i = 0; i < 3; i++) {

	mov	r10, r14
	npad	13
$LL12@main:

; 39   : 	for(uint64_t exp = exponent, g = base; exp;) {

	mov	rcx, QWORD PTR [r10]
	mov	r9d, 1
	mov	r8, rbx
	test	rcx, rcx
	je	SHORT $LN4@main
$LL16@main:

; 40   : 		if(exp & 1) { ans = (ans * g) % p; }

	test	cl, 1
	je	SHORT $LN19@main
	mov	rax, r8
	xor	edx, edx
	imul	rax, r9
	div	rsi
	mov	r9, rdx
$LN19@main:

; 41   : 		exp = (exp >> 1);
; 42   : 		g = (g * g) % p;

	imul	r8, r8
	xor	edx, edx
	mov	rax, r8
	div	rsi
	shr	rcx, 1
	mov	r8, rdx
	jne	SHORT $LL16@main

; 50   : 		if(r <= 1) return false;

	cmp	r9, 1
	jbe	SHORT $LN4@main

; 43   : 	}
; 44   : 	return ans;
; 45   : }
; 46   : bool is_p2826_gen(const int32_t a) {
; 47   : 	for(int i = 0; i < 3; i++) {

	add	r10, 8
	cmp	r10, rbp
	jl	SHORT $LL12@main

; 57   : 			printf("Generator: %i\n", i);

	mov	edx, edi
	lea	rcx, OFFSET FLAT:??_C@_0P@DMHBINJL@Generator?3?5?$CFi?6@
	call	printf
$LN4@main:

; 51   : 	}
; 52   : 	return true;
; 53   : }
; 54   : int search_p2826_gens(const int Num) {
; 55   : 	for(int i = 0; i < Num; i++) {

	inc	edi
	inc	rbx
	cmp	edi, 10
	jl	SHORT $LL6@main

; 68   : int main(void) { return search_p2826_gens(10); }

	mov	rbx, QWORD PTR [rsp+48]
	xor	eax, eax
	mov	rbp, QWORD PTR [rsp+56]
	mov	rsi, QWORD PTR [rsp+64]
	mov	rdi, QWORD PTR [rsp+72]
	add	rsp, 32					; 00000020H
	pop	r14
	ret	0
main	ENDP
_TEXT	ENDS
END
