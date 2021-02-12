Targets=primeq.exe

All: $(Targets)

CFLAGS=/D_CRT_SECURE_NO_WARNINGS=1 /O2 /nologo /GL /W4 /MD /I./ /openmp /arch:AVX2 /std:c17 /Zi /Qvec-report:1

CXXFLAGS=/D_CRT_SECURE_NO_WARNINGS=1 /O2 /nologo /GL /W4 /MD /EHsc /I./ /fp:fast /arch:AVX2 /Zc:__cplusplus /openmp /openmp:experimental /Zi

pf4.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

p3013.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

p3122.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

adc.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

primeq.exe: $*.c $*.h
	$(CC) $*.c $(CFLAGS) /Fe:$@
