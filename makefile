Targets=adc.exe pf4.exe p3122.exe bitrev.exe

All: $(Targets)

CFLAGS=/D_CRT_SECURE_NO_WARNINGS=1 /O2 /nologo /GL /W4 /MD /I./ /arch:AVX2 /Zi /Qvec-report:1

bitrev.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

pf4.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

p3013.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

p3122.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

adc.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@