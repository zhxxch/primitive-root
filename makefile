Targets=p3013.exe adc.exe pf4.exe

All: $(Targets)

CFLAGS=/D_CRT_SECURE_NO_WARNINGS=1 /O2 /nologo /GL /W4 /MD /I./ /Zi /Qvec-report:1

pf4.exe: $*.c
	$(CC) $** $(CFLAGS) /FAs /Fe:$@

p3013.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@

adc.exe: $*.c
	$(CC) $** $(CFLAGS) /Fe:$@