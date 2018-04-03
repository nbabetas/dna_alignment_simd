CC=			gcc
CFLAGS=		-g -Wall -Wno-unused-function -O2 -mavx2
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
LOBJS=		utils.o malloc_wrap.o 

PROG=		ksw
INCLUDES=	
LIBS=		-lz 
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ksw:libksw.a $(AOBJS) ksw.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) ksw.o -o $@ -L. -lksw $(LIBS)


libksw.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a



# DO NOT DELETE THIS LINE -- make depend depends on it.

ksw.o: ksw.h malloc_wrap.h kseq.h
malloc_wrap.o: malloc_wrap.h
utils.o: utils.h ksort.h malloc_wrap.h kseq.h
