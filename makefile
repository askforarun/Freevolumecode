CC=gcc
CFLAGS=-lm
widom:	
	$(CC)-o \
widom \
widom.c \
getnatoms.c \
getljcoeff.c \
gettimestep.c \
getboxdim.c \
Atom.c \
cubecentres.c \
water.c -$(CFLAGS)

