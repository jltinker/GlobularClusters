hd = $(HOME)/cosmo/lib
LIB = -lm -L${hd} -lcutil


CC = gcc
CFLAGS = -O2


OBJS1 = bad_mcmc.o gasdev.o ran1.o gaussj.o jacobi.o
bad_mcmc:	$(OBJS1)
	$(CC) -o $@ $(OBJS1) $(LIB)
	cp -f $@ $(HOME)/exec/$@

