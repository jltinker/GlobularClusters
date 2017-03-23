hd = $(HOME)/cosmo/lib
LIB = -lm -L${hd} -lcutil


CC = gcc
CFLAGS = -O2


OBJS01 = gc_evolve.o qromo.o polint.o midpnt.o 
gc_evolve: $(OBJS01)
	$(CC) -o $@ $(OBJS01) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS02 = gc_population.o qromo.o polint.o midpnt.o zbrent.o ran1.o gasdev.o hdr.o input_params.o
gc_population: $(OBJS02)
	$(CC) -o $@ $(OBJS02) $(LIB)
	cp -f $@ $(HOME)/exec/$@

