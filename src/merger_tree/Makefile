hd= /Users/tinker


CC= gcc
CLINK= gcc 
CFLAGS=  $(FFTW_CFLAGS)
CLIB= -L${HOME}/cosmo/lib -lcutil -lm 
XDIR= $(HOME)/exec


OBJS1=	main.o ran1.o gasdev.o zbrent.o erfinv.o get_redshift.o
merger_tree:	$(OBJS1)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1) $(CLIB)
	cp $@ $(HOME)/exec/

clean:
	rm -f *.o
