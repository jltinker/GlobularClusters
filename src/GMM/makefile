#
#  Makefile for GMM code
#
#  Written by Oleg Gnedin on February 19, 2010
#
#-------------------  macro definitions  -------------------------------
FC       = f77
CC       = c++
OPTS     = -O
FFLAGS   = -O
CFLAGS   = -O
LIB      = -lm
LDFLAGS  = -lcc_dynamic
#---------------------  implicit rules  --------------------------------
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ----- GMM -----
OBJ = gmm.o

gmm: $(OBJ) makefile
	$(CC) $(OPTS) -o $@ $(OBJ) $(LIB)

$(OBJ): makefile


# ----- Dip -----
OBJ2 = dip.o

dip: $(OBJ2) makefile
	$(FC) $(OPTS) -o $@ $(OBJ2) $(LIB)

$(OBJ2): makefile
