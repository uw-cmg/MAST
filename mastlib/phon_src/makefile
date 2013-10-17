include ../make.inc
# this is needed on a Mac OS X
#LIB = /System/Library/Frameworks/vecLib.framework/vecLib
LIB     = -L../lapack -llapack -L../blas -lblas

SOURCE = nrtype.o nrutil.o diagonalize.o data.o phon.o cubicsym.o \
	hexsym.o equiva.o get_displ.o get_wig.o set_forces.o smooth.o \
	sgama.o sgama1.o numeric_force.o dyn_mat.o rshell.o build_supercell.o \
	set_lattice.o hpsort.o generate_punti.o  \
	reader.o drdatab.o dlexlib.o
# zheev.o dsyev.o

phon:  $(SOURCE)
	$(FC) $(FREE) $(OFLAG) $(SOURCE) $(LINK) $(LIB) -o phon

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) $(OFLAG) $(FREE) $(INCS) -c $*.f
