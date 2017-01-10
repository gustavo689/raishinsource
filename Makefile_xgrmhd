NAME = xgrmhd
#FC = gfortran
FC = mpif90
#FC = ifort

FFLAGS = -O2 -Wall
#FFLAGS = -O2 -lmpi

OBJS = pram.o main.o bnd.o calcha.o calflx.o calsf.o caluu.o cdtcfl.o correction.o coord.o ct.o func.o grid.o hll.o mdgrmhd.o mpisub.o physfunc.o rec.o recov.o restar.o rkt.o
LD = $(FC)
LFLAGS = 
#LFLAGS = -lmpi
LIBS = 
COMMON_MOD = pram.f90 

.SUFFIXES:.o .f90

$(NAME):$(OBJS)
	$(LD) $(LFLAGS) $(OBJS) -o $(NAME).exe $(LIBS)

.f90.o:
	$(FC) $(FFLAGS) -c $<


clean:
	rm -f core *.o *.exe *.mod

$(OBJS): $(COMMON_MOD)