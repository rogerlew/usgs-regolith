REG = regolith

OBJT90 = 
OBJT95 = regolini.o aspect.o dbisect.o derose.o gauss_aprx.o grid_count.o grid_size.o hsolve.o \
         laplacian.o lcsd_depth.o lrsc_depth.o nasd_depth.o ndsd_depth.o nsd_a_depth.o nsd_depth.o \
         regolith_main.o rgbanner.o ssizgrd.o wetness_ndx.o xyslope.o
OBJT77 = irdgrd.o isvgrd.o srdgrd.o ssvgrd.o


LIBS =

CC = gcc -O3
CCFLAGS = 
FC = f95 -w -O3 
FFLAGS =
F90 = f95 
F90FLAGS = -w -O0
LDFLAGS = 

all: $(REG)
#-----------------------------------------------------------------

$(REG): $(OBJT95) $(OBJT77) 
	$(F90) $(CCLIBS) $(LDFLAGS) -o $@ $(OBJT95) $(OBJT77) $(CCFLAGS) $(LIBS)

#-----------------------------------------------------------------
clean:
	rm -f $(REG)
	rm -rf $(OBJT95) $(OBJT77)
	rm -rf *.mod *.exe *.stackdump

.SUFFIXES: $(SUFFIXES) .f90 .f .c .f95

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<

.c.o:
	$(CC) $(CCINCLUDE) -c -w $<

.f95.o:
	$(F90) $(F90FLAGS) -c $<

