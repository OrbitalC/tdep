include Makefile.inc
CODE = self_energy
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)selfenergy.o\
$(OBJECT_PATH)thermal_conductivity.o\
$(OBJECT_PATH)density_of_state.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHmpi) $(incLPATHhdf)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHmpi) $(incIPATHhdf)
LIBS = -lolle -lflap $(blaslapackLIBS) $(incLIBSmpi) $(incLIBShdf)

#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o:\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)selfenergy.o\
$(OBJECT_PATH)thermal_conductivity.o\
$(OBJECT_PATH)density_of_state.o
	$(F90) $(OPT) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)density_of_state.o: $(OBJECT_PATH)selfenergy.o
	$(F90) $(OPT) $(F90FLAGS) -c density_of_state.f90 $(LIBS) -o $@
$(OBJECT_PATH)thermal_conductivity.o: $(OBJECT_PATH)selfenergy.o
	$(F90) $(OPT) $(F90FLAGS) -c thermal_conductivity.f90 $(LIBS) -o $@
$(OBJECT_PATH)selfenergy.o:
	$(F90) $(OPT) $(F90FLAGS) -c selfenergy.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(OPT) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
