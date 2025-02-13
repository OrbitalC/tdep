include Makefile.inc
CODE=anharmonic_free_energy
PROG=../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)lo_thermodynamic_helpers.o $(OBJECT_PATH)options.o $(OBJECT_PATH)thirdorder.o $(OBJECT_PATH)fourthorder.o $(OBJECT_PATH)epot.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi)

#OPT = -O0 -fbacktrace -fcheck=all -Wall
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(LIBS)

$(OBJECT_PATH)main.o: main.f90 $(OBJECT_PATH)options.o $(OBJECT_PATH)lo_thermodynamic_helpers.o $(OBJECT_PATH)thirdorder.o $(OBJECT_PATH)fourthorder.o $(OBJECT_PATH)epot.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)thirdorder.o: $(OBJECT_PATH)options.o $(OBJECT_PATH)lo_thermodynamic_helpers.o
	$(F90) $(F90FLAGS) -c thirdorder.f90 $(LIBS) -o $@
$(OBJECT_PATH)fourthorder.o: $(OBJECT_PATH)options.o
	$(F90) $(F90FLAGS) -c fourthorder.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
$(OBJECT_PATH)epot.o: epot.f90 $(OBJECT_PATH)lo_thermodynamic_helpers.o
	$(F90) $(F90FLAGS) -c epot.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_thermodynamic_helpers.o: lo_thermodynamic_helpers.f90
	$(F90) $(F90FLAGS) -c lo_thermodynamic_helpers.f90 $(LIBS) -o $@
