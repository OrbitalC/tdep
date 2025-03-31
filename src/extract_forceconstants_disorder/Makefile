include Makefile.inc
CODE = extract_forceconstants_disorder
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)type_forceconstant_disorder_thirdorder.o\
$(OBJECT_PATH)fit_helper.o

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

$(OBJECT_PATH)main.o: \
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)type_forceconstant_disorder_thirdorder.o\
$(OBJECT_PATH)fit_helper.o
	$(F90) $(OPT) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_forceconstant_disorder_thirdorder.o:
	$(F90) $(OPT) $(F90FLAGS) -c type_forceconstant_disorder_thirdorder.f90 $(LIBS) -o $@
$(OBJECT_PATH)fit_helper.o:\
	fit_helper.f90\
	$(OBJECT_PATH)type_forceconstant_disorder_thirdorder.o
	$(F90) $(OPT) $(F90FLAGS) -c fit_helper.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(OPT) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
