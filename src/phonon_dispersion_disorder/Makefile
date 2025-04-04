include Makefile.inc
CODE = phonon_dispersion_disorder
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)type_phonons_disorder.o\
$(OBJECT_PATH)type_structurefactor_dispersions.o\
$(OBJECT_PATH)type_projected_phonons.o\
$(OBJECT_PATH)type_projected_phonons_bandstructure.o\
$(OBJECT_PATH)type_forceconstant_disorder_secondorder_dynamicalmatrix.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHmpi) $(incLPATHhdf)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHmpi) $(incIPATHhdf)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBSmpi) $(incLIBShdf)

F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o: \
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)type_phonons_disorder.o\
$(OBJECT_PATH)type_structurefactor_dispersions.o\
$(OBJECT_PATH)type_projected_phonons.o\
$(OBJECT_PATH)type_projected_phonons_bandstructure.o\
$(OBJECT_PATH)type_forceconstant_disorder_secondorder_dynamicalmatrix.o
	$(F90) $(OPT) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_phonons_disorder.o:\
	$(OBJECT_PATH)type_forceconstant_disorder_secondorder_dynamicalmatrix.o
	$(F90) $(OPT) $(F90FLAGS) -c type_phonons_disorder.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_structurefactor_dispersions.o:\
	$(OBJECT_PATH)type_forceconstant_disorder_secondorder_dynamicalmatrix.o\
	$(OBJECT_PATH)type_phonons_disorder.o
	$(F90) $(OPT) $(F90FLAGS) -c type_structurefactor_dispersions.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_projected_phonons.o:\
	$(OBJECT_PATH)type_forceconstant_disorder_secondorder_dynamicalmatrix.o\
	$(OBJECT_PATH)type_phonons_disorder.o
	$(F90) $(OPT) $(F90FLAGS) -c type_projected_phonons.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_projected_phonons_bandstructure.o:\
	$(OBJECT_PATH)type_projected_phonons.o
	$(F90) $(OPT) $(F90FLAGS) -c type_projected_phonons_bandstructure.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_forceconstant_disorder_secondorder_dynamicalmatrix.o:
	$(F90) $(OPT) $(F90FLAGS) -c type_forceconstant_disorder_secondorder_dynamicalmatrix.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(OPT) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
