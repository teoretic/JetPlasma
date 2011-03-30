# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include mktemplate


.DEFAULT:
	-touch $@
all: plasma
alloc.o: ./alloc.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./alloc.f90
boundc.o: ./boundc.f90 names.o trot.o limiters.o flux.o initials.o gravity.o geometry.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./boundc.f90
flux.o: ./flux.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./flux.f90
for77test.o: ./for77test.f impl.inc
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -I.	./for77test.f
geometry.o: ./geometry.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./geometry.f90
gravity.o: ./gravity.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./gravity.f90
initials.o: ./initials.f90 names.o trot.o radiationTR.o solution.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./initials.f90
limiters.o: ./limiters.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./limiters.f90
mesh.o: ./mesh.f90 names.o multies.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./mesh.f90
multies.o: ./multies.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./multies.f90
names.o: ./names.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./names.f90
plasma.o: ./plasma.f90 alloc.o trot.o time.o mesh.o solution.o flux.o boundc.o names.o multies.o solvers.o geometry.o radiationTRAB.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plasma.f90
pressure_correct.o: ./pressure_correct.f90 names.o trot.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./pressure_correct.f90
radiationTR.o: ./radiationTR.f90 names.o multies.o for77test.o mesh.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./radiationTR.f90
radiationTRAB.o: ./radiationTRAB.f90 radiationTR.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./radiationTRAB.f90
solution.o: ./solution.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./solution.f90
solvers.o: ./solvers.f90 geometry.o gravity.o names.o limiters.o flux.o trot.o multies.o boundc.o solution.o pressure_correct.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./solvers.f90
time.o: ./time.f90 names.o trot.o solution.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./time.f90
trot.o: ./trot.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./trot.f90
SRC = ./alloc.f90 ./trot.f90 ./limiters.f90 ./names.f90 ./mesh.f90 ./solution.f90 ./radiationTRAB.f90 ./pressure_correct.f90 ./for77test.f ./plasma.f90 ./initials.f90 ./boundc.f90 ./geometry.f90 ./multies.f90 ./solvers.f90 ./gravity.f90 ./radiationTR.f90 ./flux.f90 ./time.f90 impl.inc
OBJ = alloc.o trot.o limiters.o names.o mesh.o solution.o radiationTRAB.o pressure_correct.o for77test.o plasma.o initials.o boundc.o geometry.o multies.o solvers.o gravity.o radiationTR.o flux.o time.o
clean: neat
	-rm -f .cppdefs $(OBJ) plasma
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
plasma: $(OBJ) 
	$(LD) $(OBJ) -o plasma  $(LDFLAGS)
