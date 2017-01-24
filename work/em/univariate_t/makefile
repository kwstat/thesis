# makefile 

CC = g++

#For AXP / CXX
#CFLAGS = -O2 -fprm d

#For AXP / GCC
#CFLAGS = -O2 -mfp-rounding-mode=d

#For MIPS
CFLAGS = -O2 

# Include directories
INCLUDE = -I/home/kwright/include

# Link directories
LINKDIRS = -L/home/kwright/lib/@sys

# Link libraries
LINKLIBS = -lProfilOpt -lProfilNLI -lProfilMISC \
				-lProfil -lProfilLIN -lBias -lm


uni_t: uni_t.cc
	$(CC) -o uni_t $(CFLAGS) $(INCLUDE) uni_t.cc $(LINKDIRS) $(LINKLIBS)

uni_t2: uni_t2.cc
	$(CC) -o uni_t2 $(CFLAGS) $(INCLUDE) uni_t2.cc $(LINKDIRS) $(LINKLIBS)

mult_t.closed: mult_t.closed.cc
	$(CC) -o mult_t.closed $(CFLAGS) $(INCLUDE) mult_t.closed.cc $(LINKDIRS) $(LINKLIBS)

scalar: scalar.cc
	$(CC) -o scalar $(CFLAGS) $(INCLUDE) scalar.cc $(LINKDIRS) $(LINKLIBS)

real: real.cc
	cxx -o real $(CFLAGS)   real.cc $(LINKDIRS) $(LINKLIBS)


# ---------------------------------------------------------------------------

EXE = uni2
OBJ = uni2.o gridlist.o
$(EXE): $(OBJ)
		$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LINKDIRS) $(LINKLIBS)
uni2.o: uni2.cc gridlist.o
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o uni2.o uni2.cc
gridlist.o: gridlist.cc gridlist.h
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o gridlist.o gridlist.cc

clean:
		rm *.o


