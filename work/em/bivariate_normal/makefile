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


norm: norm.cc
	$(CC) -o norm $(CFLAGS) $(INCLUDE) norm.cc $(LINKDIRS) $(LINKLIBS)


real: real.cc
	cxx -o real $(CFLAGS)   real.cc $(LINKDIRS) $(LINKLIBS)


EXE = binorm
OBJ = binorm.o gridlist.o
$(EXE): $(OBJ)
		$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LINKDIRS) $(LINKLIBS)
binorm.o: binorm.cc
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o binorm.o binorm.cc
gridlist.o: gridlist.cc gridlist.h
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o gridlist.o gridlist.cc


