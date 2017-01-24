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



EXE = genetic
OBJ = genetic.o gridlist.o
$(EXE): $(OBJ)
		$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LINKDIRS) $(LINKLIBS)
genetic.o: genetic.cc
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o genetic.o genetic.cc
gridlist.o: gridlist.cc gridlist.h
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o gridlist.o gridlist.cc


