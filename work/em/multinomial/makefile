pkg = odrpack

SHLIB = $(pkg)$(SHLIB_EXT)

all: $(SHLIB)

clean:
	@-rm -f *.o *.a

---------------------------------------------------------------------------
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


simple: simple.cc
	$(CC) -o simple $(CFLAGS) $(INCLUDE) simple.cc $(LINKDIRS) $(LINKLIBS)

mult: mult.cc
	$(CC) -o mult $(CFLAGS) $(INCLUDE) mult.cc $(LINKDIRS) $(LINKLIBS)

Ex2: Ex2.cc 
				$(CC) $(CFLAGS) $(INCLUDE) -o Ex2 Ex2.cc $(LINK)

multInterval: multInterval.cc
	$(CC) -o multInterval $(CFLAGS) $(INCLUDE) multInterval.cc $(LINKDIRS) $(LINKLIBS)


EXE = mult2
OBJ = mult2.o gridlist.o
$(EXE): $(OBJ)
		$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LINKDIRS) $(LINKLIBS)
mult2.o: mult2.cc
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o mult2.o mult2.cc
gridlist.o: gridlist.cc gridlist.h
		$(CC) $(CFLAGS) $(INCLUDE) -I. -c -o gridlist.o gridlist.cc


