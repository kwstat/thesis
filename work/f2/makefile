# makefile 

CC = g++
CFLAGS = -O2
INCLUDE = -I/home/kwright/include
LINK = -L/home/kwright/lib/dec -lProfilOpt -lProfilNLI -lProfilMISC \
       -lProfil -lProfilLIN -lBias -lm

# Here we teach the makefile how to handle .cc extensions by default
.SUFFIXES: .o .cc
.cc:
				$(CC) $(CFLAGS) $(INCLUDE) -o $@ $< $(LINK)

# Now for specific cases
f2: f2.cc 
				$(CC) $(CFLAGS) $(INCLUDE) -o f2 f2.cc $(LINK)

