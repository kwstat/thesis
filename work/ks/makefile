# makefile 

CC = g++
CFLAGS = -O2 -mfp-rounding-mode=d
INCLUDE = -I/home/kwright/include
LINK = -L/home/kwright/lib/axp -lProfilOpt -lProfilNLI -lProfilMISC \
       -lProfil -lProfilLIN -lBias -lm

# Here we teach the makefile how to handle .cc extensions by default
.SUFFIXES: .o .cc
.cc:
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $< $(LINK)

# Now for specific cases
ks: ks.cc 
	$(CC) $(CFLAGS) $(INCLUDE) -o ks ks.cc $(LINK)

