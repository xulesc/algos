prg_name := psc
prg_srcs := $(wildcard *.c)
prg_objs := ${program_srcs:.c=.o}
CFLAGS := -ffast-math -O3 $(profile)
LFLAGS := -lm -static $(profile)

## For the phdNode 32 bit setup we want to compile the programs
## using the 32 bit v 3.4 gcc setup at custom location in the
## installation
ARCH := $(shell getconf LONG_BIT)
cc_32 := gcc -mtune=pentium -march=pentium
#cc_32 := /opt/i386-gnu/gcc/gcc-3.4.6/bin/gcc
cc_64 := gcc
cc := $(cc_$(ARCH))

all: PSC

PSC: psc.c
	$(cc) -o psc psc.c $(CFLAGS) $(LFLAGS) 

clean:
	rm -f psc


                
