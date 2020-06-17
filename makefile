# Simple makefile to build SQDFT code on mps at GA Tech
# This make file should be located outside the src folder

all:sqdft # This rule compiles all the rules and is the default rule to compile when we type "make"

#LIBS = -I${MKLROOT}/include -L${MKLROOT}/lib/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread  # external libraries required

# ---------compiler to be used----------
CC = mpic++ #g++ #mpicc 

# ---------flags to be passed while compiling--------
  #  -I gives info about the location of header files
  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings
  #  -O3 turns on highest level of optimization
CFLAGS = -g -funroll-loops -I./inc -g -Wall #-Wno-unused #-funroll-all-loops
SOURCEC = ./src/main.cpp ./src/readfiles.cpp ./src/deallocate.cpp ./src/spline.cpp ./src/initialize.cpp ./src/poisson.cpp ./src/anderson.cpp ./src/energy.cpp ./src/scf.cpp ./src/sq.cpp ./src/forces.cpp ./src/nonlocal.cpp ./src/md.cpp
SOURCEH = ./inc/ds_sq.h ./inc/func_sq.h ./inc/headers.h 
OBJSC = ./src/main.o ./src/readfiles.o ./src/deallocate.o ./src/spline.o ./src/initialize.o ./src/poisson.o ./src/anderson.o ./src/energy.o ./src/scf.o ./src/sq.o ./src/forces.o ./src/nonlocal.o ./src/md.o


%.o: %.cpp ${SOURCEH}
	${CC} -c -o $@ $< ${CFLAGS} 

sqdft:${OBJSC} 
	${CC} -o ./lib/sqdft $^ ${CFLAGS} ${LIBS}

.PHONY: clean
clean:
	rm -f ./src/*.o *~ core core*
