# Compiler commands
FC = gfortran
CC = g++

# Compiler flags
FFLAGS = -c
CFLAGS = -c

# Linker flags
LDFLAGS = -lgfortran -lGLU -lglut -lGL

# Executable name
EXECUTABLE = my_program

# Source files
FORTRAN_SOURCE = ynogk.f90
CPP_SOURCE = Classes2.cpp

# Object files
FORTRAN_OBJECT = $(FORTRAN_SOURCE:.f90=.o)
CPP_OBJECT = $(CPP_SOURCE:.cpp=.o)

all: $(EXECUTABLE)

$(EXECUTABLE): $(CPP_OBJECT) $(FORTRAN_OBJECT)
	$(CC) $(CPP_OBJECT) $(FORTRAN_OBJECT) -o $@ $(LDFLAGS)

$(FORTRAN_OBJECT): $(FORTRAN_SOURCE)
	$(FC) $(FFLAGS) $< -o $@

$(CPP_OBJECT): $(CPP_SOURCE)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(FORTRAN_OBJECT) $(CPP_OBJECT) $(EXECUTABLE)

