# Compiler commands
FC = gfortran
CC = g++

# Compiler flags
FFLAGS = -c
CFLAGS = -c

# Linker flags
LDFLAGS = -lgfortran -lGLU -lglut -lGL
GLDFLAGS = -lGLU -lglut -lGL

# Executable name
EXECUTABLE = my_program
RND_EXECUTABLE = my_renderer

# Source files
FORTRAN_SOURCE = ynogk.f90
CPP_SOURCE = Classes.cpp
RND_SOURCE = Renderer.cpp

# Object files
FORTRAN_OBJECT = $(FORTRAN_SOURCE:.f90=.o)
CPP_OBJECT = $(CPP_SOURCE:.cpp=.o)
RND_OBJECT = $(RND_SOURCE:.cpp=.o)


all: $(EXECUTABLE)

$(EXECUTABLE): $(CPP_OBJECT) $(FORTRAN_OBJECT)
	$(CC) $(CPP_OBJECT) $(FORTRAN_OBJECT) -o $@ $(LDFLAGS)

$(FORTRAN_OBJECT): $(FORTRAN_SOURCE)
	$(FC) $(FFLAGS) $< -o $@

$(CPP_OBJECT): $(CPP_SOURCE)
	$(CC) $(CFLAGS) $< -o $@

renderer: $(RND_EXECUTABLE)
$(RND_EXECUTABLE): $(RND_OBJECT)
	$(CC) $(RND_OBJECT) -o $@ $(GLDFLAGS)

$(RND_OBJECT): $(RND_SOURCE)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(FORTRAN_OBJECT) $(CPP_OBJECT) $(EXECUTABLE)

