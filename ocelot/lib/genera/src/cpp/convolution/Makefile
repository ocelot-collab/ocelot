
TO_BUILD = convolution

all: $(TO_BUILD)


COMPILER ?= g++
OPENMP_FLAGS ?= -fopenmp
CFLAGS = -fPIC -O3 $(OPENMP_FLAGS) -Wall

#COMPILER = g++-4.8
#CFLAGS = -fPIC -O3 -m64 -fopenmp -Wall


SOURCE = convolution.cpp

OBJ = convolution.o

#$(OBJ): $(SOURCE)
#	$(COMPILER) -c $(CFLAGS) $(SOURCE)

$(TO_BUILD): $(SOURCE)
	$(COMPILER) -shared $(CFLAGS) -o $(TO_BUILD).so $(SOURCE)


.PHONY : clean

clean:
	rm *.o *.a *.so

