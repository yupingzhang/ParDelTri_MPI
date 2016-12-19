OBJS := main.o sparsegrid.o triangulator.o

EXECUTABLE := triangulator

##CC := g++
##CC := mpiicc
CXX := mpiicpc

#CFLAGS := -O0 -std=gnu++11 -stdlib=libc++ -Wall
#CFLAGS := -O3 -mkl
CFLAGS := -mkl

LDFLAGS := -lm

%.o: %.c %.h
	$(CXX) $(CFLAGS) -c $< -o $@

$(EXECUTABLE): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@

clean:
	-rm -f *.o