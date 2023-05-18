CXX = g++

CXXFLAGS = -g -Wall -std=gnu++17

task2: matrix.o cl_matrix.o
	$(CXX) matrix.o cl_matrix.o -o task2
	rm matrix.o cl_matrix.o

matrix.o: matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp

cl_matrix.o: cl_matrix.cpp cl_matrix.h
	$(CXX) $(CXXFLAGS) -c cl_matrix.cpp

clean:
	rm *.o task2