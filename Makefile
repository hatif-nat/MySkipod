FILES=test.cpp
APP_NAME=t
FLAGS=-fopenmp -o

all: omp mpi
omp: test.cpp
	g++ test.cpp $(FLAGS) $(APP_NAME)
mpi: test_mpi.cpp
	mpic++ test_mpi.cpp -o  test_mpi
omp: 
clean:
	rm -rf $(APP_NAME)
	rm -rf test_mpi
