FILES=../omp/src/test_omp.cpp
APP_NAME=../t
FLAGS=-fopenmp -o

all: omp mpi
omp: 
	g++ ../omp/src/test_omp.cpp -fopenmp -o ../test_omp
mpi: 
	mpic++ ../mpi/src/test_mpi.cpp -o  ../test_mpi
runmpi: 
	mpirun -np 12 ../test_mpi
clean:
	rm -rf ../test_mpi
	rm -rf ../test_omp
