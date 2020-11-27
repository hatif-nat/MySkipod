all: omp mpi
omp: 
	g++ ../omp/src/test_omp.cpp -fopenmp -o ../test_omp
mpi: 
	mpicxx ../mpi/src/test_mpi.cpp -o  ../test_mpi

clean:
	rm -rf ../test_mpi
	rm -rf ../test_omp
