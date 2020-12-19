all : omp mpi
omp:
	xlc++ ../omp/src/test_omp.cpp -qsmp=omp -o ../test_omp
mpi:
	mpicxx -o test_mpi -std=c++11 mpi.cpp

clean:
	rm -rf ../test_mpi
	rm -rf ../test_omp
	rm -rf ../.tmp
	rm -f  ../*.err ../*.out