all : omp mpi
omp:
	xlc++ ../omp/src/test_omp.cpp -qsmp=omp -o ../test_omp
mpi:
	mpicxx ../mpi/src/test_mpi.cpp -o  ../test_mpi

clean:
	rm -rf ../test_mpi
	rm -rf ../test_omp
	rm -rf ../.tmp
	rm -f  ../*.err ../*.out