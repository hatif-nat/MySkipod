#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

double **allocateMatrix(int n) {
    double *line = (double *)malloc(n*n*sizeof(double));
    double **matrix= (double **)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++)
        matrix[i] = &(line[n*i]);

    return matrix;
}

double** initMatrix(int n) {
    double **matrix = allocateMatrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = (double)(rand() % 10);
        }
    }
    return matrix;
}

int main(int argc, char **argv) {

    int id, allProc;

    double time_start, time_end;
    double **matrix;
    double determinant = 1;
    int size = atoi(argv[1]);


    MPI_Status status;

    MPI_Init(&argc, &argv);
    // id данного процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // количество процессов
    MPI_Comm_size(MPI_COMM_WORLD, &allProc);

    matrix = initMatrix(size);

    // если это главный процесс, то он инициализирует матрицу
    if (id == 0) {
        cout << size << ", " << allProc << ", ";
        time_start = MPI_Wtime();
    }

    // все процессы получают исходную матрицу и единичную матрицу
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&(matrix[0][0]), size*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double lead[size];
    int block = size / (allProc - 1),
            start_row = block * id,
            end_row = block * (id + 1);
    if (id == allProc - 1)
        end_row = size;

    for (int i = 0; i < size; i++) {

        if (id == i / block) {
            for (int j = 0; j < size ; j++) {
                lead[j] = matrix[i][j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(lead, size , MPI_DOUBLE, i / block, MPI_COMM_WORLD);

        for (int j = start_row; j < end_row; j++) {
            if (j == i) {
                double d = matrix[i][i];
                if (i < end_row && i >= start_row) {
                    MPI_Send(d, 1, MPI_Double, id, 2, MPI_COMM_WORLD); //отправялем часть определителя
                }
                if (!id) {
                    MPI_Recv(d, n * 2, MPI_FLOAT, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    determinant *= d;
                }
                for (int k = 0; k < size ; k++)
                    matrix[j][k] /= d;
                continue;
            }
            double d = matrix[j][i] / lead[i];
            for (int k = 0; k < size ; k++) {
                matrix[j][k] -= d * lead[k];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (id == 0) {

        time_end = MPI_Wtime();
        std::ofstream fout(argv[2], ios::app);
        fout.is_open();
        fout << allProc << "," << argv[1] << "," << time_end - time_start << '\n';
        fout.close();
    }
    MPI_Finalize();
    return 0;
}
