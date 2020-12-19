#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

double **allocateMatrix(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double** initMatrix(size_t n) {
    double **m = allocateMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = (double)(rand() % 100);
        }
    }
    return m;
}

int main(int argc, char **argv) {

    int id, p;
    double a, c;

    /*
        id - номер текущего процесса
        p - количесво процессов
        c, a - коэф-ты, которые нужно будет умножать элементы
    */
    double time_start, time_end;
    double **m;
    int request = 0;

    bool test = true;
    int matrix_size = 3;
    if (argc > 1) {
        if ((matrix_size = atoi(argv[1])) > 0) {
            test = false;
        }
    }

    MPI_Status status;

    MPI_Init(&argc, &argv);
    // id данного процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // количество процессов
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    m = initMatrix(matrix_size);

    // если это главный процесс, то он инициализирует матрицу
    if (id == 0) {
        cout << matrix_size << ", " << p << ", ";
        time_start = MPI_Wtime();
    }

    // все процессы получают исходную матрицу и единичную матрицу
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&(m[0][0]), matrix_size*matrix_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double lead[matrix_size];
    int block = matrix_size / (p - 1),
            start_row = block * id,
            end_row = block * (id + 1);
    if (id == p - 1)
        end_row = matrix_size;

    for (int i = 0; i < matrix_size; i++) {

        if (id == i / block) {
            for (int j = 0; j < matrix_size ; j++) {
                lead[j] = m[i][j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(lead, matrix_size , MPI_DOUBLE, i / block, MPI_COMM_WORLD);

        for (int j = start_row; j < end_row; j++) {
            if (j == i) {
                double d = m[i][i];
                for (int k = 0; k < matrix_size ; k++)
                    m[j][k] /= d;
                continue;
            }
            double d = m[j][i] / lead[i];
            for (int k = 0; k < matrix_size ; k++) {
                m[j][k] -= d * lead[k];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (id == 0) {

        time_end = MPI_Wtime();
        cout << time_end - time_start << endl;
        std::ofstream fout(argv[2], ios::app);
        fout.is_open();
        fout << p << "," << argv[1] << "," << time_end - time_start << '\n';
        fout.close();
    }
    MPI_Finalize();
    return 0;
}
