#include <iostream>
#include <random>
#include <iomanip>
#include <mpi.h>
#include <malloc.h>

using namespace std;

enum MatrixType {
    IDENTITY,
    RANDOM,
    ZERO
};

enum errorCode {
    ZERO_DET
};

double **alloc_2d(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double** initMatrix(size_t n, MatrixType type = RANDOM) {
    double **m = alloc_2d(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            switch (type) {
                case RANDOM:
                    m[i][j] = (double)(rand() % 100);
                    break;
                case IDENTITY:
                    m[i][j] = (i == j) ? 1 : 0;
                    break;
                default:
                    m[i][j] = 0;
                    break;
            }
        }
    }
    return m;
}

void printMatrix(double **m, size_t n) {
    cout << "{\n";
    for (int i = 0; i < n; i++) {
        cout << "    {";
        for (int j = 0; j < n; j++) {
            cout << ((m[i][j] >= 0.0 && m[i][j] < 10.0) ? " " : "") << fixed << setprecision(3) << m[i][j] << (j != n - 1 ? ", " : "");
        }
        cout << "}" << (i != n - 1 ? "," : "") << endl;
    }
    cout << '}';
}

int main(int argc, char **argv) {

    int id, p;
    double a, c;

    /*
        id - номер текущего процесса
        p - количесво процессов
        c, a - коэф-ты, которые нужно будет умножать элементы
    */
    time_t time_start, time_end;
    double **m;

    bool test = true;
    int matrix_size = 3;
    if (argc > 1) {
        if ((matrix_size = stoi(argv[1])) > 0) {
            test = false;
        }
    }

    double cur_row[matrix_size];

    MPI_Status status;

    MPI_Init(&argc, &argv);
    // id данного процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // количество процессов
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    m = initMatrix(matrix_size, ZERO);

    // если это главный процесс, то он инициализирует матрицу
    if (id == 0) {
        //cout << "\nWorkers: " << p << endl;

        if (test) {
            m[0][0] = 1, m[0][1] = 2, m[0][2] = 1,
            m[1][0] = 3, m[1][1] = 5, m[1][2] = 2,
            m[2][0] = 2, m[2][1] = 3, m[2][2] = 3;
        } else {
            m = initMatrix(matrix_size, RANDOM);
        }

        if (matrix_size <= 10) {
            //cout << "\nInput matrix (generated randomly):\n\n";
            //printMatrix(m, matrix_size);
        }

        time_start = MPI_Wtime();
    }

    // все процессы получают исходную матрицу и единичную матрицу
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&(m[0][0]), matrix_size*matrix_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < matrix_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (id) {
            a = -1.0 / m[i][i];
            for (int j = id - 1; j < matrix_size; j += p - 1) {
                c = m[j][i] * a;
                if (j != i) {
                    for (int k = 0; k < matrix_size; k++) {
                        m[j][k] += c * m[i][k];
                    }
                }

                MPI_Send(m[j], matrix_size, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
            }
        } else {
            for (int j = 0; j < matrix_size; j++) {
                MPI_Recv(m[j], matrix_size, MPI_DOUBLE, j % (p - 1) + 1, 123, MPI_COMM_WORLD, &status);
            }
            double u = m[i][i];
            for (int k = 0; k < matrix_size; k++) {
                m[i][k] /= u;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&(m[0][0]), matrix_size*matrix_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (id == 0) {
        if (matrix_size <= 10 ) {
            //cout << "\n\nResult:\n\n";
            //printMatrix(I, matrix_size);
        }
        if (id == 0) {
            time_end = MPI_Wtime();
            //cout << "\nMPI execution time: " << time_end - time_start << endl;
            std::ofstream fout(argv[2], ios::app);
            fout.is_open();
            fout << p << "," << argv[1] << "," << time << '\n';
            fout.close();
        }
    }
    MPI_Finalize();
    return 0;
}