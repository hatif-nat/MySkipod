#include <iostream>
#include <random>
#include <iomanip>
#include <mpi.h>
#include <malloc.h>
#include <fstream>
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

double** initMatrix1(size_t n) {
    double **m = alloc_2d(n, n * 2);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n * 2; j++) {
            m[i][j] = (j < n) ? rand() % 100 : ((j - n == i) ? 1 : 0);
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
    double time_start, time_end;
    double **m, **I;
    int request = 0;

    bool test = true;
    int matrix_size = 3;
    if (argc > 1) {
        if ((matrix_size = stoi(argv[1])) > 0) {
            test = false;
        }
    }

    MPI_Status status;

    MPI_Init(&argc, &argv);
    // id данного процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // количество процессов
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    m = initMatrix1(matrix_size);
    I = initMatrix(matrix_size, IDENTITY);

    // если это главный процесс, то он инициализирует матрицу
    if (id == 0) {
        cout << matrix_size << ", " << p << ", ";

        if (test) {
            m[0][0] = 1, m[0][1] = 2, m[0][2] = 1, m[0][3] = 1, m[0][4] = 0, m[0][5] = 0,
            m[1][0] = 3, m[1][1] = 5, m[1][2] = 2, m[1][3] = 0, m[1][4] = 1, m[1][5] = 0,
            m[2][0] = 2, m[2][1] = 3, m[2][2] = 3, m[2][3] = 0, m[2][4] = 0, m[2][5] = 1;
        } else {
            m = initMatrix1(matrix_size);
        }

        if (matrix_size <= 10) {
            cout << "\nInput matrix (generated randomly):\n\n";
            printMatrix(m, matrix_size);
        }

        time_start = MPI_Wtime();
    }

    // все процессы получают исходную матрицу и единичную матрицу
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&(m[0][0]), matrix_size*matrix_size*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double lead[matrix_size * 2];
    int block = matrix_size / (p - 1),
            start_row = block * id,
            end_row = block * (id + 1);
    if (id == p - 1)
        end_row = matrix_size;

    for (int i = 0; i < matrix_size; i++) {

        if (id == i / block) {
            for (int j = 0; j < matrix_size * 2; j++) {
                lead[j] = m[i][j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(lead, matrix_size * 2, MPI_DOUBLE, i / block, MPI_COMM_WORLD);

        for (int j = start_row; j < end_row; j++) {
            if (j == i) {
                double d = m[i][i];
                for (int k = 0; k < matrix_size * 2; k++)
                    m[j][k] /= d;
                continue;
            }
            double d = m[j][i] / lead[i];
            for (int k = 0; k < matrix_size * 2; k++) {
                m[j][k] -= d * lead[k];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (matrix_size <= 10) {
        for (int i = 0; i < matrix_size; i++) {
            if (i / block == id) {
                for (int j = matrix_size; j < matrix_size * 2; j++) {
                    cout << ((m[i][j] >= 0.0 && m[i][j] < 10.0) ? " " : "") << fixed << setprecision(3) << m[i][j] << (j != matrix_size - 1 ? ", " : "");
                }
                cout << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) {
        if (matrix_size <= 10 ) {
            // cout << "\n\nResult:\n\n";
            // printMatrix(I, matrix_size);
        }
        if (id == 0) {
            time_end = MPI_Wtime();
            cout << time_end - time_start << endl;
            std::ofstream fout(argv[2], ios::app);
            fout.is_open();
            fout << p << "," << argv[1] << "," << time_end - time_start << '\n';
            fout.close();
        }
    }
    MPI_Finalize();
    return 0;
}


/*
#include <iostream>
#include <random>
#include <iomanip>
#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <tr1/regex>

using namespace std;

double **alloc_2d(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double** initMatrix(size_t n) {
    double **m = alloc_2d(n, n);
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
    /*double time_start, time_end;
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
}*/
/*#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

fill_matrix(int myrank, int n, int numproc, float ** A, float * leading_row){
    if (!myrank)
    {
        for (i=0; i<n; i++)
        {
            k = i / (n/numproc);
            if (k>=numproc) k = numproc-1;
            if (!k) {
                for (j = 0; j < n; j++) {
                    A[i][j] = (rand()) % 10;
                }
            }
            else
            {
                for (j=0; j<n; j++) {
                    leading_row[j] = (rand()) % 10;
                }
                MPI_Send(leading_row, n, MPI_FLOAT, k, 1, MPI_COMM_WORLD);
            }
        }
    }
}
found_non_zero_string(int myrank,int  source, float ** A, float * leading_row,
int i, int start_row, int last_row, int nrow, int numproc, int kk)
{
    int j = 0;
    source = i / (n/numproc); //процесс-источник (содержащий лидирующую строку)
    if (source >= numproc) source = numproc-1;

    // Выбираем строку с ненулевым ведущим элементом
    if (myrank == source)
    {
        int k = 0;
        if (A[i-start_row][i]==0)
        {
            for (k=i-start_row+1; k<nrow; k++)
            {
                if (A[k][i] != 0)
                {
                    // меняем строки местами
                    for (j=0; j<n; j++)
                    {
                        float t=A[k][j];
                        A[k][j]=A[i-start_row][j];
                        A[i-start_row][j]=t;
                    }
                    break;
                }
            }
        }
        if (k == nrow)
        {
            for (j=myrank+1; j<numproc && k; j++)
            {
                MPI_Send(&j, 1, MPI_INT, j, 3, MPI_COMM_WORLD); //отправляем процессу "сигнал" о замене строчек
                MPI_Recv(leading_row, n, MPI_FLOAT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (leading_row[i] != 0)
                {
                    MPI_Send(A[i-start_row], n, MPI_FLOAT, j, 2, MPI_COMM_WORLD);
                    for (kk=0; kk<n; kk++) A[i-start_row][kk] = leading_row[kk]; //заменяем строчку
                    k=0;
                }
                else MPI_Send(leading_row, n, MPI_FLOAT, j, 2, MPI_COMM_WORLD);
            }
            kk = 0;
            for (; j<numproc; j++)
                MPI_Send(&kk, 1, MPI_INT, j, 3, MPI_COMM_WORLD);
        }
        else
        {
            kk = 0; // Не требуется строчка ни от какого из процессов
            for (j=myrank+1; j<numproc; j++)
                MPI_Send(&kk, 1, MPI_INT, j, 3, MPI_COMM_WORLD);
        }
        if (k == nrow)
        {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            return -1;
        }
    }

    else if (myrank > source)
    {
        MPI_Recv(&kk, 1, MPI_INT, source, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (kk) //нужно отправить строчку
        {
            for (j=0; j<nrow-1 && A[j][i]==0; j++) ; //находим строчку с ненулевым ведущим элементом в каждом процессе
            MPI_Send(A[j], n, MPI_FLOAT, source, 2, MPI_COMM_WORLD); //отправялем ее (или последнюю, если не было найдено)
            MPI_Recv(A[j], n, MPI_FLOAT, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //получаем обратно (возможно измененная)
        }
    }
}

int main(int argc, char *argv[])
{
    int n, error, myrank, numproc, source, kk;
    int start_row, last_row, nrow;
    float **A, *leading_row;
    double start_time, end_time;
    
    double determinant = 1;

    MPI_Init(&argc, &argv); // Инициализируем коммуникатор MPI - через него будет происходить взаимодействие группы процессов

    MPI_Comm_size(MPI_COMM_WORLD, &numproc); //число процессов в коммуникаторе
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //номер процесса в коммуникаторе (нумерация с 0)

    if (!myrank) {
        n = (argc>1) ? atoi(argv[1]) : 512;
    }

    // синхронизируем
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Рассылаем n всем процессам

    nrow = n/numproc;
    start_row = nrow * myrank;
    last_row = start_row + nrow - 1;
    if (myrank==numproc-1) {
        last_row = n-1;
        nrow = last_row - start_row + 1;
    }

    // Выделяем память под элементы рабочей(расширенной) матрицы (исходная + единичная, приписанная справа) - размерность: n * (n*2)
    // Каждый процесс выделает себе nrow строчек и хранит только их
    A = (float**) malloc(nrow * (n) * sizeof(float));

    if (A==NULL){
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        return -1;
    }
    for (int i=0; i<nrow; i++)
    {
        A[i] = (float*) malloc((n) * sizeof(float));
        if (A[i]==NULL){
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            return -1;
        }
    }

    leading_row = (float*) malloc((n) * sizeof(float));

    if (leading_row==NULL){
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        return -1;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    
    // заполняем матрицу рандомными числами
    fill_matrix(myrank, n, numproc, A, leading_row);

    // Каждый процесс получает свои nrow строчек (считанные процессом 0)
    if (myrank) {
        for (int i = 0; i < nrow; i++) {
            MPI_Recv(A[i], n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (!myrank) start_time= MPI_Wtime();
    
    for (int i=0; i<n; i++)
    {

        found_non_zero_string(myrank, source, A, leading_row, i, start_row, last_row, nrow, numproc, kk );

        // Делим всю строку на ее ведущий элемент(делаем ведущий элемент = единице) и аккумулируем в определитель
        if (i >= start_row && i <= last_row) //(myrank == source)
        {
            // отправим 0-му процессу очередное диагональное значение

            MPI_Send(A[i-start_row][i], 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            for (int j=i+1; j<n; j++) A[i-start_row][j] /= A[i-start_row][i];
            A[i-start_row][i] = 1;
            for (int j=0; j<n; j++) leading_row[j] = A[i-start_row][j];

            // Отправка лидирующей строки(содержащей ведущий элемент) 'нижележащим' процессам
            for (int k=myrank + 1; k<numproc; k++)
                MPI_Send(leading_row, n, MPI_FLOAT, k, 1, MPI_COMM_WORLD);
        }

        // С помощью ведущего элемента аннулируем все расположенные под ним элементы i-го столбца
        if (i < last_row)
        {
            if (!myrank){
                float t;
                MPI_Recv(t, n, MPI_FLOAT, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                determinant *= t;
            }
            if (myrank != source){ //сам себе не отправлял
                MPI_Recv(leading_row, n, MPI_FLOAT, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                kk = 0;
            }
            else kk=i-start_row+1;

            for (int k=kk; k<nrow; k++)
            {
                for (int j=i+1; j<n; j++)
                    A[k][j] -= leading_row[j]*A[k][i];
                A[k][i] = 0;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    if (!myrank) {
        end_time = MPI_Wtime();
        std::ofstream fout(argv[2], ios::app);
        fout.is_open();
        fout << numproc << "," << argv[1] << "," << end_time - start_time << '\n';
        fout.close();
    }

    // Освобождаем память, выделенную под матрицу
    for (int i=0; i<nrow; i++) free(A[i]);
    free(A);
    free(leading_row);

    MPI_Finalize();

    return 0;
}
*/
