#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>

using namespace std;

int main(int argc,char* argv[]) {
    int i,j,k, kk;
    int n, error, myrank, numproc, source;
    MPI_Request *request;
    int start_row, last_row, nrow;
    float **A, *leading_row;
    double start_time, end_time, start_algo;

    error = MPI_Init(&argc, &argv); // Инициализируем
    if (error != MPI_SUCCESS) {
        return -1;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numproc); //число процессов в коммуникаторе
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //номер процесса в коммуникаторе (нумерация с 0)

    if (!myrank) {
        // вводим размер матрицы из аргуметов main
        n = (argc > 1) ? atoi(argv[1]) : 512;
    }

    // Приостановка процессов до выхода ВСЕХ процессов коммуникатора в заданную точку синхронизации
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Широковещательная рассылка n всем процессам


    if (!myrank) start_time = MPI_Wtime();

    nrow = n/numproc;
    start_row = nrow * myrank;
    last_row = start_row + nrow - 1;
    if (myrank==numproc-1) {
        last_row = n-1;
        nrow = last_row - start_row + 1;
    }

    // Выделяем память под элементы матрицы размерность: n * n
    // Каждый процесс выделает себе nrow строчек и хранит только их
    A = (float**) malloc(nrow * (n) * sizeof(float));
    if (A==NULL){
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        return -1;
    }
    for (i=0; i<nrow; i++)
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

    // Заполняем исходную матрицу
    if (!myrank)
    {
        srand(125);
        for (i=0; i<n; i++)
        {
            k = i / (n/numproc);
            if (k>=numproc) k = numproc-1;
            if (!k) {
                for (j = 0; j < n; j++) {
                    A[i][j] = (rand()) % 10;
                }
            }
            else {
                for (j=0; j<n; j++) {
                    leading_row[j] = (rand()) % 10;
                }
                MPI_Send(leading_row, n, MPI_FLOAT, k, 1, MPI_COMM_WORLD);
            }
        }
    }

    // Каждый процесс получает свои nrow строчек (считанные процессом 0)
    if (myrank) {
        for (i = 0; i < nrow; i++)
            MPI_Recv(A[i], n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (!myrank) start_algo= MPI_Wtime();

    // начинаем сами вычисления
    for (i=0; i<n; i++)
    {
            source = i / (n/numproc); //процесс-источник (содержащий лидирующую строку)
            if (source >= numproc) source = numproc-1;

            // Выбираем строку с ненулевым ведущим элементом
            if (myrank == source)
            {
                k = 0;
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
                        MPI_Isend(&kk, 1, MPI_INT, j, 3, MPI_COMM_WORLD, &request);
                }
                else
                {
                    kk = 0; // Не требуется строчка ни от какого из процессов
                    for (j=myrank+1; j<numproc; j++)
                        MPI_Isend(&kk, 1, MPI_INT, j, 3, MPI_COMM_WORLD, &request);
                }
                if (k == nrow)
                {
                    printf("Обратной матрицы не существует.\n"); //определитель = 0
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

            // Делим всю строку на ее ведущий элемент(делаем ведущий элемент = единице)
            if (i >= start_row && i <= last_row) //(myrank == source)
            {
                for (j=0; j<n; j++) leading_row[j] = A[i-start_row][j];

                // Отправка лидирующей строки(содержащей ведущий элемент) 'нижележащим' процессам
                for (k=myrank + 1; k<numproc; k++)
                    MPI_Isend(leading_row, n, MPI_FLOAT, k, 1, MPI_COMM_WORLD, &request);
            }

            // С помощью ведущего элемента аннулируем все расположенные под ним элементы i-го столбца
            if (i < last_row)
            {
                if (myrank != source){ //сам себе не отправлял
                    MPI_Recv(leading_row, n, MPI_FLOAT, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    kk = 0;
                }
                else kk=i-start_row+1;

                for (k=kk; k<nrow; k++)
                {
                    for (j=i+1; j<2*n; j++)
                        A[k][j] -= leading_row[j]*A[k][i]/ leading_row[i];
                    A[k][i] = 0;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    if(!myrank) {

        end_time = MPI_Wtime();

        std::ofstream fout(argv[2], ios::app);
        fout.is_open();
        fout << numproc << "," << argv[1] << "," << end_time - start_algo << '\n';
        fout.close();
    }
    // Освобождаем память, выделенную под матрицу
    for (i=0; i<nrow; i++) free(A[i]);
    free(A);
    free(leading_row);

    MPI_Finalize();

    return 0;
}

