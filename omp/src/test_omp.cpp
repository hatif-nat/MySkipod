#include <omp.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

//export OMP_THREAD_LIMIT = 600

void to_zero_columb(long double **matrix, int size, int columb, int nThreads)
{
    long double *current_row = new long double[size];
    current_row = matrix[columb];

    #pragma omp parallel for shared(matrix, current_row) num_threads(nThreads)
    for (int i = columb + 1; i < size; i++) {
        if (omp_get_num_threads() != int(nThreads)) {
            //std::cerr << omp_get_num_threads() << std::endl;
        }
        long double mnoz = matrix[i][columb] / current_row[columb];
        for (int j = columb; j < size; j++) {
            matrix[i][j] = matrix[i][j] - mnoz * current_row[j];
        }
    }

}
long double diagonal_el(long double **matrix, int size, int columb, int nThreads)
{
    //found row with no zero element in columb's columb
    long double *current_row = matrix[columb];
    bool found = 0;
    bool sign = 1;
    for (int i = columb; i < size; i++)
    {
        if (matrix[i][i] != 0) {
            found = 1;
            if (i != columb) {
                sign = 0;
                matrix[columb] = matrix[i];
                matrix[i] = current_row;
            }
            break;
        }
    }
    if (!found)
    {
        return 0;
    }
    to_zero_columb(matrix, size, columb, nThreads);

    if (sign)
    {
        return matrix[columb][columb];
    }
    else {
        return -matrix[columb][columb];
    }
}

int main(int argc, char* argv[]) {
    int size = 2;

    if(argc != 4) return -1;

    int nThreads = atoi(argv[1]);
    size =  atoi(argv[2]);

    srand(125);
    long double **matrix;
    matrix = new long double*[size];
    for (int i = 0 ; i < size ; i++)
    {
        matrix[i] = new long double[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = (rand()) % 10;
        }
    }

    long double result = 1;

    double timerOpenMp = omp_get_wtime();

    for (int i = 0; i < size; i++) {
        result *= diagonal_el(matrix, size, i, nThreads);
    }

    timerOpenMp = omp_get_wtime() - timerOpenMp;

    std::ofstream fout(argv[3], ios::app);
    fout.is_open();
    fout << nThreads << "," << size << "," << timerOpenMp << '\n';
    fout.close();

    for (int i = 0 ; i < size ; i++) {
        //delete [] matrix[i];
    }
    delete [] matrix;

    //cout.precision(17);
    //cout << result;
    return 0;
}
