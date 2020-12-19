#include <iostream>
#include <cstdlib>

using namespace std;

void to_zero_columb(long double **matrix, int size, int columb)
{
    long double *current_row = new long double[size];
    current_row = matrix[columb];

    for(int i = columb + 1; i < size; i++)
    {
        long double mnoz = matrix[i][columb] / current_row[columb];
        for (int j = columb; j < size; j++)
        {
            matrix[i][j] = matrix[i][j] - mnoz * current_row[j];
        }
    }
}
long double diagonal_el(long double **matrix, int size, int columb) {
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
    to_zero_columb(matrix, size, columb);
    if (sign)
    {
        return matrix[columb][columb];
    }
    else {
        return -matrix[columb][columb];
    }
}

int main() {
    int size = 2;
    std::srand(125);
    long double **matrix;
    matrix = new long double*[size];
    for (int i = 0 ; i < size ; i++)
    {
        matrix[i] = new long double[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = (std::rand()) % 10;
            cout << matrix[i][j] << ' ';
        }
        cout << '\n';
    }

    long double result = 1;

    for (int i = 0; i < size; i++) {

        result *= diagonal_el(matrix, size, i);

    }
    cout.precision(17);
    cout << result;
    delete [] matrix;
    return 0;
}
