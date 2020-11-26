#include <iostream>
#include <math.h>
#include <random>
// #include <omp.h>
using namespace std;

#define LEFT_BR 0
#define RIGHT_BR 255
#define EPS 1.0e-15
#define PI 3.14159265358979323846
#define EXP 2.71828182846

const double pi_half = PI/2.0;
double f(double x) { // функция в точке x 
	double x_in_substitution = x*x;
	return sin(x_in_substitution * pi_half) + pow(EXP, -(x_in_substitution));	
}

// double simpson_integral_par(double a, double b, int n) { 
//   const double h = (b-a)/n;
//   double k1 = 0, k2 = 0;
//   #pragma omp parallel for reduction(+:k1,k2) // num_threads(4)
//   for(int i = 1; i < n; i += 2) {
//     k1 += f(a + i*h);
//     k2 += f(a + (i+1)*h);
//   }
//   return h/3*(f(a) + 4*k1 + 2*k2);
// }

// double par_calc() {
// 	int M = 4;
// 	int N = 134217728;
// 	double h = (RIGHT_BR-LEFT_BR)/M;
// 	double sum = 0;
// 	#pragma omp parallel num_threads(M) reduction(+:sum)
// 	{
// 		int thread_num = omp_get_thread_num();
// 		double local_a = LEFT_BR + h*thread_num;
// 		double local_b = local_a + h;
// 		sum = simpson_integral_par(local_a, local_b, N);
// 	}
// }

double simpson_integral(double a, double b, int n) { 
  const double h = (b-a)/n;
  double k1 = 0, k2 = 0;
  for(int i = 1; i < n; i += 2) {
    k1 += f(a + i*h);
    k2 += f(a + (i+1)*h);
  }
  return h/3*(f(a) + 4*k1 + 2*k2);
}

double calc() {
	// int M = 4;
	int N = 134217728;
	// double h = (RIGHT_BR-LEFT_BR)/M;
	double sum = 0;
	sum = simpson_integral(LEFT_BR, RIGHT_BR, N);
}

int main() {

	cout << myFunc();
	// double I = EPS + 1, I1 = 0;//I-предыдущее вычисленное значение интеграла, I1-новое, с большим N.

	// for (int N = 2; (N <= 4) || (fabs(I1 - I) > EPS); N *= 2) {
	// 	double h, sum2=0, sum4=0, sum=0;
	// 	h = (RIGHT_BR - LEFT_BR) / (2.0 * N);//Шаг интегрирования.
	// 	// printf("h : %f\n" , h);
	// 	#pragma omp parallel for reduction(+:sum2)
	// 	for (int i = 1; i <= 2 * N - 1; i += 2) {   
	// 		sum4 += f(LEFT_BR + h * i);//Значения с нечётными индексами, которые нужно умножить на 4.
	// 		sum2 += f(LEFT_BR + h * (i + 1));//Значения с чётными индексами, которые нужно умножить на 2.
	// 	}
	// 	printf("I : %e\nN : %d\n\n", I, N);
	// 	sum = f(LEFT_BR) + 4 * sum4 + 2 * sum2 - f(RIGHT_BR);//Отнимаем значение f(RIGHT_BR) так как ранее прибавили его дважды. 
	// 	I = I1;
	// 	I1 = (h / 3) * sum;
	// }

	// printf("%e", I1);


	return 0;
}