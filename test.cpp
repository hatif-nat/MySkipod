#include <iostream>
#include <math.h>
#include <random>
#include <omp.h>
#include <ctime>
using namespace std;

#define LEFT_BR 0
#define RIGHT_BR 127
#define EPS 1.0e-6
#define PI 3.14159265358979323846
#define EXP 2.71828182846

const double pi_half = PI/2.0;
double f(double x) { // функция в точке x 
	double x_in_substitution = x*x;
	return sin(x_in_substitution * pi_half) + pow(EXP, -(x_in_substitution));	
}

double simpson_integral_par(double a, double b, int n) { 
  const double h = (b-a)/n;
  double k1 = 0, k2 = 0;
  #pragma omp parallel for reduction(+:k1,k2)
  for(int i = 1; i < n; i += 2) {
    k1 += f(a + i*h);
    k2 += f(a + (i+1)*h);
  }
  return h/3*(f(a) + 4*k1 + 2*k2);
}

double par_calc() {
	unsigned int M = 1; // число ядер
	int N = 134217728; // число разбиений отрезка
	double h = (RIGHT_BR-LEFT_BR)/ (double) M; // h for each process
	double I1 = EPS + 1;
	double I = 0;
	// for (int N = 2; abs(I - I1) >= EPS; N*=2) {
		// I1 = I; 
		// I = 0;
		#pragma omp parallel num_threads(M) reduction(+:I)
		{
			int thread_num = omp_get_thread_num();
			cout << thread_num << endl;
			double local_a = LEFT_BR + h*thread_num;
			double local_b = local_a + h;
			I = simpson_integral_par(local_a, local_b, N);
		}
		
	// }

	return I;
}

double simpson(double a, double b, int n) { 
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
	double I = 0, I1 = EPS;

	// for(int N = 2; abs(I - I1) >= EPS; N*=2) {
		// I1 = I;
		I = simpson(LEFT_BR, RIGHT_BR, N); 
	// }

	return I;
}

int main() {
	clock_t start = clock();
	cout << par_calc() << '\n';
	clock_t end = clock();
	cout << "runtime: " << (end - start) / (double) CLOCKS_PER_SEC << endl;
	
	start = clock();
	cout << calc() << '\n';
	end = clock();
	cout << "runtime: " << (end - start) / (double) CLOCKS_PER_SEC << endl;

	return 0;
}