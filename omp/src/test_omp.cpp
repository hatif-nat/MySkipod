#include <iostream>
#include <math.h>
#include <omp.h>
#include <ctime>
#include <fstream>
using namespace std;

const double my_exp  = 2.71828182846;
const double pi_half = acos(-1.0) / 2.;
const float left_br = 0;    /* lower limit of integration */
const float right_br = 255; /* upper limit of integration */

typedef double(*my_func)(double);

double f(double x) { // функция в точке x 
	double x_sq = x*x;
	return sin(x_sq) + pow(my_exp, -(x_sq));	
}

double simpson_integral(my_func f, double a, double b, int n) { 
  const double h = (b-a)/n;
  double k1 = 0, k2 = 0;
  for(int i = 1; i < n; i += 2) {
    k1 += f(a + i*h);
    k2 += f(a + (i+1)*h);
  }
  return h/3*(f(a) + 4*k1 + 2*k2);
}

double par_calc(int cores, int pows) {
	int local_N = pow(2,pows) / cores; // число разбиений отрезка
	double h = (right_br-left_br)/ (double) cores; // h for each process
	double I = 0;
	#pragma omp parallel num_threads(cores) reduction(+:I)
	{
		int thread_num = omp_get_thread_num();
		double local_a = left_br + h*thread_num;
		double local_b = local_a + h;
		I = simpson_integral(f, local_a, local_b, local_N);
	}

	return I;
}


int main(int argc, char* argv[]) {

	int cores = atoi(argv[1]);
	int pows =  atoi(argv[2]);

	// printf("cores: %d\npows: %d\n", cores, pows);
	// cout << argv[3];
	
	double timerOpenMp = omp_get_wtime();
	double result = par_calc(cores, pows) ;
	timerOpenMp = omp_get_wtime() - timerOpenMp;

	std::ofstream fout(argv[3], ios::app);
    fout.is_open();
    fout << cores << "," << pows << "," << timerOpenMp << '\n';
	fout.close();
	
	return 0;
}
