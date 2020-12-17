#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>

using namespace std;

const double my_exp  = 2.71828182846;
const double pi_half = acos(-1.0) / 2.;
const float left_br = 0;  /* lower limit of integration */
const float right_br = 255; /* upper limit of integration */

typedef double(*my_func)(double);

double f(double x) { 
	double x_sq = x*x;
	return sin(x_sq) + pow(my_exp, -x_sq);	
}

double simpson(my_func f,double a, int num, double h) { 
  double k1 = 0, k2 = 0;
  for(int i = 0; i < num; i += 2) {
    k1 += f(a+i*h);
    k2 += f(a+i*h+h);
  }
  return h/3*(f(a) + 4*k1 + 2*k2);
}

int main(int argc,char* argv[]) {
  int n, p, i, j, ierr,num;
  double h, result, a, b, pi;
  double my_a, my_range;
  double startwtime, my_time = 0.0, mintime = 0.0, time = 0.0;
  int myid, source, dest, tag;
  MPI_Status status;
  double my_result;   
  pi = acos(-1.0);  /* = 3.14159... */
  a = 0.;           /* lower limit of integration */
  b = 1;            /* upper limit of integration */
  n = (argc>1) ? pow(2, atoi(argv[1])) : 512;    /* number of increment within each process */

  dest = 0;         /* define the process that computes the final result */
  tag = 123;        /* set the tag to identify this particular job */
      
/* Starts MPI processes ... */

    MPI_Init(&argc,&argv);                 /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &p);     /* get number of processes */
    h = (right_br-left_br)/n;    /* length of increment */
    num = n / p;	    /* number of intervals calculated by each process*/
    my_range = (right_br-left_br)/p;
    my_a = left_br + myid*my_range;

    startwtime = MPI_Wtime();
    my_result = simpson(f,my_a,num,h);
    my_time = MPI_Wtime() - startwtime;
 
    if(myid == 0) { 
      result = my_result;
      time = my_time;
      for (i=1;i<p;i++) {
        source = i;           /* MPI process number range is [0,p-1] */
        MPI_Recv(&my_result, 1, MPI_DOUBLE, source, tag,
                      MPI_COMM_WORLD, &status);
        MPI_Recv(&my_time, 1, MPI_DOUBLE, source, tag-1,
                      MPI_COMM_WORLD, &status);              
        result += my_result;
        time = std::max(time, my_time);
      }
      std::ofstream fout(argv[2], ios::app);
      fout.is_open();
      fout << p << "," << argv[1] << "," << time << '\n';
      fout.close();
    }
    else
      MPI_Send(&my_result, 1, MPI_DOUBLE, dest, tag,
                    MPI_COMM_WORLD);      /* send my_result to intended dest.*/
      MPI_Send(&my_time, 1, MPI_DOUBLE, dest, tag-1,
                    MPI_COMM_WORLD);      /* send my_time to intended dest.*/        
    MPI_Finalize();                       /* let MPI finish up ... */
  return 0;
}

