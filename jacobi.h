#ifndef JACOBI_H
#define JACOBI_H



typedef struct tag_jacobi{
	
	int nr, size;
	double *input; 
	double *sum;
	double *fun;

}jacobi;

double* resultJacobi(jacobi jacobi, int* test, int loop, double* jacobian);

#endif
