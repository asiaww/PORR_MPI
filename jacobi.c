#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jacobi.h"
#include <mpi.h>

double* resultJacobi(jacobi jacobi, int* test, int loop, double* jacobian){
        
	jacobi.input = (double *) calloc(jacobi.size, sizeof(double));
	if (loop == 1){
		
		int z;
		for (z = 0; z < jacobi.size; z++){
			jacobi.input[z] = 2.0;	
		}
	}
	else{
		int z;
		for (z = 0; z < jacobi.size; z++){
			jacobi.input[z] = jacobian[z];	
		}
	}
	jacobi.fun = (double *) calloc(jacobi.size, sizeof(double));
        jacobi.sum = (double *) calloc(jacobi.size, sizeof(double));
	int i;

		
	double m;
	double ex;
	if (jacobi.nr == 1) {
		ex = 3.0;
		m = 3.0;
	}
	if (jacobi.nr == 2){
		ex = 5.0;
		m = 5.0;
	}
	if (jacobi.nr == 3){ 
		ex = 4.0;
		m = 4.0;
	}
	else{
		ex = 7.0;
		m = 4.0;
	}

        for (i = 0; i < jacobi.size; i ++){
		int k;            
		for (k = 0; k < jacobi.size; k ++){
                	if(i != k){
                    		jacobi.sum[i] -= (double)(test[i*jacobi.size+k])*jacobi.input[k]-(i+1)*m;
				
                	}
            	}
        }
	
	int n;
        for (int n = 0; n < jacobi.size; n++){
            jacobi.fun[n] = pow(jacobi.sum[n]/(test[n+n*jacobi.size]*jacobi.input[n]),1.0/ex);
        }
	free(jacobi.input);	
	free(jacobi.sum);
	return jacobi.fun;
}
