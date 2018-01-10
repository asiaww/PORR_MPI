#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "testFun.h"
#include "jacobi.h"
#include "richardson.h"
#include <math.h>

int main(int argc, char* argv[]) {

//zmienne potrzebne do obliczen
    funkcja pierwsza;
    jacobi jacobi;
    richardson rich;

    int n = 10000;
    double precision = 0.000000001;
    double precision_rich = 0.000001;
    pierwsza.nr = jacobi.nr = rich.nr = 4;
    pierwsza.size = jacobi.size = rich.size = n;

    double sum;

/// Tworzenie funkcji testowej ///

    int *vector;

    if (pierwsza.nr == 1){
        vector = createTestF1(pierwsza);
        free(pierwsza.test);
    }

    if (pierwsza.nr == 2){
        vector = createTestF2(pierwsza);
        free(pierwsza.test);
    }

    if (pierwsza.nr == 3){
        vector = createTestF3(pierwsza);
        free(pierwsza.test);
    }

    if (pierwsza.nr == 4){
        vector = createTestF4(pierwsza);
        free(pierwsza.test);
    }

    /// iteracyjne wyniki dla Jacobiego i Richardsona ///
    double* jacobian;
    double* prevJ;
    double* prevR;
    double* richardson;

    double* richardson_result;

    MPI_Status status;

    // Initialize the MPI environment
	MPI_Init(&argc, &argv);

    // Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    clock_t start = clock(), diff;

    if (world_rank % 2 == 0) {
    for (int loop = 1; loop < 50; loop++)
    {
        jacobian = resultJacobi(jacobi, vector, loop, jacobian);
        int i;
        double diff = 0.0;

        if (loop > 1) {
            for (i = 0; i < n; i++){
                   double abs = fabs(jacobian[i]-prevJ[i]);
                   if (abs > diff) diff = abs;
            }

            if (diff < precision) {
                break; 
            }
        }

       prevJ = jacobian;
       free(jacobi.fun);    
   }
   }

    if (world_rank % 2 == 1) {
    for (int loop2 = 1; loop2 < 50; loop2++){
        richardson = resultRichardson(rich, vector, loop2, richardson);
        int i;
        double diff = 0.0;
        
        if (loop2 > 1) {
            for (i = 0; i < n; i++){
               double abs = fabs(richardson[i]-prevR[i]);
               if (abs > diff) diff = abs;
            }
        
            if (diff < precision_rich*10) {
                break;
            }
        }

        prevR = richardson;
        free(rich.fun);
    }
  }

    //MPI_Bcast(&jacobian, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&richardson, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);    

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

///////////////// WYSWIETLANIE WYNIKOW /////////////////
if (world_rank == 0 || world_rank == 2) {
    printf("\n\n\nWyniki dla Jacobiego\n");
    int i;
    for (i = 0; i < n; i ++){
        printf("%f\n",jacobian[i]);
    }
    printf("Czas wykonania: %d\n", msec%1000);
}

if (world_rank == 1 || world_rank == 3) {
    printf("\n\n\nWyniki dla Richardsona\n");    
    int h;
    for (h = 0; h < n; h ++){
        printf("%f\n",richardson[h]);
    }   
    printf("\n");
    printf("Czas wykonania: %d\n", msec%1000);
}
    MPI_Finalize();
}
