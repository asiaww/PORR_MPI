# PORR_PRI

1. Zainstalowac MPICH2 http://mpitutorial.com/tutorials/installing-mpich2/
2. Kompilacja 

mpicc.mpich2 -std=c99 main.c testFun.c jacobi.c richardson.c -o program -lm

3. Uruchomienie
mpiexec -np 4 ./program

-np - liczba procesow
