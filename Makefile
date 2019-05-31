gcc quicksort_OpenMP.c quicksort_OpenMP.h -o quicksort_OpenMP -lm -fopenmp
gcc PSRS_OpenMP.c -o PSRS_OpenMP -lm -fopenmp
mpicc quicksort_MPI.c -o quicksort_MPI
mpicc PSRS_MPI.c -o PSRS_MPI
