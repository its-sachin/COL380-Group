mpiexec --bind-to none ./userRead.o $3 $1
mpiexec --bind-to none ./main.o $1 $2 $4