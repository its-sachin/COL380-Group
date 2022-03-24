mpic++ -std=c++11 -O2 -pg -fopenmp -o main2.o main2.cpp 
mpiexec --bind-to none ./main2.o outFolder 5 ../../test/user.txt out.txt
# mpiexec --bind-to none ./main2.o outFolder 5 dummy/vect.txt out.txt