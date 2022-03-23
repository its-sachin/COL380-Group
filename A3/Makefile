parread:
	# rm *.out
	mpic++ -std=c++17 -O2 -pg -fopenmp parRead.cpp

binread:
	rm *.out
	mpic++ -std=c++17 bin_read.cpp

check:
	mpic++ -std=c++17 check.cpp

main:
	mpic++ -std=c++11 -O2 -pg -fopenmp main.cpp 

n = 1   
f = ""
t = 4
run:
	mpirun -np $(n) ./a.out $(f) $(t)
runmain:
	mpirun -np $(n) ./a.out $(t)


# binread
# correctness and accuracy
# openmp
# time