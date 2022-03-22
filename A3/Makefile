parread:
	# rm *.out
	mpic++ -std=c++17 parRead.cpp

binread:
	rm *.out
	mpic++ -std=c++17 bin_read.cpp

check:
	mpic++ -std=c++17 check.cpp

main:
	mpic++ -std=c++11 main.cpp

n = 1   
f = ""
run:
	mpirun -np $(n) ./a.out $(f)


# binread
# correctness and accuracy
# openmp
# time