
clean:
	rm *.o

compile:
	nvcc -arch=sm_35 -std=c++11 --expt-relaxed-constexpr main.cu -o main.o

test:
	./main.o data_image.txt query_image.txt 10 0.5 1