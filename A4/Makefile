compile:
	nvcc -arch=sm_35 -std=c++11 --expt-relaxed-constexpr main.cu -o main.o
clean:
	rm *.o
data = data_image.txt
query = query_image.txt
t1 = 10
t2 = 0.5
n = 1
test:
	./main.o $(data) $(query) $(t1) $(t2) $(n)
