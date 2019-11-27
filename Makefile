main: main.cpp
	mpic++ -std=c++11 main.cpp -o main -lOpenCL

exe:
	mpiexec -np 4 ./main