
all:
	g++ main.cpp -std=c++11 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3/ -L /usr/local/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
