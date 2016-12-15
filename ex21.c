#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DELTA (0.01)
#define EPSILON (1.0/10000000000.0)

double f(double x, double y) {
	return (3*x*x+3*y*y-1)*(x*x+y*y-5)*(x*x+y*y-3*x+2) + 1;
}

void gradF(double x, double y, double * result) {
	result[0] = 18*x*x*x*x*x - 45*x*x*x*x + 4*x*x*x*(9*y*y-10) -18*x*x*(3*y*y-8) +
	            2*x*(9*y*y*y*y-20*y*y-27) - 9*y*y*y*y + 48*y*y - 15;
	result[1] = 2*y*(9*x*x*x*x - 18*x*x*x + 2*x*x*(9*y*y - 10) - 6*x*(3*y*y-8) +
	            9*y*y*y*y - 20*y*y - 27);
}

void jacobianMatrix(double x, double x_0, double y, double y_0, double ** mat) {
	double * gradient = mat[0];
	gradF(x, y, gradient);
	mat[1][0] = 2*(x - x_0);
	mat[1][1] = 2*(y - y_0);
}

void F(double x, double x_0, double y, double y_0, double * result) {
	result[0] = f(x, y);
	result[1] = (x - x_0) * (x - x_0) + (y - y_0) * (y - y_0) - DELTA*DELTA;
}

double det2x2(double ** mat) {
	return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

double cramerSolve2x2(double ** mat, double * b, double * result) {
	double det = det2x2(mat);
	double tmp[2][2];
	tmp[0][0] = b[0];
	tmp[1][0] = b[1];
	tmp[0][1] = mat[0][1];
	tmp[1][1] = mat[1][1];
	result[0] = det2x2(tmp) / det;
	tmp[0][0] = mat[0][0];
	tmp[1][0] = mat[1][0];
	tmp[0][1] = b[0];
	tmp[1][1] = b[1];
	result[1] = det2x2(tmp) / det;
}

void newtonMethod(double * x_previous, double * x_next) {
	double mat[2][2];
	double b[2];
	double difference[2];
	jacobianMatrix(x_previous[0], x_next[0],
	               x_previous[1], x_next[1],
	               mat);
	F(x_previous[0], x_next[0],
	  x_previous[1], x_next[1],
	  b);
	cramerSolve2x2(mat, b, difference);
	x_next[0] = difference[0] + x_previous[0];
	x_next[1] = difference[1] + x_previous[1];
}

double norm2(double x, double y) {
	return sqrt(x*x + y*y);
}

void initialGuess(double x, double y, double * guess) {
	double norm = norm2(x,y);
	double tangent[2];
	gradF(x, y, tangent);
	guess[0] = tangent[0] / norm;
	guess[1] = tangent[1] / norm;
}

// TODO: Hacer cosas del gnuplot

int main() {
	// TODO: Implementar cosas pa la gráfica
	return 0;
}
