#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DELTA ((double) 0.01)
#define THRESHOLD ((double) 1e-10)

double f(double x, double y) {
	return (3*x*x+3*y*y-1)*(x*x+y*y-5)*(x*x+y*y-3*x+2) + 1;
}

void gradF(double x, double y, double result[2]) {
	result[0] = (6*x)*(x*x+y*y-5)*(x*x+y*y-3*x+2) +
                (3*x*x+3*y*y-1)*(2*x)*(x*x+y*y-3*x+2) +
                (3*x*x+3*y*y-1)*(x*x+y*y-5)*(2*x - 3);
	result[1] = (6*y)*(x*x+y*y-5)*(x*x+y*y-3*x+2) +
                (3*x*x+3*y*y-1)*(2*y)*(x*x+y*y-3*x+2) +
                (3*x*x+3*y*y-1)*(x*x+y*y-5)*(2*y);
}

void jacobianMatrix(double center_point[2], double point[2], double mat[2][2]) {
	double x, y, x_0, y_0;
	x = point[0];
	y = point[1];
	x_0 = center_point[0];
	y_0 = center_point[1];
	double gradient[2];
	gradF(x, y, gradient);
	mat[0][0] = gradient[0];
    mat[0][1] = gradient[1];
	mat[1][0] = 2*(x - x_0);
	mat[1][1] = 2*(y - y_0);
}

void F(double center_point[2], double point[2], double result[2]) {
	double x, y, x_0, y_0;
	x = point[0];
	y = point[1];
	x_0 = center_point[0];
	y_0 = center_point[1];
	result[0] = f(x, y);
	result[1] = (x - x_0) * (x - x_0) + (y - y_0) * (y - y_0) - DELTA*DELTA;
}

double det2x2(double mat[2][2]) {
	return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

void cramerSolve2x2(double mat[2][2], double b[2], double result[2]) {
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

void newtonMethod(double x_previous[2], double x_next[2], double center_point[2]) {
	double mat[2][2];
	double b[2];
	double difference[2];
	jacobianMatrix(center_point,
				   x_previous,
	               mat);
	F(center_point,
	  x_previous,
	  b);
    b[0] = -b[0];
    b[1] = -b[1];
	cramerSolve2x2(mat, b, difference);
	x_next[0] = difference[0] + x_previous[0];
	x_next[1] = difference[1] + x_previous[1];
}

double norm2(double x, double y) {
	return sqrt(x*x + y*y);
}

void initialGuess(double x, double y, double guess[2]) {
	double tangent[2];
    gradF(x, y, tangent);
    double norm = norm2(tangent[0],tangent[1]);
    guess[0] = x + DELTA*DELTA*(tangent[0] / norm);
	guess[1] = y + DELTA*DELTA*(tangent[1] / norm);
}

double vectorDifference2(double v1[2], double v2[2]) {
	return sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) +
	            (v1[1] - v2[1]) * (v1[1] - v2[1]));
}

void swapVectors(double v1[2], double v2[2]) {
	double tmp[2];
	tmp[0] = v1[0];
	tmp[1] = v1[1];
	v1[0] = v2[0];
	v1[1] = v2[1];
	v2[0] = tmp[0];
	v2[1] = tmp[1];
}

void calculateNextPoint(double center_point[2], double next_point[2]) {
	double previous_point[2];
	initialGuess(center_point[0], center_point[1], previous_point);
	newtonMethod(previous_point, next_point, center_point);
	while (fabs(f(next_point[0], next_point[1])) > THRESHOLD) {
		swapVectors(previous_point, next_point);
		newtonMethod(previous_point, next_point, center_point);
	}
}

int main() {
	double initial_point[2], center_point[2], next_point[2];
	// Obtained after running ex21aux.c
	//center_point[0] = initial_point[0] = -0.59220292723488865416; // X
    //center_point[0] = initial_point[0] = -2.23490048993232548469;
    center_point[0] = initial_point[0] = 0;
	center_point[1] = initial_point[1] = -0.60311683131456095275; // Y
    //center_point[1] = initial_point[1] = 0;
    printf("%.16f %.16f\n", initial_point[0], initial_point[1]);
	do {
		calculateNextPoint(center_point, next_point);
		printf("%.16f %.16f\n", next_point[0], next_point[1]);
		swapVectors(center_point, next_point);
	} while(vectorDifference2(center_point, initial_point) > DELTA*DELTA);
	return 0;
}
