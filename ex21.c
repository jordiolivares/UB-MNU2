double f(double x, double y) {
	return (3*x*x+3*y*y-1)*(x*x+y*y-5)*(x*x+y*y-3*x+2) + 1;
}

void gradF(double x, double y, double * result) {
	result[0] = 18*x*x*x*x*x - 45*x*x*x*x + 4*x*x*x*(9*y*y-10) -18*x*x*(3*y*y-8) +
	            2*x*(9*y*y*y*y-20*y*y-27) - 9*y*y*y*y + 48*y*y - 15;
	result[1] = 2*y*(9*x*x*x*x - 18*x*x*x + 2*x*x*(9*y*y - 10) - 6*x*(3*y*y-8) +
	            9*y*y*y*y - 20*y*y - 27);
}

int main() {
	printf("f(1,0.5) = %e \n", f(1, 0.5));
	double * grad = malloc(sizeof(double) * 2);
	gradF(1, 0.5, grad);
	printf("gradF(1, 0.5) = (%e , %e) \n", grad[0], grad[1]);
	return 0;
}
