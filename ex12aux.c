#include <math.h>
#include <stdio.h>

#define THRESHOLD ((double) 1e-20)

double f(double y) {
    return (3*y*y -1) * (y*y - 5) * (y*y + 2) + 1;
}

double derivF(double y) {
    return (6 * y) * (y*y - 5) * (y*y + 2) +
           (3*y*y -1) * (2 * y) * (y*y + 2) +
           (3*y*y -1) * (y*y - 5) * (2*y);
}

double newtonMethodIteration(double x) {
    return x - (f(x) / derivF(x));
}

int main() {
    double x_0, x_1;
    x_0 = -1;
    x_1 = newtonMethodIteration(x_0);
    while (fabs(x_0 - x_1) > THRESHOLD) {
        x_0 = x_1;
        x_1 = newtonMethodIteration(x_1);
    }
    printf("y = 0\nx = %.20f\n", x_1);
    return 0;
}