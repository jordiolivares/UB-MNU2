#include "basic_structs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define THIRD ((double)(1e0 / 3e0))
#define BETA (2e0 / 3e0)
#define TOL ((double)(1e0 / 1e12))

double norm_infty(vector *vect1, vector *vect2) {
    double max = -INFINITY;
    double tmp;
    for (int i = 0; i < vect1->num_elems; i++) {
        tmp = fabs(vect1->elements[i] - vect2->elements[i]);
        if (tmp > max) {
            max = tmp;
        }
    }
    return max;
}

void print_sol(vector *vect) {
    for (int i = 0; i < 10; i++) {
        printf("x_%d = %.12f\n", i, vect->elements[i]);
    }
    return;
}

void jacobi(vector *old, vector *destination, vector *b) {
    int n = old->num_elems;
    double *oldX = old->elements;
    double *dest = destination->elements;
    double *local_b = b->elements;
    dest[0] = THIRD * (local_b[0] - oldX[2] - oldX[n - 2]);
    dest[1] = THIRD * (local_b[1] - oldX[3] - oldX[n - 1]);
    for (int i = 2; i < n - 2; i++) {
        dest[i] = THIRD * (local_b[i] - oldX[i - 2] - oldX[i + 2]);
    }
    dest[n - 2] = THIRD * (local_b[n - 2] - oldX[n - 4] - oldX[0]);
    dest[n - 1] = THIRD * (local_b[n - 1] - oldX[n - 3] - oldX[1]);
    return;
}

void sor(vector *old, vector *destination, vector *b, double w) {
    int n = old->num_elems;
    double *oldX = old->elements;
    double *dest = destination->elements;
    double *local_b = b->elements;
    dest[0] = (1 - w) * oldX[0] + (w / 3) * (local_b[0] - oldX[2] - oldX[n - 2]);
    dest[1] = (1 - w) * oldX[1] + (w / 3) * (local_b[1] - oldX[3] - oldX[n - 1]);
    for (int i = 2; i < n - 2; i++) {
        dest[i] = (1 - w) * oldX[i] + (w / 3) * (local_b[i] - dest[i - 2] - oldX[i + 2]);
    }
    dest[n - 2] = (1 - w) * oldX[n - 2] + (w / 3) * (local_b[n - 2] - dest[n - 4] - dest[0]);
    dest[n - 1] = (1 - w) * oldX[n - 1] + (w / 3) * (local_b[n - 1] - dest[n - 3] - dest[1]);
    return;
}

void gauss_seidel(vector *old, vector *destination, vector *b) {
    sor(old, destination, b, 1);
    return;
}

vector *solve(double convergence_cutoff, solver func) {
    vector *solution = malloc(sizeof(vector));
    int n = 1000000;
    solution->num_elems = n;
    solution->elements = calloc(n, sizeof(double));
    vector *oldX = malloc(sizeof(vector));
    oldX->num_elems = n;
    oldX->elements = calloc(n, sizeof(double));
    vector *b = malloc(sizeof(vector));
    b->num_elems = n;
    b->elements = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        b->elements[i] = ((double) (i + 1)) / ((double) n);
    }
    double dif = INFINITY;
    vector *tmp;
    int iterations = 0;
    while (dif > convergence_cutoff) {
        func(oldX, solution, b);
        tmp = oldX;
        oldX = solution;
        solution = tmp;
        dif = (BETA / (1 - BETA)) * norm_infty(oldX, solution);
        iterations++;
    }
    free(oldX);
    free(b);
    return solution;
}

void sor_optim(vector *old, vector *destination, vector *b) {
    sor(old, destination, b, 1.1);
    return;
}

int main(int argc, char *argv[]) {
    printf("Jacobi:\n");
    vector *sol_jacobi = solve(TOL, jacobi);
    print_sol(sol_jacobi);
    printf("Gauss-Seidel:\n");
    vector *sol_gauss = solve(TOL, gauss_seidel);
    print_sol(sol_gauss);
    printf("SOR (w = 1.1):\n");
    vector *sol_sor = solve(TOL, sor_optim);
    print_sol(sol_sor);
    return 0;
}
