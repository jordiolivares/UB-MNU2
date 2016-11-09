#include "basic_structs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define THIRD ((double)(1e0 / 3e0))
#define BETA (2e0 / 3e0)
#define TOL ((double)(1e0 / 1e12))
#define N 1000000

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

double dot(vector *vect1, vector * vect2) {
    double result = 0;
    for (int i = 0; i < vect1->num_elems; i++) {
        result += vect1->elements[i] * vect2->elements[i];
    }
    return result;
}

void transpA(vector *vect, vector *result) {
    int n = vect->num_elems;
    double *dest = result->elements;
    double *x = vect->elements;
    dest[0] = 3 * x[0] + x[2] + x[n-2];
    dest[1] = 3 * x[1] + x[3] + x[n-1];
    for (int i = 2; i < n - 2; i++) {
        dest[i] = 3 * x[i] + x[i+2] + x[i-2];
    }
    dest[n-2] = 3 * x[n-2] + x[n-4] + x[0];
    dest[n-1] = 3 * x[n-1] + x[n-3] + x[1];
    return;
}

void gradient(vector *x, vector *b, vector *result) {
    transpA(x, result);
    for (int i = 0; i < b->num_elems; i++) {
        result->elements[i] -= b->elements[i];
    }
}
void print_sol(vector *vect) {
    for (int i = 0; i < 10; i++) {
        printf("x_%d = %.12f\n", i, vect->elements[i]);
    }
    return;
}

double steepest_descent(vector *old, vector *result, vector *b, double w) {
    vector *p_k = malloc(sizeof(vector));
    vector *transp = malloc(sizeof(vector));
    transp->num_elems = N;
    p_k->num_elems = N;
    transp->elements = malloc(sizeof(double) * N);
    p_k->elements = malloc(sizeof(double) * N);
    gradient(old, b, p_k);
    transpA(p_k, transp);
    double error = dot(p_k, p_k);
    double a_k =  error / dot(transp, p_k);
    for (int i = 0; i < result->num_elems; i++) {
        result->elements[i] = old->elements[i] - w * a_k * p_k->elements[i];
    }
    free(p_k->elements);
    free(transp->elements);
    free(p_k);
    free(transp);
    return sqrt(error);
}

double steepest_descent_optim(vector *old, vector *result, vector *b) {
    return steepest_descent(old, result, b, 0.8);
}

vector *solve(double convergence_cutoff, solver_error func) {
    vector *solution = malloc(sizeof(vector));
    solution->num_elems = N;
    solution->elements = malloc(sizeof(double) * N);
    vector *oldX = malloc(sizeof(vector));
    oldX->num_elems = N;
    oldX->elements = malloc(sizeof(double) * N);
    oldX->elements[0] = 10;
    vector *b = malloc(sizeof(vector));
    b->num_elems = N;
    b->elements = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        b->elements[i] = ((double) (i + 1)) / ((double) N);
    }
    double error = INFINITY;
    vector *tmp;
    while (error > convergence_cutoff) {
        error = func(oldX, solution, b);
        tmp = oldX;
        oldX = solution;
        solution = tmp;
    }
    free(oldX);
    free(b);
    return solution;
}

int main() {
    printf("Steepest descent:\n");
    vector *sol = solve(TOL, steepest_descent_optim);
    print_sol(sol);
    return 0;
}
