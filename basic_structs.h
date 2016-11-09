typedef struct vector vector;
struct vector {
    double *elements;
    int num_elems;
};

typedef void (*solver)(vector *, vector *, vector *);
typedef double (*solver_error)(vector *, vector *, vector *);
