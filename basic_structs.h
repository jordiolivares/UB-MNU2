typedef struct vector vector;
struct vector {
    double *elements;
    int num_elems;
};

typedef void (*solver)(vector *, vector *, vector *);
