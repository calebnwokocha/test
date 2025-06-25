#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define LEAF_SIZE 64   // threshold for switching to classical multiplication

// Allocate a zeroed square matrix of size n×n
static double* alloc_matrix(int n) {
    double *M = calloc((size_t)n * n, sizeof(double));
    if (!M) { perror("calloc"); exit(EXIT_FAILURE); }
    return M;
}

// Compute next power of two ≥ x
static int next_pow2(int x) {
    int p = 1;
    while (p < x) p <<= 1;
    return p;
}

// Classical multiplication: C += A × B, size r×r
static void classical_mul(double *A, double *B, double *C, int r) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < r; i++)
        for (int k = 0; k < r; k++) {
            double aik = A[i*r + k];
            for (int j = 0; j < r; j++)
                C[i*r + j] += aik * B[k*r + j];
        }
}

// Add/Subtract: C = A ± B over r×r
static void add_sub(double *A, double *B, double *C, int r, int sign) {
    int n2 = r * r;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n2; i++)
        C[i] = A[i] + sign * B[i];
}

// Recursive Strassen multiplication: C += A × B, size r×r
static void strassen_rec(double *A, double *B, double *C, int r) {
    if (r <= LEAF_SIZE) {
        classical_mul(A, B, C, r);
        return;
    }
    int h = r / 2;
    #define SUB(M,i,j) (M + (i)*h*r + (j)*h)
    double *A11 = SUB(A,0,0), *A12 = SUB(A,0,1), *A21 = SUB(A,1,0), *A22 = SUB(A,1,1);
    double *B11 = SUB(B,0,0), *B12 = SUB(B,0,1), *B21 = SUB(B,1,0), *B22 = SUB(B,1,1);
    double *C11 = SUB(C,0,0), *C12 = SUB(C,0,1), *C21 = SUB(C,1,0), *C22 = SUB(C,1,1);

    int szb = h * h;
    double *S1 = malloc(szb * sizeof(double)), *S2 = malloc(szb * sizeof(double)),
           *S3 = malloc(szb * sizeof(double)), *S4 = malloc(szb * sizeof(double)),
           *S5 = malloc(szb * sizeof(double)), *S6 = malloc(szb * sizeof(double)),
           *S7 = malloc(szb * sizeof(double)), *S8 = malloc(szb * sizeof(double)),
           *S9 = malloc(szb * sizeof(double)), *S10 = malloc(szb * sizeof(double));
    double *P1 = alloc_matrix(h), *P2 = alloc_matrix(h), *P3 = alloc_matrix(h),
           *P4 = alloc_matrix(h), *P5 = alloc_matrix(h), *P6 = alloc_matrix(h),
           *P7 = alloc_matrix(h);

    add_sub(B12, B22, S1, h, -1);
    add_sub(A11, A12, S2, h, +1);
    add_sub(A21, A22, S3, h, +1);
    add_sub(B21, B11, S4, h, -1);
    add_sub(A11, A22, S5, h, +1);
    add_sub(B11, B22, S6, h, +1);
    add_sub(A12, A22, S7, h, -1);
    add_sub(B21, B22, S8, h, +1);
    add_sub(A11, A21, S9, h, -1);
    add_sub(B11, B12, S10, h, +1);

    #pragma omp task shared(P1) firstprivate(A11, S1)
        strassen_rec(A11, S1, P1, h);
    #pragma omp task shared(P2) firstprivate(S2, B22)
        strassen_rec(S2, B22, P2, h);
    #pragma omp task shared(P3) firstprivate(S3, B11)
        strassen_rec(S3, B11, P3, h);
    #pragma omp task shared(P4) firstprivate(A22, S4)
        strassen_rec(A22, S4, P4, h);
    #pragma omp task shared(P5) firstprivate(S5, S6)
        strassen_rec(S5, S6, P5, h);
    #pragma omp task shared(P6) firstprivate(S7, S8)
        strassen_rec(S7, S8, P6, h);
    #pragma omp task shared(P7) firstprivate(S9, S10)
        strassen_rec(S9, S10, P7, h);

    #pragma omp taskwait

    add_sub(P5, P4, S1, h, +1);
    add_sub(S1, P2, S2, h, -1);
    add_sub(S2, P6, C11, h, +1);
    add_sub(P1, P2, C12, h, +1);
    add_sub(P3, P4, C21, h, +1);
    add_sub(P5, P1, S1, h, +1);
    add_sub(S1, P3, S2, h, -1);
    add_sub(S2, P7, C22, h, -1);

    free(S1); free(S2); free(S3); free(S4); free(S5);
    free(S6); free(S7); free(S8); free(S9); free(S10);
    free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7);
}

// Wrapper: multiply Q (n×k) × Q^T (k×m) → K (n×m)
static void strassen_mul(double *A, int n, int k,
                         double *B, int kb, int m,
                         double *C) {
    int S = next_pow2(((n>m?n:m)>k? (n>m?n:m) : k));
    double *Ap = alloc_matrix(S);
    double *Bp = alloc_matrix(S);
    double *Cp = alloc_matrix(S);

    for (int i = 0; i < n; i++)
        memcpy(Ap + i*S, A + i*k, k * sizeof(double));
    for (int i = 0; i < k; i++)
        memcpy(Bp + i*S, B + i*m, m * sizeof(double));

    #pragma omp parallel
    #pragma omp single
    strassen_rec(Ap, Bp, Cp, S);

    for (int i = 0; i < n; i++)
        memcpy(C + i*m, Cp + i*S, m * sizeof(double));

    free(Ap); free(Bp); free(Cp);
}

int main() {
    int n, m;
    printf("Enter n (rows of Q) and m (columns of Q), separated by space: ");
    if (scanf("%d %d", &n, &m) != 2) {
        fprintf(stderr, "Invalid input for n and m.\n");
        return EXIT_FAILURE;
    }

    double *Q = malloc((size_t)n * m * sizeof(double)); // Query
    double *K = calloc((size_t)n * n, sizeof(double)); // Key
    double *V = calloc((size_t)n * m, sizeof(double)); // Value
    int *A = malloc((size_t)n * m * sizeof(int)); // Assignment
    if (!Q || !K || !V || !A) { perror("malloc"); return EXIT_FAILURE; }

    printf("Enter the %d x %d entries of Q, row by row:\n", n, m);
    for (int i = 0; i < n*m; i++) {
        if (scanf("%lf", &Q[i]) != 1) {
            fprintf(stderr, "Invalid matrix entry.\n");
            return EXIT_FAILURE;
        }
    }

    // Compute K = Q × Q^T. The transpose of Q is Q^T
    strassen_mul(Q, n, m, Q, m, n, K);

    // Compute V = K × Q
    strassen_mul(K, n, n, Q, m, m, V);

    printf("Q^T is transpose of Q. Result of V = (Q x Q^T) x Q:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf("%.6g%s", V[i*m + j], j+1<m ? " " : "");
        printf("\n");
    }

    // Compute A based on V: A[i] = 1 if V[i] > 0, else 0
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n * m; i++) {
        A[i] = (V[i] > 0.0) ? 1 : 0;
    }

    printf("Matrix A (1 if corresponding V > 0, else 0):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf("%d%s", A[i*m + j], j+1<m ? " " : "");
        printf("\n");
    }

    free(Q); free(K); free(V); free(A);
    return 0;
}
