#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <float.h>

#define LEAF_SIZE 64    // threshold for switching to classical multiplication
#define PAD -1.0        // sentinel for padding slots

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
static void classical_mul(const double *A, const double *B, double *C, int r) {
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < r; i++) {
        for (int k = 0; k < r; k++) {
            double aik = A[i*r + k];
            for (int j = 0; j < r; j++)
                C[i*r + j] += aik * B[k*r + j];
        }
    }
}

// Add/Subtract: C = A ± B over r×r
static void add_sub(const double *A, const double *B, double *C, int r, int sign) {
    int n2 = r * r;
    #pragma omp parallel for schedule(dynamic)
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

// Multiply square matrices A and B (n×n) into C (n×n)
static void strassen_mul_square(const double *A, const double *B, int n, double *C) {
    int S = next_pow2(n);
    double *Ap = alloc_matrix(S), *Bp = alloc_matrix(S), *Cp = alloc_matrix(S);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) memcpy(Ap + i*S, A + i*n, n * sizeof(double));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) memcpy(Bp + i*S, B + i*n, n * sizeof(double));

    #pragma omp parallel
    #pragma omp single nowait
    strassen_rec(Ap, Bp, Cp, S);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) memcpy(C + i*n, Cp + i*S, n * sizeof(double));

    free(Ap); free(Bp); free(Cp);
}

// Encode text into Q (n×n) with sentinel padding
static void encode_with_sentinel(const unsigned char *text, double *Q, size_t n) {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
        size_t base = i * n;
        Q[base] = (double)text[i];
        for (size_t j = 1; j < n; j++) {
            Q[base + j] = PAD;
        }
    }
}

// Decode with sentinel (for verification)
static unsigned char* decode_with_sentinel(const double *Q, size_t n) {
    unsigned char *text = malloc(n);
    if (!text) return NULL;
    for (size_t i = 0; i < n; i++) {
        double v = Q[i * n];
        if (v == PAD || v < 0.0 || v > 255.0 || v != floor(v)) {
            free(text);
            return NULL;
        }
        text[i] = (unsigned char)v;
    }
    return text;
}

int main() {
    omp_set_num_threads(omp_get_max_threads());
    while (1) {
        int n;
        printf("Enter n (text length and matrix dimension): "); fflush(stdout);
        if (scanf("%d", &n) != 1 || n <= 0) {
            fprintf(stderr, "Invalid input for n.\n");
            return EXIT_FAILURE;
        }

        // Read exact-length text
        unsigned char *text = malloc(n + 1);
        if (!text) { perror("malloc"); return EXIT_FAILURE; }
        printf("Enter text (%d chars): ", n);
        scanf(" ");
        fgets((char*)text, n+1, stdin);
        if (strlen((char*)text) < (size_t)n) {
            fprintf(stderr, "Text too short.\n");
            free(text);
            return EXIT_FAILURE;
        }

        size_t size = n;
        double *Q = malloc(size * size * sizeof(double)); // Query
        double *K = calloc(size * size, sizeof(double)); // Key
        double *V = calloc(size * size, sizeof(double)); // Value
        int    *A = malloc(size * size * sizeof(int));
        if (!Q || !K || !V || !A) { perror("malloc"); return EXIT_FAILURE; }

        // Encode text into square Q
        encode_with_sentinel(text, Q, size);
        free(text);

        // Compute K = Q × Q and V = K × Q
        strassen_mul_square(Q, Q, size, K);
        strassen_mul_square(K, Q, size, V);

        // Print V
        printf("Result of V = Q × Q × Q:\n");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++)
                printf("%.6g%s", V[i*size + j], j+1<size ? " " : "");
            printf("\n");
        }

        // Build A (1 if V>0, else 0)
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < size*size; i++)
            A[i] = (V[i] > 0.0) ? 1 : 0;

        printf("Matrix A (1 if V>0 else 0):\n");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++)
                printf("%d%s", A[i*size + j], j+1<size ? " " : "");
            printf("\n");
        }

        free(Q); free(K); free(V); free(A);
        printf("\n");
    }
    return 0;
}
