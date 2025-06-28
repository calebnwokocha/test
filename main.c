#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <time.h>

/*// Portable print for size_t: assumes C99+ or casts to unsigned long if needed.
#define FMT_ZU "%zu"*/

typedef struct {
    int *orig;   // original literals
    int  sz;     // number of original literals
} Clause;

// Phase 1: read clauses, skip any lit == 0
static Clause* read_clauses(int *out_M, int *out_max_id, bool **out_seen) {
    Clause *C = NULL;
    int cap = 0, M = 0, max_id = 0;
    bool *seen = calloc(1024, sizeof(bool));
    char line[2048];

    printf("Enter clauses (blank line to finish):\n");
    while (fgets(line, sizeof(line), stdin) && *line != '\n') {
        if (M == cap) {
            cap = cap ? cap * 2 : 4;
            C   = realloc(C, cap * sizeof(Clause));
        }
        // parse into a temporary vector
        int *temp = NULL, tcap = 0, tsz = 0;
        for (char *tok = strtok(line, " \t\n"); tok; tok = strtok(NULL, " \t\n")) {
            int lit = atoi(tok);
            if (lit == 0) continue;        // ← skip zeros entirely
            if (tsz == tcap) {
                tcap = tcap ? tcap*2 : 4;
                temp = realloc(temp, tcap * sizeof(int));
            }
            temp[tsz++] = lit;
            int v = abs(lit);
            if (v > max_id) max_id = v;
            if (v < 1024)  seen[v] = true;
        }
        C[M].orig = temp;
        C[M].sz   = tsz;
        M++;
    }

    *out_M      = M;
    *out_max_id = max_id;
    *out_seen   = seen;
    return C;
}

// Phase 2: build sorted list of distinct vars
static int* build_vars(bool *seen, int max_id, int *out_n) {
    int n = 0;
    for (int v = 1; v <= max_id; v++)
        if (v < 1024 && seen[v]) n++;
    int *vars = malloc(n * sizeof(int));
    for (int v = 1, i = 0; v <= max_id; v++)
        if (v < 1024 && seen[v]) vars[i++] = v;
    *out_n = n;
    return vars;
}

// Phase 3: scan, pad+“sort”, build forbidden in one pass
static void process_clause(const Clause *c,
                           int *vars, int *var_to_idx, int n,
                           int *out_sorted, unsigned char *out_forbid_Q,
                           bool *early_unsat_flag, int clause_idx)
{
    bool seen_pos[n], seen_neg[n];
    memset(seen_pos, 0, n * sizeof(bool));
    memset(seen_neg, 0, n * sizeof(bool));

    // detect p vs. ¬p
    for (int i = 0; i < c->sz; i++) {
        int lit = c->orig[i];
        int v   = abs(lit);
        int idx = var_to_idx[v];
        if (lit > 0) {
            if (seen_neg[idx]) {
                #pragma omp critical
                {
                    printf("UNSAT (clause %d contains both %d and -%d)\n\n",
                           clause_idx+1, v, v);
                    *early_unsat_flag = true;
                }
                return;
            }
            seen_pos[idx] = true;
        } else {
            if (seen_pos[idx]) {
                #pragma omp critical
                {
                    printf("UNSAT (clause %d contains both %d and -%d)\n\n",
                           clause_idx+1, v, v);
                    *early_unsat_flag = true;
                }
                return;
            }
            seen_neg[idx] = true;
        }
    }

    // build padded+sorted & forbidden
    for (int j = 0; j < n; j++) {
        if (seen_neg[j]) {
            out_sorted[j]  = -vars[j];
            out_forbid_Q[clause_idx * n + j] = 1;  // Store forbidden vector here
        } else {
            // either seen_pos[j] or padded positive
            out_sorted[j]  =  vars[j];
            out_forbid_Q[clause_idx * n + j] = 0;
        }
    }
}

// Phase 4: print one clause
static void print_clause(const int *sorted, const unsigned char *forbid_Q, int n, int idx) {
    printf("Clause %d (padded and sorted):", idx+1);
    for (int j = 0; j < n; j++)
        printf(" %d", sorted[j]);
    printf("\nClause %d forbidden a[1..%d]:", idx+1, n);
    for (int j = 0; j < n; j++)
        printf(" %d", forbid_Q[idx * n + j]);
    printf("\n\n");
}

/*// Computes minimal period via KMP failure function (linear time).
static size_t minimal_period(const unsigned char *s, size_t m) {
    size_t *pi = malloc(m * sizeof *pi), k = 0;
    if (!pi) { perror("malloc"); exit(EXIT_FAILURE); }
    pi[0] = 0;
    for (size_t i = 1; i < m; i++) {
        while (k > 0 && s[k] != s[i]) k = pi[k - 1];
        if (s[k] == s[i]) k++;
        pi[i] = k;
    }
    size_t p = m - pi[m - 1];
    free(pi);
    return (m % p == 0 ? p : m);
}


// Fisher–Yates shuffle in-place
static void fisher_yates(size_t *a, size_t m) {
    for (size_t i = m - 1; i > 0; i--) {
        size_t j = rand() % (i + 1);
        size_t tmp = a[i]; a[i] = a[j]; a[j] = tmp;
    }
}

// Test: no fixed point, and each row changes at least one entry
static bool valid_permutation(const unsigned char *Q, const size_t *o,
                              size_t n, size_t m)
{
    // No fixed points
    for (size_t j = 0; j < m; j++)
        if (o[j] == j) return false;

    // Each row changed somewhere
    for (size_t i = 0; i < n; i++) {
        const unsigned char *Qi = Q + i * m;
        bool changed = false;
        for (size_t j = 0; j < m; j++)
            if (Qi[j] != Qi[o[j]]) {
                changed = true;
                break;
            }
        if (!changed) return false;
    }
    return true;
}

// Build a derangement σ ; expected few retries until success
static bool build_safe_derangement(size_t *o,
                                   const unsigned char *Q,
                                   size_t n, size_t m)
{
    srand((unsigned)time(NULL));
    // Initialize to identity permutation
    for (size_t j = 0; j < m; j++)
        o[j] = j;

    // Retry until valid permutation found
    const int maxRetries = 10;
    for (int attempt = 0; attempt < maxRetries; attempt++) {
        fisher_yates(o, m);
        if (valid_permutation(Q, o, n, m))
            return true;
    }
    return false;
}

// Computes V = Q ★ σ (deranged) with full safety checks.
bool compute_V_safe(const unsigned char *Q, unsigned char *V,
                    size_t n, size_t m, size_t t) {
    if (m < 2) {
        fputs("Error: need m ≥ 2\n", stderr);
        return false;
    }
    // Pre-scan: check periods
    for (size_t i = 0; i < n; i++) {
        size_t p = minimal_period(Q + i * m, m);
        if (p == 1) {
            fprintf(stderr, "Warning: constant row at " FMT_ZU ", will force change\n", i);
        }
    }
    // Build derangement
    size_t *sigma = malloc(m * sizeof *sigma);
    if (!sigma) { perror("malloc"); return false; }
    if (!build_safe_derangement(sigma, Q, n, m)) {
        fputs("Error: cannot build derangement that changes all rows\n", stderr);
        free(sigma);
        return false;
    }
    // Compute V
    for (size_t i = 0; i < n; i++) {
        const unsigned char *Qi = Q + i * m;
        unsigned char *Vi = V + i * t;
        for (size_t j = 0; j < t; j++)
            Vi[j] = Qi[sigma[j % m]];
    }
    free(sigma);
    return true;
}*/

// Function to print a 2D array
void print_2d_array(unsigned char *array, size_t rows, size_t cols) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            printf("%d ", array[i * cols + j]);
        }
        printf("\n");
    }
}

int main(void) {
    while (1) {
        int M, max_id, n;
        bool *seen;
        Clause *C = read_clauses(&M, &max_id, &seen);

        int *vars       = build_vars(seen, max_id, &n);
        int *var_to_idx = calloc(max_id+1, sizeof(int));
        for (int i = 0; i < n; i++)
            var_to_idx[vars[i]] = i;

        printf("Detected %d distinct variable(s) and %d clause(s)\n\n", n, M);

        // pre-allocate shared scratch
        int  (*sorted)[n] = malloc(M * sizeof *sorted);
        unsigned char *forbid_Q  = malloc(M * n * sizeof(unsigned char));
        bool early_unsat = false;

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < M; i++) {
            if (!early_unsat) {
                process_clause(&C[i],
                               vars, var_to_idx, n,
                               sorted[i], forbid_Q,
                               &early_unsat, i);
            }
        }

        if (!early_unsat) {
            for (int i = 0; i < M; i++)
                print_clause(sorted[i], forbid_Q, n, i);
        }

/*        // Compute V using compute_V_safe
        size_t t = n; // Assuming t equals n for this example
        unsigned char *V = malloc(M * t * sizeof(unsigned char));
        if (compute_V_safe(forbid_Q, V, M, n, t)) {
            printf("Matrix V (resulting matrix):\n");
            print_2d_array(V, M, t);
        } else {
            printf("Failed to compute matrix V.\n");
        }
*/

        // cleanup
        for (int i = 0; i < M; i++)
            free(C[i].orig);
        free(C);
        free(seen);
        free(vars);
        free(var_to_idx);
        free(sorted);
        free(forbid_Q);
        //free(V);
    }
    return 0;
}
