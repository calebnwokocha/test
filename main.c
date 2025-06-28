#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <omp.h>

/* Clause structure */
typedef struct {
    int *orig;   // original literals
    int  sz;     // number of original literals
} Clause;

/* Read clauses from stdin, skip lit == 0 */
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
        int *temp = NULL, tcap = 0, tsz = 0;
        for (char *tok = strtok(line, " \t\n"); tok; tok = strtok(NULL, " \t\n")) {
            int lit = atoi(tok);
            if (lit == 0) continue;
            if (tsz == tcap) {
                tcap = tcap ? tcap * 2 : 4;
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

/* Build sorted list of distinct vars */
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

/* Process a clause: pad, sort, build forbidden vector */
static void process_clause(const Clause *c,
                           int *vars, int *var_to_idx, int n,
                           int *out_sorted, unsigned char *out_forbidden,
                           bool *early_unsat_flag, int clause_idx)
{
    bool seen_pos[n], seen_neg[n];
    memset(seen_pos, 0, n * sizeof(bool));
    memset(seen_neg, 0, n * sizeof(bool));

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

    for (int j = 0; j < n; j++) {
        if (seen_neg[j]) {
            out_sorted[j]             = -vars[j];
            out_forbidden[clause_idx*n + j] = 1;
        } else {
            out_sorted[j]             =  vars[j];
            out_forbidden[clause_idx*n + j] = 0;
        }
    }
}

/* Print padded clause and its forbidden vector */
static void print_clause(const int *sorted, const unsigned char *forbidden, int n, int idx) {
    printf("Clause %d (padded and sorted):", idx+1);
    for (int j = 0; j < n; j++) printf(" %d", sorted[j]);
    printf("\nClause %d forbidden a[1..%d]:", idx+1, n);
    for (int j = 0; j < n; j++) printf(" %d", forbidden[idx*n + j]);
    printf("\n\n");
}

/* Print binary representation */
static inline void print_binary(uint32_t v, size_t cols) {
    for (size_t b = 0; b < cols; ++b) {
        putchar((v >> (cols - 1 - b)) & 1 ? '1' : '0');
        if (b + 1 < cols) putchar(' ');
    }
    putchar('\n');
}

/* Find and print first n missing binary numbers */
void print_missing_omp(const unsigned char *arr, size_t rows, size_t cols, size_t want)
{
    if (cols == 0 || cols > 32) return;
    uint32_t U = 1U << cols;
    size_t B = (U + 63) >> 6;
    uint64_t *present = malloc(B * sizeof *present);
    if (!present) { perror("malloc"); return; }
    memset(present, 0, B * sizeof *present);

    #pragma omp parallel
    {
        uint64_t *local = malloc(B * sizeof *local);
        if (!local) abort();
        memset(local, 0, B * sizeof *local);

        #pragma omp for nowait schedule(static)
        for (size_t i = 0; i < rows; ++i) {
            uint32_t v = 0;
            for (size_t j = 0; j < cols; ++j)
                v = (v << 1) | (arr[i*cols + j] & 1);
            local[v >> 6] |= 1ULL << (v & 63);
        }

        #pragma omp critical
        for (size_t k = 0; k < B; ++k)
            present[k] |= local[k];

        free(local);
    }

    size_t printed = 0;
    for (uint32_t v = 0; v < U && printed < want; ++v) {
        if ((present[v >> 6] & (1ULL << (v & 63))) == 0) {
            print_binary(v, cols);
            ++printed;
        }
    }
    printf("\n");
    free(present);
}

int main(void) {
    while (1) {
        int M, max_id, n;
        bool *seen;
        Clause *C = read_clauses(&M, &max_id, &seen);

        int *vars = build_vars(seen, max_id, &n);
        int *var_to_idx = calloc(max_id + 1, sizeof(int));
        for (int i = 0; i < n; i++) var_to_idx[vars[i]] = i;

        printf("Detected %d distinct variable(s) and %d clause(s)\n\n", n, M);

        int (*sorted)[n] = malloc(M * sizeof *sorted);
        unsigned char *forbidden = malloc((size_t)M * n * sizeof *forbidden);
        bool early_unsat = false;

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < M; i++) {
            if (!early_unsat)
                process_clause(&C[i], vars, var_to_idx, n,
                               sorted[i], forbidden, &early_unsat, i);
        }

        if (!early_unsat) {
            for (int i = 0; i < M; i++) print_clause(sorted[i], forbidden, n, i);

            size_t total = 1U << n;
            size_t want = total > (size_t)M ? total - (size_t)M : 0;
            printf("First %lu missing assignments (out of %lu possible):\n", (unsigned long)want, (unsigned long)total);
            print_missing_omp(forbidden, M, n, want);
        }

        for (int i = 0; i < M; i++) free(C[i].orig);
        free(C); free(seen); free(vars); free(var_to_idx);
        free(sorted); free(forbidden);
    }
    return 0;
}
