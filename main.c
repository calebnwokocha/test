#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <omp.h>


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
                           int *out_sorted, unsigned char *out_forbidden,
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
            out_forbidden[clause_idx * n + j] = 1;  // Store forbidden vector here
        } else {
            // either seen_pos[j] or padded positive
            out_sorted[j]  =  vars[j];
            out_forbidden[clause_idx * n + j] = 0;
        }
    }
}

// Phase 4: print one clause
static void print_clause(const int *sorted, const unsigned char *forbidden, int n, int idx) {
    printf("Clause %d (padded and sorted):", idx+1);
    for (int j = 0; j < n; j++)
        printf(" %d", sorted[j]);
    printf("\nClause %d forbidden a[1..%d]:", idx+1, n);
    for (int j = 0; j < n; j++)
        printf(" %d", forbidden[idx * n + j]);
    printf("\n\n");
}

// Convert row into limbs
static void row_to_limbs(const unsigned char *row, uint64_t *limb, int C, int LIMBS) {
    memset(limb, 0, LIMBS * sizeof *limb);
    for (int j = 0; j < C; j++) {
        int limb_i = j / 64;
        int off    = j % 64;            // pack bit 0 → LSB, bit 1 → next, ...
        limb[limb_i] |= (uint64_t)(row[j] & 1) << off;
    }
}

// Convert limbs back to row
static void limbs_to_row(const uint64_t *limb, unsigned char *row, int C) {
    for (int j = 0; j < C; j++) {
        int limb_i = j / 64;
        int off    = j % 64;
                        //printf("no error");
        row[j] = (limb[limb_i] >> off) & 1;
                //printf("no error");
    }
}

// Lexicographic compare
static int limbs_cmp(const void *p, const void *q, void *arg) {
    int LIMBS = *(int *)arg;
    const uint64_t *a = p, *b = q;
    for (int i = LIMBS - 1; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return  1;
    }
    return 0;
}

// Increment limbs; return 1 if overflow beyond C bits
static int limbs_inc(uint64_t *a, int LIMBS) {
    for (int i = LIMBS - 1; i >= 0; i--) {
        uint64_t old = a[i];
        a[i] = old + 1;
        if (a[i] > old) return 0;
    }
    return 1;
}

// Sort with context argument
static void sort_r(void *base, size_t nmemb, size_t size,
    int (*cmp)(const void *, const void *, void *), void *arg)
{
    char *a = base, *tmp = malloc(size);
    for(size_t gap = nmemb/2; gap > 0; gap /= 2) {
        for(size_t i = gap; i < nmemb; i++) {
            memcpy(tmp, a + i*size, size);
            size_t j = i;
            while(j >= gap && cmp(a + (j-gap)*size, tmp, arg) > 0) {
                memcpy(a + j*size, a + (j-gap)*size, size);
                j -= gap;
            }
            memcpy(a + j*size, tmp, size);
        }
    }
    free(tmp);
}

static void process(int M, int C, int LIMBS, const unsigned char *forbidden) {
    // Pack rows
    uint64_t *A = malloc((size_t)M * LIMBS * sizeof *A);
    //#pragma omp parallel for
    for (int i = 0; i < M; i++)
        row_to_limbs(forbidden + (size_t)i * C, A + (size_t)i * LIMBS, C, LIMBS);

    // Sort
    sort_r(A, M, LIMBS * sizeof(uint64_t), limbs_cmp, &LIMBS);

    // Compute offsets
    uint64_t *offs = malloc((size_t)M * sizeof *offs);
    uint64_t G = 0;
    for (int i = 0; i + 1 < M; i++) {
        uint64_t *a = A + (size_t)i * LIMBS;
        uint64_t *b = A + (size_t)(i + 1) * LIMBS;
        // pick the most significant limb (top 64 bits of your C‑bit vector):
        uint64_t ai = a[LIMBS - 1];
        uint64_t bi = b[LIMBS - 1];
        uint64_t gap = (bi > ai) ? bi - ai - 1 : 0;
        offs[i] = G;
        G += gap;
    }
    offs[M - 1] = G;

    size_t total = G * (size_t)C;
    // Allocate output
    unsigned char *out = malloc(total);
    if (!out) {
        fprintf(stderr, "ERROR: cannot allocate %" PRIu64 "*%d bytes = %zu\n",
                G, C, total);
        exit(1);
    }

// Generate gaps
    //#pragma omp parallel for schedule(dynamic)    // re‑enable if thread‑safe
    for (int i = 0; i + 1 < M; i++) {
        uint64_t base    = offs[i];
        size_t   offset  = base * (size_t)C;
        uint64_t ai      = A[i * LIMBS + (LIMBS - 1)];
        uint64_t bi      = A[(i+1) * LIMBS + (LIMBS - 1)];
        uint64_t gap     = (bi > ai) ? bi - ai - 1 : 0;

    // bounds‐check one last time if you like:
        if (offset + gap * (size_t)C > total) {
            fprintf(stderr, "OUT OF BOUNDS: i=%d gap=%" PRIu64 "\n", i, gap);
            exit(1);
        }

       uint64_t *cur = malloc(LIMBS * sizeof *cur);
        memcpy(cur, A + (size_t)i * LIMBS, LIMBS * sizeof *cur);

        for (uint64_t k = 0; k < gap; k++) {
            limbs_inc(cur, LIMBS);
            limbs_to_row(cur, out + offset * C, C);
            offset++;
        }

        free(cur);
    }

    // Print
    for (uint64_t i = 0, G0 = offs[M - 1]; i < G0; i++) {
        for (int j = 0; j < C; j++) putchar(out[i * C + j] ? '1' : '0');
        putchar('\n');
    }

    free(A);
    free(offs);
    free(out);
}

int main(void) {
    while (1) {
        int M, max_id, n;
        bool *seen;
        Clause *clause = read_clauses(&M, &max_id, &seen);

        int *vars       = build_vars(seen, max_id, &n);
        int *var_to_idx = calloc(max_id+1, sizeof(int));
        for (int i = 0; i < n; i++)
            var_to_idx[vars[i]] = i;

        printf("Detected %d distinct variable(s) and %d clause(s)\n\n", n, M);

        // pre-allocate shared scratch
        int  (*sorted)[n] = malloc(M * sizeof *sorted);
        unsigned char *forbidden = malloc(M * n * sizeof(unsigned char));
        bool early_unsat = false;

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < M; i++) {
            if (!early_unsat) {
                process_clause(&clause[i],
                               vars, var_to_idx, n,
                               sorted[i], forbidden,
                               &early_unsat, i);
            }
        }

        if (!early_unsat) {
            for (int i = 0; i < M; i++)
                print_clause(sorted[i], forbidden, n, i);
        }

        int LIMBS = (n + 63) / 64;
        process(M, n, LIMBS, forbidden);

        // cleanup
        for (int i = 0; i < M; i++)
            free(clause[i].orig);
        free(clause);
        free(seen);
        free(vars);
        free(var_to_idx);
        free(sorted);
        free(forbidden);
    }
    return 0;
}
