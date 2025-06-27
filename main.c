#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
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
                           int *out_sorted, bool *out_forbid,
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
            out_forbid[j]  = true;
        } else {
            // either seen_pos[j] or padded positive
            out_sorted[j]  =  vars[j];
            out_forbid[j]  = false;
        }
    }
}

// Phase 4: print one clause
static void print_clause(const int *sorted, const bool *forbid, int n, int idx) {
    printf("Clause %d (padded and sorted):", idx+1);
    for (int j = 0; j < n; j++)
        printf(" %d", sorted[j]);
    printf("\nClause %d forbidden a[1..%d]:", idx+1, n);
    for (int j = 0; j < n; j++)
        printf(" %d", forbid[j] ? 1 : 0);
    printf("\n\n");
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
        bool (*forbid)[n] = malloc(M * sizeof *forbid);
        bool early_unsat = false;

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < M; i++) {
            if (!early_unsat) {
                process_clause(&C[i],
                               vars, var_to_idx, n,
                               sorted[i], forbid[i],
                               &early_unsat, i);
            }
        }

        if (!early_unsat) {
            for (int i = 0; i < M; i++)
                print_clause(sorted[i], forbid[i], n, i);
        }

        // cleanup
        for (int i = 0; i < M; i++)
            free(C[i].orig);
        free(C);
        free(seen);
        free(vars);
        free(var_to_idx);
        free(sorted);
        free(forbid);
    }
    return 0;
}
