#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

typedef struct {
    int *lit;
    int sz;
} Clause;

int cmp_abs(const void *a, const void *b) {
    int x = abs(*(int*)a), y = abs(*(int*)b);
    return (x > y) - (x < y);
}

Clause* read_clauses(int **vars_out, int *n_out, int *m_out) {
    Clause *C = NULL;
    int M = 0, cap = 0, max_var = 0;
    int *var_seen = calloc(1024, sizeof(int)); // dynamic tracking

    printf("Enter clauses (blank line to finish):\n");
    while (1) {
        char line[1024];
        if (!fgets(line, sizeof(line), stdin) || line[0]=='\n') break;

        if (M+1 > cap) {
            cap = cap ? cap*2 : 4;
            C = realloc(C, cap * sizeof *C);
        }

        Clause *c = &C[M];
        c->lit = malloc(8 * sizeof(int));
        c->sz = 0;
        char *tok = strtok(line, " \t\n");
        while (tok) {
            int lit = atoi(tok);
            if (c->sz % 8 == 7) // expand every 8
                c->lit = realloc(c->lit, (c->sz+8)*sizeof(int));
            c->lit[c->sz++] = lit;

            int av = abs(lit);
            if (av > max_var) max_var = av;
            if (av >= 0 && av < 1024) var_seen[av] = 1;

            tok = strtok(NULL, " \t\n");
        }
        M++;
    }

    // build distinct variable list
    int total = 0;
    for (int v = 1; v <= max_var; v++)
        if (v < 1024 && var_seen[v]) total++;
    int *vars = malloc(total * sizeof(int));
    for (int v = 1, i = 0; v <= max_var; v++)
        if (v < 1024 && var_seen[v]) vars[i++] = v;
    free(var_seen);

    *vars_out = vars;
    *n_out = total;
    *m_out = M;
    return C;
}

void pad_clause(Clause *c, int n, int *vars) {
    bool seen[n];
    memset(seen, 0, n);
    for (int i = 0; i < c->sz; i++)
        for (int j = 0; j < n; j++)
            if (abs(c->lit[i]) == vars[j]) seen[j] = true;

    for (int j = 0; j < n; j++) {
        if (!seen[j]) {
            c->lit[c->sz++] = vars[j];
        }
    }
}

bool* build_forbidden(Clause *c, int n, int *vars) {
    bool *a = calloc(n+1, sizeof(bool));
    for (int i = 0; i < c->sz; i++) {
        int lit = c->lit[i];
        for (int j = 0; j < n; j++) {
            if (abs(lit) == vars[j]) {
                a[j+1] = (lit < 0);
                break;
            }
        }
    }
    return a;
}

// After reading and before padding/sorting, do:
bool check_unsat_clause(int *lits, int sz) {
    // Build a small hash or boolean array indexed by the actual literal values
    // For simplicity, use dynamic range based on max absolute.
    for (int i = 0; i < sz; i++) {
        for (int j = i + 1; j < sz; j++) {
            if (lits[i] + lits[j] == 0) return true; // found x and -x
        }
    }
    return false;
}

int main(void) {
    int *vars, n, m;
    Clause *C = read_clauses(&vars, &n, &m);
    printf("Detected %d distinct variables, %d clauses\n", n, m);

    bool **forbidden = malloc(m * sizeof(bool*));
    for (int i = 0; i < m; i++) {
        if (check_unsat_clause(C[i].lit, C[i].sz)) {
            printf("UNSAT (clause %d contains a variable and its negation)\n", i + 1);
            return 0; // terminate early
        }
        pad_clause(&C[i], n, vars);
        qsort(C[i].lit, C[i].sz, sizeof(int), cmp_abs);
        forbidden[i] = build_forbidden(&C[i], n, vars);
    }

    for (int i = 0; i < m; i++) {
        // Print the padded and sorted clause
        printf("Clause %d padded and sorted literals:", i + 1);
        for (int j = 0; j < C[i].sz; j++) {
            printf(" %d", C[i].lit[j]);
        }
        printf("\n");

        // Print the forbidden assignment vector
        printf("Clause %d forbidden a[1..%d]:", i + 1, n);
        for (int j = 1; j <= n; j++) {
            printf(" %d", forbidden[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < m; i++) {
        free(C[i].lit);
        free(forbidden[i]);
    }
    free(forbidden);
    free(C);
    free(vars);
    return 0;
}
