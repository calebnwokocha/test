#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

typedef struct {
    int *lit;
    int sz;
} Clause;

int cmp_int(const void *a, const void *b) {
    return (*(int*)a) - (*(int*)b);
}

// a) Read clauses and determine n, m, and sz
Clause* read_clauses(int *n, int *m) {
    int max_var = 0;
    int max_clause_size = 0;
    Clause *C = malloc(sizeof(*C));
    *m = 0;

    printf("Enter clauses (press Enter on an empty line to finish):\n");
    while (true) {
        C[*m].lit = malloc(10 * sizeof(int)); // initial allocation
        printf("Clause %d: enter literals (space-separated):\n", *m + 1);

        char line[1024];
        if (!fgets(line, sizeof(line), stdin) || line[0] == '\n') {
            break; // end input on empty line
        }

        int sz = 0;
        char *token = strtok(line, " \n");
        while (token) {
            C[*m].lit[sz++] = atoi(token);
            if (abs(C[*m].lit[sz - 1]) > max_var) {
                max_var = abs(C[*m].lit[sz - 1]);
            }
            token = strtok(NULL, " \n");
        }
        C[*m].sz = sz;
        if (sz > max_clause_size) {
            max_clause_size = sz;
        }
        (*m)++;
        C = realloc(C, (*m + 1) * sizeof(*C)); // reallocate for next clause
    }

    *n = max_var;
    printf("Detected n = %d variables, m = %d clauses, max clause size = %d\n", *n, *m, max_clause_size);
    return C;
}

// b) Ensure all variables appear (positively if missing)
void pad_clause(Clause *c, int n) {
    bool present[n + 1];
    for (int v = 1; v <= n; v++) present[v] = false;
    for (int i = 0; i < c->sz; i++)
        present[abs(c->lit[i])] = true;
    for (int v = 1; v <= n; v++)
        if (!present[v])
            c->lit[c->sz++] = v;
}

// c) Sort literals
void sort_clause(Clause *c) {
    qsort(c->lit, c->sz, sizeof(int), cmp_int);
}

// d) Build forbidden vector
bool* build_forbidden(Clause *c, int n) {
    bool *a = calloc(n + 1, sizeof(bool));
    for (int j = 0; j < c->sz; j++) {
        int lit = c->lit[j];
        a[abs(lit)] = (lit < 0);
    }
    return a;
}

int main(void) {
    int n, m;
    Clause *C = read_clauses(&n, &m);

    bool **forbidden = malloc(m * sizeof(bool*));
    for (int i = 0; i < m; i++) {
        pad_clause(&C[i], n);
        sort_clause(&C[i]);
        forbidden[i] = build_forbidden(&C[i], n);
    }

    for (int i = 0; i < m; i++) {
        printf("Clause %d forbidden a[1..%d]:", i + 1, n);
        for (int v = 1; v <= n; v++)
            printf(" %d", forbidden[i][v]);
        printf("\n");
    }

    for (int i = 0; i < m; i++) {
        free(C[i].lit);
        free(forbidden[i]);
    }
    free(C);
    free(forbidden);
    return 0;
}
