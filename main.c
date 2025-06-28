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

// Multiply decimal‐string `s` by 2 in place.
//   s is e.g. "12345" → becomes "24690"
static void dec_mul2(char *s) {
    int carry = 0;
    for (char *p = s; *p; p++) {
        int d = (*p - '0')*2 + carry;
        *p  = '0' + (d % 10);
        carry = d / 10;
    }
    if (carry) {                          // overflow a digit
        size_t len = strlen(s);
        memmove(s+1, s, len+1);
        s[0] = '0' + carry;
    }
}

// Add a single bit 0 or 1 to decimal‐string `s`
//   s is "24690", bit=1 → becomes "24691"
static void dec_add1(char *s, int bit) {
    int carry = bit;
    for (char *p = s + strlen(s) - 1; p >= s && carry; p--) {
        int d = (*p - '0') + carry;
        *p  = '0' + (d % 10);
        carry = d / 10;
    }
    if (carry) {
        size_t len = strlen(s);
        memmove(s+1, s, len+1);
        s[0] = '0' + carry;
    }
}

// Convert a C‑bit binary row[0..C-1] into a decimal string.
// Returns a malloc’d char* which must be free’d by caller.
static char* binrow_to_dec(const unsigned char *row, int C) {
    // worst‐case decimal digits ≈ C*log10(2) + 2
    int bufsz = (int)(C * 0.30103) + 3;
    char *s = calloc(bufsz,1);
    strcpy(s, "0");
    for (int j = 0; j < C; j++) {
        dec_mul2(s);
        if (row[j]) dec_add1(s, 1);
    }
    return s;
}

// Compute dst = b − a − 1, where b > a.
// All strings are decimal, no leading zeros (except "0").
// Returns a malloc’d string in dst (must free).
static char* dec_gap(const char *b, const char *a) {
    // First compute temp = b − a
    int lb = strlen(b), la = strlen(a);
    int maxl = lb > la ? lb : la;
    char *temp = calloc(maxl+2,1);
    int borrow = 0;
    // right‑align both
    for (int i = 0; i < maxl; i++) {
        int db = (i < lb ? b[lb-1-i]-'0' : 0);
        int da = (i < la ? a[la-1-i]-'0' : 0);
        int d  = db - da - borrow;
        if (d < 0) { d += 10; borrow = 1; }
        else       borrow = 0;
        temp[maxl-1-i] = '0' + d;
    }
    // strip leading zeros
    char *p = temp;
    while (*p == '0' && *(p+1)) p++;
    memmove(temp, p, strlen(p)+1);

    // now subtract 1
    borrow = 1;
    for (int i = strlen(temp)-1; i >= 0 && borrow; i--) {
        int d = (temp[i]-'0') - borrow;
        if (d < 0) { d += 10; borrow = 1; }
        else       borrow = 0;
        temp[i] = '0' + d;
    }
    // strip leading zero if any
    p = temp;
    while (*p == '0' && *(p+1)) p++;
    memmove(temp, p, strlen(p)+1);
    return temp;
}

// Compare two decimal strings a and b (no leading zeros)
// returns −1,0,+1
static int dec_cmp(const char *a, const char *b) {
    int la = strlen(a), lb = strlen(b);
    if (la < lb) return -1;
    if (la > lb) return  1;
    return strcmp(a,b);
}

// Add 1 to decimal string in place:
static void dec_incr(char *s) { dec_add1(s, 1); }

// In‐place divide decimal string by 2, returning remainder 0 or 1
static int dec_div2(char *s) {
    int carry = 0;
    for (char *p = s; *p; p++) {
        int d = carry*10 + (*p - '0');
        *p = '0' + (d / 2);
        carry = d % 2;
    }
    // strip leading zeros
    char *q = s;
    while (*q == '0' && *(q+1)) q++;
    if (q != s) memmove(s, q, strlen(q)+1);
    return carry;
}

// Convert decimal string `s` into a C‑bit row[0..C-1] (array of 0/1 unsigned char).
//   `row[0]` becomes the most‐significant bit (leftmost).
static void dec_to_binrow(const char *s_in, unsigned char *row, int C) {
    // make a mutable copy
    char *s = strdup(s_in);
    // extract bits from LSB→MSB into a temp buffer
    unsigned char *tmp = malloc(C);
    for (int j = C-1; j >= 0; j--) {
        tmp[j] = dec_div2(s);
    }
    free(s);
    // copy into row
    memcpy(row, tmp, C);
    free(tmp);
}

// wrapper for qsort: compare two char* by numeric value
static int dec_rows_cmp(const void *pa, const void *pb) {
    const char *a = *(const char**)pa;
    const char *b = *(const char**)pb;
    return dec_cmp(a, b);
}

// comparator for qsort on uint64_t
static int uint64_cmp(const void *pa, const void *pb) {
    uint64_t a = *(const uint64_t*)pa;
    uint64_t b = *(const uint64_t*)pb;
    if (a < b) return -1;
    if (a > b) return  1;
    return 0;
}

static void process(int M, int C, const unsigned char *forbidden) {
    if (C <= 64) {
        // ───── Fast 64‑bit path ─────
        uint64_t *vals = malloc(M * sizeof *vals);
        for (int i = 0; i < M; i++) {
            uint64_t v = 0;
            for (int j = 0; j < C; j++)
                v = (v << 1) | (forbidden[i*C + j] & 1);
            vals[i] = v;
        }
        qsort(vals, M, sizeof *vals, uint64_cmp);

        // enumerate gaps
        for (int i = 0; i + 1 < M; i++) {
            uint64_t lo = vals[i], hi = vals[i+1];
            for (uint64_t x = lo + 1; x < hi; x++) {
                for (int j = 0; j < C; j++)
                    putchar(((x >> (C-1-j)) & 1) ? '1' : '0');
                putchar('\n');
            }
        }
        free(vals);

    } else {
        // ── Decimal‑string path (handles any C) ──
        // binrow_to_dec, dec_gap, dec_cmp, dec_incr, dec_to_binrow, dec_rows_cmp
        char **dec_rows = malloc(M * sizeof(char*));
        for (int i = 0; i < M; i++)
            dec_rows[i] = binrow_to_dec(forbidden + i*C, C);

        qsort(dec_rows, M, sizeof(char*), dec_rows_cmp);

        for (int i = 0; i + 1 < M; i++) {
            char *low  = dec_rows[i];
            char *high = dec_rows[i+1];

            // compute gap = high − low − 1
            char *gapstr = dec_gap(high, low);

            // enumerate low+1 … high−1
            char *cur = strdup(low);
            dec_incr(cur);
            while (dec_cmp(cur, high) < 0) {
                unsigned char row[C];
                dec_to_binrow(cur, row, C);
                for (int j = 0; j < C; j++)
                    putchar(row[j] ? '1' : '0');
                putchar('\n');
                dec_incr(cur);
            }

            free(cur);
            free(gapstr);
        }
        for (int i = 0; i < M; i++)
            free(dec_rows[i]);
        free(dec_rows);
    }
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

        process(M, n, forbidden);

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
