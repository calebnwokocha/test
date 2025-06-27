#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#define MAX_VARS 64
#define MAX_CLAUSES 10000
#define HT_SIZE 200003

typedef struct {
    int v[3];
    int sign[3];
} Clause;

static Clause clauses[MAX_CLAUSES];
static int N, M;
static uint64_t ht[HT_SIZE];
static uint8_t used_ht[HT_SIZE];
static uint8_t assignment[MAX_VARS];

// Simple hash: multiplication + modulo (prime table)
static inline size_t hsh(uint64_t x) {
    return (x * 11400714819323198485ULL) % HT_SIZE;
}

static void forbid(uint64_t x) {
    size_t i = hsh(x);
    while (used_ht[i] && ht[i] != x) i = (i + 1) % HT_SIZE;
    used_ht[i] = 1;
    ht[i] = x;
}

static int is_forbidden(uint64_t x) {
    size_t i = hsh(x);
    while (used_ht[i]) {
        if (ht[i] == x) return 1;
        i = (i + 1) % HT_SIZE;
    }
    return 0;
}

// Read input: N variables, M clauses, each clause with 3 literals
static int read_instance(void) {
    if (scanf("%d %d", &N, &M) != 2) return 0;
    if (N < 1 || N > MAX_VARS) {
        fprintf(stderr, "Error: variable count must be 1–%d.\n", MAX_VARS);
        return 0;
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < 3; j++) {
            int lit;
            if (scanf("%d", &lit) != 1) return 0;
            clauses[i].v[j] = abs(lit) - 1;
            clauses[i].sign[j] = (lit > 0) ? +1 : -1;
        }
    }
    return 1;
}

// Compute all forbidden full assignments
static void build_forbidden_table(void) {
    for (int i = 0; i < M; i++) {
        uint64_t mask = 0;
        for (int j = 0; j < 3; j++) {
            int idx = clauses[i].v[j];
            if (clauses[i].sign[j] == -1)
                mask |= 1ULL << idx;
        }
        forbid(mask);
    }
}

// Sample a random full assignment not in the forbidden set
static void sample_assignment(void) {
    srand((unsigned)time(NULL));
    while (1) {
        uint64_t bits = ((uint64_t)rand() << 32) ^ rand();
        bits &= (N == 64 ? (uint64_t)-1 : ((1ULL << N) - 1));
        if (is_forbidden(bits)) continue;
        for (int i = 0; i < N; i++)
            assignment[i] = (bits >> i) & 1;
        return;
    }
}

int main(void) {
    if (!read_instance()) return 1;
    build_forbidden_table();
    sample_assignment();
    for (int i = 0; i < N; i++) printf("%u ", assignment[i]);
    printf("\n");
    return 0;
}
