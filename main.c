#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

//
// Data structures and hashing utilities
//

// Clause of three literals in {-n..-1} ∪ {1..n}
typedef struct { int lit[3]; } Clause;

// Open‐address hash table for 64-bit keys
typedef struct {
    uint64_t *keys;
    uint8_t  *used;
    size_t    mask;
} OAHTable;

// Initialize to capacity = next power of two ≥ cap
static void oah_init(OAHTable *t, size_t cap){
    size_t s=1; while(s<cap) s<<=1;
    t->keys = calloc(s, sizeof *t->keys);
    t->used = calloc(s, sizeof *t->used);
    t->mask = s-1;
}

// Simple 64-bit mix hash
static inline size_t oah_hash(OAHTable *t, uint64_t x){
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    return (size_t)x & t->mask;
}

// Insert x; return 1 if already present, 0 if newly added
static int oah_insert(OAHTable *t, uint64_t x){
    size_t i = oah_hash(t,x);
    while(t->used[i]){
        if(t->keys[i]==x) return 1;
        i = (i+1) & t->mask;
    }
    t->used[i]=1; t->keys[i]=x;
    return 0;
}

//
// Global problem data
//
static Clause *clauses;
static int     n_vars, n_clauses;
static OAHTable forbid_tbl, seen_tbl;
static char   *assignment;

// Read n, m, then m clauses of three ints each
static int read_instance(){
    if(scanf("%d %d",&n_vars,&n_clauses)!=2) return 0;
    clauses   = malloc(sizeof(Clause)*n_clauses);
    assignment= malloc(n_vars);
    for(int i=0;i<n_clauses;i++){
        if(scanf("%d %d %d",
           &clauses[i].lit[0],
           &clauses[i].lit[1],
           &clauses[i].lit[2])!=3) return 0;
    }
    return 1;
}

// Build forbidden mask table: one mask per clause
static void build_forbid(){
    oah_init(&forbid_tbl, n_clauses*2);
    for(int i=0;i<n_clauses;i++){
        uint64_t mask=0;
        for(int j=0;j<3;j++){
            int v = abs(clauses[i].lit[j]) - 1;
            // literal false when bit = opposite of sign:
            if(clauses[i].lit[j]<0) mask |= (1ULL<<v);
        }
        oah_insert(&forbid_tbl, mask);
    }
}

// Evaluate clause satisfaction
static int sat_clause(const Clause *c){
    for(int j=0;j<3;j++){
        int v = abs(c->lit[j]) - 1;
        int val = assignment[v];
        if((c->lit[j]>0 && val) || (c->lit[j]<0 && !val))
            return 1;
    }
    return 0;
}

// Main solver combining all requirements
static int solve(){
    // allocate seen table for up to min(2^n, 1e6) samples
    double total = n_vars<63 ? pow(2.0,n_vars) : 1e6;
    oah_init(&seen_tbl, (size_t)fmin(total,1e6));

    int max_restarts = 3*n_vars, max_flips = 3*n_vars;
    uint64_t full_mask = n_vars==64 ? ~0ULL : ((1ULL<<n_vars)-1);

    for(int r=0;r<max_restarts;r++){
        // sample until fresh and not forbidden
        uint64_t x;
        do{
            x = ((uint64_t)rand()<<32) ^ rand();
            x &= full_mask;
        } while(oah_insert(&seen_tbl,x) || oah_insert(&forbid_tbl,x));
        // initialize assignment
        for(int i=0;i<n_vars;i++)
            assignment[i] = (x>>i)&1;

        // local flips per Schöning
        for(int f=0;f<max_flips;f++){
            int uns=-1;
            for(int i=0;i<n_clauses;i++){
                if(!sat_clause(&clauses[i])){ uns=i; break; }
            }
            if(uns<0) return 1;  // found solution

            // flip a random literal's variable
            int j = rand()%3;
            int v = abs(clauses[uns].lit[j]) - 1;
            assignment[v] ^= 1;
            x ^= (1ULL<<v);

            // if seen or forbidden, skip further flips
            if(oah_insert(&seen_tbl,x) || oah_insert(&forbid_tbl,x)){
                if(seen_tbl.mask+1 >= full_mask+1) return 0;
                break;
            }
        }
    }
    return 0;
}

int main(){
    srand(time(NULL));
    if(!read_instance()){ fprintf(stderr,"ERROR\n"); return 1; }
    build_forbid();
    if(solve()){
        for(int i=0;i<n_vars;i++) printf("%d ",assignment[i]);
        printf("\n");
    } else {
        printf("UNSAT\n");
    }
    return 0;
}
