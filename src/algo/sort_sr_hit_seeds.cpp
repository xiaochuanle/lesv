#include "sort_sr_hit_seeds.h"

#include <boost/algorithm/sorting/spread_sort.hpp>

using namespace boost;

/////////////////

struct DDFKmerMatch_SoffLessThan {
    inline bool operator()(const DDFKmerMatch &x, const DDFKmerMatch &y) const { return x.soff < y.soff; }
};

struct DDFKmerMatch_SoffRightShift {
        inline idx operator()(const DDFKmerMatch &x, const unsigned offset) { return x.soff >> offset; }
};

extern "C"
void
sort_ddfkm_soff_lt(size_t n, DDFKmerMatch* a)
{
    integer_sort(a, a + n, DDFKmerMatch_SoffRightShift(), DDFKmerMatch_SoffLessThan());
}

////////////////

struct DDFKmerMatch_ContextLessThan {
    inline bool operator()(const DDFKmerMatch &x, const DDFKmerMatch &y) const { return x.context < y.context; }
};

struct DDFKmerMatch_ContextRightShift {
        inline int operator()(const DDFKmerMatch &x, const unsigned offset) { return x.context >> offset; }
};

extern "C"
void
sort_ddfkm_context_lt(size_t n, DDFKmerMatch* a)
{
    integer_sort(a, a + n, DDFKmerMatch_ContextRightShift(), DDFKmerMatch_ContextLessThan());
}

////////////////

struct DDFKmer_HashLessThan {
    inline bool operator()(const DDFKmer& x, const DDFKmer& y) const { return x.hash < y.hash; }
};

struct DDFKmer_HashRightShift {
    inline u64 operator()(const DDFKmer &x, const unsigned offset) { return x.hash >> offset; }
};

extern "C"
void
sort_ddfk_hash_lt(size_t n, DDFKmer* a)
{
    integer_sort(a, a + n, DDFKmer_HashRightShift(), DDFKmer_HashLessThan());
}


/////////////////////

struct ChainSeed_SoffLessThan {
    inline bool operator()(const ChainSeed& x, const ChainSeed& y) const { return x.soff < y.soff; }
};

struct ChainSeed_SoffRightShift {
    inline idx operator()(const ChainSeed &x, const unsigned offset) { return x.soff >> offset; }
};

extern "C"
void
sort_chain_seed_soff_lt(size_t n, ChainSeed* a)
{
    integer_sort(a, a + n, ChainSeed_SoffRightShift(), ChainSeed_SoffLessThan());
}