#include "lookup_table.h"

#include "../corelib/ksort.h"
#include "hash_list_bucket_sort.h"

typedef struct {
    u64 hash;
    i64 offset;
} KmerHashAndOffset;

static u64
calc_num_kmers(const text_t* text,
	const int kmer_size,
	const int window_size)
{
	u64 nkm = 0;
	const int num_seqs = seqdb_num_seqs(text);
	for (int i = 0; i < num_seqs; ++i) {
		u64 size = seqdb_seq_size(text, i);
		if (size < kmer_size) continue;
		u64 n = (size - kmer_size) / window_size + 1;
		nkm += n;
	}	
	return nkm;
}

static KmerHashAndOffset*
get_khao_array(const text_t* db,
    const int kmer_size,
    const int window_size,
    u64* khao_count)
{
    hbn_timing_begin(__FUNCTION__);

    //hbn_assert(db->unpacked_seq != NULL);
    u64 num_kmers = calc_num_kmers(db, kmer_size, window_size);
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)calloc(num_kmers, sizeof(KmerHashAndOffset));
    const int kIntersect = kmer_size > window_size;
    const int kStride = kmer_size - window_size;
    const u64 kIntersectMask = kIntersect ? ((U64_ONE << (kStride<<1)) - 1) : 0;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1);
    const int num_subjects = seqdb_num_seqs(db);
    size_t cnt = 0;

    for (int i = 0; i < num_subjects; ++i) {
        const u64 subject_size = seqdb_seq_size(db, i);
        if (subject_size < kmer_size) continue;
        const u64 start = seqdb_seq_offset(db, i);

        if (!kIntersect) {
            for (u64 j = 0; j <= subject_size - kmer_size; j += window_size) {
                u64 hash = 0;
                for (int k = 0; k < kmer_size; ++k) {
                    const u64 pos = start + j + k;
                    //const u8 c = db->unpacked_seq[pos];
					const u8 c = _get_pac(db->packed_seq, pos);
                    hash = (hash << 2) | c;
                }
                hbn_assert(hash < kMaxHashValue);
                KmerHashAndOffset khao = { hash, start + j };
                khao_array[cnt++] = khao;
            }
        } else {
            u64 hash = 0;
            for (int j = 0; j < kmer_size; ++j) {
                const u64 pos = start + j;
                //const u8 c = db->unpacked_seq[pos];
				const u8 c = _get_pac(db->packed_seq, pos);
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            KmerHashAndOffset khao = { hash, start };
            khao_array[cnt++] = khao;
            for (u64 j = window_size; j <= subject_size - kmer_size; j += window_size) {
                hash &= kIntersectMask;
                for (int k = kStride; k < kmer_size; ++k) {
                    const u64 pos = start + j + k;
                    //const u8 c = db->unpacked_seq[pos];
					const u8 c = _get_pac(db->packed_seq, pos);
                    hash = (hash << 2) | c;
                }
                hbn_assert(hash < kMaxHashValue);
                KmerHashAndOffset khao = { hash, start + j };
                khao_array[cnt++] = khao;
            }
        }
    }
    hbn_assert(cnt == num_kmers);
    *khao_count = num_kmers;
    hbn_timing_end(__FUNCTION__);
    return khao_array;
}

static u64
remove_repetitive_kmers(KmerHashAndOffset* khao_array,
    const u64 khao_count,
    const int max_kmer_occ)
{
    u64 distinct_kmers = 0;
    u64 removed_distinct_kmers = 0;
    u64 removed_kemrs = 0;

    u64 i = 0;
    while (i < khao_count) {
        ++distinct_kmers;
        u64 j = i + 1;
        while (j < khao_count && khao_array[i].hash == khao_array[j].hash) ++j;
        u64 n = j - i;
        if (n > max_kmer_occ) {
            ++removed_distinct_kmers;
            removed_kemrs += n;
            for (u64 k = i; k < j; ++k) khao_array[k].offset = U64_MAX;
        }
        i = j;
    }

    {
        char buf1[64], buf2[64], buf3[64];
        u64_to_string_comma(khao_count, buf1);
        u64_to_string_comma(removed_kemrs, buf2);
        double perc = 100.0 * removed_kemrs / khao_count;
        double_to_string(perc, buf3);
        HBN_LOG("Total kmers: %s, %s (%s%%) are filtered out.", buf1, buf2, buf3);

        u64_to_string_comma(distinct_kmers, buf1);
        u64_to_string_comma(removed_distinct_kmers, buf2);
        perc = 100.0 * removed_distinct_kmers / distinct_kmers;
        double_to_string(perc, buf3);
        HBN_LOG("Distinct kmers: %s, %s (%s%%) are fltered out.", buf1, buf2, buf3);
    }

    i = 0;
    for (u64 k = 0; k < khao_count; ++k) {
        if (khao_array[k].offset != U64_MAX) khao_array[i++] = khao_array[k];
    }
    return i;
}

static LookupTable*
build_lktbl_array_from_khao_array(KmerHashAndOffset* khao_array, u64 khao_count)
{
	u64* offset_array = (u64*)calloc(khao_count, sizeof(u64));
    u64 i = 0;
	u64 distinct_kmers = 0;
    while (i < khao_count) {
        u64 j = i + 1;
        while (j < khao_count && khao_array[i].hash == khao_array[j].hash) ++j;
		for (u64 k = i; k < j; ++k) offset_array[k] = khao_array[k].offset;
		++distinct_kmers;
        i = j;
    }

    LktblKmerInfo* lki_array = (LktblKmerInfo*)calloc(distinct_kmers, sizeof(LktblKmerInfo));
	u64 lki_idx = 0;
	i = 0;
	while (i < khao_count) {
		u64 j = i + 1;
		while (j < khao_count && khao_array[i].hash == khao_array[j].hash) ++j;
		LktblKmerInfo lki;
		lki.hash = khao_array[i].hash;
		lki.offset_list_idx = i;
		lki.cnt = j - i;
		lki_array[lki_idx++] = lki;
		i = j;
	}
	hbn_assert(lki_idx == distinct_kmers);

	LookupTable* lktbl = (LookupTable*)calloc(1, sizeof(LookupTable));
	lktbl->offset_list = offset_array;
	lktbl->offset_list_size = khao_count;
	lktbl->lki_array = lki_array;
	lktbl->lki_count = distinct_kmers;
	return lktbl;
}

static u64 
hash_extractor(void* list, const u64 i)
{
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)(list);
    return khao_array[i].hash;
}

static u64 
offset_extractor(void* list, const u64 i)
{
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)(list);
    return khao_array[i].offset;
}

static void
set_khao_array_item_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
    KmerHashAndOffset* src_khao_array = (KmerHashAndOffset*)(src);
    KmerHashAndOffset* dst_khao_array = (KmerHashAndOffset*)(dst);
    dst_khao_array[dst_idx] = src_khao_array[src_idx];
}

LookupTable*
build_lookup_table(const text_t* db,
    const int kmer_size,
    const int window_size,
    const int max_kmer_occ,
    const int num_threads)
{
    u64 khao_count = 0;
    KmerHashAndOffset* khao_array = get_khao_array(db, kmer_size, window_size, &khao_count);
    radix_sort(khao_array, 
        sizeof(KmerHashAndOffset), 
        khao_count, 
        num_threads,
        offset_extractor,
        hash_extractor,
        set_khao_array_item_value);
    HBN_LOG("verifying hash order of %zu elements", khao_count);
    for (u64 i = 0; i < khao_count - 1; ++i) {
        hbn_assert(khao_array[i].hash <= khao_array[i+1].hash, "i = %zu, %zu --- %zu", 
            i, khao_array[i].hash, khao_array[i+1].hash);
    }
    khao_count = remove_repetitive_kmers(khao_array, khao_count, max_kmer_occ);
	LookupTable* lktbl = build_lktbl_array_from_khao_array(khao_array,khao_count);
	free(khao_array);
	return lktbl;
}

LookupTable*
LookupTableFree(LookupTable* lktbl)
{
	free(lktbl->offset_list);
	free(lktbl->lki_array);
	free(lktbl);
	return NULL;
}