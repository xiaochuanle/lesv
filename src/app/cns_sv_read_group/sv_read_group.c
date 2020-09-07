#include "sv_read_group.h"

#include "../necat2sv/sv_signature.h"
#include "../../corelib/khash.h"
#include "../../corelib/ksort.h"

KHASH_SET_INIT_INT(int_set);

static void
s_truncate_back_space_in_line(char* line)
{
    char* p = line;
    while (1) {
        int c = *p;
        if (c == '\0') break;
        if (c == '\r') {
            *p = '\0';
            break;
        }
        if (c == '\n') {
            *p = '\0';
            break;
        }
        ++p;
    }
}

typedef struct {
    int i;
    int len;
    int score;
} SvSigInfo;

#define SvSigInfo_score_gt(a, b) ((a).score > (b).score)
KSORT_INIT(SvSigInfo_score_gt, SvSigInfo, SvSigInfo_score_gt);

static void
s_filter_outlier_svsig(const int group_id, vec_sv_sig* svsig_list)
{
    const int kMaxSvsigCnt = 50;
    if (kv_size(*svsig_list) <= kMaxSvsigCnt) return;
    int sigc = kv_size(*svsig_list);
    SvSignature* siga = (SvSignature*)calloc(sigc, sizeof(SvSignature));
    memcpy(siga, kv_data(*svsig_list), sigc * sizeof(SvSignature));
    SvSigInfo* sigia = (SvSigInfo*)calloc(sigc, sizeof(SvSigInfo));
    for (int i = 0; i < sigc; ++i) {
        if (siga[i].type == eGapAlignIns) {
            sigia[i].len = siga[i].qto - siga[i].qfrom;
        } else {
            sigia[i].len = siga[i].sto - siga[i].sfrom;
        }
        sigia[i].i = i;
        sigia[i].score = 0;
    }

    for (int i = 0; i < sigc; ++i) {
        for (int j = i + 1; j < sigc; ++j) {
            int li = sigia[i].len;
            int lj = sigia[j].len;
            hbn_assert(li > 1);
            hbn_assert(lj > 1);
            int ml = hbn_max(li, lj);
            int ms = hbn_min(li, lj);
            if (ml - ms <= ml * 0.2) {
                sigia[i].score++;
                sigia[j].score++;
            }
        }
    }

    ks_introsort_SvSigInfo_score_gt(sigc, sigia);
#if 0
    HBN_LOG("group %d, number of sigs: %d", group_id, sigc);
    for (int i = 0; i < sigc; ++i) {
        int x = sigia[i].i;
        dump_svsig(fprintf, stderr, siga[x], NULL);
        fprintf(stderr, "score = %d,\tlen = %d\n", sigia[i].score, sigia[i].len);
    }
#endif

    kv_clear(*svsig_list);
    for (int i = 0; i < kMaxSvsigCnt; ++i) {
        int x = sigia[i].i;
        kv_push(SvSignature, *svsig_list, siga[x]);
    }
    sfree(siga);
    sfree(sigia);
}

static int
s_load_next_sv_read_group(FILE* in, vec_sv_sig* svsig_list, int* group_id_)
{
    int group_id, cnt;
    char line[HBN_MAX_PATH_LEN];
    if (!fgets(line, HBN_MAX_PATH_LEN, in)) return 0;
    HBN_SCANF(sscanf, line, 2, "%d%d", &group_id, &cnt);
    SvSignature svsig;
    kv_clear(*svsig_list);
    khash_t(int_set)* added_ids = kh_init(int_set);
    for (int i = 0; i < cnt; ++i) {
        if (!fgets(line, HBN_MAX_PATH_LEN, in)) {
            HBN_ERR("unexpected end of file in reading group %d, only %d of %d signatures are read",
                group_id, i+1, cnt);            
        }
        s_truncate_back_space_in_line(line);
        sread_sv_signature(line, &svsig);
        khiter_t pos = kh_get(int_set, added_ids, svsig.qid);
        if (pos == kh_end(added_ids)) {
            kv_push(SvSignature, *svsig_list, svsig);
            int r = 0;
            kh_put(int_set, added_ids, svsig.qid, &r);
            hbn_assert(r == 1);
        }
    }
    *group_id_ = group_id;
    kh_destroy(int_set, added_ids);
    s_filter_outlier_svsig(group_id, svsig_list);
    return 1;
}

SvReadGroup*
SvReadGroupLoad(const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    const int subject_id,
    FILE* sv_read_group_file,
    FILE* pac_file)
{
    kv_dinit(vec_sv_sig, svsig_list);
    int group_id = -1;
    if (!s_load_next_sv_read_group(sv_read_group_file, &svsig_list, &group_id)) return NULL;

    SvReadGroup* group = (SvReadGroup*)calloc(1, sizeof(SvReadGroup));
    group->group_id = group_id;
    group->subject_id = subject_id;
    group->sv_type = kv_front(svsig_list).type;
    group->raw_read_info_array = raw_read_info_array;
    group->raw_read_name_array = raw_read_name_array;
    group->sv_read_group_size = kv_size(svsig_list);
    group->svri_array = (SvReadInfo*)calloc(group->sv_read_group_size, sizeof(SvReadInfo));
    kv_init(group->fwd_sv_read);
    kv_init(group->rev_sv_read);
    size_t res = 0;
    size_t max_res = 0;
    for (size_t i = 0; i < group->sv_read_group_size; ++i) {
        SvSignature svsig = kv_A(svsig_list, i);
        int id = svsig.qid;
        int size = raw_read_info_array[id].seq_size;
        SvReadInfo* svri = group->svri_array + i;
        svri->global_id = id;
        svri->local_id = i;
        svri->name = raw_read_name_array + raw_read_info_array[id].hdr_offset;
        svri->size = size;
        svri->offset_in_sv_read = res;
        svri->raw_seq_size = size;
        svri->raw_seq_from = 0;
        svri->raw_seq_to = size;
        svri->fsfrom = svsig.fsfrom;
        svri->fsto = svsig.fsto;
        svri->fsqdir = svsig.qdir;
        res += size;
        max_res = hbn_max(size, max_res);
    }
    kv_resize(u8, group->fwd_sv_read, res);
    kv_resize(u8, group->rev_sv_read, res);
    kv_fill(group->fwd_sv_read, 0);
    kv_fill(group->rev_sv_read, 0);
    
    size_t n = (max_res + 3) / 4;
    u8* pac = (u8*)calloc(n, sizeof(u8));
    for (int i = 0; i < group->sv_read_group_size; ++i) {
        int global_id = group->svri_array[i].global_id;
        hbn_assert(group->svri_array[i].local_id == i);
        hbn_assert(group->svri_array[i].size == raw_read_info_array[global_id].seq_size);
        size_t offset = raw_read_info_array[global_id].seq_offset;
        hbn_assert((offset%4) == 0);
        offset /= 4;
        size_t size = group->svri_array[i].size;
        n = (size+3)/4;
        fseek(pac_file, offset, SEEK_SET);
        hbn_fread(pac, 1, n, pac_file);

        u8* fwd = kv_data(group->fwd_sv_read) + group->svri_array[i].offset_in_sv_read;
        u8* rev = kv_data(group->rev_sv_read) + group->svri_array[i].offset_in_sv_read + size - 1;
        for (size_t p = 0; p < size; ++p) {
            u8 c = _get_pac(pac, p);
            u8 rc = 3 - c;
            *fwd =  c; ++fwd;
            *rev = rc; --rev;
        }
        memset(pac, 0, n);
    }
    kv_destroy(svsig_list);
    free(pac);
    return group;
}

SvReadGroup*
SvReadGroupNew(const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    const int group_id,
    const int subject_id,
    EGapAlignOpType sv_type,
    const int num_reads)
{
    SvReadGroup* group = (SvReadGroup*)calloc(1, sizeof(SvReadGroup));
    group->group_id = group_id;
    group->subject_id = subject_id;
    group->sv_type = sv_type;
    group->raw_read_info_array = raw_read_info_array;
    group->raw_read_name_array = raw_read_name_array;
    group->sv_read_group_size = num_reads;
    group->svri_array = (SvReadInfo*)calloc(num_reads, sizeof(SvReadInfo));
    kv_init(group->fwd_sv_read);
    kv_init(group->rev_sv_read);
    return group;
}

SvReadGroup*
SvReadGroupFree(SvReadGroup* group)
{
    if (!group) return NULL;
    sfree(group->svri_array);
    kv_destroy(group->fwd_sv_read);
    kv_destroy(group->rev_sv_read);
    sfree(group);
    return NULL;
}

const u8* 
SvReadGroup_ExtractSeq(SvReadGroup* group, const int local_read_id, const int strand)
{
    hbn_assert(local_read_id < group->sv_read_group_size);
    size_t offset = group->svri_array[local_read_id].offset_in_sv_read;
    const u8* seq = (strand == FWD) ? kv_data(group->fwd_sv_read) : kv_data(group->rev_sv_read);
    return seq + offset;
}

int
SvReadGroup_ReadSize(SvReadGroup* group, const int local_read_id)
{
    hbn_assert(local_read_id < group->sv_read_group_size);
    return group->svri_array[local_read_id].size;
}

int 
SvReadGroup_GlobalReadId(SvReadGroup* group, const int local_read_id)
{
    hbn_assert(local_read_id < group->sv_read_group_size);
    return group->svri_array[local_read_id].global_id;    
}

const char*
SvReadGroup_ReadName(SvReadGroup* group, const int local_read_id)
{
    hbn_assert(local_read_id < group->sv_read_group_size);
    const char* name = group->svri_array[local_read_id].name;
    return name;
}