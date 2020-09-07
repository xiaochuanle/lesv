#include "cns_one_group.h"

#include <pthread.h>

#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/refine_align.h"
#include "../../algo/cns_extend_chain_seed_list.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/cns_read_header.h"
#include "../necat2sv/sv_read_group_file_name.h"
#include "map_results.h"

#define MAX_CNS_COV 15

static BOOL
check_ovlp_mapping_range(const int qb, const int qe, const int qs,
						 const int sb, const int se, const int ss,
						 double ratio)
{
	const int oq = qe - qb;
	const int qqs = qs * ratio;
	const int os = se - sb;
	const int qss = ss * ratio;
	return oq >= qqs || os >= qss;
}

static BOOL
extend_chain_seed_array_simple(HbnTracebackData* traceback_data,
    Ksw2Data* ksw_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const HbnProgramOptions* opts,
    const double min_perc_identity,
    const int read_id,
    const u8* read,
    const int read_length,
    const int subject_id,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* ident_perc_)
{
    int r = 0;
    r = hbn_traceback(traceback_data,
            read,
            read_length,
            subject,
            subject_length,
            chain_seed_array,
            chain_seed_count,
            opts->ovlp_cov_res,
            min_perc_identity,
            FALSE);
    //HbnTracebackDataDump(fprintf, stderr, traceback_data);
    if (!r) return FALSE;
    int qoff = traceback_data->qoff;
    int qend = traceback_data->qend;
    int soff = traceback_data->soff;
    int send = traceback_data->send;
    double ident_perc = traceback_data->ident_perc;
    if (!check_ovlp_mapping_range(qoff, qend, read_length,
            soff, send, subject_length, opts->ovlp_cov_perc / 100.0)) return FALSE;

    const char* qas = traceback_data->qas;
    const char* sas = traceback_data->sas;
    int aln_size = traceback_data->qae - qas;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            0, read, qoff, qend, traceback_data->qas,
            0, subject, soff, send, sas,
            aln_size, TRUE);

    const char* new_qas = NULL;
    const char* new_qae = NULL;
    const char* new_sas = NULL;
    const char* new_sae = NULL;

    refine_align(traceback_data->ksw,
            read,
            qoff,
            qend,
            read_length,
            subject,
            soff,
            send,
            subject_length,
            qas,
            sas,
            aln_size,
            &traceback_data->ext_qabuf,
            &traceback_data->ext_sabuf);
            qas = ks_s(traceback_data->ext_qabuf);
            sas = ks_s(traceback_data->ext_sabuf);
            aln_size = ks_size(traceback_data->ext_qabuf);

    truncate_align_bad_ends(qas, sas, aln_size, &qoff, &qend, &soff, &send, 
        &new_qas, &new_qae, &new_sas, &new_sae);
    aln_size = new_qae - new_qas;
    double eff_perc_identity = calc_effective_ident_perc(new_qas, new_sas, aln_size);
    if (eff_perc_identity < min_perc_identity) return FALSE;

    ks_clear(*qaln);
    kputsn(new_qas, aln_size, qaln);
    ks_clear(*saln);
    kputsn(new_sas, aln_size, saln);
    hbn_assert(ks_size(*qaln) == ks_size(*saln));
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        read_id, read, qoff, qend, ks_s(*qaln),
        subject_id, subject, soff, send, 
        ks_s(*saln), ks_size(*qaln), TRUE);
    *qoff_ = qoff;
    *qend_ = qend;
    *soff_ = soff;
    *send_ = send;
    *ident_perc_ = calc_ident_perc(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), NULL, NULL);
    //HBN_LOG("ident_perc = %g", *ident_perc_);
    return TRUE;       
}

static BOOL 
extend_chain_seed_array_ksw(HbnTracebackData* traceback_data,
    Ksw2Data* ksw2_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const HbnProgramOptions* opts,
    const double min_perc_identity,
    const int read_id,
    const u8* read,
    const int read_length,
    const int subject_id,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* ident_perc_)
{
    int r = 0;
    r = hbn_traceback(traceback_data,
            read,
            read_length,
            subject,
            subject_length,
            chain_seed_array,
            chain_seed_count,
            opts->ovlp_cov_res,
            min_perc_identity,
            FALSE);
    //HBN_LOG("perc_identity = %g", traceback_data->ident_perc);
    if (!r) return FALSE;
    //HbnTracebackDataDump(fprintf, stderr, traceback_data);
    int qoff = traceback_data->qoff;
    int qend = traceback_data->qend;
    int soff = traceback_data->soff;
    int send = traceback_data->send;
    double ident_perc = traceback_data->ident_perc;
    if (!check_ovlp_mapping_range(qoff, qend, read_length,
            soff, send, subject_length, opts->ovlp_cov_perc / 100.0)) return FALSE;
    
    int distance = hbn_max(qend - qoff, send - soff);
    distance = 1.2 * distance * (100.0 - ident_perc) / 100.0;
    if (distance == 0) distance = (send + qend - qoff - soff) * 0.05;
    //HBN_LOG("distance = %d", distance);

    r = nw_ksw2_extd2(ksw2_data,
            read_id,
            read,
            qoff,
            qend,
            read_length,
            subject_id,
            subject,
            soff,
            send,
            subject_length,
            opts->ovlp_cov_res,
            0.0,
            distance,
            &qoff,
            &qend,
            &soff,
            &send,
            &ident_perc,
            &traceback_data->qabuf,
            &traceback_data->sabuf);
    if (!r) {
        int distance1 = hbn_max(qend - qoff, send - soff);
        qoff = traceback_data->qoff;
        qend = traceback_data->qend;
        soff = traceback_data->soff;
        send = traceback_data->send;
        ident_perc = traceback_data->ident_perc;
        r = nw_ksw2_extd2(ksw2_data,
                read_id,
                read,
                qoff,
                qend,
                read_length,
                subject_id,
                subject,
                soff,
                send,
                subject_length,
                opts->ovlp_cov_res,
                0.0,
                distance1,
                &qoff,
                &qend,
                &soff,
                &send,
                &ident_perc,
                &traceback_data->qabuf,
                &traceback_data->sabuf);
        if (r) {
            HBN_LOG("[%d, %d] x [%d, %d] dist = %d is rescued with dist = %d, edlib_perc = %g, ksw_perc = %g",
                qoff, qend, soff, send, distance, distance1, traceback_data->ident_perc, ident_perc);
        }
    }
    if (!r) return FALSE;   

    const char* qas = NULL;
    const char* qae = NULL;
    const char* sas = NULL;
    const char* sae = NULL;
    r = truncate_align_bad_ends(ks_s(traceback_data->qabuf),
            ks_s(traceback_data->sabuf),
            ks_size(traceback_data->qabuf),
            &qoff,
            &qend,
            &soff,
            &send,
            &qas,
            &qae,
            &sas,
            &sae);
    if (!r) {
        HBN_LOG("truncate fail");
        return FALSE;
    }
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            0, read, qoff, qend, qas,
            0, subject, soff, send, sas,
            qae - qas, TRUE);
    double eff_perc_identity = calc_effective_ident_perc(qas, sas, sae - sas);
    if (eff_perc_identity < min_perc_identity) {
        HBN_LOG("eff_perc_identity = %g is too low", eff_perc_identity);
        return FALSE;
    }
    ks_clear(*qaln);
    kputsn(qas, qae - qas, qaln);
    ks_clear(*saln);
    kputsn(sas, sae - sas, saln);
    hbn_assert(ks_size(*qaln) == ks_size(*saln));
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        read_id, read, qoff, qend, ks_s(*qaln),
        subject_id, subject, soff, send, 
        ks_s(*saln), ks_size(*qaln), TRUE);
    *qoff_ = qoff;
    *qend_ = qend;
    *soff_ = soff;
    *send_ = send;
    *ident_perc_ = calc_ident_perc(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), NULL, NULL);
    return TRUE;
}

static BOOL
meap_consensus_one_segment(FCCnsData* cns_data,
    int *sfrom,
    int *sto,
    const int subject_size,
    const int min_size,
    kstring_t* cns_seq)
{
    build_backbone(kv_data(cns_data->tag_list),
        kv_size(cns_data->tag_list),
        subject_size,
        cns_data->dci_alloc,
        cns_data->li_alloc,
        &cns_data->item_list,
        &cns_data->cov_list);
    
    consensus_backbone_segment(kv_data(cns_data->item_list),
        *sfrom,
        *sto,
        kv_data(cns_data->cov_list),
        cns_seq,
        sfrom,
        sto,
        INDEL_COV_FACTOR);

    if (ks_size(*cns_seq) < min_size) return FALSE;
    return TRUE;
}

static void
correct_one_sv_read(CnsThreadData* data, 
    SvReadGroup* svr_group,
    const int local_read_id,
    SvReadGroup* cns_group,
    const int cns_iter)
{
    hbn_assert(local_read_id < svr_group->sv_read_group_size);
    hbn_assert(svr_group->svri_array[local_read_id].local_id == local_read_id);
    const int subject_id = svr_group->svri_array[local_read_id].global_id;
    const char* subject_name = data->raw_read_info_array[subject_id].hdr_offset + data->raw_read_name_array;
    const int subject_length = svr_group->svri_array[local_read_id].size;
    if (subject_length == 0) return;
    const u8* fwd_subject = SvReadGroup_ExtractSeq(svr_group, local_read_id, FWD);
    const u8* rev_subject = SvReadGroup_ExtractSeq(svr_group, local_read_id, REV);
    //HBN_LOG("correcting %d:%d:%s", local_read_id, subject_id, subject_name);
    InitHitFindData* hit_finder = data->hit_finder;
    InitHitFindData_Init(hit_finder, subject_id, subject_name,
        fwd_subject, rev_subject, subject_length);
    u8* cov_stats = (u8*)calloc(subject_length, sizeof(u8));
    kstring_t* qaln = &data->qaln;
    kstring_t* saln = &data->saln;
    vec_chain_seed* chain_seed_list = &data->chain_seed_list;
    int num_extended_can = 0;
    int num_added_aln = 0;
    FCCnsDataClear(data->cns_data);

    for (int i = 0; i < svr_group->sv_read_group_size; ++i) {
        if (i == local_read_id) continue;
        const int query_id = i;
        const char* query_name = SvReadGroup_ReadName(svr_group, i);
        const u8* fwd_query = SvReadGroup_ExtractSeq(svr_group, i, FWD);
        const u8* rev_query = SvReadGroup_ExtractSeq(svr_group, i, REV);
        const int query_length = SvReadGroup_ReadSize(svr_group, i);
        if (query_length == 0) continue;
        InitHitFindData_AddQuery(hit_finder, query_id, query_name,
            fwd_query, rev_query, query_length);
        InitHitFindData_FindHits(hit_finder, F_R);
        if (kv_empty(hit_finder->hit_list)) {
		    //HBN_LOG("hit find fail, score = %d", cns_hit_array[i].score);
		    continue;
	    }
        HbnInitHit* hit = kv_data(hit_finder->hit_list);
        InitHitFindData_SetupCnsAlignInfoFromHit(hit_finder, hit, chain_seed_list);
        ChainSeed* csa = kv_data(*chain_seed_list);
        const int csc = kv_size(*chain_seed_list);
        if (csc == 0) {
		    //HBN_LOG("cs cnt = 0");
		    continue;
	    }
        const u8* query = (hit->sdir == FWD) ? fwd_query : rev_query;
        validate_mem(HBN_LOG_ARGS_DEFAULT, query, fwd_subject, csa, csc);
        int sb = csa[0].soff;
        int se = csa[csc-1].soff + csa[csc-1].length;
        hbn_assert(sb >= 0);
        hbn_assert(se <= subject_length, "name = %d, i = %d, csc = %d, sb = %d, se = %d, sl = %d, hit_score = %d",
            subject_name, i, csc, sb, se, subject_length, hit->score);
        if (i >= MAX_CNS_COV && subject_subseq_cov_is_full(cov_stats, sb, se, MAX_CNS_COV)) {
            //HBN_LOG("subseq is full");
            continue;
        }
        ++num_extended_can;

        int qb, qe;
        double perc_identity;
        int r;
        if (cns_iter == 1) {
            r = extend_chain_seed_array_simple(data->tbck_data,
                    data->ksw_data,
                    csa,
                    csc,
                    data->opts,
                    data->opts->cns1_perc_identity,
                    query_id,
                    query,
                    query_length,
                    subject_id,
                    fwd_subject,
                    subject_length,
                    qaln,
                    saln,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity);
        } else {
            r = extend_chain_seed_array_ksw(data->tbck_data,
                    data->ksw_data,
                    csa,
                    csc,
                    data->opts,
                    data->opts->cns2_perc_identity,
                    query_id,
                    query,
                    query_length,
                    subject_id,
                    fwd_subject,
                    subject_length,
                    qaln,
                    saln,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity);
        }
        if (!r) {
            //print_init_hit_range(hit);
            //for (size_t x = 0; x < kv_size(hit_finder->hit_list) && x < 5; ++x) {
            //    hit = &kv_A(hit_finder->hit_list, x);
            //    fprintf(stderr, "%d\t", x);
            //    print_init_hit_range(hit);
            //}
            //HBN_LOG("extension fail, score = %d, csc = %d", hit->score, csc);
            continue;
        }
        //HBN_LOG("add ovlp [%s, %d, %d, %d] x [%d, %d, %d] %g", 
        //    query_name, qb, qe, query_length,
        //    sb, se, subject_length, perc_identity);
        //s_find_sv_signature(query_name, ks_s(*qaln), ks_s(*saln), ks_size(*qaln),
        //    0, qb, qe, sb, se);
        for (int x = sb; x < se; ++x) ++cov_stats[x];
        ++num_added_aln;
        make_align_tags_from_ovlp(ks_s(*qaln),
            ks_s(*saln),
            ks_size(*qaln),
            qb,
            qe,
            sb,
            se,
            DEFAULT_CNS_WEIGHT,
            &data->cns_data->tag_list);

        if (num_added_aln >= MAX_CNS_COV 
            && 
            subject_subseq_cov_is_full(cov_stats, 0, subject_length, MAX_CNS_COV)) {
            //HBN_LOG("subject cover is full");
            break;
        }
    }

    int from = 0, to = 0;
    int i = 0;
    int max_size = 0;
    while (i < subject_length) {
        while (i < subject_length && cov_stats[i] < data->opts->min_cov) ++i;
        if (i >= subject_length) break;
        int j = i + 1;
        while (j < subject_length && cov_stats[j] >= data->opts->min_cov) ++j;
        if (j - i >= data->opts->min_size) {
            if (j - i > max_size) {
                max_size = j - i;
                from = i;
                to = j;
            }
        }
        i = j;
    }
    sfree(cov_stats);
    if (max_size < data->opts->min_size) return;
    //HBN_LOG("from = %d, to = %d", from, to);

    kstring_t* cns_subseq = qaln;
    ks_clear(*cns_subseq);
    if (!meap_consensus_one_segment(data->cns_data,
        &from,
        &to,
        subject_length,
        data->opts->min_size,
        cns_subseq)) return;

    SvReadInfo* svri = cns_group->svri_array + local_read_id;
    SvReadInfo* src = svr_group->svri_array + local_read_id;
    svri->global_id = src->global_id;
    svri->local_id = src->local_id;
    svri->name = src->name;
    svri->size = from + ks_size(*cns_subseq) + (subject_length - to);
    svri->raw_seq_size = subject_length;
    svri->raw_seq_from = from;
    svri->raw_seq_to = svri->raw_seq_from + ks_size(*cns_subseq);
    svri->fsfrom = src->fsfrom;
    svri->fsto = src->fsto;
    svri->fsqdir = src->fsqdir;

    hbn_assert(data->out_lock != NULL);
    pthread_mutex_lock(data->out_lock);
    svri->offset_in_sv_read = kv_size(cns_group->fwd_sv_read);
    /// forward consensus sequence
    if (from) {
        int nitem = from;
        const u8* ptr = fwd_subject;
        kv_push_v(u8, cns_group->fwd_sv_read, ptr, nitem);
    }
    if (1)
    {
        int nitem = ks_size(*cns_subseq);
        const u8* ptr = (u8*)(ks_s(*cns_subseq));
        kv_push_v(u8, cns_group->fwd_sv_read, ptr, nitem);
    }
    if (subject_length > to) {
        int nitem = subject_length - to;
        const u8* ptr = fwd_subject + to;
        kv_push_v(u8, cns_group->fwd_sv_read, ptr, nitem);
    }

    /// reverse consensus sequence
    size_t n = kv_size(cns_group->fwd_sv_read);
    while (n > svri->offset_in_sv_read) {
        --n;
        u8 c = kv_A(cns_group->fwd_sv_read, n);
        c = 3 - c;
        kv_push(u8, cns_group->rev_sv_read, c);
    }
    pthread_mutex_unlock(data->out_lock);
}

/////////////////////// mt

static void
s_dump_sv_read_info(const int group_id,
    SvReadInfo* svri,
    const char* read_name,
    const int read_lid,
    const int read_gid,
    const u8* read,
    const char* subject_name,
    const int subject_id,
    const u8* subject,
    const int total_ssize,
    const int qb,
    const int qe,
    const int sb,
    const int se,
    const char* qaln,
    const char* saln,
    const int aln_size,
    FILE* fasta_out,
    FILE* sam_out,
    pthread_mutex_t* out_lock)
{
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        read_lid, read, qb, qe, qaln,
        subject_id, subject, sb, se, saln, aln_size, TRUE);
    //HBN_LOG("[%d, %d, %d, %d, %d] x [%d, %d]", read_lid, svri->fsqdir, qb, qe, svri->size, sb, se);
    int cns_qb = qb;
    int cns_qe = qe;
    if (svri->fsqdir == FWD) {
        cns_qb = svri->raw_seq_from;
        cns_qe = svri->raw_seq_to;
    } else {
        cns_qb = qb;
        cns_qe = qe;
    }

    int aln_from = 0;
    int qi = qb;
    if (cns_qb > qi) {
        while (qi < cns_qb) {
            if (qaln[aln_from] != GAP_CHAR) ++qi;
            ++aln_from;
        }
    }
    int aln_to = aln_size;
    qi = qe;
    if (qi > cns_qe) {
        while (qi > cns_qe) {
            if (qaln[aln_to - 1] != GAP_CHAR) --qi;
            --aln_to;
        }
    }
    if (aln_from >= aln_to) return;

    int qif = qb, qie = 0;
    int sif = sb, sie = 0;
    int ai = 0;
    while (ai < aln_from) {
        if (qaln[ai] != GAP_CHAR) ++qif;
        if (saln[ai] != GAP_CHAR) ++sif;
        ++ai;
    }
    hbn_assert(ai == aln_from);

    qie = qif;
    sie = sif;
    while (ai < aln_to) {
        if (qaln[ai] != GAP_CHAR) ++qie;
        if (saln[ai] != GAP_CHAR) ++sie;
        ++ai;
    }
    hbn_assert(ai == aln_to);

    //HBN_LOG("qif = %d, qie = %d, sif = %d, sie = %d, aln_from = %d, aln_to = %d",
    //    qif, qie, sif, sie, aln_from, aln_to);

    const char* q = qaln + aln_from;
    const char* s = saln + aln_from;
    const int as = aln_to - aln_from;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        read_lid, read, qif, qie, q,
        subject_id, subject, sif, sie, s, as, TRUE);

    int score = 0;
    double perc_identity = calc_ident_perc(q, s, as, NULL, &score);
    double eff_perc_identity = calc_effective_ident_perc(q, s, as);

#if 0
    HBN_LOG("[%d:%d:%s, %d, %d, %d] x [%d, %d], perc_identity = %g, eff_perc_identity = %g", 
        read_lid, read_gid, read_name, qif, qie, svri->size,
        sif, sie, 
        perc_identity, eff_perc_identity);
#endif
    
    ks_dinit(fasta);
    ks_dinit(assembled_name);
    char gid_str[64];
    char sid_str[64];
    u64_to_fixed_width_string_r(group_id, gid_str, HBN_DIGIT_WIDTH);
    u64_to_fixed_width_string_r(subject_id, sid_str, HBN_DIGIT_WIDTH);
    ksprintf(&assembled_name, "%s_s%s_g%s", read_name, sid_str, gid_str);
    ksprintf(&fasta, ">%s\n", ks_s(assembled_name));
    for (int i = qif; i < qie; ++i) {
        int c = read[i];
        c = DECODE_RESIDUE(c);
        kputc(c, &fasta);
    }
    kputc('\n', &fasta);

    ks_dinit(sam);
    ks_dinit(rg_info);
    make_sam_rg_info(subject_id, &rg_info);
    print_one_sam_result(ks_s(assembled_name), svri->fsqdir, 0, qie - qif, qie - qif,
        subject_name, svri->fsfrom + sif, svri->fsfrom + sie, total_ssize,
        q, s, as, score, perc_identity, eff_perc_identity, 
        TRUE, ks_s(rg_info), &sam);

    pthread_mutex_lock(out_lock);
    if (fasta_out) hbn_fwrite(ks_s(fasta), 1, ks_size(fasta), fasta_out);
    if (sam_out) hbn_fwrite(ks_s(sam), 1, ks_size(sam), sam_out);
    pthread_mutex_unlock(out_lock);

    ks_destroy(assembled_name);
    ks_destroy(rg_info);
    ks_destroy(sam);
    ks_destroy(fasta);
}

static void*
map_sv_read_func(void* params)
{
    CnsThreadData* data = (CnsThreadData*)(params);
    SvReadGroup* group = data->cns_group;
    const CSeqDB* db = data->db;
    const int subject_id = group->subject_id;
    const char* subject_name = seqdb_seq_name(db, subject_id);
    kv_dinit(vec_u8, subject_v);

    while (1) {
        int i = -1;
        pthread_mutex_lock(data->svr_id_lock);
        i = *data->next_svr_id;
        ++(*data->next_svr_id);
        pthread_mutex_unlock(data->svr_id_lock);
        if (i >= group->sv_read_group_size) break;
        if (group->svri_array[i].size == 0) continue;

        const int read_lid = group->svri_array[i].local_id;
        const int read_gid = group->svri_array[i].global_id;
        const u8* read = SvReadGroup_ExtractSeq(group, i, group->svri_array[i].fsqdir);
        const int read_length = group->svri_array[i].size;
        const char* read_name = group->svri_array[i].name;
        seqdb_extract_subsequence(db, subject_id, group->svri_array[i].fsfrom,
            group->svri_array[i].fsto, FWD, &subject_v);
        const u8* subject = kv_data(subject_v);
        const int subject_length = kv_size(subject_v);
        int qb, qe, sb, se;
        double perc_identity;
        int distance = hbn_max(read_length, subject_length);
        distance = distance * 0.2;
        int r = nw_ksw2_extd2(data->ksw_data,
                    read_lid,
                    read,
                    0,
                    read_length,
                    read_length,
                    subject_id,
                    subject,
                    0,
                    subject_length,
                    subject_length,
                    0,
                    0,
                    distance,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity,
                    &data->qaln,
                    &data->saln);
        if (!r) {
            distance = hbn_max(read_length, subject_length);
            r = nw_ksw2_extd2(data->ksw_data,
                    read_lid,
                    read,
                    0,
                    read_length,
                    read_length,
                    subject_id,
                    subject,
                    0,
                    subject_length,
                    subject_length,
                    0,
                    0,
                    distance,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity,
                    &data->qaln,
                    &data->saln);            
        }
        if (!r) {
            HBN_LOG("ksw fail for %d", read_lid);
            continue;
        }

        s_dump_sv_read_info(group->group_id,
            group->svri_array + i,
            read_name,
            read_lid,
            read_gid,
            read,
            subject_name,
            subject_id,
            subject,
            seqdb_seq_size(db, subject_id),
            qb,
            qe,
            sb,
            se,
            ks_s(data->qaln),
            ks_s(data->saln),
            ks_size(data->qaln),
            data->fasta_out,
            data->sam_out,
            data->out_lock);
    }
    kv_destroy(subject_v);
    return NULL;
}

void
map_sv_read_group_mt(CnsThreadData** data, SvReadGroup* group)
{
    const int num_threads = data[0]->opts->num_threads;
    pthread_t jobs[num_threads];
    *data[0]->next_svr_id = 0;
    for (int i = 0; i < num_threads; ++i) {
        data[i]->raw_group = NULL;
        data[i]->cns_group = group;
        pthread_create(jobs + i, NULL, map_sv_read_func, data[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], NULL);
        data[i]->raw_group = NULL;
        data[i]->cns_group = NULL;
    }    
}

static void*
s_svr_group_cns_func(void* params)
{
    CnsThreadData* data = (CnsThreadData*)(params);
    SvReadGroup* raw_group = data->raw_group;
    SvReadGroup* cns_group = data->cns_group;
    hbn_assert(raw_group);
    hbn_assert(cns_group);
    while (1) {
        int svr_id = 0;
        pthread_mutex_lock(data->svr_id_lock);
        svr_id = *data->next_svr_id;
        ++(*data->next_svr_id);
        pthread_mutex_unlock(data->svr_id_lock);
        if (svr_id >= raw_group->sv_read_group_size) break;
        if (raw_group->svri_array[svr_id].size == 0) continue;
        hbn_assert(data->cns_iter == 1 || data->cns_iter == 2);
        correct_one_sv_read(data, raw_group, svr_id, cns_group, data->cns_iter);
    }
    return NULL;
}

static void
s_cns_one_group_mt(CnsThreadData** data,
    SvReadGroup* raw_group,
    SvReadGroup* cns_group,
    const int cns_iter)
{
    const int num_threads = data[0]->opts->num_threads;
    pthread_t jobs[num_threads];
    *data[0]->next_svr_id = 0;
    for (int i = 0; i < num_threads; ++i) {
        data[i]->raw_group = raw_group;
        data[i]->cns_group = cns_group;
        data[i]->cns_iter = cns_iter;
        pthread_create(jobs + i, NULL, s_svr_group_cns_func, data[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], NULL);
        data[i]->cns_iter = -1;
    }
}

int 
cns_one_group_mt(CnsThreadData** data, 
    FILE* sv_read_group_file,
    const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* packed_read_file,
    const int subject_id)
{
    SvReadGroup* svr_group =  SvReadGroupLoad(raw_read_info_array,
                                    raw_read_name_array,
                                    subject_id,
                                    sv_read_group_file,
                                    packed_read_file);
    if (!svr_group) return 0;
#if 0
    if (svr_group->group_id != 111) {
        svr_group = SvReadGroupFree(svr_group);
        return 1;
    }
#endif
    HBN_LOG("correcting group %d, %d reads", svr_group->group_id, svr_group->sv_read_group_size);
    char job_name[256];
    SvReadGroup* cns1_group = NULL;
    SvReadGroup* cns2_group = NULL;

    {
        sprintf(job_name, "correction 1");
        hbn_timing_begin(job_name);
        cns1_group = SvReadGroupNew(raw_read_info_array, 
                                raw_read_name_array, 
                                svr_group->group_id, 
                                svr_group->subject_id, 
                                svr_group->sv_type,
                                svr_group->sv_read_group_size);
        s_cns_one_group_mt(data, svr_group, cns1_group, 1);
        hbn_timing_end(job_name);
    }

    {
        sprintf(job_name, "correction 2");
        hbn_timing_begin(job_name);
        cns2_group = SvReadGroupNew(raw_read_info_array, 
                                raw_read_name_array,
                                cns1_group->group_id, 
                                cns1_group->subject_id, 
                                cns1_group->sv_type,
                                cns1_group->sv_read_group_size);
        s_cns_one_group_mt(data, cns1_group, cns2_group, 2);
        hbn_timing_end(job_name);
    }

    //calc_ovlp_identity_in_svr_group(data, cns2_group);
	//HBN_LOG("************************* identity of cns1:");
    //map_sv_read_group_mt(data, cns1_group);
	//HBN_LOG("************************* identity of cns2:");
    if (0)
    {
        sprintf(job_name, "map corrected reads");
        hbn_timing_begin(job_name);
        map_sv_read_group_mt(data, cns2_group);
        hbn_timing_end(job_name);
    }
    //HBN_LOG("************************** pm identity for read 0");
    //s_calc_ovlp_identity_in_svr_group_for_one_read(data, cns2_group, 0);

    SvReadGroupFree(svr_group);
    SvReadGroupFree(cns1_group);
    SvReadGroupFree(cns2_group);
    //exit(0);
    return 1;    
}

//////// st

static void
s_dump_corrected_sv_group(CnsThreadData* data)
{
    SvReadGroup* group = data->cns_group;
    const int group_id = group->group_id;
    const int subject_id = group->subject_id;
    ks_dinit(fasta);
    ks_dinit(assembled_name);

    for (int i = 0; i < group->sv_read_group_size; ++i) {
        if (group->svri_array[i].size == 0) continue;

        const int read_lid = group->svri_array[i].local_id;
        const int read_gid = group->svri_array[i].global_id;
        const u8* read = SvReadGroup_ExtractSeq(group, i, FWD);
        const int read_length = group->svri_array[i].size;
        const char* read_name = group->svri_array[i].name;

        make_sv_read_header(read_name, group->svri_array[i].fsqdir, subject_id,
            group_id,
            group->svri_array[i].fsfrom, group->svri_array[i].fsto, &fasta);
        make_cns_read_header(ks_s(fasta), 0, 0, 
            group->svri_array[i].raw_seq_from,
            group->svri_array[i].raw_seq_to,
            group->svri_array[i].size,
            group->svri_array[i].raw_seq_size,
            &assembled_name);

        ks_clear(fasta);
        ksprintf(&fasta, ">%s\n", ks_s(assembled_name));
        for (int p = 0; p < read_length; ++p) {
            int c = read[p];
            c = DECODE_RESIDUE(c);
            kputc(c, &fasta);
        }
        kputc('\n', &fasta);
        hbn_assert(data->out_lock);
        pthread_mutex_lock(data->out_lock);
        hbn_fwrite(ks_s(fasta), 1, ks_size(fasta), data->fasta_out);
        pthread_mutex_unlock(data->out_lock);
    }

    ks_destroy(assembled_name);
    ks_destroy(fasta);
}

static void
s_cns_one_group_st(CnsThreadData* data,
    SvReadGroup* raw_group,
    SvReadGroup* cns_group,
    const int cns_iter)
{
    int next_svr_id = 0;
    pthread_mutex_t svr_id_lock = PTHREAD_MUTEX_INITIALIZER;
    data->next_svr_id = &next_svr_id;
    data->svr_id_lock = &svr_id_lock;
    data->raw_group = raw_group;
    data->cns_group = cns_group;
    data->cns_iter = cns_iter;
    s_svr_group_cns_func(data);
    data->next_svr_id = NULL;
    data->svr_id_lock = NULL;
    data->cns_iter = -1;
}

void
map_sv_read_group_st(CnsThreadData* data, SvReadGroup* group)
{
    int next_svr_id = 0;
    pthread_mutex_t svr_id_lock = PTHREAD_MUTEX_INITIALIZER;
    data->next_svr_id = &next_svr_id;
    data->svr_id_lock = &svr_id_lock;
    data->raw_group = NULL;
    data->cns_group = group;
    map_sv_read_func(data);
    data->next_svr_id = NULL;
    data->svr_id_lock = NULL;
    data->raw_group = NULL;
    data->cns_group = NULL;
}

int 
cns_one_group_st_worker(CnsThreadData* data, 
    FILE* sv_read_group_file,
    pthread_mutex_t* sv_read_group_load_lock,
    const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* packed_read_file,
    const int subject_id)
{
    pthread_mutex_lock(sv_read_group_load_lock);
    SvReadGroup* svr_group =  SvReadGroupLoad(raw_read_info_array,
                                    raw_read_name_array,
                                    subject_id,
                                    sv_read_group_file,
                                    packed_read_file);
    pthread_mutex_unlock(sv_read_group_load_lock);
    if (!svr_group) return 0;
    char job_name[256];
    sprintf(job_name, "correct group %d (%d reads)", svr_group->group_id, svr_group->sv_read_group_size);
    hbn_timing_begin(job_name);

    SvReadGroup* cns1_group = SvReadGroupNew(raw_read_info_array, 
                                raw_read_name_array, 
                                svr_group->group_id, 
                                svr_group->subject_id, 
                                svr_group->sv_type,
                                svr_group->sv_read_group_size);
    s_cns_one_group_st(data, svr_group, cns1_group, 1);

    SvReadGroup* cns2_group = SvReadGroupNew(raw_read_info_array, 
                                raw_read_name_array,
                                cns1_group->group_id, 
                                cns1_group->subject_id, 
                                cns1_group->sv_type,
                                cns1_group->sv_read_group_size);
    s_cns_one_group_st(data, cns1_group, cns2_group, 2);

    s_dump_corrected_sv_group(data);

    //map_sv_read_group_st(data, cns2_group);

    SvReadGroupFree(svr_group);
    SvReadGroupFree(cns1_group);
    SvReadGroupFree(cns2_group);
    hbn_timing_end(job_name);
    return 1;    
}

typedef struct {
    CnsThreadData* data;
    FILE* sv_read_group_file;
    pthread_mutex_t* sv_read_group_load_lock;
    FILE* packed_read_file;
    int subject_id;
} SvReadGroupCnsStData;

static void*
s_cns_one_group_st_func(void* params)
{
    SvReadGroupCnsStData* data = (SvReadGroupCnsStData*)(params);
    while (1) {
        int r = cns_one_group_st_worker(data->data, 
                    data->sv_read_group_file,
                    data->sv_read_group_load_lock,
                    data->data->raw_read_info_array,
                    data->data->raw_read_name_array,
                    data->packed_read_file,
                    data->subject_id);
        if (!r) break;
        //break;
    }
    return NULL;
}

void
cns_suject_svr_groups_st(int argc, char* argv[],
    HbnProgramOptions* opts,
    CSeqDB* db,
    CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* pac_file,
    const int sid)
{
    char path[HBN_MAX_PATH_LEN];
    make_subject_sv_read_group_path(opts->svr_group_dir, sid, path);
    HBN_LOG("sv_read_group_path: %s", path);
    hbn_dfopen(sv_group_file, path, "r");
    pthread_mutex_t sv_group_file_load_lock = PTHREAD_MUTEX_INITIALIZER;
    make_subject_corrected_normal_sv_read_path(opts->wrk_dir, sid, path);
    hbn_dfopen(fasta_out, path, "w");
    FILE* sam_out = NULL;

    if (sam_out)
    {
        ks_dinit(rg_info);
        make_sam_rg_info(sid, &rg_info);
        print_sam_prolog(sam_out, kSamVersion, HBN_PACKAGE_VERSION, ks_s(rg_info), 
            DEFAULT_SAMPLE, db, sid, argc, argv);
        ks_destroy(rg_info);
    }
    pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
    const int num_threads = opts->num_threads;
    CnsThreadData* data[num_threads];
    SvReadGroupCnsStData* st_data[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        data[i] = CnsThreadDataNew(opts,
                        db,
                        raw_read_info_array,
                        raw_read_name_array,
                        fasta_out,
                        sam_out,
                        &out_lock,
                        NULL,
                        NULL);
        st_data[i] = (SvReadGroupCnsStData*)calloc(1, sizeof(SvReadGroupCnsStData));
        st_data[i]->data = data[i];
        st_data[i]->sv_read_group_file = sv_group_file;
        st_data[i]->sv_read_group_load_lock = &sv_group_file_load_lock;
        st_data[i]->packed_read_file = pac_file;
        st_data[i]->subject_id = sid;
    }

    pthread_t jobs[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobs + i, NULL, s_cns_one_group_st_func, st_data[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], NULL);
    }

    for (int i = 0; i < num_threads; ++i) {
        data[i] = CnsThreadDataFree(data[i]);
        sfree(st_data[i]);
    }
    hbn_fclose(sv_group_file);
    hbn_fclose(fasta_out);
    if (sam_out) hbn_fclose(sam_out);
}