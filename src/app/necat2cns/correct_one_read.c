#include "correct_one_read.h"

#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/refine_align.h"
#include "../../algo/cns_extend_chain_seed_list.h"
#include "../../corelib/cns_read_header.h"

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

static BOOL 
extend_chain_seed_array(HbnTracebackData* traceback_data,
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
    if (opts->use_ksw) {
        return cns_extend_chain_seed_array_with_ksw(traceback_data,
                    ksw2_data,
                    chain_seed_array,
                    chain_seed_count,
                    opts->ovlp_cov_res,
                    opts->ovlp_cov_perc,
                    opts->perc_identity,
                    read_id,
                    read,
                    read_length,
                    subject_id,
                    subject,
                    subject_length,
                    qaln,
                    saln,
                    qoff_,
                    qend_,
                    soff_,
                    send_,
                    ident_perc_);
    } else {
        return cns_extend_chain_seed_array(traceback_data,
                    ksw2_data,
                    chain_seed_array,
                    chain_seed_count,
                    opts->ovlp_cov_res,
                    opts->ovlp_cov_perc,
                    opts->perc_identity,
                    read_id,
                    read,
                    read_length,
                    subject_id,
                    subject,
                    subject_length,
                    qaln,
                    saln,
                    qoff_,
                    qend_,
                    soff_,
                    send_,
                    ident_perc_);        
    }
}

static void
s_do_correct_without_replace_mode(CnsThreadData* data,
    RawReadCnsInfo* cns_info,
    const char* subject_name,
    const int subject_length,
    const u8* cov_stats,
    const int num_extended_can,
    const int num_added_aln,
    kstring_t* cns_subseq,
    kstring_t* cns_fasta)
{
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
    if (max_size < data->opts->min_size) return;
    if (max_size < subject_length * 0.8) {
        //HBN_LOG("consensus subsequence is too short [%d, %d, %d]", from, to, subject_length);
        return;
    }
    //HBN_LOG("from = %d, to = %d", from, to);

    ks_clear(*cns_subseq);
    if (!meap_consensus_one_segment(data->cns_data,
        &from,
        &to,
        subject_length,
        data->opts->min_size,
        cns_subseq)) return;

    for (size_t p = 0; p < ks_size(*cns_subseq); ++p) {
        int c = ks_A(*cns_subseq, p);
        c = DECODE_RESIDUE(c);
        ks_A(*cns_subseq, p) = c;
    }    

    cns_info->cns_from = from;
    cns_info->cns_to = to;
    cns_info->cns_read_size = ks_size(*cns_subseq);
    cns_info->raw_read_size = subject_length;
    make_cns_read_header(subject_name,
        num_extended_can,
        num_added_aln,
        0,
        ks_size(*cns_subseq),
        ks_size(*cns_subseq),
        subject_length,
        &data->tbck_data->qabuf);
    ks_clear(*cns_fasta);
    ksprintf(cns_fasta, ">%s\n", ks_s(data->tbck_data->qabuf));
    kputsn(ks_s(*cns_subseq), ks_size(*cns_subseq), cns_fasta);
    kputc('\n', cns_fasta);
    pthread_mutex_lock(data->out_lock);
    cns_info->cns_fasta_offset = ks_size(*data->cns_fasta_array);
    cns_info->cns_fasta_size = ks_size(*cns_fasta);
    kputsn(ks_s(*cns_fasta), ks_size(*cns_fasta), data->cns_fasta_array);
    pthread_mutex_unlock(data->out_lock);
}

static void
s_do_correct_with_replace_mode(CnsThreadData* data,
    RawReadCnsInfo* cns_info,
    const char* subject_name,
    const int subject_length,
    const u8* fwd_subject,
    const u8* cov_stats,
    const int num_extended_can,
    const int num_added_aln,
    kstring_t* cns_subseq,
    kstring_t* cns_fasta)
{
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
    if (max_size < data->opts->min_size) return;
    //HBN_LOG("from = %d, to = %d", from, to);

    ks_clear(*cns_fasta);
    if (!meap_consensus_one_segment(data->cns_data,
        &from,
        &to,
        subject_length,
        data->opts->min_size,
        cns_fasta)) return;  

    for (size_t p = 0; p < ks_size(*cns_fasta); ++p) {
        int c = ks_A(*cns_fasta, p);
        hbn_assert(c >= 0 && c < 4);
    }    

    ks_clear(*cns_subseq);  
    int cns_subseq_offset = 0;
    int cns_subseq_end = ks_size(*cns_fasta);
    if (from) {
        cns_subseq_offset = from;
        int nitem = from;
        const char* ptr = (const char*)fwd_subject;
        for (int p = 0; p < nitem; ++p) {
            int c = ptr[p];
            hbn_assert(c >= 0 && c < 4);
        }
        kputsn(ptr, nitem, cns_subseq);
        //HBN_LOG("size = %zu", ks_size(*cns_subseq));
        for (size_t p = 0; p < ks_size(*cns_subseq); ++p) {
            int c = ks_A(*cns_subseq, p);
            hbn_assert(c >= 0 && c < 4, "p = %zu, c = %d, size = %zu", p, c, ks_size(*cns_subseq));
        }    
    }
    if (1) {
        int nitem = ks_size(*cns_fasta);
        const char* ptr = ks_s(*cns_fasta);
        for (int p = 0; p < nitem; ++p) {
            int c = ptr[p];
            hbn_assert(c >= 0 && c < 4);
        }
        kputsn(ptr, nitem, cns_subseq);
        //HBN_LOG("size = %zu", ks_size(*cns_subseq));
        for (size_t p = 0; p < ks_size(*cns_subseq); ++p) {
            int c = ks_A(*cns_subseq, p);
            hbn_assert(c >= 0 && c < 4, "p = %zu, c = %d, size = %zu", p, c, ks_size(*cns_subseq));
        }   
    }
    if (subject_length > to) {
        cns_subseq_end = ks_size(*cns_subseq);
        int nitem = subject_length - to;
        const char* ptr = (const char*)(fwd_subject + to);
        for (int p = 0; p < nitem; ++p) {
            int c = ptr[p];
            hbn_assert(c >= 0 && c < 4);
        }
        kputsn(ptr, nitem, cns_subseq);
       // HBN_LOG("size = %zu", ks_size(*cns_subseq));
        for (size_t p = 0; p < ks_size(*cns_subseq); ++p) {
            int c = ks_A(*cns_subseq, p);
            hbn_assert(c >= 0 && c < 4, "p = %zu, c = %d, size = %zu", p, c, ks_size(*cns_subseq));
        }   
    }
    for (size_t p = 0; p < ks_size(*cns_subseq); ++p) {
        int c = ks_A(*cns_subseq, p);
        hbn_assert(c >= 0 && c < 4, "p = %zu, c = %d, size = %zu", p, c, ks_size(*cns_subseq));
        c = DECODE_RESIDUE(c);
        ks_A(*cns_subseq, p) = c;
    }    

    cns_info->cns_from = 0;
    cns_info->cns_to = subject_length;
    cns_info->cns_read_size = ks_size(*cns_subseq);
    cns_info->raw_read_size = subject_length;
    make_cns_read_header(subject_name,
        num_extended_can,
        num_added_aln,
        cns_subseq_offset,
        cns_subseq_end,
        ks_size(*cns_subseq),
        subject_length,
        &data->tbck_data->qabuf);

    ks_clear(*cns_fasta);
    ksprintf(cns_fasta, ">%s\n", ks_s(data->tbck_data->qabuf));
    kputsn(ks_s(*cns_subseq), ks_size(*cns_subseq), cns_fasta);
    kputc('\n', cns_fasta);
    pthread_mutex_lock(data->out_lock);
    cns_info->cns_fasta_offset = ks_size(*data->cns_fasta_array);
    cns_info->cns_fasta_size = ks_size(*cns_fasta);
    kputsn(ks_s(*cns_fasta), ks_size(*cns_fasta), data->cns_fasta_array);
    pthread_mutex_unlock(data->out_lock);
}

static BOOL
s_is_target_query(const char* query_name)
{
    const char* target1 = "80f02c09-1ecd-4568-a256-93bffd841e80";
    const char* target2 = "20e4c353-6b85-4e8e-86e3-63984260199d";
    const char* pos = strstr(query_name, target1);
    if (pos != NULL) return TRUE;
    pos = strstr(query_name, target2);
    if (pos != NULL) return TRUE;
    return FALSE;
}

static void
s_find_sv_signature(
    const char* query_name,
    const char* qas,
    const char* sas,
    const int aln_size,
    const int subject_offset,
    int qb,
    int qe,
    int sb,
    int se)
{
    const int E = 40;
    int i = 0;
    int qi = qb;
    int si = subject_offset + sb;
    while (i < aln_size) {
        if (qas[i] != GAP_CHAR && sas[i] != GAP_CHAR) {
            ++i;
            ++qi;
            ++si;
            continue;
        }

        if (qas[i] == GAP_CHAR) {
            hbn_assert(sas[i] != GAP_CHAR);
            int j = i + 1;
            while (j < aln_size && qas[j] == GAP_CHAR) {
                hbn_assert(sas[j] != GAP_CHAR);
                ++j;
            }
            int n = j - i;
            if (n >= E) {
                HBN_LOG("**** find DEL pos at [%d] x [%d, %d), len = %d", qi, si, si + n, n);
            }
            si += n;
            i = j;
            continue;
        }

        if (sas[i] == GAP_CHAR) {
            hbn_assert(qas[i] != GAP_CHAR);
            int j = i + 1;
            while (j < aln_size && sas[j] == GAP_CHAR) {
                hbn_assert(qas[j] != GAP_CHAR);
                ++j;
            }
            int n = j - i;
            if (n >= E) {
                HBN_LOG("**** find INS pos at [%d, %d] x [%d], length = %d", qi, qi + n, si, n);
            }
            qi += n;
            i = j;
            continue;
        }
    }
    hbn_assert(qi == qe);
    hbn_assert(si == subject_offset + se);
}

void
correct_one_read(CnsThreadData* data, int cns_info_id)
{
    hbn_assert(cns_info_id < data->cns_info_count);
    RawReadCnsInfo* cns_info = data->cns_info_array + cns_info_id;
    const int subject_id = cns_info->oid;
    hbn_assert(subject_id < data->raw_reads->dbinfo.num_seqs);
    const char* subject_name = RawReadReader_ReadName(data->raw_reads, subject_id);
    const int subject_length = RawReadReader_ReadSize(data->raw_reads, subject_id);
    if (subject_length < data->opts->min_size) return;
    //if (!s_is_target_query(subject_name)) return;
	//const char* target = "SRR6702603.53971";
    //if (strcmp(subject_name, target)) return;
    const u8* fwd_subject = RawReadReader_ExtractRead(data->raw_reads, subject_id, FWD, &data->fwd_subject);
    const u8* rev_subject = RawReadReader_ExtractRead(data->raw_reads, subject_id, REV, &data->rev_subject);
    InitHitFindData* hit_finder = data->hit_finder;
    InitHitFindData_Init(hit_finder, subject_id, subject_name, fwd_subject, rev_subject, subject_length);
    FCCnsDataClear(data->cns_data);
    HbnConsensusInitHit* cns_hit_array = data->hit_array + cns_info->can_from;
    int cns_hit_count =  cns_info->can_to - cns_info->can_from;
    //HBN_LOG("can_from = %zu, can_to = %zu", cns_info->can_from, cns_info->can_to);
    if (cns_hit_count < data->opts->min_cov) return;
    ks_introsort_cns_hit_score_gt(cns_hit_count, cns_hit_array);
    kv_resize(u8, data->cov_stats, subject_length);
    kv_zero(u8, data->cov_stats);
    u8* cov_stats = kv_data(data->cov_stats);
    kstring_t* qaln = &data->qaux;
    kstring_t* saln = &data->saux;
    vec_chain_seed* chain_seed_list = &data->chain_seed_list;
    int num_extended_can = 0;
    int num_added_aln = 0;
    HBN_LOG("correct subject %d:%s:%d:%d", subject_id, subject_name, subject_length, cns_hit_count);

    for (int i = 0; i < cns_hit_count && i < MAX_CNS_SR; ++i) {
        const int query_id = cns_hit_array[i].qid;
        const char* query_name = RawReadReader_ReadName(data->raw_reads, query_id);
        //HBN_LOG("align sr %d:%s, score = %d", i, query_name, cns_hit_array[i].score);
        RawReadReader_ExtractRead(data->raw_reads, query_id, FWD, &data->fwd_read);
        RawReadReader_ExtractRead(data->raw_reads, query_id, REV, &data->rev_read);
        const u8* fwd_query = kv_data(data->fwd_read);
        const u8* rev_query = kv_data(data->rev_read);
        int query_length = RawReadReader_ReadSize(data->raw_reads, query_id);
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
#if 0
        if (!check_chain_seed_array_cov(csa, csc, query_length, subject_length, 
            data->opts->ovlp_cov_res, data->opts->ovlp_cov_perc)) {
		    //print_init_hit_range(hit);
		    //HBN_LOG("hit too short");
            continue;
        }
#endif
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
        //HBN_LOG("extend hit range [%d, %d] x [%d, %d]", csa[0].qoff, csa[csc-1].qoff + csa[csc-1].length, csa[0].soff, csa[csc-1].soff + csa[csc-1].length);
        r = extend_chain_seed_array(data->tbck_data,
                    data->ksw_data,
                    csa,
                    csc,
                    data->opts,
                    data->opts->perc_identity,
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
    //exit(0);

    if (data->opts->use_replace_mode) {
        s_do_correct_with_replace_mode(data,
            cns_info,
            subject_name,
            subject_length,
            fwd_subject,
            cov_stats,
            num_extended_can,
            num_added_aln,
            qaln,
            saln);
    } else {
        s_do_correct_without_replace_mode(data,
            cns_info,
            subject_name,
            subject_length,
            cov_stats,
            num_extended_can,
            num_added_aln,
            qaln,
            saln);
    }
//exit(0);
}
