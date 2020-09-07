#include "align_subseqs.h"

#include "../../algo/refine_align.h"
#include "../../algo/hbn_traceback_aux.h"

SubseqAlignResult*
SubseqAlignResultNew()
{
    SubseqAlignResult* result = (SubseqAlignResult*)calloc(1, sizeof(SubseqAlignResult));
    ks_init(result->qaux);
    ks_init(result->saux);
    return result;
}

SubseqAlignResult*
SubseqAlignResultFree(SubseqAlignResult* result)
{
    ks_destroy(result->qaux);
    ks_destroy(result->saux);
    free(result);
    return NULL;
}

static void
s_extract_align_info(RawReadReader* query_reader,
    RawReadReader* subject_reader,
    int query_gid,
    int query_dir,
    int query_from,
    int query_to,
    int subject_gid,
    int subject_from,
    int subject_to,
    vec_u8* query_v,
    const u8** query,
    int* qsubseq_size,
    const char** query_name,
    vec_u8* subject_v,
    const u8** subject,
    int* ssubseq_size,
    const char** subject_name)
{
    RawReadReader_ExtractSubRead(query_reader, query_gid, query_from, query_to, query_dir, query_v);
    *query = kv_data(*query_v);
    *qsubseq_size = kv_size(*query_v);
    *query_name = RawReadReader_ReadName(query_reader, query_gid);

    RawReadReader_ExtractSubRead(subject_reader, subject_gid, subject_from, subject_to, FWD, subject_v);
    *subject = kv_data(*subject_v);
    *ssubseq_size = kv_size(*subject_v);
    *subject_name = RawReadReader_ReadName(subject_reader, subject_gid);
}

BOOL
align_and_refine_subseq_with_edlib(RawReadReader* query_reader,
    RawReadReader* subject_reader,
    HbnTracebackData* tbck_data,
    const int max_dist,
    int query_gid, 
    int query_dir,
    int query_from,
    int query_to,
    int subject_gid,
    int subject_from,
    int subject_to,
    SubseqAlignResult* result)
{
//if (query_gid != 71 || subject_gid != 3) return FALSE;
//HBN_LOG("[%d, %d, %d, %d] x [%d, %d, %d, %d]", query_gid, query_dir, query_from, query_to, subject_gid, FWD, subject_from, subject_to);
    BOOL status = TRUE;
    kv_dinit(vec_u8, query_v);
    kv_dinit(vec_u8, subject_v);
    kv_dinit(vec_chain_seed, hit_seed_list);
    const u8* query = NULL;
    const char* query_name = NULL;
    int qsubseq_size = 0;
    const u8* subject = NULL;
    const char* subject_name = NULL;
    int ssubseq_size = 0;
    const char* qas = NULL;
    const char* qae = NULL;
    const char* sas = NULL;
    const char* sae = NULL;
    int qb, qe, sb, se;
    double perc_identity, eff_perc_identity;
    int dist, score;
    s_extract_align_info(query_reader, subject_reader, query_gid, query_dir, query_from, query_to,
        subject_gid, subject_from, subject_to, 
        &query_v, &query, &qsubseq_size, &query_name,
        &subject_v, &subject, &ssubseq_size, &subject_name);

    int r = edlib_nw(tbck_data->edlib,
                query,
                qsubseq_size,
                subject,
                ssubseq_size,
                -1,
                &tbck_data->qabuf,
                &tbck_data->sabuf); 
    if (!r) {
        HBN_LOG("[%d:%s, %d, %d] x [%d:%s, %d, %d] edlib fail", 
            query_gid, query_name, query_from, query_to, 
            subject_gid, subject_name, subject_from, subject_to);
        status = FALSE;
        goto clean_align_and_refine_subseq_with_edlib;       
    } 

    int aln_size = ks_size(tbck_data->qabuf);
    qas = ks_s(tbck_data->qabuf);
    qae = qas + aln_size;
    sas = ks_s(tbck_data->sabuf);
    sae = sas + aln_size;
    qb = 0;
    qe = qsubseq_size;
    sb = 0;
    se = ssubseq_size;
    r = truncate_align_bad_ends(qas, sas, aln_size, &qb, &qe, &sb, &se, &qas, &qae, &sas, &sae);
    if (!r) {
        HBN_LOG("[%d:%s, %d, %d] x [%d:%s, %d, %d] truncate fail", 
            query_gid, query_name, query_from, query_to, 
            subject_gid, subject_name, subject_from, subject_to);
        status = FALSE;
        goto clean_align_and_refine_subseq_with_edlib;        
    }
    aln_size = qae - qas;
    int x = (qe) - (qb);
    int y = (se) - (sb);
    if (x < qsubseq_size * 0.8 && y < ssubseq_size * 0.8) {
        HBN_LOG("align too short: (%d:%s) [%d, %d] => [%d, %d] and (%d:%s) [%d, %d] => [%d, %d]",
            query_gid, query_name, query_from, query_to, (qb) + query_from, (qe) + query_from,
            subject_gid, subject_name, subject_from, subject_to, (sb) + subject_from, (se) + subject_from);        
        status = FALSE;
        goto clean_align_and_refine_subseq_with_edlib; 
    }    

    r = refine_align(tbck_data->ksw,
            query,
            qb,
            qe,
            qsubseq_size,
            subject,
            sb,
            se,
            ssubseq_size,
            qas,
            sas,
            aln_size,
            &result->qaux,
            &result->saux);
    if (!r) {
        HBN_LOG("[%d:%s, %d, %d] x [%d:%s, %d, %d] refine fail", 
            query_gid, query_name, query_from, query_to, 
            subject_gid, subject_name, subject_from, subject_to);
        status = FALSE;
        goto clean_align_and_refine_subseq_with_edlib;
    }
    hbn_assert(ks_size(result->qaux) == ks_size(result->saux));
    aln_size = ks_size(result->qaux);
    qas = ks_s(result->qaux);
    qae = qas + aln_size;
    sas = ks_s(result->saux);
    sae = sas + aln_size;
	//dump_align_string(qas, sas, aln_size, stderr);
    qb = (qb) + query_from;
    qe = (qe) + query_from;
    sb = (sb) + subject_from;
    se = (se) + subject_from;
    perc_identity = calc_ident_perc(qas, sas, aln_size, &dist, &score);
    eff_perc_identity = calc_effective_ident_perc(qas, sas, aln_size);

clean_align_and_refine_subseq_with_edlib:
    if (status) {
        result->qb = qb;
        result->qe = qe;
        result->sb = sb;
        result->se = se;
        result->qas = qas;
        result->qae = qae;
        result->sas = sas;
        result->sae = sae;
        result->dist = dist;
        result->score = score;
        result->perc_identity = perc_identity;
        result->eff_perc_identity = eff_perc_identity;
    }
    kv_destroy(query_v);
    kv_destroy(subject_v);
    kv_destroy(hit_seed_list);
    return status;
}

BOOL
align_and_refine_subseq_with_ksw(RawReadReader* query_reader,
    RawReadReader* subject_reader,
    HbnTracebackData* tbck_data,
    const int max_dist,
    int query_gid, 
    int query_dir,
    int query_from,
    int query_to,
    int subject_gid,
    int subject_from,
    int subject_to,
    SubseqAlignResult* result)
{
    BOOL status = TRUE;
    kv_dinit(vec_u8, query_v);
    kv_dinit(vec_u8, subject_v);
    kv_dinit(vec_chain_seed, hit_seed_list);
    const u8* query = NULL;
    const char* query_name = NULL;
    int qsubseq_size = 0;
    const u8* subject = NULL;
    const char* subject_name = NULL;
    int ssubseq_size = 0;
    const char* qas = NULL;
    const char* qae = NULL;
    const char* sas = NULL;
    const char* sae = NULL;
    int qb, qe, sb, se;
    double perc_identity, eff_perc_identity;
    int dist, score;
    s_extract_align_info(query_reader, subject_reader, query_gid, query_dir, query_from, query_to,
        subject_gid, subject_from, subject_to, 
        &query_v, &query, &qsubseq_size, &query_name,
        &subject_v, &subject, &ssubseq_size, &subject_name);
    int max_subseq_size = hbn_max(qsubseq_size, ssubseq_size);

    qb = 0;
    qe = 0;
    sb = 0;
    se = 0;
    int band_width = (max_dist < 0) ? max_subseq_size : max_dist;
    int r = nw_ksw2_extd2(tbck_data->ksw,
                query_gid,
                query,
                0,
                qsubseq_size,
                qsubseq_size,
                subject_gid,
                subject,
                0,
                ssubseq_size,
                ssubseq_size,
                0,
                0,
                band_width,
                &qb,
                &qe,
                &sb,
                &se,
                &perc_identity,
                &result->qaux,
                &result->saux);
    if ((r == FALSE) && (band_width < KSW_MAX_BW) &&  (band_width < max_subseq_size)) {
        band_width = max_subseq_size;
        r = nw_ksw2_extd2(tbck_data->ksw,
                query_gid,
                query,
                0,
                qsubseq_size,
                qsubseq_size,
                subject_gid,
                subject,
                0,
                ssubseq_size,
                ssubseq_size,
                0,
                0,
                band_width,
                &qb,
                &qe,
                &sb,
                &se,
                &perc_identity,
                &result->qaux,
                &result->saux);
	if (r) {
		HBN_LOG("rescue: [%d:%s, %d, %d, %d, %d] x [%d:%s, %d, %d, %d, %d]", query_gid, query_name, query_dir, qb, qe, qsubseq_size,
			subject_gid, subject_name, FWD, sb, se, ssubseq_size);
	}
    }
    if (!r) {
        HBN_LOG("[%d:%s, %d, %d] x [%d:%s, %d, %d] ksw fail, band_width = %d", 
            query_gid, query_name, query_from, query_to, 
            subject_gid, subject_name, subject_from, subject_to, band_width);
        status = FALSE;
        goto clean_align_and_refine_subseq_with_ksw;        
    }

    int aln_size = ks_size(result->qaux);
    qas = ks_s(result->qaux);
    qae = qas + aln_size;
    sas = ks_s(result->saux);
    sae = sas + aln_size;
    r = truncate_align_bad_ends(qas, sas, aln_size, &qb, &qe, &sb, &se, &qas, &qae, &sas, &sae);
    if (!r) {
        HBN_LOG("[%d:%s, %d, %d] x [%d:%s, %d, %d] truncate fail", 
            query_gid, query_name, query_from, query_to, 
            subject_gid, subject_name, subject_from, subject_to);
        status = FALSE;
        goto clean_align_and_refine_subseq_with_ksw;        
    }
    aln_size = qae - qas;
    int x = (qe) - (qb);
    int y = (se) - (sb);
    if (x < qsubseq_size * 0.8 && y < ssubseq_size * 0.8) {
        HBN_LOG("align too short: (%d:%s) [%d, %d] => [%d, %d] and (%d:%s) [%d, %d] => [%d, %d]",
            query_gid, query_name, query_from, query_to, (qb) + query_from, (qe) + query_from,
            subject_gid, subject_name, subject_from, subject_to, (sb) + subject_from, (se) + subject_from);        
        status = FALSE;
        goto clean_align_and_refine_subseq_with_ksw;
    }    
    qb = (qb) + query_from;
    qe = (qe) + query_from;
    sb = (sb) + subject_from;
    se = (se) + subject_from;

    //dump_align_string(qas, sas, aln_size, stderr);
    perc_identity = calc_ident_perc(qas, sas, aln_size, &dist, &score);
    eff_perc_identity = calc_effective_ident_perc(qas, sas, aln_size);

clean_align_and_refine_subseq_with_ksw:
    if (status) {
        result->qb = qb;
        result->qe = qe;
        result->sb = sb;
        result->se = se;
        result->qas = qas;
        result->qae = qae;
        result->sas = sas;
        result->sae = sae;
        result->dist = dist;
        result->score = score;
        result->perc_identity = perc_identity;
        result->eff_perc_identity = eff_perc_identity;
    }
    kv_destroy(query_v);
    kv_destroy(subject_v);
    kv_destroy(hit_seed_list);
    return status;
}
