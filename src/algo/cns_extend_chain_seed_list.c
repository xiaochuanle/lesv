#include "cns_extend_chain_seed_list.h"

#include "refine_align.h"

BOOL
check_ovlp_position(const int qb, const int qe, const int qs,
    const int sb, const int se, const int ss)
{
    const int E = 300;
    if (qb <= E && qs - qe <= E) return TRUE;
    if (sb <= E && ss - se <= E) return TRUE;
    if (qs - qe <= E && sb <= E) return TRUE;
    if (qb <= E && ss - se <= E) return TRUE;
    return FALSE;
}

BOOL
check_ovlp_cov(const int qb, const int qe, const int qs,
						const int sb, const int se, const int ss,
                        int res_coverage,
						double perc_coverage)
{
    const double ratio = perc_coverage / 100.0;
	const int oq = qe - qb;
	const int qqs = qs * ratio;
	const int os = se - sb;
	const int qss = ss * ratio;
	return (oq >= res_coverage) || (os >= res_coverage) || (oq >= qqs) || (os >= qss);
}

BOOL
check_chain_seed_array_cov(ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int read_length,
    const int subject_length,
    const int res_coverage,
    const double perc_coverage)
{
    ChainSeed cs = chain_seed_array[0];
    int qb = cs.qoff;
    int sb = cs.soff;
    cs = chain_seed_array[chain_seed_count-1];
    int qe = cs.qoff + cs.length;
    int se = cs.soff + cs.length;
    double ratio = perc_coverage / 100.0 - 0.15;
    if (ratio < .0001) ratio = .0001;
    int r = (qe - qb >= res_coverage * 0.85)
            ||
            (se - sb >= res_coverage * 0.85)
            ||
            (qe - qb >= read_length * ratio) 
            || 
            (se - sb >= subject_length * ratio);
    return r;
}

BOOL
subject_subseq_cov_is_full(u8* cov_stats, int soff, int send, const int kMaxCov)
{
    int n = 0;
    for (int i = soff; i < send; ++i) 
        if (cov_stats[i] >= kMaxCov) ++n;
    if (send - soff >= n + 200) return FALSE;
    return TRUE;
}

void 
normalize_gaps(const char* qstr, 
    const char* tstr, 
    const int aln_size, 
    kstring_t* qnorm, 
    kstring_t* tnorm, 
    const BOOL push)
{
    ks_clear(*qnorm);
    ks_clear(*tnorm);

#ifndef NDEBUG
    int qcnt = 0, tcnt = 0;
    for (int i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != GAP_CHAR) ++qcnt;
        if (tc != GAP_CHAR) ++tcnt;
    }
#endif

    // convert mismatches to indels
    for (int i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != tc && qc != GAP_CHAR && tc != GAP_CHAR) {
            kputc(GAP_CHAR, qnorm);
            kputc(qc, qnorm);
            kputc(tc, tnorm);
            kputc(GAP_CHAR, tnorm);
        } else { 
            kputc(qc, qnorm);
            kputc(tc, tnorm);
        }
    }

    // push gaps to the right, but not pass the end
    if (push)
    {
        int qlen = ks_size(*qnorm);
        char* qn = ks_s(*qnorm);
        int tlen = ks_size(*tnorm);
        char* tn = ks_s(*tnorm);
        for (int i = 0; i < qlen - 1; ++i)
        {
            // push target gaps
            if (tn[i] == GAP_CHAR)
            {
                int j = i;
                while (1)
                {
                    const char c = tn[++j];
                    if (c != GAP_CHAR || j > qlen - 1)
                    {
                        if (c == qn[i]) { tn[i] = c; tn[j] = GAP_CHAR; }
                        break;
                    }
                }
            }
            // push query gaps
            if (qn[i] == GAP_CHAR)
            {
                int j = i;
                while (1)
                {
                    const char c = qn[++j];
                    if (c != GAP_CHAR || j > tlen - 1)
                    {
                        if (c == tn[i]) { qn[i] = c; qn[j] = GAP_CHAR; }
                        break;
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    int qcnt2 = 0, tcnt2 = 0;
    for (size_t i = 0; i < ks_size(*qnorm); ++i) 
        if (ks_A(*qnorm, i) != GAP_CHAR) ++qcnt2;
    for (size_t i = 0; i < ks_size(*tnorm); ++i)
        if (ks_A(*tnorm, i) != GAP_CHAR) ++tcnt2;
    hbn_assert(qcnt == qcnt2);
    hbn_assert(tcnt == tcnt2);
#endif
}

BOOL
cns_extend_chain_seed_array(HbnTracebackData* traceback_data,
    Ksw2Data* ksw_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int min_cov_res,
    const double min_perc_coverage,
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
            1,
            min_perc_identity,
            FALSE);
    //HbnTracebackDataDump(fprintf, stderr, traceback_data);
    if (!r) {
        return FALSE;
    }

    int qoff = traceback_data->qoff;
    int qend = traceback_data->qend;
    int soff = traceback_data->soff;
    int send = traceback_data->send;
    if (!check_ovlp_position(qoff, qend, read_length, soff, send, subject_length)) {
        if (0)
        HBN_LOG("invalid ovlp position  [%d, %d, %d, %d] x [%d, %d, %d] %g",
          read_id, qoff, qend, read_length,
            soff, send, subject_length, traceback_data->ident_perc);
        return FALSE;
    }
    r = check_ovlp_cov(qoff, qend, read_length,
            soff, send, subject_length, 
            min_cov_res,
            min_perc_coverage);
    if (!r) {
        if (0)
        HBN_LOG("cov test fail  [%d, %d, %d, %d] x [%d, %d, %d] %g",
          read_id, qoff, qend, read_length,
            soff, send, subject_length, traceback_data->ident_perc);
        return r;
    }

    const char* qas = traceback_data->qas;
    const char* sas = traceback_data->sas;
    int aln_size = traceback_data->qae - qas;
    const char* qae = qas + aln_size;
    const char* sae = sas + aln_size;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            0, read, qoff, qend, traceback_data->qas,
            0, subject, soff, send, sas,
            aln_size, TRUE);

    const char* new_qas = NULL;
    const char* new_qae = NULL;
    const char* new_sas = NULL;
    const char* new_sae = NULL;

    if (!truncate_align_bad_ends(qas, sas, aln_size,
            &qoff, &qend, &soff, &send, 
            &new_qas, &new_qae, &new_sas, &new_sae)) {
        //HBN_LOG("truncate fail");
        return FALSE;
    }
    qas = new_qas;
    qae = new_qae;
    sas = new_sas;
    sae = new_sae;
    aln_size = qae - qas;

#if 0
    normalize_gaps(qas,
        sas,
        aln_size,
        qaln,
        saln,
        TRUE);
#else
    ks_clear(*qaln);
    kputsn(qas, aln_size, qaln);
    ks_clear(*saln);
    kputsn(sas, aln_size, saln);
#endif 

    *qoff_ = qoff;
    *qend_ = qend;
    *soff_ = soff;
    *send_ = send;
    *ident_perc_ = calc_ident_perc(qas, sas, aln_size, NULL, NULL);
    return TRUE;
}

BOOL
cns_extend_chain_seed_array_with_refine(HbnTracebackData* traceback_data,
    Ksw2Data* ksw_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int min_cov_res,
    const double min_perc_coverage,
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
            1,
            min_perc_identity,
            TRUE);
    if (0) HbnTracebackDataDump(fprintf, stderr, traceback_data);
    if (!r) {
        //HBN_LOG("csa extension fail, min_cov_res = %d, min_perc_identity = %g", min_cov_res, min_perc_identity);
        return FALSE;
    }
    int qoff = traceback_data->qoff;
    int qend = traceback_data->qend;
    int soff = traceback_data->soff;
    int send = traceback_data->send;
    if (!check_ovlp_position(qoff, qend, read_length, soff, send, subject_length)) {
        if (0)
        HBN_LOG("invalid ovlp position  [%d, %d, %d, %d] x [%d, %d, %d] %g",
          read_id, qoff, qend, read_length,
            soff, send, subject_length, traceback_data->ident_perc);
        return FALSE;
    }
    r = check_ovlp_cov(qoff, qend, read_length,
            soff, send, subject_length, 
            min_cov_res,
            min_perc_coverage);
    if (!r) {
        if (0)
        HBN_LOG("cov test fail  [%d, %d, %d, %d] x [%d, %d, %d] %g",
          read_id, qoff, qend, read_length,
            soff, send, subject_length, traceback_data->ident_perc);
        return r;
    }

    const char* qas = traceback_data->qas;
    const char* sas = traceback_data->sas;
    int aln_size = traceback_data->qae - qas;
    const char* qae = qas + aln_size;
    const char* sae = sas + aln_size;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            0, read, qoff, qend, traceback_data->qas,
            0, subject, soff, send, sas,
            aln_size, TRUE);

    const char* new_qas = NULL;
    const char* new_qae = NULL;
    const char* new_sas = NULL;
    const char* new_sae = NULL;

    if (!truncate_align_bad_ends(qas, sas, aln_size,
            &qoff, &qend, &soff, &send, 
            &new_qas, &new_qae, &new_sas, &new_sae)) {
        //HBN_LOG("truncate fail");
        return FALSE;
    }
    qas = new_qas;
    qae = new_qae;
    sas = new_sas;
    sae = new_sae;
    aln_size = qae - qas;

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
    if (eff_perc_identity < min_perc_identity) {
        HBN_LOG("eff_perc_identity = %g, min = %g", eff_perc_identity, min_perc_identity);
        return FALSE;
    }

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

BOOL 
cns_extend_chain_seed_array_with_ksw(HbnTracebackData* traceback_data,
    Ksw2Data* ksw2_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int min_cov_res,
    const double min_perc_coverage,
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
            min_cov_res,
            min_perc_identity,
            FALSE);
    //HbnTracebackDataDump(fprintf, stderr, traceback_data);
    //HBN_LOG("perc_identity = %g", traceback_data->ident_perc);
    if (!r) return FALSE;
    int qoff = traceback_data->qoff;
    int qend = traceback_data->qend;
    int soff = traceback_data->soff;
    int send = traceback_data->send;
    r = check_ovlp_cov(qoff, qend, read_length,
            soff, send, subject_length, 
            min_cov_res, min_perc_coverage);
    if (!r) return r;
    
    int distance = hbn_max(qend - qoff, send - soff);
    distance = 1.2 * distance * (100.0 - traceback_data->ident_perc) / 100.0;
    if (distance == 0) distance = (send + qend - qoff - soff) * 0.05;
    //HBN_LOG("distance = %d", distance);

    double ident_perc = -1.0;
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
            min_cov_res,
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
                min_cov_res,
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