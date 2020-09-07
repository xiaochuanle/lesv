#include "resolve_broken_ovlp.h"

BOOL 
is_correct_ovlp_positions(const int qb, const int qe, const int qs,
    const int sb, const int se, const int ss)
{
    const int E = 300;
    if (qb <= E && qe + E >= qs) return TRUE;
    if (sb <= E && se + E >= ss) return TRUE;
    if (qs - qe <= E && sb <= E) return TRUE;
    if (ss - se <= E && qb <= E) return TRUE;
    return FALSE;
}

BOOL
is_overhang_ovlp(const int qb, const int qe, const int qs, 
    const int sb, const int se, const int ss)
{
    if (is_correct_ovlp_positions(qb, qe, qs, sb, se, ss)) return FALSE;
    const int E = 300;
    int r = (qb <= E) || (qs - qe <= E) || (sb <= E) || (ss - se <= E);
    return r;
}

int
s_extend_hits_for_broken_map(HbnTracebackData* tbck_data, 
    SupportReadHitFindData* hit_finder,
    HbnInitHit* hit_array, 
    int hit_count,
    vec_chain_seed* chain_seed_list,
    const char* read_name,
    const int read_length,
    const int subject_id,
    const char* subject_name,
    const u8* subject,
    const int subject_length,
    const int min_cov_res,
    const double min_perc_coverage,
    const double min_perc_identity,
    int* qdir_,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln,
    vec_m4* m4_list)
{
    kv_clear(*m4_list);
    M4Record m4; memset(&m4, 0, sizeof(M4Record));
    for (int i = 0; i < hit_count; ++i) {
        HbnInitHit* hit = hit_array + i;
        int query_id = 0;
        int query_dir = -1;
        const char* query_name = NULL;
        int query_length = 0;
        const u8* query = NULL;
        ChainSeed* csa = NULL;
        int csc = 0;
        SupportReadHitFindData_SetupCnsAlignInfoFromHit(hit_finder,
            hit, chain_seed_list,
            &query_id, &query_dir, &query_length, &query_name, &query);
        csa = kv_data(*chain_seed_list);
        csc = kv_size(*chain_seed_list);
        if (csc == 0) continue;
        validate_mem(HBN_LOG_ARGS_DEFAULT, query, subject, csa, csc);
        int sb = csa[0].soff;
        int se = csa[csc-1].soff + csa[csc-1].length;
        hbn_assert(sb >= 0);
        int qb, qe;
        double perc_identity;
        double eff_perc_identity;
        int r;
        r = cns_extend_chain_seed_array_with_ksw(tbck_data,
                    tbck_data->ksw,
                    csa,
                    csc,
                    min_cov_res,
                    min_perc_coverage,
                    0.0,
                    query_id,
                    query,
                    query_length,
                    subject_id,
                    subject,
                    subject_length,
                    qaln,
                    saln,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity);
        if (!r) continue;
        //dump_align_string(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), stderr);
        //print_init_hit_range(hit, 1, query_length, subject_length);
        //HBN_LOG("ksw map [%d:%s, %d, %d, %d] x [%d:%s, %d, %d, %d] %g", query_id, read_name, qb, qe, query_length,
        //    subject_id, subject_name, sb, se, subject_length, perc_identity);
        eff_perc_identity = calc_effective_ident_perc(ks_s(*qaln), ks_s(*saln), ks_size(*qaln));
        if (eff_perc_identity < min_perc_identity) continue;
        if (is_correct_ovlp_positions(qb, qe, query_length, sb, se, subject_length)) {
            *qdir_ = hit->qdir;
            *qoff_ = qb;
            *qend_ = qe;
            *soff_ = sb;
            *send_ = se;
            *perc_identity_ = perc_identity;
            *eff_perc_identity_ = eff_perc_identity;
            return TRUE;
        }
        m4.qdir = hit->qdir;
        m4.qoff = qb;
        m4.qend = qe;
        m4.qsize = query_length;
        m4.soff = sb;
        m4.send = se;
        m4.ssize = subject_length;
        m4.ident_perc = eff_perc_identity;
        kv_push(M4Record, *m4_list, m4);
    }
    return FALSE;
}

static BOOL 
two_m4s_are_dual(M4Record* m1, 
    M4Record* m2,
    int* qdir_,
    int* qfrom_,
    int* qto_,
    int* sfrom_,
    int* sto_)
{
    if (m1->qdir != m2->qdir) return FALSE;
    const int E = 50;

    /// case 1: soff, look for send or qend
    ///           m1q  _______*___________*____________
    ///           m1s         *___________*____________________________
    if (m1->soff <= E) {
        /// case 1-1: 
        ///       m2q             ______*________________*_____
        ///       m2s          _________*________________*
        if (m2->ssize - m2->send <= E) {
            if (m2->qend > m1->qoff) {
                *qdir_ = m1->qdir;
                *qfrom_ = m1->qoff;
                *qto_ = m2->qend;
                *sfrom_ = m1->soff;
                *sto_ = m2->send;
                return TRUE;
            }
        }
        /// case 1-2:
        ///       m2q          _______*___________________*
        ///       m2s              ___*___________________*____
        if (m2->qsize - m2->qend <= E) {
            if (m2->send > m1->soff) {
                *qdir_ = m1->qdir;
                *qfrom_ = m1->qoff;
                *qto_ = m2->qend;
                *sfrom_ = m1->soff;
                *sto_ = m2->send;
                return TRUE;
            }
        }
    }

    /// case 2: send, look for soff or qoff
    ///           m1q               -----*---------------*----
    ///           m1s          __________*_______________*
    if (m1->ssize - m1->send <= E) {
        /// case 2-1
        ///       m2q    _____*__________________*_________    
        ///       m2s         *__________________*_______
        if (m2->soff <= E) {
            if (m2->qoff < m1->qend) {
                *qdir_ = m1->qdir;
                *qfrom_ = m2->qoff;
                *qto_ = m1->qend;
                *sfrom_ = m2->soff;
                *sto_ = m2->send;
                return TRUE;
            }
        }
        /// case 2-2
        ///       m2q       *___________________*_______
        ///       m2s  _____*___________________*___
        if (m2->qoff <= E) {
            if (m2->soff < m1->send) {
                *qdir_ = m1->qdir;
                *qfrom_ = m2->qoff;
                *qto_ = m1->qend;
                *sfrom_ = m2->soff;
                *sto_ = m1->send;
                return TRUE;
            }
        }
    }

    /// case 3: qoff, look for qend or send
    ///           m1q            *___________*____________
    ///           m1s      ______*___________*________
    if (m1->qoff <= E) {
        /// case 3-1: qend
        ///       m2q           ______*_______________________*
        ///       m2s     ____________*_______________________*_______
        if (m2->qsize - m2->qend <= E) {
            if (m1->soff < m2->send) {
                *qdir_ = m1->qdir;
                *qfrom_ = m1->qoff;
                *qto_ = m2->qend;
                *sfrom_ = m1->soff;
                *sto_ = m2->send;
                return TRUE;
            }
        }
        /// case 3-2: send
        ///       m2q            _____*_____________________*______ 
        ///       m2s        _________*_____________________*
        if (m2->ssize - m2->send <= E) {
            if (m1->qoff < m2->qend) {
                *qdir_ = m1->qdir;
                *qfrom_ = m1->qoff;
                *qto_ = m2->qend;
                *sfrom_ = m1->soff;
                *sto_ = m2->send;
                return TRUE;
            }
        }
    }

    /// case 4: qend, look for qoff or soff
    ///           m1q            -------*---------------*
    ///           m1s                ___*_______________*_______
    if (m1->qsize - m1->qend <= E) {
        /// case 4-1: qoff
        ///       m2q       *______________________*______
        ///       m2s   ____*______________________*_________
        if (m2->qoff <= E) {
            if (m2->soff < m1->send) {
                *qdir_ = m1->qdir;
                *qfrom_ = m2->qoff;
                *qto_ = m1->qend;
                *sfrom_ = m2->soff;
                *sto_ = m1->send;
                return TRUE;
            }
        }
        /// case 4-2: soff
        ///       m2q    ______*_________________*___
        ///       m2s          *_________________*_______
        if (m2->soff <= E) {
            if (m2->qoff < m1->qend) {
                *qdir_ = m1->qdir;
                *qfrom_ = m2->qoff;
                *qto_ = m1->qend;
                *sfrom_ = m2->soff;
                *sto_ = m1->send;
                return TRUE;
            }
        }
    }

    return FALSE;
}

static BOOL
s_chain_dual_m4s(HbnTracebackData* tbck_data, 
    SupportReadHitFindData* hit_finder, 
    M4Record* m1, M4Record* m2,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const u8* subject,
    const int subject_length,
    const double min_perc_identity,
    int qdir,
    int qoff,
    int qend,
    int soff,
    int send,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln)
{
    //HBN_LOG("chaining");
    //DUMP_M4_RECORD(fprintf, stderr, *m1);
    //DUMP_M4_RECORD(fprintf, stderr, *m2);
    int qb, qe, sb, se;
    double perc_identity;
    double eff_perc_identity;

    const u8* query = (qdir == FWD) ? fwd_read : rev_read;
    int band_width = hbn_max(qend - qoff, send - soff);
    band_width *= 0.2;
    int r = nw_ksw2_extd2(tbck_data->ksw,
                0,
                query,
                qoff,
                qend,
                qend - qoff,
                0,
                subject,
                soff,
                send,
                send - soff,
                0,
                0,
                band_width,
                &qb,
                &qe,
                &sb,
                &se,
                &perc_identity,
                &tbck_data->qabuf,
                &tbck_data->sabuf);
    if (!r) {
        band_width = hbn_max(qend - qoff, send - soff);
        r = nw_ksw2_extd2(tbck_data->ksw,
                0,
                query,
                qoff,
                qend,
                qend - qoff,
                0,
                subject,
                soff,
                send,
                send - soff,
                0,
                0,
                band_width,
                &qb,
                &qe,
                &sb,
                &se,
                &perc_identity,
                &tbck_data->qabuf,
                &tbck_data->sabuf);        
    }
    if (!r) return FALSE;

    int aln_size = ks_size(tbck_data->qabuf);
    const char* qas = ks_s(tbck_data->qabuf);
    const char* qae = qas + aln_size;
    const char* sas = ks_s(tbck_data->sabuf);
    const char* sae = sas + aln_size;
    r = truncate_align_bad_ends(qas, sas, aln_size, &qb, &qe, &sb, &se, &qas, &qae, &sas, &sae);
    if (!r) return r;

    perc_identity = calc_ident_perc(qas, sas, aln_size, NULL, NULL);
    eff_perc_identity = calc_effective_ident_perc(qas, sas, aln_size);
    if (eff_perc_identity < min_perc_identity) return FALSE;

    if (!is_correct_ovlp_positions(qb, qe, read_length, sb, se, subject_length)) return FALSE;

    *qoff_ = qb;
    *qend_ = qe;
    *soff_ = sb;
    *send_ = se;
    *perc_identity_ = perc_identity;
    *eff_perc_identity_ = eff_perc_identity;
    ks_clear(*qaln); kputsn(qas, qae - qas, qaln);
    ks_clear(*saln); kputsn(sas, sae - sas, saln);
    return TRUE;
}

static int
s_find_dual_m4(HbnTracebackData* tbck_data, 
    SupportReadHitFindData* hit_finder, 
    M4Record* m4a, 
    int m4c,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const u8* subject,
    const int subject_length,
    const double min_perc_identity,
    int* qdir_,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln)
{
    //for (int i = 0; i < m4c; ++i) DUMP_M4_RECORD(fprintf, stderr, m4a[i]);
     
    for (int i = 0; i < m4c - 1; ++i) {
        M4Record* mi = m4a + i;
        for (int j = i + 1; j < m4c; ++j) {
            M4Record* mj = m4a + j;
            if (two_m4s_are_dual(mi, mj, qdir_, qoff_, qend_, soff_, send_)) {
                if (s_chain_dual_m4s(tbck_data, hit_finder, mi, mj, fwd_read, rev_read, read_length,
                        subject, subject_length,
                        min_perc_identity, *qdir_, *qoff_, *qend_,
                        *soff_, *send_, 
                        qoff_, qend_, soff_, send_,
                        perc_identity_, eff_perc_identity_,
                        qaln, saln)) {
                    return TRUE;
                }
            }
            if (two_m4s_are_dual(mj, mi, qdir_, qoff_, qend_, soff_, send_)) {
                if (s_chain_dual_m4s(tbck_data, hit_finder, mi, mj, fwd_read, rev_read, read_length,
                        subject, subject_length,
                        min_perc_identity, *qdir_, *qoff_, *qend_,
                        *soff_, *send_, 
                        qoff_, qend_, soff_, send_,
                        perc_identity_, eff_perc_identity_,
                        qaln, saln)) {
                    return TRUE;
                }                
            }
        }
    }
    return FALSE;
}

int
resolve_broken_ovlp(HbnTracebackData* tbck_data, 
    SupportReadHitFindData* hit_finder, 
    vec_chain_seed* chain_seed_list,
    HbnInitHit* hit_array, 
    int hit_count,
    const char* read_name,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const int subject_id,
    const char* subject_name,
    const u8* subject,
    const int subject_length,
    const int min_cov_res,
    const double min_perc_coverage,
    const double min_perc_identity,
    int* qdir_,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln)
{
    kv_dinit(vec_m4, m4_list);
    if (s_extend_hits_for_broken_map(tbck_data, hit_finder, hit_array, hit_count,
            chain_seed_list,
            read_name, read_length,
            subject_id, subject_name, subject, subject_length,
            min_cov_res, min_perc_coverage, min_perc_identity,
            qdir_, qoff_, qend_, soff_, send_, 
            perc_identity_, eff_perc_identity_,
            qaln, saln, &m4_list)) {
        return TRUE;
    }
    int r = s_find_dual_m4(tbck_data,
                hit_finder,
                kv_data(m4_list),
                kv_size(m4_list),
                fwd_read,
                rev_read,
                read_length,
                subject,
                subject_length,
                min_perc_identity,
                qdir_,
                qoff_,
                qend_,
                soff_,
                send_,
                perc_identity_,
                eff_perc_identity_,
                qaln,
                saln);
    return r;
}