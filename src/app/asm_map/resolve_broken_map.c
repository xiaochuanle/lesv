#include "resovle_broken_map.h"

#include "../../algo/hbn_traceback_aux.h"
#include "../../ncbi_blast/setup/hsp2string.h"

static BOOL 
two_hsps_are_dual(BlastHSP* h1,
    BlastHSP* h2,
    int* qdir_,
    int* qfrom_,
    int* qto_,
    int* sfrom_,
    int* sto_)
{
    int qdir1 = h1->hbn_query.strand;
    int qdir2 = h2->hbn_query.strand;
    if (qdir1 != qdir2) return FALSE;

    int qb1 = h1->hbn_query.offset;
    int qe1 = h1->hbn_query.end;
    int qb2 = h2->hbn_query.offset;
    int qe2 = h2->hbn_query.end;
    hbn_assert(h1->hbn_query.seq_size == h2->hbn_query.seq_size);
    int qsize = h1->hbn_query.seq_size;
    int sb1 = h1->hbn_subject.offset;
    int se1 = h1->hbn_subject.end;
    int sb2 = h2->hbn_subject.offset;
    int se2 = h2->hbn_subject.end;
    hbn_assert(h1->hbn_subject.strand == FWD);
    hbn_assert(h2->hbn_subject.strand == FWD);
    hbn_assert(h1->hbn_subject.seq_size == h2->hbn_subject.seq_size);
    int ssize = h1->hbn_subject.seq_size;
    
    const int E = 50;

    /// case 1: soff, look for send or qend
    ///           m1q  _______*___________*____________
    ///           m1s         *___________*____________________________
    if (sb1 <= E) {
        /// case 1-1: 
        ///       m2q             ______*________________*_____
        ///       m2s          _________*________________*
        if (ssize - se2 <= E) {
            if (qe2 > qb1) {
                *qdir_ = qdir1;
                *qfrom_ = qb1;
                *qto_ = qe2;
                *sfrom_ = sb1;
                *sto_ = se2;
                return TRUE;
            }
        }
        /// case 1-2:
        ///       m2q          _______*___________________*
        ///       m2s              ___*___________________*____
        if (qsize - qe2 <= E) {
            if (se2 > sb1) {
                *qdir_ = qdir1;
                *qfrom_ = qb1;
                *qto_ = qe2;
                *sfrom_ = sb1;
                *sto_ = se2;
                return TRUE;
            }
        }
    }

    /// case 2: send, look for soff or qoff
    ///           m1q               -----*---------------*----
    ///           m1s          __________*_______________*
    if (ssize - se1 <= E) {
        /// case 2-1
        ///       m2q    _____*__________________*_________    
        ///       m2s         *__________________*_______
        if (sb2 <= E) {
            if (qb2 < qe1) {
                *qdir_ = qdir1;
                *qfrom_ = qb2;
                *qto_ = qe1;
                *sfrom_ = sb2;
                *sto_ = se1;
                return TRUE;
            }
        }
        /// case 2-2
        ///       m2q       *___________________*_______
        ///       m2s  _____*___________________*___
        if (qb2 <= E) {
            if (sb2 < se1) {
                *qdir_ = qdir1;
                *qfrom_ = qb2;
                *qto_ = qe1;
                *sfrom_ = sb2;
                *sto_ = se1;
                return TRUE;
            }
        }
    }

    /// case 3: qoff, look for qend or send
    ///           m1q            *___________*____________
    ///           m1s      ______*___________*________
    if (qb1 <= E) {
        /// case 3-1: qend
        ///       m2q           ______*_______________________*
        ///       m2s     ____________*_______________________*_______
        if (qsize - qe2 <= E) {
            if (sb1 < se2) {
                *qdir_ = qdir1;
                *qfrom_ = qb1;
                *qto_ = qe2;
                *sfrom_ = sb1;
                *sto_ = se2;
                return TRUE;
            }
        }
        /// case 3-2: send
        ///       m2q            _____*_____________________*______ 
        ///       m2s        _________*_____________________*
        if (ssize - se2 <= E) {
            if (qb1 < qe2) {
                *qdir_ = qdir1;
                *qfrom_ = qb1;
                *qto_ = qe2;
                *sfrom_ = sb1;
                *sto_ = se2;
                return TRUE;
            }
        }
    }

    /// case 4: qend, look for qoff or soff
    ///           m1q            -------*---------------*
    ///           m1s                ___*_______________*_______
    if (qsize - qe1 <= E) {
        /// case 4-1: qoff
        ///       m2q       *______________________*______
        ///       m2s   ____*______________________*_________
        if (qb2 <= E) {
            if (sb2 < se1) {
                *qdir_ = qdir1;
                *qfrom_ = qb2;
                *qto_ = qe1;
                *sfrom_ = sb2;
                *sto_ = se1;
                return TRUE;
            }
        }
        /// case 4-2: soff
        ///       m2q    ______*_________________*___
        ///       m2s          *_________________*_______
        if (sb2 <= E) {
            if (qb2 < qe1) {
                *qdir_ = qdir1;
                *qfrom_ = qb2;
                *qto_ = qe1;
                *sfrom_ = sb2;
                *sto_ = se1;
                return TRUE;
            }
        }
    }

    return FALSE;
}

static BOOL 
is_correct_ovlp(const int qb, const int qe, const int qs, 
    const int sb, const int se, const int ss)
{
    const int E = 100;
    if (qb <= E && qs - qe <= E) return TRUE;
    if (sb <= E && ss - se <= E) return TRUE;
    if (ss - se <= E && qb <= E) return TRUE;
    if (qs - qe <= E && sb <= E) return TRUE;
    return FALSE;
}

static BOOL
s_chain_dual_hsps(HbnTracebackData* tbck_data, 
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
    int qb, qe, sb, se;
    double perc_identity;
    double eff_perc_identity;

    const int E = 100;
    if (qoff <= E) qoff = 0;
    if (read_length - qend <= E) qend = read_length;
    if (soff <= E) soff = 0;
    if (subject_length - send <= E) send = subject_length;

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
    if (!is_correct_ovlp(qb, qe, read_length, sb, se, subject_length)) return FALSE;

    int aln_size = ks_size(tbck_data->qabuf);
    const char* qas = ks_s(tbck_data->qabuf);
    const char* qae = qas + aln_size;
    const char* sas = ks_s(tbck_data->sabuf);
    const char* sae = sas + aln_size;
    //r = truncate_align_bad_ends(qas, sas, aln_size, &qb, &qe, &sb, &se, &qas, &qae, &sas, &sae);
    //if (!r) return r;

    perc_identity = calc_ident_perc(qas, sas, aln_size, NULL, NULL);
    eff_perc_identity = calc_effective_ident_perc(qas, sas, aln_size);
    if (eff_perc_identity < min_perc_identity) return FALSE;

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
s_find_dual_hsps(HbnTracebackData* data, 
    BlastHSP* hsp_array,
    int hsp_count,
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
    int* score_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln)
{
    BOOL dump_chain_hsps = 0;
    for (int i = 0; i < hsp_count - 1; ++i) {
        BlastHSP* hi = hsp_array + i;
        for (int j = i + 1; j < hsp_count; ++j) {
            BlastHSP* hj = hsp_array + j;
            if (two_hsps_are_dual(hi, hj, qdir_, qoff_, qend_, soff_, send_)) {
                if (s_chain_dual_hsps(data, fwd_read, rev_read, read_length,
                        subject, subject_length,
                        min_perc_identity, *qdir_, *qoff_, *qend_,
                        *soff_, *send_, 
                        qoff_, qend_, soff_, send_,
                        perc_identity_, eff_perc_identity_,
                        qaln, saln)) {
                    *score_ = hi->hsp_info.raw_score + hj->hsp_info.raw_score;
                    if (dump_chain_hsps) {
                        HBN_LOG("chain two hsps:");
                        dump_blasthsp(fprintf, stderr, *hi);
                        dump_blasthsp(fprintf, stderr, *hj);
                        HBN_LOG("into [%d, %d, %d] x [%d, %d, %d], score = %d, %g, %g", *qoff_, *qend_, read_length,
                            *soff_, *send_, subject_length, *score_, *perc_identity_, *eff_perc_identity_);
                    }
                    return TRUE;
                }
            }
            if (two_hsps_are_dual(hj, hi, qdir_, qoff_, qend_, soff_, send_)) {
                if (s_chain_dual_hsps(data, fwd_read, rev_read, read_length,
                        subject, subject_length,
                        min_perc_identity, *qdir_, *qoff_, *qend_,
                        *soff_, *send_, 
                        qoff_, qend_, soff_, send_,
                        perc_identity_, eff_perc_identity_,
                        qaln, saln)) {
                    *score_ = hi->hsp_info.raw_score + hj->hsp_info.raw_score;
                    if (dump_chain_hsps) {
                        HBN_LOG("chain two hsps:");
                        dump_blasthsp(fprintf, stderr, *hi);
                        dump_blasthsp(fprintf, stderr, *hj);
                        HBN_LOG("into [%d, %d, %d] x [%d, %d, %d], score = %d, %g, %g", *qoff_, *qend_, read_length,
                            *soff_, *send_, subject_length, *score_, *perc_identity_, *eff_perc_identity_);
                    }
                    return TRUE;
                }                
            }
        }
    }
    return FALSE;
}

int
resolve_broken_map(HbnTracebackData* tbck_data, 
    BlastHSP* hsp_array,
    int hsp_cnt,
    const char* read_name,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const int subject_id,
    const char* subject_name,
    const u8* subject,
    const int subject_length,
    const double min_perc_identity,
    int* qdir_,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    int* score_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln)
{
    int r = s_find_dual_hsps(tbck_data,
                hsp_array,
                hsp_cnt,
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
                score_,
                perc_identity_,
                eff_perc_identity_,
                qaln,
                saln);
    return r;
}

