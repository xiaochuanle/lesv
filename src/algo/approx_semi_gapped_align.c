#include "approx_semi_gapped_align.h"

#include "ksw2.h"
#include "hbn_traceback_aux.h"
#include "../corelib/fasta.h"

static const int kMatLen = 8;
static const int kMaxDelta = 5000;

ApproxSemeGappedAlignData*
ApproxSemeGappedAlignDataNew()
{
    ApproxSemeGappedAlignData* data = (ApproxSemeGappedAlignData*)calloc(1, sizeof(ApproxSemeGappedAlignData));

    data->ksw = Ksw2DataNew();
    data->ksw->reward = 2;
    data->ksw->penalty = 8;
    data->ksw->go = 4;
    data->ksw->ge = 2;
    data->ksw->go1 = 24;
    data->ksw->ge1 = 1;
    data->ksw->ambi_penalty = 1;
    data->ksw->zdrop = 30000;
    data->ksw->band_width = 30000;
    data->ksw->end_bonus = -1;
    ksw_gen_simple_mat(5, data->ksw->mat, data->ksw->reward, data->ksw->penalty, data->ksw->ambi_penalty);

    data->qid = -1;
    data->sid = -1;
    data->fwd_query = NULL;
    data->query_length = 0;
    data->fwd_subject = NULL;
    data->subject_length = 0;
    kv_init(data->rev_query);
    kv_init(data->rev_subject);
    ks_init(data->qaux);
    ks_init(data->saux);
    ks_init(data->fwd_qaln);
    ks_init(data->fwd_saln);
    ks_init(data->rev_qaln);
    ks_init(data->rev_saln);
    ks_init(data->bp_qaln);
    ks_init(data->bp_saln);
    ks_init(data->qaln);
    ks_init(data->saln);
    return data;
}

ApproxSemeGappedAlignData*
ApproxSemeGappedAlignDataFree(ApproxSemeGappedAlignData* data)
{
    data->ksw = Ksw2DataFree(data->ksw);
    kv_destroy(data->rev_query);
    kv_destroy(data->rev_subject);
    ks_destroy(data->qaux);
    ks_destroy(data->saux);
    ks_destroy(data->fwd_qaln);
    ks_destroy(data->fwd_saln);
    ks_destroy(data->rev_qaln);
    ks_destroy(data->rev_saln);
    ks_destroy(data->bp_qaln);
    ks_destroy(data->bp_saln);
    ks_destroy(data->qaln);
    ks_destroy(data->saln);
    free(data);
    return NULL;
}

void
ApproxSemeGappedAlignData_Init(ApproxSemeGappedAlignData* data,
    const u8* query, const int query_length,
    const u8* subject, const int subject_length)
{
    data->fwd_query = query;
    data->query_length = query_length;
    data->fwd_subject = subject;
    data->subject_length = subject_length;
    kv_resize(u8, data->rev_query, query_length);
    for (int p = query_length - 1; p >= 0; --p) {
        kv_A(data->rev_query, p) = query[query_length - 1 - p];
    }
    kv_resize(u8, data->rev_subject, subject_length);
    for (int p = subject_length - 1; p >= 0; --p) {
        kv_A(data->rev_subject, p) = subject[subject_length - 1 - p];
    }
    ks_clear(data->fwd_qaln);
    ks_clear(data->fwd_saln);
    ks_clear(data->rev_qaln);
    ks_clear(data->rev_saln);
    ks_clear(data->bp_qaln);
    ks_clear(data->bp_saln);
    ks_clear(data->qaln);
    ks_clear(data->saln);
    data->f_qoff = 0;
    data->f_qend = -1;
    data->f_soff = 0;
    data->f_send = -1;
    data->r_qoff = 0;
    data->r_qend = -1;
    data->r_soff = 0;
    data->r_send = -1;
    data->perc_identity = -1.0;
    data->eff_perc_identity = -1.0;
    data->qoff = 0;
    data->qend = -1;
    data->soff = 0;
    data->send = -1;
}

static BOOL
calc_next_subseq_size(int qleft, int tleft, int* qsize, int* tsize)
{
	const int block_size = 50000;
	const int max_block_size = 60000;
	BOOL last_block = 0;
	if (qleft <= max_block_size || tleft <= max_block_size) {
		last_block = 1;
		*qsize = qleft;
		*tsize = tleft;
	} else {
		*qsize = block_size;
		*tsize = block_size;
	}
	return last_block;
}

static BOOL 
get_next_subseqs(int* qfrom, int* qto, const int qsize,
    int* sfrom, int* sto, int ssize)
{
    int qf = *qfrom;
    int ql = qsize - qf;
    int sf = *sfrom;
    int sl = ssize - sf;
    int qs = 0;
    int ss = 0;
    BOOL last_subseq = FALSE;

    last_subseq = calc_next_subseq_size(ql, sl, &qs, &ss);

    *qfrom = qf;
    *qto = qf + qs;
    *sfrom = sf;
    *sto = sf + ss;
    return last_subseq;
}

static void
ksw_right_extend(ApproxSemeGappedAlignData* data, 
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    int* qend,
    int* send,
    kstring_t* qaln,
    kstring_t* saln,
    const BOOL truncate_bad_ends)
{
    Ksw2Data* ksw = data->ksw;
    const int flag = KSW_EZ_EXTZ_ONLY;
    ksw_extz_t ez; 
    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_reset_extz(&ez);
    ksw_extd2_sse(ksw->km, query_length, query,
        subject_length, subject,
        5, ksw->mat, ksw->go, ksw->ge, ksw->go1, ksw->ge1, 
        ksw->band_width, ksw->zdrop, ksw->end_bonus, flag, &ez);

    *qend = 0;
    *send = 0;
    ks_clear(*qaln);
    ks_clear(*saln);
    int qi = 0;
    int si = 0;
    for (int k = 0; k < ez.n_cigar; ++k) {
        int op_num = ez.cigar[k]>>4;
        int op_type = "MIDN"[ez.cigar[k]&0xf];
        int c = 0;
        switch (op_type)
        {
        case 'M':
            for (int t = 0; t < op_num; ++t, ++qi, ++si) {
                c = query[qi];
                hbn_assert(c >= 0 && c < 4);
                c = DECODE_RESIDUE(c);
                kputc(c, qaln);
                c = subject[si];
                hbn_assert(c >= 0 && c < 4);
                c = DECODE_RESIDUE(c);
                kputc(c, saln);
            }
            break;
        case 'I':
            for (int t = 0; t < op_num; ++t, ++qi) {
                c = query[qi];
                hbn_assert(c >= 0 && c < 4);
                c = DECODE_RESIDUE(c);
                kputc(c, qaln);
                kputc(GAP_CHAR, saln);
            }
            break;
        case 'D':
            for (int t = 0; t < op_num; ++t, ++si) {
                kputc(GAP_CHAR, qaln);
                c = subject[si];
                hbn_assert(c >= 0 && c < 4);
                c = DECODE_RESIDUE(c);
                kputc(c, saln);
            }
            break;            
        default:
            HBN_LOG("invalid op_type: %d", op_type);
            break;
        }
    }
    kfree(data->km, ez.cigar);

    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        data->qid, query, 0, qi, ks_s(*qaln),
        data->sid, subject, 0, si, ks_s(*saln),
        ks_size(*qaln), TRUE);

    if (!truncate_bad_ends) {
        *qend = qi;
        *send = si;
        return;
    }

    int m = 0;
    int ai = ks_size(*qaln);
    while (ai && m < kMatLen) {
        char qc = ks_A(*qaln, ai-1);
        char sc = ks_A(*saln, ai-1);
        m = (qc == sc) ? (m+1) : 0;
        --ai;
        if (qc != GAP_CHAR) --qi;
        if (sc != GAP_CHAR) --si;
    }

    if (m < kMatLen) {
        qi = 0;
        si = 0;
        ks_clear(*qaln);
        ks_clear(*saln);
        while (qi < query_length && si < subject_length) {
            char qc = DECODE_RESIDUE(query[qi]);
            char sc = DECODE_RESIDUE(subject[si]);
            if (qc != sc) break;
            kputc(qc, qaln);
            kputc(sc, saln);
            ++qi;
            ++si;
        }
        *qend = qi;
        *send = si;
        return;
    }

    ai += kMatLen;
    qi += kMatLen;
    si += kMatLen;
    qaln->l = ai;
    saln->l = ai;
    *qend = qi;
    *send = si;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        data->qid, query, 0, qi, ks_s(*qaln),
        data->sid, subject, 0, si, ks_s(*saln),
        ks_size(*qaln), TRUE);
//dump_align_string(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), stderr);
}

static BOOL
s_find_left_break_point(ApproxSemeGappedAlignData* data, int* qoff, int* soff)
{
    int q_max = hbn_min(data->f_qend, data->query_length - data->r_qend);
    int s_max = hbn_min(data->f_send, data->subject_length - data->r_send);
    int qi = 0;
    int si = 0;
    int ai = 0;
    const char* qaln = ks_s(data->fwd_qaln);
    const char* saln = ks_s(data->fwd_saln);
    hbn_assert(ks_size(data->fwd_qaln) == ks_size(data->fwd_saln));
    int aln_size = ks_size(data->fwd_qaln);
    while (ai < aln_size && qi < q_max && si < s_max) {
        if (qaln[ai] != GAP_CHAR) ++qi;
        if (saln[ai] != GAP_CHAR) ++si;
        ++ai;
    }
    hbn_assert(ai <= aln_size);
    hbn_assert(qi == q_max || si == s_max);
    //HBN_LOG("******* qi = %d, q_max = %d, si = %d, s_max = %d",
    //    qi, q_max, si, s_max);
    
    if (qi <= kMaxDelta || si <= kMaxDelta) {
        *qoff = 0;
        *soff = 0;
        return TRUE;
    }

    int qd = 0, sd = 0;
    while (qd < kMaxDelta && sd < kMaxDelta) {
        --ai;
        if (qaln[ai] != GAP_CHAR) ++qd;
        if (saln[ai] != GAP_CHAR) ++sd;
    }
    hbn_assert(qd == kMaxDelta || sd == kMaxDelta);
    qi -= qd;
    si -= sd;

    int m = 0;
    while (qi < q_max && si < s_max) {
        hbn_assert(ai < aln_size);
        m = (qaln[ai] == saln[ai]) ? (m+1) : 0;
        if (qaln[ai] != GAP_CHAR) ++qi;
        if (saln[ai] != GAP_CHAR) ++si;
        ++ai;
        if (m == kMatLen) break;
    }

    if (m < kMatLen) return FALSE;
    qi -= kMatLen;
    si -= kMatLen;
    ai -= kMatLen;
    ++qi;
    ++si;
    ++ai;
    *qoff = qi;
    *soff = si;
    return TRUE;
}

static void
s_dump_seqs(ApproxSemeGappedAlignData* data)
{
    hbn_dfopen(out, "x.fasta", "w");
    fprintf(out, ">%d\n", data->qid);
    for (int i = 0; i < data->query_length; ++i) {
        int c = data->fwd_query[i];
        hbn_assert(c >= 0 && c < 4);
        c = DECODE_RESIDUE(c);
        fprintf(out, "%c", c);
    }
    fprintf(out, "\n");

    fprintf(out, ">%d\n", data->sid);
    for (int i = 0; i < data->subject_length; ++i) {
        int c = data->fwd_subject[i];
        hbn_assert(c >= 0 && c < 4);
        c = DECODE_RESIDUE(c);
        fprintf(out, "%c", c);
    }
    fprintf(out, "\n");    

    hbn_fclose(out);
}

static BOOL
s_find_right_break_point(ApproxSemeGappedAlignData* data, int* qoff, int* soff)
{
    int q_max = hbn_min(data->r_qend, data->query_length - data->f_qend);
    int s_max = hbn_min(data->r_send, data->subject_length - data->f_send);
    int qi = 0;
    int si = 0;
    int ai = 0;
    const char* qaln = ks_s(data->rev_qaln);
    const char* saln = ks_s(data->rev_saln);
    hbn_assert(ks_size(data->rev_qaln) == ks_size(data->rev_saln));
    const int aln_size = ks_size(data->rev_qaln);
    while (ai < aln_size && qi < q_max && si < s_max) {
        if (qaln[ai] != GAP_CHAR) ++qi;
        if (saln[ai] != GAP_CHAR) ++si;
        ++ai;
    }
    hbn_assert(ai <= aln_size);
    if (qi < q_max && si < s_max) s_dump_seqs(data);
    hbn_assert(qi == q_max || si == s_max, 
	"qid = %d, sid = %d, qi = %d, q_max = %d, si = %d, s_max = %d, ai = %d, aln_size = %d, r_qoff = %d, r_qend = %d",
	data->qid, data->sid, qi, q_max, si, s_max, ai, aln_size, data->r_qoff, data->r_qend);
    //HBN_LOG("******* qi = %d, q_max = %d, si = %d, s_max = %d",
    //    data->query_length - qi, data->query_length - q_max, data->subject_length - si, data->subject_length - s_max);

    {
        int x = 0, y = 0, z = 0;
        for (z = 0; z < ai; ++z) {
            if (qaln[z] != GAP_CHAR) ++x;
            if (saln[z] != GAP_CHAR) ++y;
        }
        hbn_assert(x == qi);
        hbn_assert(y == si);
        hbn_assert(z == ai);
    }

    if (qi <= kMaxDelta || si <= kMaxDelta) {
        *qoff = data->query_length;
        *soff = data->subject_length;
        return TRUE;
    }

    int qd = 0, sd = 0;
    while (qd < kMaxDelta && sd < kMaxDelta) {
        --ai;
        hbn_assert(ai >= 0);
        if (qaln[ai] != GAP_CHAR) ++qd;
        if (saln[ai] != GAP_CHAR) ++sd;
    }
    hbn_assert(qd == kMaxDelta || sd == kMaxDelta);
    qi -= qd;
    si -= sd;

    {
        int x = 0, y = 0, z = 0;
        for (z = 0; z < ai; ++z) {
            if (qaln[z] != GAP_CHAR) ++x;
            if (saln[z] != GAP_CHAR) ++y;
        }
        hbn_assert(x == qi);
        hbn_assert(y == si);
        hbn_assert(z == ai);
    }

    int m = 0;
    while (qi < q_max && si < s_max) {
        hbn_assert(ai < aln_size);
        hbn_assert(qaln[ai] != GAP_CHAR || saln[ai] != GAP_CHAR);
        m = (qaln[ai] == saln[ai]) ? (m+1) : 0;
        if (qaln[ai] != GAP_CHAR) ++qi;
        if (saln[ai] != GAP_CHAR) ++si;
        ++ai;
        if (m == kMatLen) break;
    }

    {
        int x = 0, y = 0, z = 0;
        for (z = 0; z < ai; ++z) {
            if (qaln[z] != GAP_CHAR) ++x;
            if (saln[z] != GAP_CHAR) ++y;
        }
        hbn_assert(x == qi);
        hbn_assert(y == si);
        hbn_assert(z == ai);
    }

    if (m < kMatLen) return FALSE;
    qi -= kMatLen;
    si -= kMatLen;
    ai -= kMatLen;

    {
        int x = 0, y = 0, z = 0;
        for (z = 0; z < ai; ++z) {
            if (qaln[z] != GAP_CHAR) ++x;
            if (saln[z] != GAP_CHAR) ++y;
        }
        hbn_assert(x == qi);
        hbn_assert(y == si);
        hbn_assert(z == ai);

        for (int x = 0; x < kMatLen; ++x) hbn_assert(qaln[ai+x] == saln[ai+x]);
    }

    ++ai;
    ++qi;
    ++si;
    *qoff = data->query_length - qi;
    *soff = data->subject_length - si;
    return TRUE;    
}

static void
s_ApproxSemeGappedAlignData_RightExtend(ApproxSemeGappedAlignData* data, int fwd_extend)
{
    const u8* query = NULL;
    const u8* subject = NULL;
    kstring_t* qaln = NULL;
    kstring_t* saln = NULL;
    int* qbeg = NULL;
    int* qend = NULL;
    int* sbeg = NULL;
    int* send = NULL;
    if (fwd_extend) {
        query = data->fwd_query;
        subject = data->fwd_subject;
        qaln = &data->fwd_qaln;
        saln = &data->fwd_saln;
        qbeg = &data->f_qoff;
        qend = &data->f_qend;
        sbeg = &data->f_soff;
        send = &data->f_send;
    } else {
        query = kv_data(data->rev_query);
        subject = kv_data(data->rev_subject);
        qaln = &data->rev_qaln;
        saln = &data->rev_saln;
        qbeg = &data->r_qoff;
        qend = &data->r_qend;
        sbeg = &data->r_soff;
        send = &data->r_send;        
    }

    const int qsize = data->query_length;
    const int ssize = data->subject_length;
    int qfrom = 0, qto, sfrom = 0, sto;
    
    while (1) {
        BOOL last_subseq = get_next_subseqs(&qfrom, &qto, qsize, &sfrom, &sto, ssize);
	HBN_LOG("extend [%d, %d] x [%d, %d], last_subseq = %d", qfrom, qto, sfrom, sto, last_subseq);
        const u8* q = query + qfrom;
        const int ql = qto - qfrom;
        const u8* s = subject + sfrom;
        const int sl = sto - sfrom;
        int qe = 0;
        int se = 0;
        ksw_right_extend(data, q, ql, s, sl, &qe, &se, &data->qaux, &data->saux, !last_subseq);
	HBN_LOG("qe = %d, ql = %d, se = %d, sl = %d", qe, ql, se, sl);

        if (last_subseq || (ql - qe > 2000 && sl - se > 2000)) {
            kputsn(ks_s(data->qaux), ks_size(data->qaux), qaln);
            kputsn(ks_s(data->saux), ks_size(data->saux), saln);
            *qbeg = 0;
            *qend = qfrom + qe;
            *sbeg = 0;
            *send = sfrom + se;
            break;
        } else {
            kputsn(ks_s(data->qaux), ks_size(data->qaux) - kMatLen, qaln);
            kputsn(ks_s(data->saux), ks_size(data->saux) - kMatLen, saln);
            qfrom += (qe - kMatLen);
            sfrom += (se - kMatLen);
        }
    }

    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        data->qid, query, *qbeg, *qend, ks_s(*qaln),
        data->sid, subject, *sbeg, *send, ks_s(*saln),
        ks_size(*qaln), TRUE);
}

static BOOL
s_fill_large_indel(ApproxSemeGappedAlignData* data,
    int qleft, int qright,
    int sleft, int sright)
{
    int qb, qe, sb, se;
    double perc_identity;
    int ql = qright - qleft;
    int sl = sright - sleft;
    int r = nw_ksw2_extd2(data->ksw,
                0,
                data->fwd_query,
                qleft,
                qright,
                data->query_length,
                0,
                data->fwd_subject,
                sleft,
                sright,
                data->subject_length,
                0,
                0,
                hbn_max(ql, sl),
                &qb,
                &qe,
                &sb,
                &se,
                &perc_identity,
                &data->bp_qaln,
                &data->bp_saln);
    if (!r) {
        HBN_LOG("fill large indel fail");
        return r;
    }
    //dump_align_string(ks_s(data->bp_qaln), ks_s(data->bp_saln), ks_size(data->bp_qaln), stderr);
    double eff_perc_identity = calc_effective_ident_perc(ks_s(data->bp_qaln), 
                                    ks_s(data->bp_saln), ks_size(data->bp_qaln));
    //HBN_LOG("ql = %d, sl = %d, perc_identity = %g, eff_perc_identity = %g", 
    //    ql, sl, perc_identity, eff_perc_identity);
    return TRUE;
}

BOOL 
ApproxSemeGappedAlignData_SemiGappedAlign(ApproxSemeGappedAlignData* data,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length)
{
    ApproxSemeGappedAlignData_Init(data, query, query_length, subject, subject_length);
    s_ApproxSemeGappedAlignData_RightExtend(data, TRUE);
    const int kMaxOverhang = 2000;
    if (query_length - data->f_qend <= kMaxOverhang && subject_length - data->f_send <= kMaxOverhang) {
        data->qoff = data->f_qoff;
        data->qend = data->f_qend;
        data->soff = data->f_soff;
        data->send = data->f_send;
        kputsn(ks_s(data->fwd_qaln), ks_size(data->fwd_qaln), &data->qaln);
        kputsn(ks_s(data->fwd_saln), ks_size(data->fwd_saln), &data->saln);
        data->perc_identity = calc_ident_perc(ks_s(data->qaln),
                                ks_s(data->saln),
                                ks_size(data->qaln),
                                NULL, NULL);
        data->eff_perc_identity = calc_effective_ident_perc(ks_s(data->qaln),
                                    ks_s(data->saln),
                                    ks_size(data->qaln));
        validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            data->qid, data->fwd_query, data->qoff, data->qend, ks_s(data->qaln),
            data->sid, data->fwd_subject, data->soff, data->send, ks_s(data->saln),
            ks_size(data->qaln), TRUE);
        //HBN_LOG("fwd extension successs");
        return TRUE;
    }
    s_ApproxSemeGappedAlignData_RightExtend(data, FALSE);
    if (query_length - data->r_qend <= kMaxOverhang && subject_length - data->r_send <= kMaxOverhang) {
        data->qoff = data->query_length - data->r_qend;
        data->qend = data->query_length - data->r_qoff;
        data->soff = data->subject_length - data->r_send;
        data->send = data->subject_length - data->r_soff;
        for (size_t p = ks_size(data->rev_qaln); p; --p) {
            int qc = ks_A(data->rev_qaln, p-1);
            kputc(qc, &data->qaln);
            int sc = ks_A(data->rev_saln, p-1);
            kputc(sc, &data->saln);
        }
        data->perc_identity = calc_ident_perc(ks_s(data->qaln),
                                ks_s(data->saln),
                                ks_size(data->qaln),
                                NULL, NULL);
        data->eff_perc_identity = calc_effective_ident_perc(ks_s(data->qaln),
                                    ks_s(data->saln),
                                    ks_size(data->qaln)); 
        validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            data->qid, data->fwd_query, data->qoff, data->qend, ks_s(data->qaln),
            data->sid, data->fwd_subject, data->soff, data->send, ks_s(data->saln),
            ks_size(data->qaln), TRUE);     
        //HBN_LOG("rev extension success");  
        return TRUE;
    }
    //ApproxSemeGappedAlignData_DumpAlignPos(fprintf, stderr, data);

    int qbl = 0, qbr = 0, sbl = 0, sbr = 0;
    if (!s_find_left_break_point(data, &qbl, &sbl)) return FALSE;
    if (!s_find_right_break_point(data, &qbr, &sbr)) return FALSE;
    //HBN_LOG("qbl = %d, qbr = %d, sbl = %d, sbr = %d", qbl, qbr, sbl, sbr);
    int qx = qbr - qbl;
    int sx = sbr - sbl;
    if (qx > 60000 || sx > 60000) {
	ApproxSemeGappedAlignData_DumpAlignPos(fprintf, stderr, data);
        HBN_LOG("indel region [%d, %d, %d] x [%d, %d, %d] is too long",
            qbl, qbr, qx, sbl, sbr, sx);
    }
    if (qx < 0 || qx > 60000 || sx < 0 || sx > 60000 || (abs(qx - sx) > 40000)) return FALSE;
    if (!s_fill_large_indel(data, qbl, qbr, sbl, sbr)) {
        HBN_LOG("fill large indel fail");
        return FALSE;
    }
    data->qoff = 0;
    data->qend = data->query_length;
    data->soff = 0;
    data->send = data->subject_length;
    int ai = 0;
    int qi = 0;
    int si = 0;
    while (qi < qbl) {
        int c = ks_A(data->fwd_qaln, ai);
        if (c != GAP_CHAR) ++qi;
        c = ks_A(data->fwd_saln, ai);
        if (c != GAP_CHAR) ++si;
        ++ai;
    }
    hbn_assert(qi == qbl && si == sbl);
    kputsn(ks_s(data->fwd_qaln), ai, &data->qaln);
    kputsn(ks_s(data->fwd_saln), ai, &data->saln);
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            data->qid, data->fwd_query, 0, qbl, ks_s(data->qaln),
            data->sid, data->fwd_subject, 0, sbl, ks_s(data->saln),
            ks_size(data->qaln), TRUE);   
    kputsn(ks_s(data->bp_qaln), ks_size(data->bp_qaln), &data->qaln);
    kputsn(ks_s(data->bp_saln), ks_size(data->bp_saln), &data->saln);
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            data->qid, data->fwd_query, 0, qbr, ks_s(data->qaln),
            data->sid, data->fwd_subject, 0, sbr, ks_s(data->saln),
            ks_size(data->qaln), TRUE); 
    qi = data->query_length;
    si = data->subject_length;
    //HBN_LOG("rqoff = %d, rqend = %d, qi = %d, rsoff = %d, rsend = %d, si = %d",
    //    data->r_qoff, data->r_qend, qi, data->r_soff, data->r_send, si);
    ai = 0;
    hbn_assert(ks_size(data->rev_qaln) == ks_size(data->rev_saln));
    while (qi > qbr) {
        int c = ks_A(data->rev_qaln, ai);
        if (c != GAP_CHAR) --qi;
        c = ks_A(data->rev_saln, ai);
        if (c != GAP_CHAR) --si;
        ++ai;
    }
    hbn_assert(qi == qbr && si == sbr, "qi = %d, qbr = %d, si = %d, sbr = %d, ai = %d, al = %d",
        qi, qbr, si, sbr, ai, ks_size(data->rev_qaln));
    while (ai) {
        --ai;
        int c = ks_A(data->rev_qaln, ai);
        kputc(c, &data->qaln);
        c = ks_A(data->rev_saln, ai);
        kputc(c, &data->saln);
    }
    data->perc_identity = calc_ident_perc(ks_s(data->qaln),
                                ks_s(data->saln),
                                ks_size(data->qaln),
                                NULL, NULL);
    data->eff_perc_identity = calc_effective_ident_perc(ks_s(data->qaln),
                                    ks_s(data->saln),
                                    ks_size(data->qaln));            
    return TRUE;
}

void
ApproxSemeGappedAlignData_UnitTest(const char* fasta_path)
{
    HbnFastaReader* fasta = HbnFastaReaderNew(fasta_path);
    kv_dinit(vec_u8, query);
    kv_dinit(vec_u8, subject);
    ApproxSemeGappedAlignData* data = ApproxSemeGappedAlignDataNew();
    while (!HbnLineReaderAtEof(fasta->line_reader)) {
        HbnFastaReaderReadOneSeq(fasta);
        kv_clear(query);
        for (size_t p = 0; p < ks_size(fasta->sequence); ++p) {
            int c = ks_A(fasta->sequence, p);
            c = nst_nt16_table[c];
	    kv_push(u8, query, c);
        }
        HbnFastaReaderReadOneSeq(fasta);
        kv_clear(subject);
        for (size_t p = 0; p < ks_size(fasta->sequence); ++p) {
            int c = ks_A(fasta->sequence, p);
            c = nst_nt16_table[c];
	    kv_push(u8, subject, c);
        }        

        HBN_LOG("mapping %zu .v.s. %zu", kv_size(query), kv_size(subject));
        ApproxSemeGappedAlignData_SemiGappedAlign(data, 
            kv_data(query), kv_size(query),
            kv_data(subject), kv_size(subject));
        ApproxSemeGappedAlignData_Dump(fprintf, stderr, data);
        //break;
    }
    kv_destroy(query);
    kv_destroy(subject);
    ApproxSemeGappedAlignDataFree(data);
}
