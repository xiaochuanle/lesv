#include "refine_align.h"

#include "hbn_traceback_aux.h"

static void
make_trace_points(const char* qaln,
    const char* saln,
    const int align_size,
    const int qoff,
    const int soff,
    vec_chain_seed* chain_seed_list)
{
    kv_clear(*chain_seed_list);
    int i = 0;
    const int E = 20;
    while (i < align_size) {
        while (i < align_size && qaln[i] != saln[i]) ++i;
        if (i >= align_size) break;
        int j = i + 1;
        while (j < align_size && qaln[j] == saln[j]) ++j;
        if (j - i >= E) {
            ChainSeed cs;
            cs.qoff = i;
            cs.length = j - i;
            kv_push(ChainSeed, *chain_seed_list, cs);
        } else if (i == 0) {
            ChainSeed cs;
            cs.qoff = 0;
            cs.length = j - i;
            kv_push(ChainSeed, *chain_seed_list, cs);
        } else if (j == align_size) {
            ChainSeed cs;
            cs.qoff = i;
            cs.length = j - i;
            kv_push(ChainSeed, *chain_seed_list, cs);
        }
        i = j;
    }
    
    i = 0;
    int qidx = 0;
    int sidx = 0;
    for (size_t k = 0; k < kv_size(*chain_seed_list); ++k) {
        ChainSeed* cs = &kv_A(*chain_seed_list, k);
        while (i < cs->qoff) {
            hbn_assert(qaln[i] != GAP_CHAR || saln[i] != GAP_CHAR);
            if (qaln[i] != GAP_CHAR) ++qidx;
            if (saln[i] != GAP_CHAR) ++sidx;
            ++i;
        }
        cs->qoff = qoff + qidx;
        cs->soff = soff + sidx;
        //HBN_LOG("%zu, [%d, %d, %d]", k, cs->qoff, cs->soff, cs->length);
    }
}

static void
apped_match_subseq(const u8* query,
    const int qfrom,
    const int qto,
    const u8* subject,
    const int sfrom,
    const int sto,
    kstring_t* qaln,
    kstring_t* saln)
{
    //HBN_LOG("qf = %d, qt = %d, sf = %d, st = %d", qfrom, qto, sfrom, sto);
    hbn_assert(qto - qfrom == sto - sfrom);
    int n = qto - qfrom;
    for (int i = 0; i < n; ++i) {
        int c1 = query[qfrom + i];
        c1 = DECODE_RESIDUE(c1);
        kputc(c1, qaln);
        int c2 = subject[sfrom + i];
        c2 = DECODE_RESIDUE(c2);
        kputc(c2, saln);
        hbn_assert(c1 == c2, "qfrom = %d, qto = %d, sfrom = %d, sto = %d, qi = %d, si = %d, c1 = %c, c2 = %c", 
            qfrom, qto, sfrom, sto,
            qfrom + i, sfrom + i, c1, c2);
    }
}

static void
dump_sub_align(const char* qaln,
    const char* saln,
    int qstart,
    int sstart,
    int qfrom,
    int qto,
    int sfrom,
    int sto)
{
    int qi = qstart;
    int i = 0;
    while (qi < qfrom) {
        if (qaln[i] != GAP_CHAR) ++qi;
        ++i;
    }
    while (qi < qto) {
        fprintf(stderr, "%c", qaln[i]);
        if (qaln[i] != GAP_CHAR) ++qi;
        ++i;
    }
    fprintf(stderr, "\n");

    int si = sstart;
    i = 0;
    while (si < sfrom) {
        if (saln[i] != GAP_CHAR) ++si;
        ++i;
    }
    while (si < sto) {
        fprintf(stderr, "%c", saln[i]);
        if (saln[i] != GAP_CHAR) ++si;
        ++i;
    }
    fprintf(stderr, "\n");
}

int fill_chain_gaps(Ksw2Data* data,
        const int qid,
        const u8* query,
        const int qfrom,
        const int qto,
        const int qsize,
        const int sid,
        const u8* subject,
        const int sfrom,
        const int sto,
        const int ssize,
        const int min_align_size,
        const double min_ident_perc,
        int* qoff,
        int* qend,
        int* soff,
        int* send,
        double* ident_perc,
        kstring_t* qaln,
        kstring_t* saln)
{
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int flag = 0;
    int max_distance = hbn_max(qto - qfrom, sto - sfrom);
if (max_distance > 40000) return 0;
    max_distance = hbn_min(max_distance, 40000);
    int max_band_width = -1;
    int zdrop = -1;
    const u8* qsubseq = query + qfrom;
    const int qsubseq_size = qto - qfrom;
    const u8* ssubseq = subject + sfrom;
    const int ssubseq_size = sto - sfrom;
    ksw_extd2_sse(data->km, qsubseq_size, qsubseq, ssubseq_size, ssubseq, 5, data->mat,
        data->go, data->ge, data->go1, data->ge1, max_band_width, zdrop, data->end_bonus, flag, &ez);
    if (ez.n_cigar == 0) return 0;

    int qi = 0;
    int si = 0;
    ks_clear(*qaln);
    ks_clear(*saln);
    for (int k = 0; k < ez.n_cigar; ++k) {
        int op_num = ez.cigar[k]>>4;
        int op_type = "MIDN"[ez.cigar[k]&0xf];
        int c = 0;
        switch (op_type)
        {
        case 'M':
            for (int t = 0; t < op_num; ++t, ++qi, ++si) {
                c = qsubseq[qi];
                c = DECODE_RESIDUE(c);
                kputc(c, qaln);
                c = ssubseq[si];
                c = DECODE_RESIDUE(c);
                kputc(c, saln);
            }
            break;
        case 'I':
            for (int t = 0; t < op_num; ++t, ++qi) {
                c = qsubseq[qi];
                c = DECODE_RESIDUE(c);
                kputc(c, qaln);
                kputc(GAP_CHAR, saln);
            }
            break;
        case 'D':
            for (int t = 0; t < op_num; ++t, ++si) {
                kputc(GAP_CHAR, qaln);
                c = ssubseq[si];
                c = DECODE_RESIDUE(c);
                kputc(c, saln);
            }
            break;            
        default:
            HBN_LOG("invalid op_type: %d", op_type);
            break;
        }
    }
    if (qi < qsubseq_size || si < ssubseq_size) {
        HBN_LOG("[%d, %d, %d, %d] x [%d, %d, %d, %d] %g%% sw fail, dist = %d", qid, *qoff, *qend, qsize, sid, *soff, *send, ssize, *ident_perc, max_distance);
        return 0;
    }
    hbn_assert(qi == qsubseq_size);
    hbn_assert(si == ssubseq_size);
    kfree(data->km, ez.cigar);

    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, qsubseq, 0, qsubseq_size, ks_s(*qaln),
        0, ssubseq, 0, ssubseq_size, ks_s(*saln),
        ks_size(*qaln), TRUE);
    *ident_perc = calc_ident_perc(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), NULL, NULL);
    if (*ident_perc < min_ident_perc) return 0;
    *qoff = qfrom;
    *qend = qto;
    *soff = sfrom;
    *send = sto;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, *qoff, *qend, ks_s(*qaln),
        0, subject, *soff, *send, ks_s(*saln),
        ks_size(*qaln), TRUE);

    return 1;    
}

int
refine_align(Ksw2Data* ksw,
    const u8* query,
    const int query_align_from,
    const int query_align_to,
    const int query_length,
    const u8* subject,
    const int subject_align_from,
    const int subject_align_to,
    const int subject_length,
    const char* org_qaln,
    const char* org_saln,
    const int org_align_size,
    kstring_t* qaln,
    kstring_t* saln)
{
    ks_clear(*qaln);
    ks_clear(*saln);
    kv_dinit(vec_chain_seed, chain_seed_list);
    make_trace_points(org_qaln, org_saln, org_align_size, query_align_from, subject_align_from, &chain_seed_list);
    ChainSeed* seed_array = kv_data(chain_seed_list);
    int seed_count = kv_size(chain_seed_list);
    const int E = 10;
    const int E2 = E * 2;
    int qfrom = 0, qto = 0;
    int sfrom = 0, sto = 0;
    int r = 1;
    if (seed_array[0].length >= E2) {
        qfrom = seed_array[0].qoff;
        qto = seed_array[0].qoff + seed_array[0].length;
        sfrom = seed_array[0].soff;
        sto = seed_array[0].soff + seed_array[0].length;
        apped_match_subseq(query, qfrom, qto - E, subject, sfrom, sto - E, qaln, saln);
        qfrom = qto - E;
        sfrom = sto - E;
    } else {
        qfrom = seed_array[0].qoff;
        sfrom = seed_array[0].soff;
    }

    int qoff, qend, soff, send;
    int distance;
    double ident_perc;
    ks_dinit(qaux);
    ks_dinit(saux);
    for (int i = 0; i < seed_count - 1; ++i) {
        ChainSeed sj = seed_array[i+1];
        //HBN_LOG("%d: (%d, %d, %d)", i, sj.qoff, sj.soff, sj.length);
        if (sj.length >= E2) {
            qto = sj.qoff + E;
            sto = sj.soff + E;
            //run_nw(query + qfrom, qto - qfrom, subject + sfrom, sto - sfrom, data, &data->ext_qabuf, &data->ext_sabuf, qae, sae);
            distance = hbn_max(qto - qfrom, sto - sfrom);
            r = fill_chain_gaps(ksw,
                0,
                query,
                qfrom,
                qto,
                query_length,
                0,
                subject,
                sfrom,
                sto,
                subject_length,
                1,
                0.0,
                &qoff,
                &qend,
                &soff,
                &send,
                &ident_perc,
                &qaux,
                &saux);
	    if (!r) return r;
#if 0
HBN_LOG("old align [%d, %d] x [%d, %d]", qfrom, qto, sfrom, sto);
dump_sub_align(org_qaln, org_saln, query_align_from, subject_align_from, qfrom, qto, sfrom, sto);
HBN_LOG("new align");
hbn_fwrite(ks_s(qaux), 1, ks_size(qaux), stderr);
fprintf(stderr, "\n");
hbn_fwrite(ks_s(saux), 1, ks_size(saux), stderr);
fprintf(stderr, "\n");
#endif
            kputsn(ks_s(qaux), ks_size(qaux), qaln);
            kputsn(ks_s(saux), ks_size(saux), saln);
            hbn_assert(r);
            qfrom = qto;
            qto = sj.qoff + sj.length - E;
            sfrom = sto;
            sto = sj.soff + sj.length - E;
            apped_match_subseq(query, qfrom, qto, subject, sfrom, sto, qaln, saln);
            qfrom = qto;
            sfrom = sto;
        } else {
            qto = sj.qoff + sj.length / 2;
            sto = sj.soff + sj.length / 2;
            //run_nw(query + qfrom, qto - qfrom, subject + sfrom, sto - sfrom, data, &data->ext_qabuf, &data->ext_sabuf, qae, sae);
            distance = hbn_max(qto - qfrom, sto - sfrom);
            r = fill_chain_gaps(ksw,
                0,
                query,
                qfrom,
                qto,
                query_length,
                0,
                subject,
                sfrom,
                sto,
                subject_length,
                1,
                0.0,
                &qoff,
                &qend,
                &soff,
                &send,
                &ident_perc,
                &qaux,
                &saux);
	    if (!r) return r;
#if 0
HBN_LOG("old align [%d, %d] x [%d, %d]", qfrom, qto, sfrom, sto);
dump_sub_align(org_qaln, org_saln, query_align_from, subject_align_from, qfrom, qto, sfrom, sto);
HBN_LOG("new align");
hbn_fwrite(ks_s(qaux), 1, ks_size(qaux), stderr);
fprintf(stderr, "\n");
hbn_fwrite(ks_s(saux), 1, ks_size(saux), stderr);
fprintf(stderr, "\n");
#endif
            kputsn(ks_s(qaux), ks_size(qaux), qaln);
            kputsn(ks_s(saux), ks_size(saux), saln);
            hbn_assert(r);
            qfrom = qto;
            sfrom = sto;
        }
    }
    ks_destroy(qaux);
    ks_destroy(saux);

    ChainSeed se = seed_array[seed_count-1];
    qto = se.qoff + se.length;
    sto = se.soff + se.length;
    apped_match_subseq(query, qfrom, qto, subject, sfrom, sto, qaln, saln);

#if 0
HBN_LOG("final, align_size = %zu", ks_size(*qaln));
hbn_fwrite(ks_s(*qaln), 1, 30, stderr);
fprintf(stderr, "\n");
hbn_fwrite(ks_s(*saln), 1, 30, stderr);
fprintf(stderr, "\n");
#endif

    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, query_align_from, query_align_to, ks_s(*qaln),
        0, subject, subject_align_from, subject_align_to, ks_s(*saln), 
        ks_size(*qaln), TRUE);

    kv_destroy(chain_seed_list);
    return 1;
}
