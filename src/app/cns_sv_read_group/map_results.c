#include "map_results.h"

const char* kSamVersion = "1.6";

void
make_sam_rg_info(const int subject_id, kstring_t* rg)
{
    char sid_str[64];
    u64_to_fixed_width_string_r(subject_id, sid_str, HBN_DIGIT_WIDTH);
    ksprintf(rg, "subject%s", sid_str);
}

void
print_sam_prolog(FILE* out, 
    const char* sam_version, 
    const char* prog_version, 
    const char* rg_info,
    const char* rg_sample,
    const CSeqDB* db,
    const int subject_id,
    int argc, 
    char* argv[])
{
    fprintf(out, "@HD\t");
    fprintf(out, "VN:%s\t", sam_version);
    fprintf(out, "SO:unknown\t");
    fprintf(out, "GO:query\n");

    fprintf(out, "@SQ\tSN:%s\tLN:%zu\n", seqdb_seq_name(db, subject_id), seqdb_seq_size(db, subject_id));

    fprintf(out, "@PG\t");
    fprintf(out, "ID:%s\t", argv[0]);
    fprintf(out, "VN:%s\t", prog_version);
    fprintf(out, "CL:");
    for (int i = 0; i < argc - 1; ++i) fprintf(out, "%s ", argv[i]);
    fprintf(out, "%s", argv[argc-1]);
    fprintf(out, "\t");
    fprintf(out, "PN:%s\n", argv[0]);
    if (rg_info) {
        fprintf(out, "@RG\tID:%s", rg_info);
        if (rg_sample) fprintf(out, "\tSM:%s", rg_sample);
        fprintf(out, "\n");
    } else {
        if (rg_sample) fprintf(out, "@RG\tSM:%s\n", rg_sample);
    }
}

static void
print_cigar_core(const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    const char* q = qaln;
    const char* s = saln;
    int i = 0;
    while (i < aln_size) {
        char type = 'N';
        int cnt = 0;
        if (q[i] == GAP_CHAR) {  // delete from subject (gap in query)
            type = 'D';
            while (i < aln_size && q[i] == GAP_CHAR) {
                hbn_assert(s[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else if (s[i] == GAP_CHAR) { // insert into subject (gap in subject)
            type = 'I';
            while (i < aln_size && s[i] == GAP_CHAR) {
                hbn_assert(q[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else { // substitution
            type = 'M';
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (i < aln_size && q[i] != GAP_CHAR && s[i] != GAP_CHAR) {
                ++i;
                ++cnt;
            }
        }
        ksprintf(out, "%d%c", cnt, type);
    }
}

void
print_sam_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    if (qoff) ksprintf(out, "%dS", qoff);
    print_cigar_core(qaln, saln, aln_size, out);
    if (qend < qsize) ksprintf(out, "%dS", qsize - qend);
}

void
print_paf_cigar(const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    ksprintf(out, "cg:Z:");
    print_cigar_core(qaln, saln, aln_size, out);
}

void
print_md(const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    ksprintf(out, "MD:Z:");
    const char* q = qaln;
    const char* s = saln;
    int i = 0, l_md = 0;
    while (i < aln_size) {
        int j = i;
        if (s[i] == GAP_CHAR) { // insert into subject (gap in  query)
            while (j < aln_size && s[j] == GAP_CHAR) {
                hbn_assert(q[j] != GAP_CHAR);
                ++j;
            }
        } else if (q[i] == GAP_CHAR) { // delete from subject (gap in subject)
            while (j < aln_size && q[j] == GAP_CHAR) {
                hbn_assert(s[j] != GAP_CHAR);
                ++j;
            }
            ksprintf(out, "%d^", l_md);
            kputsn(&s[i], j - i, out);
            l_md = 0;
        } else { // substitution
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (j < aln_size && q[j] != GAP_CHAR && s[j] != GAP_CHAR) ++j;
            for (int k = i; k < j; ++k) {
                if (q[k] != s[k]) {
                    ksprintf(out, "%d%c", l_md, s[k]);
                    l_md = 0;
                } else {
                    ++l_md;
                }
            }
        }
        i = j;
    }   
    if (l_md > 0) ksprintf(out, "%d", l_md); 
}

#if 0
                                         1 | string | Query sequence name                                     |
                                      |  2 |  int   | Query sequence length                                   |
                                      |  3 |  int   | Query start coordinate (0-based)                        |
                                      |  4 |  int   | Query end coordinate (0-based)                          |
                                      |  5 |  char  | `+' if query/target on the same strand; `-' if opposite |
                                      |  6 | string | Target sequence name                                    |
                                      |  7 |  int   | Target sequence length                                  |
                                      |  8 |  int   | Target start coordinate on the original strand          |
                                      |  9 |  int   | Target end coordinate on the original strand            |
                                      | 10 |  int   | Number of matching bases in the mapping                 |
                                      | 11 |  int   | Number bases, including gaps, in the mapping            |
                                      | 12 |  int   | Mapping quality (0-255 with 255 for missing)  
#endif 

void
print_one_paf_result(
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_score,
    const BOOL dump_cigar,
    const BOOL dump_md,
    kstring_t* out)
{
    const char tab = '\t';
    ksprintf(out, "%s", qname); /// 1) query name
    kputc(tab, out);
    ksprintf(out, "%d", qsize); /// 2) query length
    kputc(tab, out);
    ksprintf(out, "%d", qoff); /// 3) query start coordinate (0-based)
    kputc(tab, out);
    ksprintf(out, "%d", qend); /// 4) query end coordinate (0-based)
    kputc(tab, out);
    /// 5) '+' if query and subject on the same strand; '-' if opposite
    /// the subject sequence is always on forward strand
    if (qdir == FWD) {
        kputc('+', out);
    } else {
        kputc('-', out);
    }
    kputc(tab, out);
    ksprintf(out, "%s", sname); /// 6) subject name
    kputc(tab, out);
    ksprintf(out, "%d", ssize); /// 7) subject length
    kputc(tab,out);
    ksprintf(out, "%d", soff); /// 8) subject start coordinate (0-based)
    kputc(tab, out);
    ksprintf(out, "%d", send); /// 9) subject end coordinate (0-based)
    kputc(tab, out);
    
    int num_ident = 0;
    const char* q = qaln;
    const char* s = saln;
    for (int i = 0; i < aln_size; ++i) if (q[i] == s[i]) ++num_ident;
    ksprintf(out, "%d", num_ident); /// 10) number of matching bases in the alignment
    kputc(tab, out);
    ksprintf(out, "%d", aln_size); /// 11) number of bases in the alignment (including gaps and mismatch bases)
    kputc(tab, out);
    ksprintf(out, "%d", 60); /// 12) mapq
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", map_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", aln_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", map_score); ///  dp score
    if (dump_cigar) {
        kputc(tab, out);
        print_paf_cigar(qaln, saln, aln_size, out);
    }
    if (dump_md) {
        kputc(tab, out);
        print_md(qaln, saln, aln_size, out);
    }
    kputc('\n', out);
}

void
print_one_sam_result(
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_score,
    const double perc_identity,
    const double eff_perc_identity,
    const BOOL dump_md,
    const char* rg_sample,
    kstring_t* out)
{
    const char tab = '\t';
    int flag = 0;
    if (qdir ==  REV) flag |= 0x10; // reverse query strand
    ksprintf(out, "%s", qname); /// 1) query name
    kputc(tab, out);
    ksprintf(out, "%d", flag); /// 2) flag
    kputc(tab, out);
    ksprintf(out, "%s", sname); /// 3) subject name
    kputc(tab, out);
    ksprintf(out, "%d", soff + 1); /// 4) left most subject position (1-based)
    kputc(tab, out);
    ksprintf(out, "%d", 60); /// 5) mapq
    kputc(tab, out);
    print_sam_cigar(qoff, qend, qsize, qaln, saln, aln_size, out); /// 6) cigar
    kputc(tab, out);
    ksprintf(out, "*"); /// 7) rnext
    kputc(tab, out);
    ksprintf(out, "*"); /// 8) pnext
    kputc(tab, out);
    ksprintf(out, "%d", 0); /// 9) subject length
    kputc(tab, out);
    /// 10) aligned subsequence
    int num_ident = 0;
    const char* q = qaln;
    const char* s = saln;
    for (int i = 0; i < aln_size; ++i) {
        if (q[i] == s[i]) ++num_ident;
        if (q[i] != GAP_CHAR) kputc(q[i], out);
    }
    kputc(tab, out);
    ksprintf(out, "*"); /// 11) quality score
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", map_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", aln_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", map_score); ///  dp score
    if (dump_md) {
        kputc(tab, out);
        print_md(qaln, saln, aln_size, out);
    }
    if (rg_sample) {
        kputc(tab, out);
        ksprintf(out, "RG:Z:%s", rg_sample);
    }
    /// qs
    kputc(tab, out);
    ksprintf(out, "qs:i:%d", qoff);
    /// qe
    kputc(tab, out);
    ksprintf(out, "qe:i:%d", qend);
    // ql
    kputc(tab, out);
    ksprintf(out, "ql:i:%d", qsize);
    /// identity
    kputc(tab, out);
    ksprintf(out, "mc:f:%g", perc_identity);
    /// effective identity
    kputc(tab, out);
    ksprintf(out, "ec:f:%g", eff_perc_identity);
    kputc('\n', out);
}

void
print_one_m4_result(
    const char* qname,
    const int qid,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int sid,
    const int soff,
    const int send,
    const int ssize,
    const int map_score,
    const double perc_identity,
    const BOOL binary,
    kstring_t* line,
    kstring_t* out)
{
    M4Record m4;
    m4.qid = qid;
    m4.qdir = qdir;
    m4.qoff = qoff;
    m4.qend = qend;
    m4.qsize = qsize;
    m4.sid = sid;
    m4.sdir = FWD;
    m4.soff = soff;
    m4.send = send;
    m4.ssize = ssize;
    m4.ident_perc = perc_identity;
    m4.score = map_score;
    
    if (binary) {
        const char* input = (const char*)(&m4);
        const int input_len = sizeof(M4Record);
        kputsn(input, input_len, out);
    } else {
        ks_clear(*line);
        DUMP_M4_RECORD(ksprintf, line, m4);
        kputsn(ks_s(*line), ks_size(*line), out);
    }
}