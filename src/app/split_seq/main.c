#include "../../corelib/line_reader.h"
#include "../../corelib/khash.h"
#include "../../corelib/fasta.h"
#include "../../corelib/seq_tag_report.h"
#include "../../corelib/string2hsp.h"
#include "../../corelib/ksort.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/db_format.h"
#include "../../ncbi_blast/setup/hsp2string.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>

static int sSegLen = 0;
static int sOvlpLen = 0;
static int sMinLastSegLen = 0;

static void
print_usage(const char* pn)
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s seg_len ovlp_len min_last_seg_len input\n", pn);
}

static void
split_one_subject(const char* subject_name, const char* subject, int subject_length)
{
    int from = 0;
    while (from < subject_length) {
        int to = from + sSegLen;
        if (to > subject_length) to = subject_length;
        int left = subject_length - to;
        if (left < sMinLastSegLen) to = subject_length;
        const char* s = subject + from;
        int l = to - from;
        if (from == 0 && to == subject_length) fprintf(stdout, ">%s\n", subject_name);
        else fprintf(stdout, ">%s_%d_%d\n", subject_name, from, to);
        hbn_fwrite(s, 1, l, stdout);
        fprintf(stdout, "\n");
        from = (to < subject_length) ? (to - sOvlpLen) : subject_length;
    }
}

static void
split_one_file(const char* file_name)
{
    HBN_LOG("spliting sequences in %s", file_name);
    int seqs = 0;
    size_t res = 0;
    HbnFastaReader* fasta = HbnFastaReaderNew(file_name);
    while (!HbnLineReaderAtEof(fasta->line_reader)) {
        HbnFastaReaderReadOneSeq(fasta);
        ++seqs;
        res += ks_size(fasta->sequence);
        kputc('\0', &fasta->name);
        split_one_subject(ks_s(fasta->name), ks_s(fasta->sequence), ks_size(fasta->sequence));
    }
    HbnFastaReaderFree(fasta);
    char seq_str[64], res_str[64];
    u64_to_string_comma(seqs, seq_str);
    u64_to_string_datasize(res, res_str);
    HBN_LOG("%s sequences (%s)", seq_str, res_str);
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        print_usage(argv[0]);
        return 1;
    }
    sSegLen = atoi(argv[1]);
    sOvlpLen = atoi(argv[2]);
    sMinLastSegLen = atoi(argv[3]);
    const char* input = argv[4];
    hbn_assert(sOvlpLen < sSegLen);

    EDbFormat fmt = hbn_guess_db_format(input);
    if (fmt == eDbFormatEmptyFile) {
        HBN_LOG("file %s is empty", input);
        return 0;
    }

    if (fmt == eDbFormatUnknown) {
        HbnLineReader* line_reader = HbnLineReaderNew(input);
        while (!HbnLineReaderAtEof(line_reader)) {
            HbnLineReaderReadOneLine(line_reader);
            kstring_t* line = &line_reader->line;
            if (ks_empty(*line)) continue;
            kputc('\0', line);
            if (truncate_both_end_spaces(ks_s(*line)) == 0) continue;
            split_one_file(ks_s(*line));
        }
        HbnLineReaderFree(line_reader);
    } else {
        split_one_file(input);
    }   

    return 0;
}