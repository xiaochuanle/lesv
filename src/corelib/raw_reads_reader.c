#include "raw_reads_reader.h"

#include "../../corelib/cstr_util.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"

RawReadReader*
RawReadReaderNew(const char* db_dir, const char* db_title, const BOOL use_batch_mode)
{
    RawReadReader* reader = (RawReadReader*)calloc(1, sizeof(RawReadReader));
    reader->db_dir = db_dir;
    reader->dbinfo = seqdb_load_volume_info(db_dir, db_title, 0);
    hbn_assert(reader->dbinfo.hdr_offset_from == 0);
    hbn_assert(reader->dbinfo.seq_offset_from == 0);
    hbn_assert(reader->dbinfo.seq_start_id == 0);
    reader->raw_read_name_array = load_seq_headers(db_dir, db_title, 0, reader->dbinfo.hdr_offset_to);
    reader->raw_read_info_array = load_seq_infos(db_dir, db_title, 0, reader->dbinfo.num_seqs);
    reader->raw_read_offset_array = NULL;
    reader->use_batch_mode = use_batch_mode;
    reader->raw_read_offset_array = (size_t*)malloc(sizeof(size_t) * reader->dbinfo.num_seqs);
    if (!use_batch_mode) {
        reader->packed_raw_seq_array = seqdb_load_pac(db_dir, db_title, 0, reader->dbinfo.seq_offset_to);
        for (int i = 0; i < reader->dbinfo.num_seqs; ++i) {
            size_t offset = reader->raw_read_info_array[i].seq_offset;
            hbn_assert((offset % 4) == 0);
            reader->raw_read_offset_array[i] = offset + 1;
        }
    } else {
        char path[HBN_MAX_PATH_LEN];
        make_packed_seq_path(db_dir, db_title, path);
        hbn_fopen(reader->packed_raw_seq_stream, path, "rb");
        memset(reader->raw_read_offset_array, 0, sizeof(size_t) * reader->dbinfo.num_seqs);
    }
    return reader;
}

RawReadReader*
RawReadReaderFree(RawReadReader* reader)
{
    if (!reader) return NULL;
    if (reader->raw_read_name_array) sfree(reader->raw_read_name_array);
    if (reader->raw_read_info_array) sfree(reader->raw_read_info_array);
    if (reader->raw_read_offset_array) sfree(reader->raw_read_offset_array);
    if (reader->packed_raw_seq_array) sfree(reader->packed_raw_seq_array);
    if (reader->packed_raw_seq_stream) hbn_fclose(reader->packed_raw_seq_stream);
    sfree(reader);
    return NULL;
}

u8*
RawReadReader_ExtractRead(RawReadReader* reader, int id, int strand, vec_u8* seqv)
{
    hbn_assert(id < reader->dbinfo.num_seqs, "id = %d, num_seqs = %d", id, reader->dbinfo.num_seqs);
    hbn_assert(strand == FWD || strand == REV);
    hbn_assert(reader->raw_read_offset_array[id]);
    size_t res_from = reader->raw_read_offset_array[id] - 1;
    size_t res_cnt = reader->raw_read_info_array[id].seq_size;
    size_t res_to = res_from + res_cnt;
    hbn_assert((res_from % 4) == 0);

    kv_resize(u8, *seqv, res_cnt);
    int pos = 0;
    if (strand == FWD) {
        for (size_t i = res_from; i < res_to; ++i) {
            u8 c = _get_pac(reader->packed_raw_seq_array, i);
            kv_A(*seqv, pos) = c;
            ++pos;
        }
    } else {
        size_t i = res_to;
        while (i > res_from) {
            --i;
            u8 c = _get_pac(reader->packed_raw_seq_array, i);
            c = 3 - c;
            kv_A(*seqv, pos) = c;
            ++pos;
        }
    }

    return kv_data(*seqv);
}

u8*
RawReadReader_ExtractSubRead(RawReadReader* reader, int id, int from, int to, int strand, vec_u8* seqv)
{
    hbn_assert(id < reader->dbinfo.num_seqs);
    hbn_assert(strand == FWD || strand == REV);
    hbn_assert(to > from);
    hbn_assert(to <= reader->raw_read_info_array[id].seq_size);
    if (strand == REV) {
        int size = reader->raw_read_info_array[id].seq_size;
        int x = size - to;
        int y = size - from;
        from = x;
        to = y;
    }
    size_t res_from, res_to, res_cnt;
    hbn_assert(reader->raw_read_offset_array[id] > 0, "id = %d", id);
    res_from = reader->raw_read_offset_array[id] - 1 + from;
    res_cnt = to - from;
    res_to = res_from + res_cnt;

    kv_resize(u8, *seqv, res_cnt);
    int pos = 0;
    if (strand == FWD) {
        for (size_t i = res_from; i < res_to; ++i) {
            u8 c = _get_pac(reader->packed_raw_seq_array, i);
            kv_A(*seqv, pos) = c;
            ++pos;
        }
    } else {
        size_t i = res_to;
        while (i > res_from) {
            --i;
            u8 c = _get_pac(reader->packed_raw_seq_array, i);
            kv_A(*seqv, pos) = 3 - c;
            ++pos;
        }
    }
    return kv_data(*seqv);
}

const char* 
RawReadReader_ReadName(RawReadReader* reader, int id)
{
    hbn_assert(id < reader->dbinfo.num_seqs);

    return
    reader->raw_read_name_array
    +
    reader->raw_read_info_array[id].hdr_offset;
}

int
RawReadReader_ReadSize(RawReadReader* reader, int id)
{
    hbn_assert(id < reader->dbinfo.num_seqs);
    return reader->raw_read_info_array[id].seq_size;
}

void
RawReadReader_LoadRawReadFromCnsHitArray(const HbnConsensusInitHit* hit_array,
    const size_t hit_count,
    RawReadReader* reader)
{
    if (!reader->use_batch_mode) return;
    if (reader->packed_raw_seq_array) sfree(reader->packed_raw_seq_array);
    hbn_assert(reader->raw_read_offset_array);
    memset(reader->raw_read_offset_array, 0, sizeof(size_t) * reader->dbinfo.num_seqs);

    for (int i = 0; i < hit_count; ++i) {
        int qid = hit_array[i].qid;
        int sid = hit_array[i].sid;
        reader->raw_read_offset_array[qid] = 1;
        reader->raw_read_offset_array[sid] = 1;
    }

    size_t num_res = 0;
    for (int i = 0; i < reader->dbinfo.num_seqs; ++i) {
        if (!reader->raw_read_offset_array[i]) continue;
        size_t s = reader->raw_read_info_array[i].seq_size;
        s = (s + 3) >> 2;
        s <<= 2;
        num_res += s;
    }
    hbn_assert((num_res % 4) == 0);

    size_t pac_bytes = num_res >> 2;
    reader->packed_raw_seq_array = (u8*)calloc(pac_bytes, 1);
    size_t pac_byte_idx = 0;
    int loaded_seqs = 0;
    size_t loaded_res = 0;
    for (int i = 0; i < reader->dbinfo.num_seqs; ++i) {
        if (!reader->raw_read_offset_array[i]) continue;
        size_t s = reader->raw_read_info_array[i].seq_offset;
        hbn_assert((s % 4) == 0);
        s >>= 2;
        size_t n = reader->raw_read_info_array[i].seq_size;
        ++loaded_seqs;
        loaded_res += n;
        n = (n + 3) >> 2;
        reader->raw_read_offset_array[i] = (pac_byte_idx<<2)|1;
        fseek(reader->packed_raw_seq_stream, s, SEEK_SET);
        u8* p = reader->packed_raw_seq_array + pac_byte_idx;
        hbn_fread(p, 1, n, reader->packed_raw_seq_stream);
        pac_byte_idx += n;
    }
    hbn_assert(pac_byte_idx == pac_bytes);

    char seq_str[64];
    char res_str[64];
    u64_to_string_comma(loaded_seqs, seq_str);
    u64_to_string_datasize(loaded_res, res_str);
    HBN_LOG("Load %s reads (%s)", seq_str, res_str);
}

void
RawReadReader_SetReadFlag(RawReadReader* reader, const int id)
{
    hbn_assert(id < reader->dbinfo.num_seqs);
    reader->raw_read_offset_array[id] = 1;
}

void
RawReadReader_LoadFlaggedReads(RawReadReader* reader)
{
    if (!reader->use_batch_mode) return;
    if (reader->packed_raw_seq_array) sfree(reader->packed_raw_seq_array);
    hbn_assert(reader->raw_read_offset_array);

    size_t num_res = 0;
    for (int i = 0; i < reader->dbinfo.num_seqs; ++i) {
        //HBN_LOG("i = %d, offset = %zu", i, reader->raw_read_offset_array[i]);
        if (!reader->raw_read_offset_array[i]) continue;
        size_t s = reader->raw_read_info_array[i].seq_size;
        s = (s + 3) >> 2;
        s <<= 2;
        num_res += s;
    }
    hbn_assert((num_res % 4) == 0);

    HBN_LOG("num_res = %zu", num_res);
    size_t pac_bytes = num_res >> 2;
    reader->packed_raw_seq_array = (u8*)calloc(pac_bytes, 1);
    size_t pac_byte_idx = 0;
    int loaded_seqs = 0;
    size_t loaded_res = 0;
    for (int i = 0; i < reader->dbinfo.num_seqs; ++i) {
        if (!reader->raw_read_offset_array[i]) continue;
        size_t s = reader->raw_read_info_array[i].seq_offset;
        hbn_assert((s % 4) == 0);
        s >>= 2;
        size_t n = reader->raw_read_info_array[i].seq_size;
        ++loaded_seqs;
        loaded_res += n;
        n = (n + 3) >> 2;
        reader->raw_read_offset_array[i] = (pac_byte_idx<<2)|1;
        fseek(reader->packed_raw_seq_stream, s, SEEK_SET);
        u8* p = reader->packed_raw_seq_array + pac_byte_idx;
        hbn_fread(p, 1, n, reader->packed_raw_seq_stream);
        pac_byte_idx += n;
    }
    hbn_assert(pac_byte_idx == pac_bytes, "pac_byte_idx = %zu, pac_bytes = %zu",
        pac_byte_idx, pac_bytes);

    char seq_str[64];
    char res_str[64];
    u64_to_string_comma(loaded_seqs, seq_str);
    u64_to_string_datasize(loaded_res, res_str);
    HBN_LOG("Load %s reads (%s)", seq_str, res_str);        
}