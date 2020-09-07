#ifndef __RAW_READS_READER_H
#define __RAW_READS_READER_H

#include "../../corelib/seqdb.h"
#include "../../corelib/gapped_candidate.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char* db_dir;
    const char* db_title;
    CSeqDBInfo dbinfo;
    const char* raw_read_name_array;
    const CSeqInfo* raw_read_info_array;
    size_t* raw_read_offset_array;
    FILE* packed_raw_seq_stream;
    u8* packed_raw_seq_array;
    BOOL use_batch_mode;
} RawReadReader;

RawReadReader*
RawReadReaderNew(const char* db_dir, const char* db_title, const BOOL use_batch_mode);

RawReadReader*
RawReadReaderFree(RawReadReader* reader);

u8*
RawReadReader_ExtractRead(RawReadReader* reader, int id, int strand, vec_u8* seqv);

u8*
RawReadReader_ExtractSubRead(RawReadReader* reader, int id, int from, int to, int strand, vec_u8* seqv);

const char* 
RawReadReader_ReadName(RawReadReader* reader, int id);

int
RawReadReader_ReadSize(RawReadReader* reader, int id);

void
RawReadReader_LoadRawReadFromCnsHitArray(const HbnConsensusInitHit* hit_array,
    const size_t hit_count,
    RawReadReader* reader);

void
RawReadReader_LoadFlaggedReads(RawReadReader* reader);

void
RawReadReader_SetReadFlag(RawReadReader* reader, const int id);

#ifdef __cplusplus
}
#endif

#endif // __RAW_READS_READER_H