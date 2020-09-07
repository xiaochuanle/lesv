#ifndef __DETECT_DUPLICATE_NAME_H
#define __DETECT_DUPLICATE_NAME_H

#include "seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

void detect_duplicate_sequence_names(const char* db_dir, const char* db_title, const char* fasta_input_path);

#ifdef __cplusplus
}
#endif

#endif // __DETECT_DUPLICATE_NAME_H