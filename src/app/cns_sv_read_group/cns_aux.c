#include "cns_aux.h"


CnsThreadData*
CnsThreadDataNew(HbnProgramOptions* opts, 
    CSeqDB* db,
    const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* fasta_out,
    FILE* sam_out,
    pthread_mutex_t* out_lock,
    int* svr_id,
    pthread_mutex_t* svr_id_lock)
{
    CnsThreadData* data = (CnsThreadData*)calloc(1, sizeof(CnsThreadData));
    data->opts = opts;
    data->db = db;
    data->raw_read_info_array = raw_read_info_array;
    data->raw_read_name_array = raw_read_name_array;
    data->hit_finder = InitHitFindDataNew(opts->memsc_kmer_size,
                            opts->memsc_kmer_window,
                            opts->memsc_score);
    data->tbck_data = HbnTracebackDataNew();
    data->ksw_data = Ksw2DataNew();
    ksw2_extd2_set_params(data->ksw_data);
    data->cns_data = FCCnsDataNew();
    data->fasta_out = fasta_out;
    data->sam_out = sam_out;
    data->out_lock = out_lock;
    data->next_svr_id = svr_id;
    data->svr_id_lock = svr_id_lock;
    ks_init(data->qaln);
    ks_init(data->saln);
    kv_init(data->chain_seed_list);
    return data;
}

CnsThreadData*
CnsThreadDataFree(CnsThreadData* data)
{
    if (!data) return NULL;
    data->hit_finder = InitHitFindDataFree(data->hit_finder);
    data->tbck_data = HbnTracebackDataFree(data->tbck_data);
    data->cns_data = FCCnsDataFree(data->cns_data);
    Ksw2DataFree(data->ksw_data);
    ks_destroy(data->qaln);
    ks_destroy(data->saln);
    kv_destroy(data->chain_seed_list);
    sfree(data);
    return NULL;
}