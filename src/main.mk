ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhbn.a

SOURCES      := \
	./corelib/build_db.c \
	./corelib/cmd_arg.c \
	./corelib/cstr_util.c \
	./corelib/cns_read_header.c \
	./corelib/db_format.c \
	./corelib/detect_duplicate_name.cpp \
	./corelib/fasta.c \
	./corelib/gapped_candidate.c \
	./corelib/hbn_aux.c \
	./corelib/hbn_format.c \
	./corelib/hbn_hit.c \
	./corelib/hbn_package_version.c \
	./corelib/kstring.c \
	./corelib/line_reader.c \
	./corelib/m4_record.c \
	./corelib/name2id_map.c \
	./corelib/partition_mt.c \
	./corelib/seqdb_summary.c \
	./corelib/seqdb.c \
	./corelib/seq_tag.c \
	./corelib/seq_tag_report.cpp \
	./corelib/small_object_alloc.c \
	./corelib/string2hsp.c \
	./corelib/raw_reads_reader.c \
	./algo/fccns/fccns_align_tag.c \
	./algo/fccns/fccns_aux.c \
	./algo/fccns/fccns.c \
	./algo/align.c \
	./algo/approx_semi_gapped_align.c \
	./algo/chain_dp.c \
	./algo/cns_extend_chain_seed_list.c \
	./algo/dalign.c \
	./algo/diff_gapalign.cpp \
	./algo/edlib.cpp \
	./algo/edlib_wrapper.c \
	./algo/hash_list_bucket_sort.c \
	./algo/hbn_traceback.c \
	./algo/hbn_traceback_aux.c \
	./algo/init_hit_finder.c \
	./algo/kalloc.c \
	./algo/ksw2_extd2_sse.c \
	./algo/ksw2_extz2_sse.c \
	./algo/ksw2_wrapper.c \
	./algo/lookup_table.c \
	./algo/refine_align.c \
	./algo/hbn_word_finder.c \
	./algo/sort_sr_hit_seeds.cpp \
	./ncbi_blast/c_ncbi_blast_aux.c \
	./ncbi_blast/ncbi_blast_aux.cpp \
	./ncbi_blast/cmdline_args/blast_args.cpp \
	./ncbi_blast/cmdline_args/cmdline_flags.cpp \
	./ncbi_blast/cmdline_args/format_flags.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_allow.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_desc.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_types.cpp \
	./ncbi_blast/str_util/ncbistr_util.cpp \
	./ncbi_blast/str_util/ncbistr.cpp \
	./ncbi_blast/str_util/str_cmp.cpp \
	./ncbi_blast/str_util/numeric_str_interconv.cpp \
	./ncbi_blast/str_util/str_util.cpp \
	./ncbi_blast/setup/blast_encoding.c \
	./ncbi_blast/setup/blast_hits.c \
	./ncbi_blast/setup/blast_message.c \
	./ncbi_blast/setup/blast_options.c \
	./ncbi_blast/setup/blast_stat.c \
	./ncbi_blast/setup/blast_parameters.c \
	./ncbi_blast/setup/blast_program.c \
	./ncbi_blast/setup/blast_types.cpp \
	./ncbi_blast/setup/boost_erf.c \
	./ncbi_blast/setup/hsp2string.cpp \
	./ncbi_blast/setup/ncbi_math.c \
	./ncbi_blast/setup/blast_query_info.c \
	./ncbi_blast/setup/blast_sequence_blk.c \
	./ncbi_blast/setup/gapinfo.c

SRC_INCDIRS  := ./third_party/spreadsortv2

SUBMAKEFILES := \
	./app/asm_map/main.mk \
	./app/map/main.mk \
	./app/hbndb/viewhbndb.mk \
	./app/hbndb/makehbndb.mk \
	./app/hbndb/hbndb2fasta.mk \
	./app/split_seq/main.mk \
	./app/necat2sv/m4x.mk \
	./app/necat2sv/find_sv_reads.mk \
	./app/necat2sv/find_sv_signature.mk \
	./app/necat2sv/make_sv_read_groups.mk \
	./app/necat2sv/map_cns_sv_read.mk \
	./app/cns_sv_read_group/main.mk \
	./scripts/main.mk \
