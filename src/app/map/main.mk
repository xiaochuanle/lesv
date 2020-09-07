ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := qx2map
SOURCES  := \
	chain_and_extend_kmer_matches.c \
	cmdline_args.cpp \
	hbn_align_one_volume.c \
	hbn_build_seqdb.c \
	hbn_extend_subseq_hit.c \
	hbn_extend_subseq_hit_diff.c \
	hbn_find_subseq_hit.c \
	hbn_job_control.c \
	hbn_options.c \
	hbn_subseq_hit.c \
	hbn_task_struct.c \
	main.c \
	mecat_results.c \

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=