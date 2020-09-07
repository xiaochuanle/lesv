ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := qx2msvrg
SOURCES  := align_subseqs.c \
            find_one_sv_group.cpp \
            sv_reads.cpp \
            sv_signature.cpp \
            make_sv_read_groups.c \
            sv_read_file_name.cpp \
            sv_signature_file_name.cpp \
            sv_read_group_file_name.cpp \

SRC_INCDIRS  := . ../../third_party/spreadsortv2

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=