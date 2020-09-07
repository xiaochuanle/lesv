ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := qx2svr
SOURCES  := align_subseqs.c \
            find_sv_reads.c \
            sv_reads.cpp \
            trf_array.cpp \
            sv_read_file_name.cpp \
            sv_signature_file_name.cpp \
            sv_read_group_file_name.cpp \

SRC_INCDIRS  := . ../../third_party/spreadsortv2

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=