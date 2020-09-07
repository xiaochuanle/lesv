ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := qx2csvrg
SOURCES  := \
				 cmdline_args.cpp \
				 cns_aux.c \
				 cns_one_group.c \
				 main.c \
				 map_results.c \
				 sv_read_group.c \
				 ../necat2sv/sv_signature.cpp \
				 ../necat2sv/sv_read_group_file_name.cpp \

SRC_INCDIRS  := . ../../third_party/spreadsortv2

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=