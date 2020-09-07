ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := qx2asmmap
SOURCES  := \
	cmdline_args.cpp \
	hbn_options.c \
	map_one_part.c \
	map_one_read.c \
	main.c \
	map_aux.c \
	resolve_broken_map.c \
	../map/mecat_results.c \

SRC_INCDIRS  := . ../../third_party/spreadsortv2

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=