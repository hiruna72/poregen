CC       = cc
CXX      = g++
#CFLAGS   += -g -Wall -O2
CFLAGS   += -g -Wall
LANGFLAG  = -x c++ -std=c++11
LDFLAGS  += $(LIBS) -lz -lm -lpthread
BUILD_DIR = build
CPPFLAGS += -I slow5lib/include/ -I $(BUILD_DIR)/htslib

HTS_VERSION = 1.9

ifeq ($(zstd),1)
LDFLAGS		+= -lzstd
endif

# change the tool name to what you want
BINARY = poregen

OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/poregen.o \
      $(BUILD_DIR)/subtool0.o \
      $(BUILD_DIR)/gmove.o \
      $(BUILD_DIR)/kmer_freq.o \
      $(BUILD_DIR)/reform.o \
      $(BUILD_DIR)/thread.o \

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test

$(BINARY): $(OBJ) slow5lib/lib/libslow5.a $(BUILD_DIR)/lib/libhts.a
	$(CXX) $(CFLAGS) $(OBJ) slow5lib/lib/libslow5.a $(BUILD_DIR)/lib/libhts.a $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/poregen.o: src/poregen.cpp src/misc.h src/error.h src/poregen.h
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/subtool0.o: src/subtool0.c src/error.h
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/gmove.o: src/gmove.cpp src/error.h src/ksort.h
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/kmer_freq.o: src/kmer_freq.cpp src/error.h
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/reform.o: src/reform.cpp src/error.h $(BUILD_DIR)/lib/libhts.a
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c
	$(CXX) $(LANGFLAG) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

# follow the main.o above and add more objects here if needed

all:$(BINARY)

slow5lib/lib/libslow5.a:
	$(MAKE) -C slow5lib zstd=$(zstd) no_simd=$(no_simd) zstd_local=$(zstd_local) lib/libslow5.a

$(BUILD_DIR)/lib/libhts.a:
	@if command -v curl; then \
		curl -o $(BUILD_DIR)/htslib.tar.bz2 -L https://github.com/samtools/htslib/releases/download/$(HTS_VERSION)/htslib-$(HTS_VERSION).tar.bz2; \
	else \
		wget -O $(BUILD_DIR)/htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/$(HTS_VERSION)/htslib-$(HTS_VERSION).tar.bz2; \
	fi
	tar -xf $(BUILD_DIR)/htslib.tar.bz2 -C $(BUILD_DIR)
	mv $(BUILD_DIR)/htslib-$(HTS_VERSION) $(BUILD_DIR)/htslib
	rm -f $(BUILD_DIR)/htslib.tar.bz2
	cd $(BUILD_DIR)/htslib && \
	./configure --prefix=`pwd`/../ --enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no && \
	make -j8 && \
	make install

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o
	make -C slow5lib clean

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/* autom4te.cache

test: $(BINARY)
	./test/test.sh

