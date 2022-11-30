CC       = cc
CXX      = g++
CPPFLAGS += -I slow5lib/include/
CFLAGS   += -g -Wall -O2
LANGFLAG  = -x c++ -std=c++11
LDFLAGS  += $(LIBS) -lz -lm -lpthread
BUILD_DIR = build

ifeq ($(zstd),1)
LDFLAGS		+= -lzstd
endif

# change the tool name to what you want
BINARY = poregen

OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/poregen.o \
      $(BUILD_DIR)/subtool0.o \
      $(BUILD_DIR)/gmove.o \
      $(BUILD_DIR)/thread.o \

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test

$(BINARY): $(OBJ) slow5lib/lib/libslow5.a
	$(CXX) $(CFLAGS) $(OBJ) slow5lib/lib/libslow5.a $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/poregen.o: src/poregen.c src/misc.h src/error.h src/poregen.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/subtool0.o: src/subtool0.c src/error.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/gmove.o: src/gmove.cpp src/error.h src/ksort.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

# follow the main.o above and add more objects here if needed

all:$(BINARY)

slow5lib/lib/libslow5.a:
	$(MAKE) -C slow5lib zstd=$(zstd) no_simd=$(no_simd) zstd_local=$(zstd_local) lib/libslow5.a

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o
	make -C slow5lib clean

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/* autom4te.cache

test: $(BINARY)
	./test/test.sh

