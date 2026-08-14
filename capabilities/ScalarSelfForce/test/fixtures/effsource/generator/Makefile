# Makefile — one-command regeneration of the EffectiveSource-C oracle fixtures
# (ESSF node effsource-oracle-fixtures, W2-D).
#
#   make fixtures                 # rebuild both capture binaries + rewrite ../<*.jl>
#   make fixtures OUTDIR=/tmp/x   # write elsewhere (used by the reproducibility test)
#
# The ref-code/EffectiveSource mirror is READ-ONLY: sources are COPIED into BUILD_DIR
# (never compiled in place; override BUILD_DIR to a scratch dir to keep the tree clean).
# Requires: cc/gcc, make, GSL via pkg-config (GSL 2.7.1 + gcc 13.3.0 at capture time).

# Regeneration-only input: the READ-ONLY upstream EffectiveSource C mirror. The default
# records the capture site; override EFFSOURCE_SRC (or GRTOOLS_EFFSOURCE_ROOT via the
# test harness) at a site with its own mirror. The committed fixtures never need this.
EFFSOURCE_SRC ?= /data/haiyangw/claude/emri-disk/ref-code/EffectiveSource
BUILD_DIR     ?= build
OUTDIR        ?= $(abspath ..)

CC         ?= gcc
CFLAGS     ?= -O2 -std=c99
GSL_CFLAGS := $(shell pkg-config --cflags gsl)
GSL_LIBS   := $(shell pkg-config --libs gsl)

all: $(BUILD_DIR)/capture_elliptic $(BUILD_DIR)/capture_comoving

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# scratch copies of the mirror sources (mirror itself is never touched)
$(BUILD_DIR)/kerr-circular.c: $(EFFSOURCE_SRC)/kerr-circular.c | $(BUILD_DIR)
	cp $< $@
$(BUILD_DIR)/kerr-circular-comoving.c: $(EFFSOURCE_SRC)/kerr-circular-comoving.c | $(BUILD_DIR)
	cp $< $@
$(BUILD_DIR)/effsource.h: $(EFFSOURCE_SRC)/effsource.h | $(BUILD_DIR)
	cp $< $@

$(BUILD_DIR)/capture_elliptic: capture.c $(BUILD_DIR)/kerr-circular.c $(BUILD_DIR)/effsource.h
	$(CC) $(CFLAGS) $(GSL_CFLAGS) -I$(BUILD_DIR) -o $@ capture.c \
	    $(BUILD_DIR)/kerr-circular.c $(GSL_LIBS) -lm

$(BUILD_DIR)/capture_comoving: capture.c $(BUILD_DIR)/kerr-circular-comoving.c $(BUILD_DIR)/effsource.h
	$(CC) $(CFLAGS) $(GSL_CFLAGS) -DEFFSOURCE_VARIANT_COMOVING -I$(BUILD_DIR) -o $@ capture.c \
	    $(BUILD_DIR)/kerr-circular-comoving.c $(GSL_LIBS) -lm

fixtures: all
	mkdir -p $(OUTDIR)
	$(BUILD_DIR)/capture_elliptic $(OUTDIR)
	$(BUILD_DIR)/capture_comoving $(OUTDIR)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all fixtures clean
