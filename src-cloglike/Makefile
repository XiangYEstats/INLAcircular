CC = cc
CFLAGS = -O2 -std=c11 -Wall -Wextra -Wpedantic
CPPFLAGS =
BUILD_DIR = build

CLOGLIKE_SOURCE = INLAcirc_cloglike.c
CLOGLIKE_OBJECT = $(BUILD_DIR)/INLAcirc_cloglike.o
TEST_BINARY = $(BUILD_DIR)/INLAcirc_cloglike_test

.PHONY: all check clean

all: $(CLOGLIKE_OBJECT)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(CLOGLIKE_OBJECT): $(CLOGLIKE_SOURCE) INLAcirc_cloglike.h \
	INLAcirc_common.h INLAcirc_cgeneric_compat.h $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -fPIC -c $(CLOGLIKE_SOURCE) -o $@

$(TEST_BINARY): $(CLOGLIKE_SOURCE) tests/INLAcirc_cloglike_test.c \
	INLAcirc_cloglike.h INLAcirc_common.h INLAcirc_cgeneric_compat.h $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. $(CLOGLIKE_SOURCE) \
		tests/INLAcirc_cloglike_test.c -lm -o $@

check: $(CLOGLIKE_OBJECT) $(TEST_BINARY)
	$(TEST_BINARY)
	@symbols="$$(nm -g --defined-only $(CLOGLIKE_OBJECT) | awk '{print $$3}')"; \
	bad="$$(printf '%s\n' "$$symbols" | grep -Ev '^INLAcirc_' || true)"; \
	test -z "$$bad" || { printf 'Unnamespaced symbols:\n%s\n' "$$bad"; exit 1; }

clean:
	rm -rf $(BUILD_DIR)
