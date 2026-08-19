CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -Wpedantic
CPPFLAGS ?=

BUILD_DIR = build
CLOGLIKE_SOURCE = INLAcirc_cloglike.c
PRODUCTION_OBJECT = $(BUILD_DIR)/INLAcirc_cloglike.o
TEST_OBJECT = $(BUILD_DIR)/INLAcirc_cloglike_test_module.o
TEST_BINARY = $(BUILD_DIR)/INLAcirc_cloglike_test
TEST_CPPFLAGS = -DINLACIRC_TEST_CGENERIC_COMPAT

.PHONY: all check clean

# The production target intentionally requires the cgeneric.h supplied by
# the INLA source tree.  That header is not distributed in this directory.
all: $(PRODUCTION_OBJECT)

$(BUILD_DIR):
	mkdir -p $@

$(PRODUCTION_OBJECT): $(CLOGLIKE_SOURCE) INLAcirc_cloglike.h \
	INLAcirc_common.h cgeneric.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -fPIC -c $(CLOGLIKE_SOURCE) -o $@

$(TEST_OBJECT): $(CLOGLIKE_SOURCE) INLAcirc_cloglike.h \
	INLAcirc_common.h tests/INLAcirc_cgeneric_test_compat.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(TEST_CPPFLAGS) $(CFLAGS) -I. -fPIC \
		-c $(CLOGLIKE_SOURCE) -o $@

$(TEST_BINARY): $(CLOGLIKE_SOURCE) tests/INLAcirc_cloglike_test.c \
	INLAcirc_cloglike.h INLAcirc_common.h \
	tests/INLAcirc_cgeneric_test_compat.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(TEST_CPPFLAGS) $(CFLAGS) -I. \
		$(CLOGLIKE_SOURCE) tests/INLAcirc_cloglike_test.c \
		-lm -o $@

check: $(TEST_OBJECT) $(TEST_BINARY)
	./$(TEST_BINARY)
	@symbols="$$(nm -g --defined-only $(TEST_OBJECT) | awk '{print $$3}')"; \
	bad="$$(printf '%s\n' "$$symbols" | grep -Ev '^INLAcirc_' || true)"; \
	test -z "$$bad" || { printf 'Unnamespaced symbols:\n%s\n' "$$bad"; exit 1; }

clean:
	rm -rf $(BUILD_DIR)
