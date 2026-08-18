CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wpedantic
CPPFLAGS =

CLOGLIKE_SOURCE = INLAcirc_cloglike.c
CLOGLIKE_OBJECT = INLAcirc_cloglike.o
TEST_BINARY = INLAcirc_cloglike_test

.PHONY: all check clean

all: $(CLOGLIKE_OBJECT)

$(CLOGLIKE_OBJECT): $(CLOGLIKE_SOURCE) INLAcirc_cloglike.h \
	INLAcirc_common.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -fPIC -c $(CLOGLIKE_SOURCE) -o $@

$(TEST_BINARY): $(CLOGLIKE_SOURCE) tests/INLAcirc_cloglike_test.c \
	INLAcirc_cloglike.h INLAcirc_common.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. $(CLOGLIKE_SOURCE) \
		tests/INLAcirc_cloglike_test.c -lm -o $@

check: $(CLOGLIKE_OBJECT) $(TEST_BINARY)
	./$(TEST_BINARY)
	@symbols="$$(nm -g --defined-only $(CLOGLIKE_OBJECT) | awk '{print $$3}')"; \
	bad="$$(printf '%s\n' "$$symbols" | grep -Ev '^INLAcirc_' || true)"; \
	test -z "$$bad" || { printf 'Unnamespaced symbols:\n%s\n' "$$bad"; exit 1; }

clean:
	rm -rf $(CLOGLIKE_OBJECT) $(TEST_BINARY)

