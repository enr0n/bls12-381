CC = gcc

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

INCLUDE = -Isrc -Iinclude
CFLAGS = -fPIC -Wall -Wextra -Werror -O2 $(INCLUDE) -g
LDFLAGS = -shared

RM = rm -f

TARGET_LIB = libbls12381.so

SRCS = $(wildcard src/*.c)
OBJS = $(SRCS:.c=.o)

.PHONY: all
all: $(TARGET_LIB)

$(TARGET_LIB): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	-$(RM) $(TARGET_LIB) $(OBJS) $(SRCS:.c=.d)

.PHONY: install
install: $(TARGET_LIB)
	install -d $(PREFIX)/lib/
	install -m 644 $(TARGET_LIB) $(PREFIX)/lib/
	install -d $(PREFIX)/include/
	install -m 644 include/BLS12_381.h $(PREFIX)/include/
