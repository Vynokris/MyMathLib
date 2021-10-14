CFLAGS=-O0 -g -Wall -Iinclude
CPPFLAGS=-MMD
LDFLAGS=
LDLIBS=

OBJS=tests/test_math.o src/my_math.o
DEPS=$(OBJS:.o=.d)

.PHONY: all clean

all: test_math

-include $(DEPS)

test_math: $(OBJS)
	gcc $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o: %.c
	gcc -c $(CFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) test_math