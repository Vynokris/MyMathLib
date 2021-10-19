CFLAGS=-O0 -g -Wall -Wextra -Wno-unused-parameter
CPPFLAGS=-MMD -Iexternals/include
LDLIBS=-lraylib -lm -ldl -lpthread
LDFLAGS=-Lexternals/libs

OBJS=tests/test_math.o
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