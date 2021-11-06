CFLAGS=-O0 -g -Wall -Wextra -Wno-unused-parameter
CPPFLAGS= -Iexternals/include
LDLIBS=-lraylib -lm -ldl -lpthread
LDFLAGS=-Lexternals/libs

OBJS=tests/test_math.o

.PHONY: all clean

all: test_math

test_math: $(OBJS)
	g++ $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o: %.cpp
	g++ -c $(CFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(OBJS) test_math && clear