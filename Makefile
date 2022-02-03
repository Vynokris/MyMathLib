CFLAGS   = -O0 -g -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable
CPPFLAGS = -Iexternals/include
LDLIBS   = -lm -lpthread

OBJS = externals/include/my_math.o tests/test_math.o

.PHONY: all clean

all: test_math

test_math: $(OBJS)
	g++ $(CFLAGS) $(CPPFLAGS) $(LDLIBS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) test_math && clear