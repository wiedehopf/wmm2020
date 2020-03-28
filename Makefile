# compile the example program using the minimal geomag.h library

DIALECT = -std=c11
CFLAGS += $(DIALECT) -O2 -g -W -D_DEFAULT_SOURCE -Wall -Werror

LIBS = -lm

# Main programs
#

%.o: %.c *.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

example: geomag.o example.c
	$(CC) -g -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f example *.o

