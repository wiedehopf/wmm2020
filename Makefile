# on Unix, compile desired program using
# 
#  make geomag   ! C version of magpoint
#

DIALECT = -std=c11
CFLAGS += $(DIALECT) -O2 -g -W -D_DEFAULT_SOURCE -Wall -Werror

LIBS = -lm

# Main programs
#

%.o: %.c *.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

geomag: geomag.o
	$(CC) -g -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f geomag

