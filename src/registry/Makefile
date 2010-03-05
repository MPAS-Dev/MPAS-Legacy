.SUFFIXES: .c .o

OBJS = parse.o dictionary.o gen_inc.o fortprintf.o

all: parse

parse: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	$(RM) *.o parse

.c.o:
	$(CC) -c $<
