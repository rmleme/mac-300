CFLAGS=-Wall -g ${INC} ${LIB}
CC=gcc

all: ep2-1 ep2-2

ep2-1.o: ep2-1.c
	$(CC) $(CFLAGS) -c $*.c

ep2-1: ep2-1.o
	${CC} ${CFLAGS} -o ep2-1 ep2-1.o -lm

ep2-2.o: ep2-2.c
	$(CC) $(CFLAGS) -c $*.c

ep2-2: ep2-2.o
	${CC} ${CFLAGS} -o ep2-2 ep2-2.o -lm

clean: 
	rm -f *.o ep2-1 ep2-2 *~