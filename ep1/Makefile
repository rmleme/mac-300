CFLAGS=-Wall -g ${INC} ${LIB}
CC=gcc

cholesky.o: cholesky.c
	$(CC) $(CFLAGS) -c $*.c

cholesky : cholesky.o
	${CC} ${CFLAGS} -o cholesky cholesky.o -lm

LU.o: LU.c
	$(CC) $(CFLAGS) -c $*.c

LU : LU.o
	${CC} ${CFLAGS} -o LU LU.o -lm
clean: 
	rm -f *.o cholesky LU *~