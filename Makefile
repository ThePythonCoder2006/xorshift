CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic
DBFLAGS = -ggdb3 -gdwarf-5
FASTFLAGS = -O3
LFLAGS = -mavx2

fast:
	$(CC) .\square_matrix.c -o mat.exe $(CFLAGS) $(FASTFLAGS) $(LFLAGS)
main:
	$(CC) .\square_matrix.c -o mat.exe $(CFLAGS) $(LFLAGS)
db:
	$(CC) .\square_matrix.c -o mat.exe $(CFLAGS) $(DBFLAGS) $(LFLAGS)