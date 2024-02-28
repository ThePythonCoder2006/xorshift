CC = gcc

N = 32# default

DEFINE_N = -DNUM="($(N)U)"

CFLAGS = -Wall -Wextra -pedantic
DBFLAGS = -ggdb3 -gdwarf-5
FASTFLAGS = -O2
LFLAGS = -mavx2

SQR_MAT_SRC = .\square_matrix.c .\constants.c

fast:
	$(CC) $(DEFINE_N) $(SQR_MAT_SRC) -o mat $(CFLAGS) $(FASTFLAGS) $(LFLAGS)
plain:
	$(CC) $(DEFINE_N) $(SQR_MAT_SRC) -o mat $(CFLAGS) $(LFLAGS)
db:
	$(CC) $(DEFINE_N) $(SQR_MAT_SRC) -o mat $(CFLAGS) $(DBFLAGS) $(LFLAGS)
pre:
	$(CC) $(DEFINE_N) $(SQR_MAT_SRC) -o mat_pre.c -E $(CFLAGS) $(LFLAGS)

lut:
	$(CC) $(DEFINE_N) .\generate_lut.c -o lut $(CFLAGS) $(LFLAGS) $(FASTFLAGS)