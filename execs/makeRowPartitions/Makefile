SRC=Framework.C
LIBPATH=../gsl/lib
INCLPATH=../gsl/include

CC=g++
CFLAGS = -g -std=c++0x
LFLAG = -lgsl -lgslcblas 

ALL: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH) -I $(INCLPATH)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o makeRowPartitions
clean:
	rm makeRowPartitions
