CC = g++
CFILES = test_mult.cpp
all: vector
	${CC} -Wall -llapack -lblas -L . -lvector ${CFILES} -o linear.out
vector:
	${CC} -c vector.cpp -o libvector.a
clean:
	rm -f *.out *.a
