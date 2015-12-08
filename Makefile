CC = mpicxx
CFILES = test_mult.cpp
ifeq ($(val), 1)
	DFLAGS = -DPRECOND
endif
all: vector
	${CC} -Wall -llapack -lblas -L . -lvector ${CFILES} -o linear.out ${DFLAGS}
vector:
	${CC} -c vector.cpp -o libvector.a
clean:
	rm -f *.out *.a
