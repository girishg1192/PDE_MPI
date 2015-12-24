CC = mpicxx
CFILES = test_mult.cpp
ifeq ($(precond), 1)
	DFLAGS = -DPRECOND
endif
ifeq ($(scale), 1)
	DFLAGS := -DSCALING
endif
all: vector
	${CC} -Wall -L . -lvector ${CFILES} -o linear.out ${DFLAGS}
weak: vector
	${CC} -Wall -L . -lvector ${CFILES} -o weak.out ${DFLAGS} -DSCALING
vector:
	${CC} -c vector.cpp -o libvector.a
clean:
	rm -f *.out *.a
