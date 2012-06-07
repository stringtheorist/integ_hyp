all: gauss_legendre_test log_test

test_utils.o: integ_hyp.o integ_hyp.h
	gcc -c -lm -lgsl -lgslcblas -I/usr/local/include test_utils.c

log_test: log_test.o integ_hyp.o integ_hyp.h test_utils.o test_utils.h
	gcc -lm -lgsl -lgslcblas -I/usr/local/include log_test.o integ_hyp.o test_utils.o -olog_test

log_test.o: log_test.c integ_hyp.o integ_hyp.h test_utils.o test_utils.h
	gcc -c -lm -lgsl -lgslcblas -I/usr/local/include log_test.c

test_funcs.o: test_funcs.c test_funcs.c
	gcc -c -lm -lgsl -lgslcblas -I/usr/local/include test_funcs.c  

gauss_legendre_test: gauss_legendre_test.o integ_hyp.o integ_hyp.h
	gcc -lm -lgsl -lgslcblas -I/usr/local/include gauss_legendre_test.o integ_hyp.o -ogauss_legendre_test

gauss_legendre_test.o: gauss_legendre_test.c integ_hyp.o integ_hyp.h
	gcc -c -lm -lgsl -lgslcblas -I/usr/local/include gauss_legendre_test.c

integ_hyp.o: integ_hyp.c integ_hyp.h
	gcc -c -lm -lgsl -lgslcblas -I/usr/local/include integ_hyp.c  

clean: 
	rm -rf *.o gauss_legendre_test
