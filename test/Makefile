.PHONY: all
all: test-fp2 test-fp6 test-fp12 test-G1 test-G2 test-pairing
	./tmp/test-fp2
	./tmp/test-fp6
	./tmp/test-fp12
	./tmp/test-G1
	./tmp/test-G2
	./tmp/test-pairing

test-fp2:
	mkdir -p tmp
	gcc -o tmp/test-fp2 -I../include -I../src fp2_test.c ../src/finite_field.c -lgmp

test-fp6:
	mkdir -p tmp
	gcc -o tmp/test-fp6 -I../include -I../src fp6_test.c ../src/finite_field.c -lgmp

test-fp12:
	mkdir -p tmp
	gcc -o tmp/test-fp12 -I../include -I../src fp12_test.c ../src/finite_field.c -lgmp

test-G1:
	mkdir -p tmp
	gcc -o tmp/test-G1 -I../include -I../src G1_test.c ../src/finite_field.c ../src/G1.c -lgmp

test-G2:
	mkdir -p tmp
	gcc -o tmp/test-G2 -I../include -I../src G2_test.c ../src/finite_field.c ../src/G2.c -lgmp

test-pairing:
	mkdir -p tmp
	gcc -o tmp/test-pairing -I../include -I../src pairing_test.c ../src/*.c -lgmp
