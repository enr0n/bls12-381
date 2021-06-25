.PHONY: all
all: test-fp2 test-fp6 test-fp12 test-G1
	./tmp/test-fp2
	./tmp/test-fp6
	./tmp/test-fp12
	./tmp/test-G1

test-fp2:
	mkdir -p tmp
	gcc -o tmp/test-fp2 -I../src fp2_test.c ../src/finite_field.c -lgmp

test-fp6:
	mkdir -p tmp
	gcc -o tmp/test-fp6 -I../src fp6_test.c ../src/finite_field.c -lgmp

test-fp12:
	mkdir -p tmp
	gcc -o tmp/test-fp12 -I../src fp12_test.c ../src/finite_field.c -lgmp

test-G1:
	mkdir -p tmp
	gcc -o tmp/test-G1 -I../src G1_test.c ../src/finite_field.c ../src/G1.c -lgmp
