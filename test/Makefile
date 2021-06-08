.PHONY: all
all: test-fp2 test-fp6 test-fp12
	./tmp/test-fp2
	./tmp/test-fp6
	./tmp/test-fp12

test-fp2:
	mkdir -p tmp
	gcc -o tmp/test-fp2 -I../src fp2_test.c ../src/*.c -lgmp

test-fp6:
	mkdir -p tmp
	gcc -o tmp/test-fp6 -I../src fp6_test.c ../src/*.c -lgmp

test-fp12:
	mkdir -p tmp
	gcc -o tmp/test-fp12 -I../src fp12_test.c ../src/*.c -lgmp

