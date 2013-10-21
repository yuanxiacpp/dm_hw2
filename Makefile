all: project2

project2: project2.c
	gcc -o project2 project2.c -lm

clean:
	rm *~ project2
