all: lab1/ lab2/ lab3/
	(cd lab1; make all)
	(cd lab2; make all)
	(cd lab3; make all)

check: lab1/ lab2/ lab3/
	(cd lab1; make check)
	(cd lab2; make check)
	(cd lab3; make check)

.PHONY: clean

clean: lab1/ lab2/ lab3/
	(cd lab1; make clean)
	(cd lab2; make clean)
	(cd lab3; make clean)
