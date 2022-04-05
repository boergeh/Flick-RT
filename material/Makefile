all:	b t

b:	
	cd gas; make b
	clang++ -std=c++20 -o test_all test_all.cpp
t:
	cd gas; make t
	./test_all
c:
	cd gas; make c
	rm -f ./test_all *~ ./water_input/*~
