all:	build test

build:	
	cd gas; make build
	$(FLICK_COMPILER) -o test_all test_all.cpp
test:
#	@cd gas; make test
	@./test_all
clean:
	cd gas; make clean
	rm -f ./test_all *~ ./water_input/*~
