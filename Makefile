all:	build test

build:	
	cd numeric; make build
	cd geometry; make build
	cd material; make build
	cd transporter; make build
	cd main; make build
test:	
	@cd numeric; make test
	@cd geometry; make test
	@cd material; make test
	@cd transporter; make test
clean:	
	cd numeric; make clean
	cd geometry; make clean
	cd material; make clean
	cd transporter; make clean
	cd main; make clean
	rm -f *~
