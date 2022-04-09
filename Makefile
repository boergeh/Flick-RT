
# Read README.md for information on compilation and running

all:	check-env build test

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

check-env:
ifndef FLICK_PATHH
	@echo
	@echo "Cannot find flick environmental variables."
	@cat Prerequisites
	@echo "Try again running make after updating your shell script."
	@echo
	@exit 1
endif

AA=$(shell echo $0)
t:
	echo $(AA)
