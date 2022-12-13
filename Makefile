
# Read README.md for information on compilation and running

all:	check-env build test

build:	
	cd numeric; make obj link
	cd geometry; make obj link
	cd polarization; make obj link
	cd component; make obj link
	cd material; make obj link
	cd coating; make obj link
	cd material/gas; make obj link
	cd material/aerosols; make obj link
	cd transporter; make obj link
	cd radiator; make obj link
	cd model; make obj link
	cd main; make obj link
test:	
	@cd numeric; make test
	@cd geometry; make test
	@cd polarization; make test
	@cd component; make test
	@cd material; make test
	@cd material/gas; make test
	@cd material/aerosols; make test
	@cd coating; make test
	@cd transporter; make test
	@cd model; make test
	@cd radiator; make test
clean:	
	cd numeric; make clean
	cd geometry; make clean
	cd polarization; make clean
	cd component; make clean	
	cd material; make clean
	cd material/gas; make clean
	cd material/aerosols; make clean
	cd coating; make clean
	cd transporter; make clean
	cd radiator; make clean
	cd model; make clean
	cd main; make clean
	rm -f *~
check-env:
ifndef FLICK_PATH
	@echo
	@echo "Cannot find Flick RT environmental variables."
	@cat Prerequisites
	@echo "Try to run make again after updating your shell script."
	@echo
	@exit 1
endif
