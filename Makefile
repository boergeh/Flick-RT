
# Read README.md for information on compilation and running


all:	check-env build test

with_python:	all with_python

build:
	cd environment; make obj link
	cd astronomy; make obj link
	cd numeric/linalg; make obj link
	cd numeric; make obj link
	cd numeric/legendre; make obj link
	cd numeric/wigner; make obj link
	cd numeric/spherical_harmonics; make obj link
	cd mie; make obj link
	cd geometry; make obj link
	cd polarization; make obj link
	cd component; make obj link
	cd material; make obj link
	cd coating; make obj link
	cd material/gas; make obj link
	cd material/aerosols; make obj link
	cd material/water; make obj link	
	cd material/ice; make obj link	
	cd material/marine_cdom; make obj link	
	cd material/marine_particles; make obj link	
	cd accurt_api; make obj link
	cd transporter; make obj link
	cd radiator; make obj link
	cd model; make obj link
	cd main; make obj link
test:
	@cd environment; make test
	@cd astronomy; make test
	@cd numeric/linalg; make test
	@cd numeric; make test
	@cd numeric/legendre; make test
	@cd numeric/wigner; make test
	@cd numeric/spherical_harmonics; make test
	@cd mie; make test
	@cd geometry; make test
	@cd polarization; make test
	@cd component; make test
	@cd material; make test
	@cd material/gas; make test
	@cd material/aerosols; make test
	@cd material/water; make test
	@cd material/ice; make test
	@cd material/marine_cdom; make test
	@cd material/marine_particles; make test
	@cd accurt_api; make test
	@cd coating; make test
	@cd transporter; make test
	@cd model; make test
	@cd radiator; make test
clean:
	cd environment; make clean
	cd astronomy; make clean	
	cd numeric/linalg; make clean
	cd numeric; make clean
	cd numeric/legendre; make clean
	cd numeric/wigner; make clean
	cd numeric/spherical_harmonics; make clean
	cd mie; make clean
	cd geometry; make clean
	cd polarization; make clean
	cd component; make clean	
	cd material; make clean
	cd material/gas; make clean
	cd material/aerosols; make clean
	cd material/water; make clean
	cd material/ice; make clean
	cd material/marine_particles; make clean
	cd material/marine_cdom; make clean
	cd accurt_api; make clean
	cd coating; make clean
	cd transporter; make clean
	cd radiator; make clean
	cd model; make clean
	cd main; make clean
	rm -f *~
with_python:
	@echo ''
	@echo 'Testing all python scripts. May take an hour ...'
	cd Example/python_plots; python3 test_all.py
	cd Example/accurt_calls/logo; python3 test_all.py
	cd Example/accurt_calls/atmosphere_ocean; python3 test_all.py
check-env:
ifndef FLICK_PATH
	./update_shell.sh
endif

