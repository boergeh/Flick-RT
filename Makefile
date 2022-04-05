all:	b t
b:	
	cd numeric; make b
	cd geometry; make b
	cd material; make b
	cd transporter; make b
	cd main; make b
t:	
	cd numeric; make t
	cd geometry; make t
	cd material; make t
	cd transporter; make t
	cd main; make t
c:	
	cd numeric; make c
	cd geometry; make c
	cd material; make c
	cd transporter; make c
	cd main; make c
	rm -f *~
