all:main.x
OBJ= Statistical.f90 	ClassParticle.f90 Dynamics.f90 Implementation.f90
N=50
Nwarmup=300
T=954
L=102


main.x:	$(OBJ) 
	gfortran $^ -O3  -o $@

graphs: main.x read.py
	./main.x $(N) $(Nwarmup) $(T) $(L)
	python3 read.py $(N) $(Nwarmup) $(T) $(L)

.PHONY:clean
clean:
	rm -f *.x *.o  
