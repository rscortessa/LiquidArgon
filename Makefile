all:main.x
OBJ= statNM.f90 libNM.f90 Dyn2.f90 NM.f90
N=50
Nwarmup=300
T=954
L=102
A= N Nwarmup T L

main.x:	$(OBJ) 
	gfortran $^ -O3  -o $@

graphs: main.x read.py
	./main.x $(N) $(Nwarmup) $(T) $(L)
	python3 read.py $(N) $(Nwarmup) $(T) $(L)

.PHONY:clean
clean:
	rm -f *.x *.o  
