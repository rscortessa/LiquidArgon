!THIS FILE HAS TO CONTAIN THE MAIN PROGRAM.
PROGRAM main
  USE dynamical
  IMPLICIT NONE
  CLASS(PARTICLE), DIMENSION(:), ALLOCATABLE :: cluster
  INTEGER :: m,u,ios,i,max
  REAL :: k,L,dt,du,Emm,mass,sigma,vo,T,dte,dv
  CHARACTER(len=10) :: name

  sigma=3.4 !Angstroms
  k=1.380649 !Boltzmann Constant
  Emm=120.0 !Energy in K
  mass=39.95*1.6747 !Mass Argon
  dt=0.01 !Quantity in picoseconds
  du=dt/sigma*(Emm*k/mass)**1/2
  T=94.4
  dte=1
  vo=(T/Emm)**(0.5)
  dv=dte/Emm
  name="LA.txt"
  u=20
  OPEN(UNIT=u,IOSTAT=ios, FILE=name,STATUS='replace',ACTION='write')
  
  L=10.229
  m=1000
  i=1
  max=1000
  ALLOCATE(cluster(m))
  
  CALL initcond(vo,cluster,L,10,10,10)
  CALL dynamics(cluster,L)
  DO i=0,max
     CALL printing(cluster,u,i*du,L)
     CALL verlet(du,cluster,L)
     CALL thermalize(cluster,vo,dv)
  END DO
  CLOSE(u)
END PROGRAM main
