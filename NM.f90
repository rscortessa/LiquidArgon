!THIS FILE HAS TO CONTAIN THE MAIN PROGRAM.
PROGRAM main
  USE class_particle
  IMPLICIT NONE
  CLASS(PARTICLE), DIMENSION(:), ALLOCATABLE :: cluster
  INTEGER :: m,u,ios,i,max
  REAL :: k,L,dt,Emm,sigma
  CHARACTER(len=10) :: name
  name="hola.txt"
  u=20
  OPEN(UNIT=u,IOSTAT=ios, FILE=name,STATUS='new',ACTION='write',POSITION="append")
  sigma=3.4
  L=10.229
  Emm=120.*1.380649/(39.95*1.6747)*10.0**(4)!Emin/m
  print*,Emm
  dt=10.**(-4)/sigma!0.1 nanosecond
  print*,dt
  !Choose the value
  dt=dt*Emm**(1/2.0)
  m=1000
  print*,dt
  i=1
  max=100
  ALLOCATE(cluster(m))
  CALL initcond(k,cluster,L)
  DO i=1,max
     CALL printing(cluster,u,i*dt)
     CALL verlet(dt,cluster,L)
  END DO
  CLOSE(u)
END PROGRAM main
