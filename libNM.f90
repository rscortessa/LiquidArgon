! THIS FILE SHOULD CONTAIN THE DYNAMICS OF THE SYSTEM TO STUDY.
! WE HAVE TO INCLUDE THE FORCES, VELOCITY, POSITION AND ACELERATION EQUATIONS
! THE THERMOSTAT HAS TO BE INCLUDED AS WELL.
MODULE class_particle
  USE stat
  IMPLICIT NONE
  TYPE :: particle
     REAL, DIMENSION(3) :: position, velocity, acceleration
     REAL:: mass, radius
  END TYPE particle

contains

  FUNCTION mimg(L,a)
    REAL, DIMENSION(3) :: mimg,a
    REAL, INTENT(IN) :: L
    INTEGER :: i

    DO i=1,3
       a(i)=a(i)-nint(a(i)/L)*L
    END DO

    mimg=a
    
  END FUNCTION mimg
  
  FUNCTION Vij(aux)
    REAL, DIMENSION(3) :: Vij,aux
    REAL :: norm
    norm=NORM2(aux)

    IF(norm <= 2.25) THEN
       Vij=24.0*(2.0/norm**14-1/norm**8)*aux
    ELSE
       Vij=0
    ENDIF

  END FUNCTION Vij

  
  SUBROUTINE  dynamics(ps,L)

    REAL :: norm
    REAL, INTENT(IN) :: L
    CLASS(PARTICLE) , DIMENSION(:), intent(inout) :: ps
    REAL,DIMENSION(3) :: aux
    INTEGER :: i,j,m

    m=size(ps)

    DO i=1,m
       ps(i)%acceleration=0
    END DO
    
    DO i=1,m
       DO j=i+1,m
          aux=ps(i)%position-ps(j)%position
          aux=mimg(L,aux)
          norm=NORM2(aux)
          ps(i)%acceleration=ps(i)%acceleration+Vij(aux)
          ps(j)%acceleration=ps(j)%acceleration-Vij(aux)
       END DO
    END DO
    
  END SUBROUTINE dynamics

  SUBROUTINE verlet(dt,ps,L)

    IMPLICIT NONE
    INTEGER ::i,m
    REAL,INTENT(IN) :: L
    REAL, DIMENSION(:,:), ALLOCATABLE :: aux
    REAL :: dt
    CLASS(PARTICLE) , DIMENSION(:), intent(inout) :: ps

    m=size(ps)
    ALLOCATE(aux(m,3))

    DO i=1,m
       ps(i)%position=ps(i)%position+ps(i)%velocity*dt+1/2.0*ps(i)%acceleration*dt**2
       ps(i)%position=mimg(L,ps(i)%position)
       aux(i,:)=ps(i)%acceleration
    END DO

    CALL dynamics(ps,L)

    DO i=1,m
       ps(i)%velocity=ps(i)%velocity+(aux(i,:)+ps(i)%acceleration)*dt/2
    END DO
    
  END SUBROUTINE verlet

  SUBROUTINE initcond(kt,ps,L)
    IMPLICIT NONE
    REAL , INTENT(IN) :: kt,L
    REAL, DIMENSION(:,:), ALLOCATABLE :: aux,aux2
    INTEGER :: i,m
    CLASS(PARTICLE), DIMENSION(:), intent(inout) :: ps
    m=size(ps)
    ALLOCATE(aux(m,3),aux2(m,3))
    !CALL gaussian(aux,0,1.0)
    CALL uniform(aux2,L,10,10,10)  
    DO i=1,m
       ps(i)%position=aux2(i,:)
       !ps(i)%velocity=aux(i,:)
       ps(i)%velocity=0.0
    END DO
  END SUBROUTINE initcond

  SUBROUTINE PRINTING(ps,u,dt)
    INTEGER , INTENT(IN):: u
    INTEGER :: m,i
    REAL, INTENT(IN) :: dt
    CLASS(PARTICLE), DIMENSION(:) :: ps
    m=size(ps)
    WRITE(u,*) "Time \t",dt
    DO i=1,m
       WRITE(u,*) ps(i)%position, ps(i)%velocity
    END DO
  END SUBROUTINE PRINTING
END MODULE class_particle
  
