! THIS FILE SHOULD CONTAIN THE DYNAMICS OF THE SYSTEM TO STUDY.
! WE HAVE TO INCLUDE THE FORCES, VELOCITY, POSITION AND ACELERATION EQUATIONS
! THE THERMOSTAT HAS TO BE INCLUDED AS WELL.
MODULE class_particle
  IMPLICIT NONE
  TYPE :: particle
     REAL, DIMENSION(3) :: position, velocity, acceleration
     REAL:: mass, radius
  END TYPE particle

contains
  
  SUBROUTINE  potential(ps,L)
    REAL :: norm
    REAL, INTENT(IN) :: L
    CLASS(PARTICLE) , DIMENSION(m), intent(inout) :: ps
    REAL,DIMENSION(3) :: potential,aux
    INTEGER :: i,j,m
    m=size(ps)
    DO i=1,m
       ps(i)%acceleration=0
    END DO
    
    DO i=1,m,1
    DO j=i+1,m,1
       aux=ps(i)%position-ps(j)%position
       DO l=1,3
          aux[l]=aux[l]-nint(aux[l]/L)*L
       END DO
       norm=NORM2(aux)
       IF (norm <= 2.25) THEN
          ps(i)%acceleration=ps(i)%acceleration+24*(2/norm**14-1/norm**8)*aux
          ps(j)%acceleration=ps(j)%acceleration-24*(2/norm**14-1/norm**8)*aux
       ENDIF
    END DO
    END DO
    
  END SUBROUTINE potential

  SUBROUTINE verlet(dt,ps,L)
    IMPLICIT NONE
    INTEGER ::i,m
    INTEGER,INTENT(IN) :: L
    REAL, DIMENSION(m) :: aux
    m=size(ps)
    DO i=1,m
       ps(i)%position=ps(i)%position+velocity*dt+1/2*acceleration*dt**2
       aux(i)=ps(i)%acceleration
    END DO
    CALL potential(ps,L)
    DO i=1,m
    ps(i)%velocity=ps(i)%velocity+(aux(i)+ps(i)%acceleration)*dt/2
    END DO
    
  END SUBROUTINE verlet
  
END MODULE class_particle
  
  PROGRAM test
    USE class_particle
    IMPLICIT NONE
    TYPE (particle) :: a
    a%radius=1.0
    PRINT*,a%radius
    
  END PROGRAM test
