! THIS FILE SHOULD CONTAIN THE DYNAMICS OF THE SYSTEM TO STUDY.
! WE HAVE TO INCLUDE THE FORCES, VELOCITY, POSITION AND ACELERATION EQUATIONS
! THE THERMOSTAT HAS TO BE INCLUDED AS WELL.
MODULE class_particle
  USE stat
  IMPLICIT NONE
  TYPE :: particle
     REAL, DIMENSION(3) :: x, v, a
   CONTAINS
     PROCEDURE :: setx => setx
     PROCEDURE :: setv => setv
     PROCEDURE :: seta => seta
     PROCEDURE :: state => state
  END TYPE particle

CONTAINS
  
  SUBROUTINE setx(p,xn,L)
    CLASS(PARTICLE) :: p
    REAL, DIMENSION(3) :: xn
    REAL, INTENT(IN) :: L
    INTEGER :: i
    DO i=1,3
       p%x(i)=xn(i)-NINT(xn(i)/L)*L
    END DO
  END SUBROUTINE setx
  
  SUBROUTINE setv(p,vn)  
    CLASS(PARTICLE) :: p
    REAL, DIMENSION(3) :: vn
    p%v=vn
  END SUBROUTINE setv

  SUBROUTINE seta(p,an)  
    CLASS(PARTICLE) :: p
    REAL, DIMENSION(3) :: an
    p%a=an
  END SUBROUTINE seta

  SUBROUTINE state(p,u)
    CLASS(PARTICLE) :: p
    INTEGER :: u
    WRITE(u,*) p%x(:),p%v(:),NORM2(p%v(:))**2
  END SUBROUTINE state

END MODULE class_particle



