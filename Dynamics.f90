MODULE dynamical
  USE class_particle
  USE stat
  IMPLICIT NONE
  
CONTAINS  
  FUNCTION mimg(L,a)
    REAL, DIMENSION(3) :: mimg,a
    REAL, INTENT(IN) :: L
    INTEGER :: i

    DO i=1,3
       a(i)=a(i)-nint(a(i)/L)*L
    END DO

    mimg=a
    
  END FUNCTION mimg
  
  FUNCTION aij(aux)
    REAL, DIMENSION(3) :: aij,aux
    REAL :: norm
    norm=NORM2(aux)
    IF(norm <= 2.25) THEN
       aij=24.0*(2.0/norm**14-1/norm**8)*aux
    ELSE
       aij=0
    ENDIF

  END FUNCTION aij

  FUNCTION ke(ps)
    CLASS(PARTICLE), DIMENSION(:), INTENT(INOUT) :: ps
    INTEGER :: m,i
    REAL :: ke
    m=size(ps)
    ke=0
    DO i=1,m
       ke=ke+NORM2(ps(i)%v)**2
    END DO
    ke=0.5*ke
  END FUNCTION ke

   FUNCTION pe(ps,L,v)
    CLASS(PARTICLE), DIMENSION(:), INTENT(INOUT) :: ps
    INTEGER :: m,i,j
    REAL :: pe,norm,L 
    REAL, DIMENSION(3) :: aux
    INTEGER, INTENT(IN) :: v
    m=size(ps)
    pe=0
    IF (v/= 0) THEN
       DO i=1,m
          DO j=i+1,m
             aux=ps(i)%x-ps(j)%x
             aux=mimg(L,aux)
             norm=NORM2(aux)
             WRITE(v,*) norm
             IF (norm<2.25) THEN
                pe=pe+4.0/norm**12-4.0/norm**6
             ENDIF
          END DO
       END DO
       
    ELSE
       DO i=1,m
          DO j=i+1,m
             aux=ps(i)%x-ps(j)%x
             aux=mimg(L,aux)
             norm=NORM2(aux)
             IF (norm<2.25) THEN
                pe=pe+4.0/norm**12-4.0/norm**6
             ENDIF
          END DO
       END DO
    ENDIF
  END FUNCTION pe
  
  SUBROUTINE  dynamics(ps,L)

    REAL, INTENT(IN) :: L
    CLASS(PARTICLE) , DIMENSION(:), INTENT(INOUT) :: ps
    REAL,DIMENSION(3) :: aux
    INTEGER :: i,j,m

    m=size(ps)
    aux=0
    
    DO i=1,m
       CALL ps(i)%seta(aux)
    END DO
    
    DO i=1,m
       DO j=i+1,m
          aux=ps(i)%x-ps(j)%x
          aux=mimg(L,aux)
          ps(i)%a=ps(i)%a+aij(aux)
          ps(j)%a=ps(j)%a-aij(aux)
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
       CALL ps(i)%setx(ps(i)%x+ps(i)%v*dt+1/2.0*ps(i)%a*dt**2,L)
       aux(i,:)=ps(i)%a
    END DO

    CALL dynamics(ps,L)

    DO i=1,m
       CALL ps(i)%setv(ps(i)%v+(aux(i,:)+ps(i)%a)*dt/2.0)
    END DO
    
  END SUBROUTINE verlet

  
  SUBROUTINE uniform(u,L,n,m,t)  
    REAL, DIMENSION(:,:), INTENT(INOUT) :: u
    REAL, INTENT(IN) :: L
    INTEGER :: i,j,k
    INTEGER, INTENT(IN) :: n,m,t
    DO i=0,n-1
       DO j=0,m-1
          DO k=1,t
             u(i*t*m+j*t+k,3)=L/t*(k-0.5)-L/2.0
             u(i*t*m+j*t+k,2)=L/m*(j-0.5)-L/2.0
             u(i*t*m+j*t+k,1)=L/n*(i-0.5)-L/2.0
          END DO
       END DO
    END DO
    
  END SUBROUTINE uniform

  
  SUBROUTINE initcond(vo,ps,L,q,r,s)
    IMPLICIT NONE
    REAL , INTENT(IN) :: vo,L
    REAL, DIMENSION(:,:), ALLOCATABLE :: aux,aux2
    INTEGER :: i,m
    INTEGER, INTENT(IN) :: q,r,s
    CLASS(PARTICLE), DIMENSION(:), intent(inout) :: ps
    m=size(ps)
    ALLOCATE(aux(m,3),aux2(m,3))
    CALL uniform(aux2,L,q,r,s)
    CALL gaussian(aux,vo)
    DO i=1,m
       CALL ps(i)%setx(aux2(i,:),L)
       CALL ps(i)%setv(aux(i,:))
    END DO
  END SUBROUTINE initcond


  SUBROUTINE thermalize(ps,sigma,ds)
    REAL, INTENT(IN) :: sigma,ds
    CLASS(PARTICLE), DIMENSION(:) :: ps  
    REAL :: sum
    INTEGER :: m,i
    
    m=size(ps)
    sum=2.0*ke(ps)
    sum=sum/(m*3.0)
    
    IF(ABS(sum-sigma**2)>ds) THEN
       DO i=1,m
          ps(i)%v=ps(i)%v*sigma/sum**(0.5)
       END DO
    END IF
    
  END SUBROUTINE thermalize

  
  SUBROUTINE printing(ps,u,dt,L,v)
    INTEGER , INTENT(IN):: u,v
    INTEGER :: m,i
    REAL, INTENT(IN) :: dt,L
    CLASS(PARTICLE), DIMENSION(:) :: ps
    m=size(ps)
    WRITE(u,*) dt,ke(ps),pe(ps,L,v) 
    DO i=1,m
       CALL ps(i)%state(u)
    END DO
  END SUBROUTINE printing
  
END MODULE dynamical
  
