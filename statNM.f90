! This file should contain the functions to analyze the data
! The pair correlation functions etc.
MODULE stat

CONTAINS
  
  SUBROUTINE initialize(o)
    INTEGER, INTENT(IN) :: o
    INTEGER :: n,i
    INTEGER,ALLOCATABLE :: seed(:)
    REAL :: r
    CALL random_seed(size=n)
    ALLOCATE(seed(n))
    seed = o
    ! putting
    ! arbitrary seed to all elements
    CALL random_seed(put=seed)
    DEALLOCATE(seed)   
  END SUBROUTINE initialize

  SUBROUTINE gaussian(u,seed,sigma)  
    REAL, DIMENSION(:,:), INTENT(INOUT) :: u
    INTEGER, INTENT(IN) :: seed
    REAL, INTENT(IN) :: sigma
    INTEGER :: i,ios,n
    REAL, DIMENSION(3) :: x
    CALL initialize(seed)
    n=size(u(:,1))
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(x)
    PRINT*, "BEFORE"
    PRINT*, u
    DO i=1,n
       u(i,1)=SQRT(-2.0*sigma**2*LOG(u(i,1)))*SIN(8*ATAN(1.D0)*x(1))
       u(i,2)=SQRT(-2.0*sigma**2*LOG(u(i,2)))*SIN(8*ATAN(1.D0)*x(2))
       u(i,3)=SQRT(-2.0*sigma**2*LOG(u(i,3)))*SIN(8*ATAN(1.D0)*x(3))
    END DO
    PRINT*, "AFTER"
    PRINT*,u
    PRINT*, "FINISH"
  END SUBROUTINE gaussian

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

  
END MODULE stat
