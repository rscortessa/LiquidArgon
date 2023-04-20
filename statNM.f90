! This file should contain the functions to analyze the data
! The pair correlation functions etc.
MODULE stat

CONTAINS

    SUBROUTINE initialize(o)
    INTEGER, INTENT(IN) :: o
    INTEGER :: n
    INTEGER,ALLOCATABLE :: seed(:)
    CALL random_seed(size=n)
    ALLOCATE(seed(n))
    seed = o
    ! putting
    ! arbitrary seed to all elements
    CALL random_seed(put=seed)
    DEALLOCATE(seed)   
  END SUBROUTINE initialize

  SUBROUTINE gaussian(u,sigma)  
    REAL, DIMENSION(:,:), INTENT(INOUT) :: u
    REAL, INTENT(IN) :: sigma
    REAL, DIMENSION(:,:), ALLOCATABLE :: x
    INTEGER :: i,j,n,m
    n=size(u(:,1))
    m=size(u(1,:))
    ALLOCATE(x(m,2))
    DO i=1,n
       DO j=1,m
          CALL RANDOM_NUMBER(x)
          u(i,j)=sigma*SQRT(-2.0*LOG(x(j,1)))*SIN(8.0*ATAN(1.0)*x(j,2))
    END DO
    END DO
  END SUBROUTINE gaussian

  
END MODULE stat
