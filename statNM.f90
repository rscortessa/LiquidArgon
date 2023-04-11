! This file should contain the functions to analyze the data
! The pair correlation functions etc.
MODULE stat

CONTAINS

    SUBROUTINE initialize(o)
    INTEGER, INTENT(IN) :: o
    INTEGER :: n,i
    INTEGER,ALLOCATABLE :: seed(:)
    CALL random_seed(size=n)
    ALLOCATE(seed(n))
    seed = o
    ! putting
    ! arbitrary seed to all elements
    CALL random_seed(put=seed)
    DEALLOCATE(seed)   
  END SUBROUTINE initialize

  
  SUBROUTINE sort(T)
    IMPLICIT NONE
    REAL ,DIMENSION(:), INTENT(INOUT)  :: T
    REAL, DIMENSION(:), ALLOCATABLE :: Taux
    LOGICAL ,DIMENSION(:), ALLOCATABLE :: sorting
    INTEGER :: length,i
    length=size(T)
    ALLOCATE(sorting(length))
    ALLOCATE(Taux(length))
    Taux=0
    sorting=.TRUE.
    DO i=1,length
       Taux(i)=MINVAL(T,sorting)
       sorting(MINLOC(T,sorting))=.FALSE.
    END DO
    T=Taux
    DEALLOCATE(Taux,sorting)
  END SUBROUTINE sort

  INTEGER FUNCTION diaconis(T) RESULT(N)
    REAL, DIMENSION (:), INTENT(INOUT) :: T
    REAL :: dx,min,max
    INTEGER :: M
    M=size(T)
    min=T(1)
    max=T(M)
    !COMPUTE THE NUMBER OF BINS AND ITS DISTANCE
    dx=(T(M*3/4)-T(M/4))*2/(M*1.0)**(1/3.0)
    N=FLOOR((max-min)/(dx))+1
  END FUNCTION diaconis

    SUBROUTINE histogram(T,hist,l,u)
    REAL, DIMENSION (:), INTENT(IN) :: T
    INTEGER, INTENT(IN) ::l,u
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: hist
    REAL :: dx,min,max,x
    INTEGER :: M,N,j,i
    !INITIALIZE VARIABLES
    ALLOCATE(hist(l,2))
    max=MAXVAL(T)
    min=MINVAL(T)
    N=l
    M=size(T)
    dx=(max-min)/(N*1.0)
    !CREATE THE HISTOGRAM
    hist=0
    x=min+dx
    i=1
    j=1
    hist(1,1)=x-0.5*dx
    DO WHILE(j<=M)
       IF(x>T(j)) THEN
          hist(i,2)=hist(i,2)+1.0
          j=j+1
       ELSE
          i=i+1
          x=x+dx
          hist(i,1)=x-0.5*dx
       ENDIF
    END DO
    hist(:,2)=hist(:,2)/(M*dx)
  END SUBROUTINE histogram


  SUBROUTINE gaussian(u,seed,sigma)  
    REAL, DIMENSION(:,:), INTENT(INOUT) :: u
    INTEGER, INTENT(IN) :: seed
    REAL, INTENT(IN) :: sigma
    REAL, DIMENSION(:,:), ALLOCATABLE :: x
    INTEGER :: i,j,n,m
    CALL initialize(seed)
    n=size(u(:,1))
    m=size(u(1,:))
    ALLOCATE(x(m,2))
    DO i=1,n
       DO j=1,m
          CALL RANDOM_NUMBER(x)
          u(i,j)=sigma*SQRT(-2.0*LOG(x(j,1)))*SIN(8*ATAN(1.D0)*x(j,2))
    END DO
    END DO
  END SUBROUTINE gaussian

  
END MODULE stat
