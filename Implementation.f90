!THIS FILE HAS TO CONTAIN THE MAIN PROGRAM.
PROGRAM main
  
  USE dynamical
  IMPLICIT NONE
  CLASS(PARTICLE), DIMENSION(:), ALLOCATABLE :: cluster
  !PARAMETERS OF THE SYSTEM
  REAL :: sigma,k,emm,mass,dt,du,dte,vo,dv,L,T
  INTEGER :: m,N,Nwarmup,stepdecay
  !AUXILIAR VARIABLES
  INTEGER :: u,v,w,ios1,ios2,ios3,i,seed,auxi,auxi2
  CHARACTER(len=100) :: name1,name2,name3
  CHARACTER(len=100) :: num1char,num2char,num3char,num4char

  IF(COMMAND_ARGUMENT_COUNT().NE.4)THEN
     WRITE(*,*)'ERROR, 4 COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
     STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,num1char)   ! The value identifies the files of the execution 
  CALL GET_COMMAND_ARGUMENT(2,num2char)   ! The value identifies the files of the execution 
  CALL GET_COMMAND_ARGUMENT(3,num3char)   ! The value identifies the files of the execution 
  CALL GET_COMMAND_ARGUMENT(4,num4char)   ! The value identifies the files of the execution 
 
  READ(num1char,*)N !Number of measurements
  READ(num2char,*)Nwarmup ! Number of steps to equilibrate the system RECOMMENDED Nwarmup=300
  READ(num3char,*) T ! Temperature dK 
  READ(num4char,*) L ! Length of the system RECOMMENDED sigma*10^-1  L=10.229

  
  !THE PARAMETERS OF THE SYSTEM ARE DEFINED!!!

   T=T/(10.0) ! Temperature [K]
   L=L/(10.0) ! Length [sigma]
   sigma=3.4 ! Length Scale [Angstroms]
   k=1.380649 ! Boltzmann Constant [ 10^{-23}m^2 kg s^{-2} K^{-1}]
   Emm=120.0 !Depth of the well in K [J/k_B]
   mass=39.95*1.6747 !Mass Argon [10^{-24} g ]
   dt=0.01 !Quantity in picoseconds [10^{-12} s]
   du=dt/sigma*(Emm*k/mass)**1/2 !Adimensional time
   dte=1  ! Variation allowed of the temperature 
   vo=(T/Emm)**(0.5) !Initial adimensional velocity
   dv=dte/Emm ! Variation of velocity allowed
   m=1000 ! Number of particles
   stepdecay=50 !Decay step of time were the correlations are negligible
  
  seed=20 !The seed for the random numbers is set
 
  CALL initialize(seed)
  ALLOCATE(cluster(m))

  !The name of the files which store the data are created:

  name1="WARMUPCHECK-N"//TRIM(num1char)//"Nw"//TRIM(num2char)//"T"//TRIM(num3char)//"L"//TRIM(num4char)//".txt" !WARMUPCHECK saves the potential and kinetic energy before the measurement 
  name2="HISTORY-N"//TRIM(num1char)//"Nw"//TRIM(num2char)//"T"//TRIM(num3char)//"L"//TRIM(num4char)//".txt" !HISTORY saves the positions,velocity of each time step
  name3="RPOSITIONS-N"//TRIM(num1char)//"Nw"//TRIM(num2char)//"T"//TRIM(num3char)//"L"//TRIM(num4char)//".txt" !RPOSITIONS saves the relative positions between particles of each time step
  
  !The units for the files are created

  w=50
  u=20
  v=30
  
  ! The system is simulated 10 times in order to obtain uncorrelated data for the
  ! pair correlation function
     
     CALL initcond(vo,cluster,L,10,10,10)
     CALL dynamics(cluster,L)

     OPEN(UNIT=w,IOSTAT=ios1, FILE=name1,STATUS='replace',ACTION='write')
 
     DO i=0,Nwarmup
        
        WRITE(w,*) i*du,ke(cluster),pe(cluster,L,0) !Stores the energy
        CALL verlet(du,cluster,L) !Performs the verlet integration
        CALL thermalize(cluster,vo,dv) !Controls the temperature

     END DO

     CLOSE(w)

     OPEN(UNIT=u,IOSTAT=ios2, FILE=name2,STATUS='replace',ACTION='write')
     OPEN(UNIT=v,IOSTAT=ios3, FILE=name3,STATUS='replace',ACTION='write')

     auxi2=0
     DO i=1,N*stepdecay+100

        !CHECK if it is neccesary to print the difference in the positions for the pair correlation function
        
        if(MOD(i,stepdecay)==0 .and. auxi2<11) THEN
           auxi2=auxi2+1
           auxi=v !The difference in positions are printed
        ELSE
           auxi=0 ! The difference in positions are not printed
        ENDIF

        CALL printing(cluster,u,i*du,L,auxi) !Prints the positions,velocity,energy of all the particles 
        CALL verlet(du,cluster,L) !Performs the Verlet integration
        
     END DO
     
     CLOSE(u)
     CLOSE(v)
  
END PROGRAM main
