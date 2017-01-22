!****************************************************************************************************************
!
!  vgwmolecular.f90: VGW Molecular module
!
!  Purpose: The module is to be used for propagating Variational Gaussian Wavepackets (VGW's) for
!           molecular systems in imaginary time beta (inverse temperature 1/T), where the potential
!           of interaction for all particle must be represented as a sum of Gaussian functions. The
!           parameters for the Guassian functions c_i, a_i, where c_i is the Gaussian coefficient and
!           a_i the exponent are read in from the file gauss.dat. See doumentation found in 
!           "vgw-molecular.pdf" for a detailed description.
!  
!      Use: The module can be used with any outside program or with the modules provided for Monte-Carlo 
!           simulation in the NVT/NPT ensembles found in the files:
!     
!           vgwmolmain.f90: The main program file that implements the othjer modules and reads in the
!                           the program parameters.
!
!             vgwmolmc.f90: The module that performs the Monte-Carlo steps, which include translations,
!                           volume changes, and parallel swapping.
!
!            vgwmolobs.f90: This module updates observables and ouputs the program data to files    
!             
!             vgwmodpt.f90: The module that performs the parallel swapping
!
!  Subroutines:
!
!            INITIALIZE_VGWM(N,TMAX,TMIN,PMIN,PMAX,TI,AT,BLX,NP,GSIZE,NODE,RCT,ENT,EG)
!                Initializes the variables and arrays.
!
!            VGWFSETUP(Q,CONT)
!                Sets up output files and Monte-Carlo environment parameters.
!
!            VGWGRIDINIT(TMAX,TMIN,PMIN,PMAX,EG)
!                Initializes the grid of inverse temperature (beta).
!
!            VGWZEROENRG(Q)     
!                Finds the energy of configuration Q to be used as zero energy.
!
!            VGWMQUENCH(Q)      
!                Main quenching subroutine. Propagtes VGW to inverse temperature 1/TMIN.
!            
!            VGWCQNCH(Q,LNP,UL) 
!                Continues the quenching of a VGW to a lower temperature for replica exchange.
!
!            PAIRS(Q)
!                Determines all pairs with rij less than rcutoff
!
!            INITDLSODE
!                Initializes 
!
!            LNPS(LOGZ)
!                Calculates the VGW density
!           
!            INVDET(A,M,DET)
!                Calculates the determinant of the inverse of A
!    
!            RHSM(NEQM,TT,YM,YPRIME)
!                Computes the right-hand-side of the molecular VGW differential equation
!
!            RHSFC(IIND,JIND,QIJ,U,UPV,G,UPM,MNUM,TYPEI,TYPEJ)
!                Computes gradient, Hessian, and energy for fully-coupled VGW's for RHSM
!    
!            RHSSP(IIND,JIND,QIJ,U,UPV,G,UPM,MNUM1,MNUM2,TYPEI,TYPEJ)
!                Computes gradient, Hessian, and energy for fully-coupled VGW's for RHSM
!
!            FILLBLOCKS(UPM)
!                Fills lower half of the matrix used to compute G
!
!            GAUSSENERGY(Q,UGAUSS)
!                Computes the classical energy of the configuration using the Gaussian potentials 
!
!   Variables:
!
!            N_ATOM : Number of atoms in the molecule
!               NEQ : Total number of equations dlsode needs to propagate 
!                     = 1 + 3*N_ATOM + N_MOLECULES*(9*ith_Mol_Size^2 + 3*ith_Mol_Size)/2
!              NMOL : Number of molecules
!           LRGMSZE : Largest molecule size
!             NREPS : Number of parallel replicas
!            NPAIRS : Number of pairs ij
!            N3ATOM : 3 x N_ATOM
!          GRIDSIZE : Size of the inverse-temperature (beta) grid
!         IARRAY(7) : Array of integer parameters fro dlsode
!           FILEIND : Output file unit number
!          VFILEIND : VMD coordinate output file unit number
!          QFILEIND : Coordinate output file unit number
!        VOLFILEIND : Volume output file unit number
!             ASIZE : The size of the data arrays
!            mynode : The current node
!      MOLDATA(:,:) : Contains the molecular data (atom type, molid)
!          IWORK(:) : Integer working array
!         NUMOBSERV : Number of observables that are computed
!          MAXGAUSS : Maximum number of Gaussians in sum
!           MAXTYPE : Maximum number of atom types
!              ATOL : Absolute error tolerance for dlsode
!  
!****************************************************************************************************************

MODULE VGWM
     IMPLICIT NONE
     INTEGER, SAVE :: N_ATOM,NEQ,N_MOL,LRGMSZE,NREPS,NPAIRS,N3ATOM,PTSIZE
     INTEGER, SAVE :: GRIDSIZE,IARRAY(7),FILEIND,VFILEIND,QFILEIND,VOLFIND,ASIZE,mynode
     INTEGER, ALLOCATABLE, SAVE :: MOLDATA(:,:),IWORK(:)
     INTEGER, PARAMETER :: NUMOBSERV=3,MAXGAUSS=20,MAXTYPE=20
     DOUBLE PRECISION, SAVE :: ATOL,TAUMAX,T,TAUI,ENOT,RCUT,RCUTSQ,VCURR,BL,BL2
     DOUBLE PRECISION, PARAMETER :: ATOMICMASS=0.020614788876D0,PI=3.1415926D0,BMAX=5.0D0,RTOL=0.0D0
     DOUBLE PRECISION, PARAMETER :: ANGTOCM=1.0D-24,NAVOG=6.02214179D23,PCONV=0.00733707D0
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: MASSARRY(:),GMAT(:),Y(:),RWORK(:),ZTOT(:),ZCURR(:),ZNEW(:)
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: BETA(:),BETALEFT(:),BETARGHT(:),VDATA(:)
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: PRESSURE(:),PRESLEFT(:),PRESRGHT(:),PARRAY(:,:) 
     DOUBLE PRECISION, SAVE :: ROECURR,ROENEW,GDATA(2*MAXGAUSS+1,MAXTYPE,MAXTYPE,2),E_ZERO,VMCONV
     CHARACTER*100 :: CNST,FNAME,VFNAME,QFNAME,VOLFNAME
     LOGICAL :: PBC,UPDATED,EXTENDED,PGRID,NRECSNT,PTINIT
     EXTERNAL DLSODE
     EXTERNAL JAC
END MODULE VGWM
  
SUBROUTINE CLEANUP_VGWM
  USE VGWM
  DEALLOCATE(MASSARRY,MOLDATA,GMAT,Y,IWORK,RWORK,ZTOT,ZCURR,ZNEW,BETA)
  DEALLOCATE(BETALEFT,BETARGHT,VDATA,PARRAY,PRESSURE,PRESLEFT,PRESRGHT)
  CLOSE(FILEIND)
  CLOSE(VFILEIND)
  CLOSE(QFILEIND)
  CLOSE(VOLFIND)
END SUBROUTINE

SUBROUTINE JAC
END SUBROUTINE


SUBROUTINE INITIALIZE_VGWM(N,TMAX,TMIN,PMIN,PMAX,TI,AT,BLX,NP,GSIZE,NODE,RCT,ENT,EG)

!**************************************************************************************
!               
!    Reads in program parameters from files 'gauss.dat' and vgwdata, initializes
!    arrays and variables.
!
!    Variables:
! 
!       N : Number of atoms
!    TMAX : Maximum temperature
!    TMIN : Minimum temperature (quenching temperature)
!    PMIN : Minimum pressure (if > -1, then constant temperature simulation)
!    PMAX : Maximum pressure 
!      TI : Initital temperature for quenching (quenching starts at TI)
!      AT : Absolute error tolerance for dlsode()
!     BLX : Maximum cube length
!      NP : Number of parallel replicas (nodes)
!   GSIZE : Number of temperautre grid points between each node
!    NODE : The node (process) number on which the code is run
!     RCT : The radius at which the pair potential is truncated.
!     ENT : Energy scaling value
!      EG : If EG = 1 Use equidistant temperature grid, otherwise use geometric grid 
!
!**************************************************************************************
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K,N,NSETS,ITYPE,MSIZE1,MSIZE2,NG,CP,CNT,ACNT,NAM,IND,NM,LSZE,NP,GSIZE,PCNT,NODE,EG
  DOUBLE PRECISION :: TMAX,TMIN,BLX,COE,AT,ALPHA,TI,TARRAY(NP),ENT,PMIN,PMAX,RCT
  CHARACTER*100 :: DUMMY

  N_ATOM=N                                                  ! Number of atoms
  N3ATOM=3*N_ATOM
  ATOL=AT                                                   ! Absolute error tolerance used by dlsode for propagation of VGW
  TAUI=MIN(TI,(0.5D0/TMIN)-TI)                              ! Initial inverse temperature
  TAUMAX=1.0D0/TMIN                                         ! Maximum inverse temperature
  NREPS=NP                                                  ! Number of parallel replicas
  GRIDSIZE=GSIZE                                            ! Number of gridpoints on temperature grid
  NRECSNT=.TRUE.                                            ! Message not received flag for parallel tempering (MPI)
  PTINIT=.TRUE.                                             ! Parallel tempering initialization
  RCUT=RCT                                                  ! Cut off radius   
  RCUTSQ=RCUT**2                                            
  NPAIRS=0.5D0*N_ATOM*(N_ATOM-1)                            ! Maximum number of pairs
  PCNT=0                                                    ! Number of atom pairs (for use with RCUT, the potential cutoff)
  mynode=NODE                                               ! node (0 = lowest - the head/master node)
  VMCONV=(ANGTOCM*NAVOG)/N_ATOM                             ! Molar volume conversion - NAVOG=Avogadros, ANGTOCM=Angstroms^3 to cm^3
  ENOT=ENT                                                  ! E_not for scaling energy
 
  IF(BLX.GT.0) THEN                                         ! Periodic boundary conditions ?
     PBC=.TRUE.                
     BL=BLX                                                 ! Box length
     BL2=BL/2.0D0
  ELSE
     PBC=.FALSE.
     BL=-1.0D0
  ENDIF
  
  ALLOCATE(MASSARRY(N3ATOM),MOLDATA(4,N_ATOM),PARRAY(2,NPAIRS))     ! Data matrix for observables, see VGWFUPDATE
  
  DO I=1,N_ATOM-1                                                   ! Inititalize pair array
     DO J=I+1,N_ATOM
        PCNT=PCNT+1
        PARRAY(1,PCNT)=I
        PARRAY(2,PCNT)=J
     ENDDO
  ENDDO
  
  OPEN(UNIT=7,FILE='vgwdata',STATUS='OLD')                                ! Read molecule data (mass, molecule, atom type)
    DO I=1,N_ATOM
       READ(7,*) MASSARRY(3*(I-1)+1),MOLDATA(1,I),MOLDATA(2,I),DUMMY
       MASSARRY(3*(I-1)+1)=1/(ATOMICMASS*MASSARRY(3*(I-1)+1))
       MASSARRY(3*(I-1)+2)=MASSARRY(3*(I-1)+1)
       MASSARRY(3*(I-1)+3)=MASSARRY(3*(I-1)+1)
    ENDDO
  CLOSE(7)

  OPEN(UNIT=7,FILE='gauss.dat',STATUS='OLD')

  READ(7,*) NSETS

  DO I=1,NSETS                                     ! Read Gaussian parameters from file into data array
     READ(7,*) ITYPE,MOLID1,MOLID2,NG              ! NG = number of gaussians
     GDATA(1,MOLID1,MOLID2,ITYPE)=NG               ! 4-Dimensional array to hold Gaussian coefficients and exponents 
     GDATA(1,MOLID2,MOLID1,ITYPE)=NG               ! 1 = intra molecular (full Gaussian), 2 = inter-molecular (single Guassian)
     DO J=1,NG
        READ(7,*) COE, ALPHA                       ! Coefficient and exponent for sum of Gaussians equation in paper
        GDATA(2*J,MOLID1,MOLID2,ITYPE)=COE
        GDATA(2*J+1,MOLID1,MOLID2,ITYPE)=ALPHA
        IF(MOLID1.NE.MOLID2) THEN
           GDATA(2*J,MOLID2,MOLID1,ITYPE)=COE
           GDATA(2*J+1,MOLID2,MOLID1,ITYPE)=ALPHA
        ENDIF
     ENDDO
  ENDDO  

  READ(7,*) CNST
  READ(7,*) E_ZERO

  CLOSE(7)

  NAM=0                                            ! Number of atoms in molecule (block size)
  NEQ=0                                            ! Total number of equations for dlsode to propagate
  IND=1
  LRGMSZE=1                                        ! Keeps track of largest molecular size
  N_MOL=0

  MOLID1=MOLDATA(1,1)

  DO I=1,N_ATOM                                    ! Determine size of G matrix
     MOLID1=MOLDATA(1,I)                           ! Molecule ID (mol. number)
     IF(MOLID1.EQ.MOLID2) THEN
        NAM=NAM+1
     ELSE
        N_MOL=N_MOL+1
        NEQ=NEQ+3*NAM*(3*NAM+1)/2
        DO J=1,NAM
           MOLDATA(1,IND)=N_MOL                    ! Molecule number
           MOLDATA(3,IND)=3*(J-1)                  ! Block offset for atom in molecule
           MOLDATA(4,IND)=NAM                      ! Molecule Size
           IF(NAM.GT.LRGMSZE) LRGMSZE=NAM          ! Largest molecule?
           IND=IND+1
        ENDDO
        NAM=1
     ENDIF
     MOLID1=MOLID2
  ENDDO
  
  NEQ=NEQ+3*NAM*(3*NAM+1)/2
  N_MOL=N_MOL+1

  DO J=1,NAM
     MOLDATA(1,IND)=N_MOL
     MOLDATA(3,IND)=3*(J-1)
     MOLDATA(4,IND)=NAM
     IF(NAM.GT.LRGMSZE) LRGMSZE=NAM
     IND=IND+1
  ENDDO
 
  ALLOCATE(GMAT(NEQ))                                     ! G matrix
  NEQ=NEQ+N3ATOM+1                                        ! Number of equations plus gamma plus coords

  IARRAY(5)=20+16*NEQ                                     ! LRW variable (dlsode rwork array size)
  IARRAY(6)=20                                            ! LIW variable (dlsode iwork array size)

  ALLOCATE(Y(NEQ),RWORK(IARRAY(5)),IWORK(IARRAY(6)))
  
  CNT=1
  IND=1
  ACNT=1

  DO I=1,N_MOL                                 ! Initialize G matrix
     NAM=MOLDATA(4,IND)
     IND=IND+NAM
     DO J=1,3*NAM
        GMAT(CNT)=MASSARRY(ACNT)
        ACNT=ACNT+1
        CNT=CNT+1
        DO K=J+1,3*NAM
           GMAT(CNT)=0.0D0
           CNT=CNT+1
        ENDDO
     ENDDO
  ENDDO

  CALL VGWGRIDINIT(TMAX,TMIN,PMIN,PMAX,EG)

  NM=N_MOL                           ! Return number of molecules
  LSZE=LRGMSZE                       ! Return largest molecule size
 
END SUBROUTINE INITIALIZE_VGWM


SUBROUTINE VGWFSETUP(Q,CONT)

!******************************************************************
!
!   Set-up of output files and Monte-Carlo simulation
!
!   Variables:
!            
!      Q : The atomic coordinates
!   CONT : If CONT = 1 continue the previous run
!
!*******************************************************************
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I,NODE,CONT,NLINES,LBLK,PBLK,LLINE,PLINE,EOF
  DOUBLE PRECISION :: Q(3,N_ATOM)
  CHARACTER*7 :: BLK
  CHARACTER*10 :: TMP
  CHARACTER*13 :: BUFF
 
  NODE=mynode+1
  FILEIND=100+NODE
  VFILEIND=800+NODE
  VOLFIND=900+NODE
  NLINES=0
  EOF=0
  PBLK=0
  LLINE=0

  CALL SYSTEM("mkdir vmdcoords")
  CALL SYSTEM("mkdir latestcoords")
  CALL SYSTEM("mkdir vcurr")

  IF(NODE.LT.10) THEN
     WRITE(FNAME,"(A6,I1)") "state.",NODE
     WRITE(VFNAME,"(A21,I1,A4)")"./vmdcoords/vmdcoords", NODE, ".xyz"
     WRITE(QFNAME,"(A23,I1)") "./latestcoords/lcoords.", NODE
     WRITE(VOLFNAME,"(A16,I1)") "./vcurr/currvol.", NODE
  ELSE
     WRITE(FNAME,"(A6,I2)") "state.",NODE
     WRITE(VFNAME,"(A21,I2,A4)")"./vmdcoords/vmdcoords", NODE, ".xyz"
     WRITE(QFNAME,"(A23,I2)") "./latestcoords/lcoords.", NODE
     WRITE(VOLFNAME,"(A16,I2)") "./vcurr/currvol.", NODE
  ENDIF
  
  IF(CONT.EQ.1) THEN
     NEQUIL=0
     OPEN(UNIT=FILEIND,FILE=FNAME,STATUS='OLD')
     OPEN(UNIT=QFILEIND,FILE=QFNAME,STATUS='OLD')
     OPEN(UNIT=VOLFIND,FILE=VOLFNAME,STATUS='OLD')
     DO WHILE(.NOT.EOF)
        READ(FILEIND,*,IOSTAT=EOF) BLK,TMP
        IF(.NOT.EOF) THEN
           NLINES=NLINES+1
        ENDIF
        IF(BLK.EQ.'<BLOCK>') THEN
           PBLK=LBLK
           PLINE=LLINE
           BLCK=PBLK
           READ(TMP,*) LBLK
           LLINE=NLINES
        ENDIF
     ENDDO

     IF(LBLK.GT.1) THEN
        REWIND(FILEIND)
        DO I=1,PLINE
           READ(FILEIND,*)
        ENDDO
        DO I=1,N_ATOM
           READ(QFILEIND,*) Q(1,I),Q(2,I),Q(3,I)
        ENDDO
        READ(FILEIND,*) BUFF, NTOTAL
        READ(FILEIND,*) BUFF, NTACCP
        READ(FILEIND,*) BUFF, NTRANS
        READ(FILEIND,*)
        READ(FILEIND,*) BUFF, STPLEN
        READ(FILEIND,*) BUFF, BL
        READ(FILEIND,*) BUFF, VSLEN, NVACCT, NVSTEP
        READ(FILEIND,*) BUFF, ENOT
        READ(FILEIND,*)
        READ(FILEIND,*) BUFF, NSWAPS(1), NSWAPS(2)
        READ(FILEIND,*) BUFF, NSWPATT(1), NSWPATT(2)       
     ELSE
        WRITE(*,*) "Only one block. Did not read previous data."
     ENDIF

     NSTEP=NSTEP+NTOTAL

     BL2=BL/2.0D0
     
     REWIND(FILEIND)
     
     DO I=1,LLINE-1
        READ(FILEIND,*)
     ENDDO
  ELSE
     OPEN(UNIT=FILEIND,FILE=FNAME)
     OPEN(UNIT=QFILEIND,FILE=QFNAME)
     OPEN(UNIT=VOLFIND,FILE=VOLFNAME)

     ! IF(ENOT.GT.0) CALL VGWZEROENRG(Q)
  ENDIF

  IF((VSLEN.GT.0).AND.(BL.GT.0).AND.(VSFREQ.GT.0)) CHNGV=.TRUE.
 
  VCURR=BL**3

  IF(BL.LT.0) THEN
     BL=-1.0D0
     VCURR=1.0D0
  ENDIF

  IF(mynode.eq.0) THEN
     OPEN(UNIT=97,FILE='initialdata')
     WRITE(97,*) "Initial Energy: ", ENOT
     DO I=1,N_ATOM
        WRITE(97,*) Q(1,I),Q(2,I),Q(3,I)
     ENDDO
     CLOSE(97)
  ENDIF

  IF(PBC) Q=Q/BL

  OPEN(UNIT=VFILEIND,FILE=VFNAME)

END SUBROUTINE VGWFSETUP
    
SUBROUTINE VGWGRIDINIT(TMAX,TMIN,PMIN,PMAX,EG)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,NODE,CNT,EG
  DOUBLE PRECISION :: TMAX,TMIN,PMIN,PMAX,CTE,T1,T2,TSTEP,RTEMP(NREPS+2),PRNDE(NREPS+2)

  NODE=mynode+1
  ASIZE=GRIDSIZE+1
  PTSIZE=ASIZE+N3ATOM+1

  ALLOCATE(BETA(ASIZE),BETALEFT(ASIZE),BETARGHT(ASIZE),ZTOT(ASIZE),ZCURR(ASIZE),ZNEW(ASIZE))
  ALLOCATE(VDATA(ASIZE),PRESSURE(ASIZE),PRESLEFT(ASIZE),PRESRGHT(ASIZE)) 
 
 IF(PMIN.GT.0) THEN
    PGRID=.TRUE.
    BETA=0.5D0/TMIN
 !*****************************************************
 !
 ! Insert code here to initialize pressure grid
 !
 !*****************************************************
 ELSE
    IF (TMIN < 0.000001D0) TMIN=0.000001D0                   ! to avoid devision by zero
  
    IF(EG.EQ.1) THEN      
       TSTEP=(TMAX-TMIN)/NREPS
       DO I=0,NREPS+1
          RTEMP(I+1)=TMIN+I*TSTEP
       ENDDO
    ELSE     
       CTE=(LOG(TMAX/TMIN))/(NREPS-1)
       CTE=EXP(CTE)
       DO I=0,NREPS+1
          RTEMP(I+1)=TMIN*CTE**I
       ENDDO
    ENDIF
       
    PGRID=.FALSE.

    T1=1.0D0/RTEMP(NODE+2)
    T2=1.0D0/RTEMP(NODE+1)
    TSTEP=(T2-T1)/GRIDSIZE
    
    DO I=0,GRIDSIZE
       BETARGHT(I+1)=T1+I*TSTEP
    ENDDO
    
    T1=T2
    T2=1.0D0/RTEMP(NODE)
    TSTEP=(T2-T1)/GRIDSIZE
    
    DO I=0,GRIDSIZE
       BETA(I+1)=T1+I*TSTEP
    ENDDO
    
    IF(NODE.NE.1) THEN
       T1=1.0D0/RTEMP(NODE)
       T2=1.0D0/RTEMP(NODE-1)
       TSTEP=(T2-T1)/GRIDSIZE
       DO I=0,GRIDSIZE
          BETALEFT(I+1)=T1+I*TSTEP
       ENDDO
    ENDIF
 
    ZTOT=0.0D0
    VDATA=0.0D0
    
    IF(BL.GT.0) THEN
       PRESSURE=PMAX                                                                      
       PRESLEFT=PMAX
       PRESRGHT=PMAX
    ELSE
       PRESSURE=0.0D0
       PRESRGHT=0.0D0
       PRESLEFT=0.0D0
    ENDIF
 ENDIF

 
END SUBROUTINE VGWGRIDINIT

SUBROUTINE VGWZEROENRG(Q)
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(N3ATOM),QT(N3ATOM),TAU(2),LOGP(2),UGAUSS,DTAU
  INTEGER :: I
  LOGICAL :: INSIDE
  EXTERNAL RHSM
 
  Y=0.0D0
  DTAU=0.0001D0

  CALL INITDLSODE

  IF(PBC) THEN
     QT=BL*Q
  ELSE
     QT=Q
  ENDIF

  IF(RCUT.GT.0) CALL PAIRS(QT)
  
  DO I=1,N3ATOM
     Y(I+1)=QT(I)                           ! COPY INITITAL COORDS TO Y VECTOR
  ENDDO
   
  DO I=1,NEQ-N3ATOM-1
     Y(N3ATOM+1+I)=GMAT(I)*TAUI
  ENDDO
  
  CALL GAUSSENERGY(Y(2),UGAUSS)
  Y(1)=-TAUI*UGAUSS
  
  T=TAUI
  TAU(1)=0.5D0*BMAX-DTAU
  TAU(2)=0.5D0*BMAX
  
  DO I=1,2
     CALL DLSODE(RHSM,NEQ,Y,T,TAU(I),IARRAY(1),RTOL,ATOL,IARRAY(2),IARRAY(4),IARRAY(3),RWORK,IARRAY(5),IWORK,IARRAY(6),JAC,IARRAY(7))
     CALL LNPS(LOGP(I))
  ENDDO
  
  ENOT=-(LOGP(2)-LOGP(1))/(2.0D0*DTAU)
  
  IF(ENOT.GT.0) WRITE(*,*) "Warning, energy of initial configuration (ENOT) is greater than 0."

END SUBROUTINE VGWZEROENRG

SUBROUTINE VGWMQUENCH(Q)
  USE VGWM
  IMPLICIT NONE
  DOUBLE PRECISION :: ROETMP,Q(N3ATOM),QT(N3ATOM),UGAUSS,BETAT
  INTEGER :: N_STEP,I,J,K,CNT

  EXTERNAL RHSM
 
  Y=0.0D0
  
  CALL INITDLSODE

  IF(PBC) THEN
     QT=BL*Q
  ELSE
     QT=Q
  ENDIF

  IF(RCUT.GT.0) CALL PAIRS(QT)

  DO I=1,N3ATOM
     Y(I+1)=QT(I)                           ! COPY INITITAL COORDS TO Y VECTOR
  ENDDO
  
  DO I=1,NEQ-N3ATOM-1
     Y(N3ATOM+1+I)=GMAT(I)*TAUI
  ENDDO

  CALL GAUSSENERGY(Y(2),UGAUSS)
  Y(1)=-TAUI*UGAUSS

  T=TAUI

  DO I=1,ASIZE
     CALL DLSODE(RHSM,NEQ,Y,T,0.5D0*BETA(I),IARRAY(1),RTOL,ATOL,IARRAY(2),IARRAY(4),IARRAY(3),RWORK,IARRAY(5),IWORK,IARRAY(6),JAC,IARRAY(7))
     CALL LNPS(ZNEW(I))
     ZNEW(I)=ZNEW(I)+BETA(I)*(ENOT-PRESSURE(I)*VCURR)!+N_ATOM*LOG(VCURR)
  ENDDO

  ROENEW=ZNEW(ASIZE)
 
END SUBROUTINE VGWMQUENCH
  
SUBROUTINE VGWCQNCH(Q,LNP,UL)
  USE VGWM
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(N3ATOM),QT(N3ATOM),LNP(ASIZE+2),UGAUSS
  INTEGER :: I,UL
  EXTERNAL RHSM

  Y=0.0D0

  CALL INITDLSODE
    
  IF(PBC) THEN
     QT=BL*Q
  ELSE
     QT=Q
  ENDIF

  IF(RCUT.GT.0) CALL PAIRS(QT)

  DO I=1,N3ATOM
     Y(I+1)=QT(I)                           ! COPY INITITAL COORDS TO Y VECTOR
  ENDDO

  DO I=1,NEQ-N3ATOM-1
     Y(N3ATOM+1+I)=GMAT(I)*TAUI
  ENDDO

  CALL GAUSSENERGY(Y(2),UGAUSS)
  Y(1)=-TAUI*UGAUSS

  T=TAUI
  
  IF(UL.EQ.1) THEN                                        ! Quench to lower temperature (left node)
     DO I=1,ASIZE
        CALL DLSODE(RHSM,NEQ,Y,T,0.5D0*BETALEFT(I),IARRAY(1),RTOL,ATOL,IARRAY(2),IARRAY(4),IARRAY(3),RWORK,IARRAY(5),IWORK,IARRAY(6),JAC,IARRAY(7)) 
        CALL LNPS(LNP(I))
        LNP(I)=LNP(I)+BETALEFT(I)*(ENOT-PRESLEFT(I)*VCURR)!+N_ATOM*LOG(VCURR)
     ENDDO
  ELSE                                                    ! Quench to higher temperature (right node)
      DO I=1,ASIZE                                   
        CALL DLSODE(RHSM,NEQ,Y,T,0.5D0*BETARGHT(I),IARRAY(1),RTOL,ATOL,IARRAY(2),IARRAY(4),IARRAY(3),RWORK,IARRAY(5),IWORK,IARRAY(6),JAC,IARRAY(7)) 
        CALL LNPS(LNP(I))
        LNP(I)=LNP(I)+BETARGHT(I)*(ENOT-PRESRGHT(I)*VCURR)!+N_ATOM*LOG(VCURR)
     ENDDO
  ENDIF
 
  LNP(ASIZE+1)=ROECURR

END SUBROUTINE VGWCQNCH

SUBROUTINE PAIRS(Q)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K
  DOUBLE PRECISION :: Q(3,N_ATOM),RSQ,QIJ(3)
  
  NPAIRS=0
  
  DO I=1,N_ATOM-1
     DO J=I+1,N_ATOM
        RSQ=0.0D0
        DO K=1,3
           QIJ(K)=Q(K,I)-Q(K,J)
           IF(PBC) THEN                                            ! Periodic boundary conditions
              IF(QIJ(K).GT.BL2) QIJ(K)=QIJ(K)-BL
              IF(QIJ(K).LT.-BL2) QIJ(K)=QIJ(K)+BL
           ENDIF
           RSQ=RSQ+QIJ(K)**2
        ENDDO
        IF(RSQ.LE.RCUTSQ) THEN
           NPAIRS=NPAIRS+1
           PARRAY(1,NPAIRS)=I
           PARRAY(2,NPAIRS)=J
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE PAIRS
  
SUBROUTINE INITDLSODE
  USE VGWM
  IMPLICIT NONE

  IWORK=0.0D0
  RWORK=0.0D0

  IARRAY(1)=1                    ! ITOL
  IARRAY(2)=1                    ! ITASK
  IARRAY(3)=1                    ! IOPT
  IARRAY(4)=1                    ! ISTATE
  IARRAY(5)=20+16*NEQ            ! LRW
  IARRAY(6)=20                   ! LIW
  IARRAY(7)=10                   ! MF

  IWORK(6)=100000
  IWORK(5)=4
  IWORK(7)=0
  IWORK(8)=0
  IWORK(9)=0
  IWORK(10)=0
  RWORK(5)=0.0D0
  RWORK(6)=0.0D0
  RWORK(7)=0.0D0
  RWORK(8)=0.0D0
  RWORK(9)=0.0D0
  RWORK(10)=0.0D0

END SUBROUTINE

SUBROUTINE LNPS(LOGZ)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K,CNT,DIM,INFO,IND
  DOUBLE PRECISION :: C(3*LRGMSZE,3*LRGMSZE),LOGZ,GAMMA,DET
  DOUBLE PRECISION :: DETL,WORK((3*LRGMSZE)**2)
  CHARACTER*1 :: UPLO="U"
 
  GAMMA=Y(1)
  CNT=2+N3ATOM
  C=0.0D0
  IND=1
  DET=0.0D0

  DO I=1,N_MOL
     DIM=3*MOLDATA(4,IND)
     IND=IND+MOLDATA(4,IND)
     DO J=1,DIM
        DO K=J,DIM
           C(J,K)=Y(CNT)
           CNT=CNT+1
        ENDDO
     ENDDO
    
     DETL=1.0D0
     CALL DPOTRF(UPLO,DIM,C,DIM,INFO)

     DO J=1,DIM
        DETL=DETL*C(J,J)
     ENDDO
    
     DET=DET+2.0D0*LOG(DETL)

  ENDDO

  LOGZ=2.0D0*GAMMA-0.5D0*DET

END SUBROUTINE 

SUBROUTINE INVDET(A,M,DET)
  IMPLICIT NONE
  DOUBLE PRECISION :: DET,A(3,3),M(3,3)

  M(1,1) = A(2,2)*A(3,3)-A(2,3)**2
  M(2,2) = A(1,1)*A(3,3)-A(1,3)**2
  M(3,3) = A(1,1)*A(2,2)-A(1,2)**2
  M(1,2) = -A(1,2)*A(3,3)+A(1,3)*A(2,3)
  M(1,3) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
  M(2,3) = -A(1,1)*A(2,3)+A(1,3)*A(1,2)
  DET = M(1,1)*A(1,1)+M(1,2)*A(1,2)+M(1,3)*A(1,3)

END SUBROUTINE INVDET  

SUBROUTINE RHSM(NEQM,TT,YM,YPRIME)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K,PCNT,L,NEQM,CNT,NDIM,NAM,IND,MASSCNT
  DOUBLE PRECISION :: QIJ(3),U,TRUXX,GUX,GUG,GU(3*LRGMSZE,3*LRGMSZE),TT 
  DOUBLE PRECISION :: UPV(N3ATOM),UPM(3*LRGMSZE,3*LRGMSZE,N_MOL),DIM
  DOUBLE PRECISION :: YM(NEQM),YPRIME(NEQM),G(3*LRGMSZE,3*LRGMSZE,N_MOL)

  CNT=2+N3ATOM
  NDIM=3*LRGMSZE
  IND=1

  DO I=1,N_MOL
     NAM=MOLDATA(4,IND)                       ! Inititalize G matrix
     IND=IND+NAM
     DO J=1,3*NAM
        G(J,J,I)=YM(CNT)
        CNT=CNT+1
        DO K=J+1,3*NAM
           G(K,J,I)=YM(CNT)                                 
           G(J,K,I)=YM(CNT)
           CNT=CNT+1
        ENDDO
     ENDDO
  ENDDO
     
  U=0.0D0
  UPV=0.0D0
  UPM=0.0D0
  TRUXX=0.0D0
   
  DO PCNT=1,NPAIRS  
     I=PARRAY(1,PCNT)
     J=PARRAY(2,PCNT)
     DO K=1,3
        QIJ(K)=YM(3*(I-1)+K+1)-YM(3*(J-1)+K+1)
        IF(PBC) THEN                                            ! Periodic boundary conditions
           IF(QIJ(K).GT.BL2) QIJ(K)=QIJ(K)-BL
           IF(QIJ(K).LT.-BL2) QIJ(K)=QIJ(K)+BL
        ENDIF
     ENDDO
     IF(MOLDATA(1,I).EQ.MOLDATA(1,J)) THEN                                             ! Same molecule ?
        CALL RHSFC(I,J,QIJ,U,UPV,G,UPM,MOLDATA(1,I),MOLDATA(2,I),MOLDATA(2,J))
     ELSE
        CALL RHSSP(I,J,QIJ,U,UPV,G,UPM,MOLDATA(1,I),MOLDATA(1,J),MOLDATA(2,I),MOLDATA(2,J))     ! Different molecule.
     ENDIF
  ENDDO

  CALL FILLBLOCKS(UPM)                  ! Fill lower diagonal of G matrix

  DO I=1,N_MOL
     DO J=1,NDIM
        TRUXX=TRUXX+UPM(J,J,I)*G(J,J,I)
        DO K=J+1,NDIM
           TRUXX=TRUXX+2.0D0*UPM(J,K,I)*G(K,J,I)
        ENDDO
     ENDDO
  ENDDO

  YPRIME(1)=-0.25D0*TRUXX-(U+E_ZERO)
  CNT=2
  IND=1

  DO I=1,N_MOL
     NDIM=3*MOLDATA(4,IND)
     DO J=1,NDIM
        GUX=0.0D0
        DO K=1,NDIM
           GUX=GUX-G(J,K,I)*UPV(3*(IND-1)+K)
        ENDDO
        YPRIME(CNT)=GUX
        CNT=CNT+1
     ENDDO
     IND=IND+MOLDATA(4,IND)
  ENDDO

  IND=1
  MASSCNT=1

  DO I=1,N_MOL
     GU=MATMUL(UPM(:,:,I),G(:,:,I))
     NDIM=3*MOLDATA(4,IND)
     IND=IND+MOLDATA(4,IND)
     DO J=1,NDIM
        DO K=J,NDIM
           GUG=0.0D0
           DO L=1,NDIM
              GUG=GUG-G(J,L,I)*GU(L,K)
           ENDDO
           IF(J.EQ.K) THEN
              GUG=GUG+MASSARRY(MASSCNT)
              MASSCNT=MASSCNT+1
           ENDIF
           YPRIME(CNT)=GUG
           CNT=CNT+1
        ENDDO
     ENDDO
  ENDDO    
  
END SUBROUTINE

SUBROUTINE RHSFC(IIND,JIND,QIJ,U,UPV,G,UPM,MNUM,TYPEI,TYPEJ)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K,IIND,JIND,IOF,JOF,MNUM,NG,TYPEI,TYPEJ
  DOUBLE PRECISION :: QIJ(3),UPV(3,N_ATOM),G(3*LRGMSZE,3*LRGMSZE,N_MOL),DETA,DETAG,ZQ(3),U,QZQ
  DOUBLE PRECISION :: UPM(3*LRGMSZE,3*LRGMSZE,N_MOL),A(3,3),M(3,3),AG(3,3),Z(3,3),EXPF,UX,UXX,AK,CK

  IOF=MOLDATA(3,IIND)
  JOF=MOLDATA(3,JIND)
 
  DO I=1,3
     DO J=I,3
        A(I,J)=G(IOF+I,IOF+J,MNUM)+G(JOF+I,JOF+J,MNUM)-G(IOF+I,JOF+J,MNUM)-G(IOF+J,JOF+I,MNUM)
     ENDDO
  ENDDO
        
  CALL INVDET(A,M,DETA)

  DETA=1.0D0/DETA
  A=M*DETA

  NG=INT(GDATA(1,TYPEI,TYPEJ,1))

  DO I=1,NG                                     ! Sum over intra-molecular Gaussians
     AG=A
     AK=GDATA(2*I+1,TYPEI,TYPEJ,1)
     CK=GDATA(2*I,TYPEI,TYPEJ,1)

     DO J=1,3
        AG(J,J)=AG(J,J)+AK
     ENDDO

     CALL INVDET(AG,M,DETAG)

     Z=-(AK**2/DETAG)*M

     DO J=1,3
        Z(J,J)=Z(J,J)+AK
     ENDDO

     Z(3,1)=Z(1,3)
     Z(2,1)=Z(1,2)
     Z(3,2)=Z(2,3)

     ZQ=MATMUL(Z,QIJ)
     QZQ=DOT_PRODUCT(QIJ,ZQ)
     ZQ=2.0D0*ZQ
     
     EXPF=SQRT(DETA/DETAG)*EXP(-QZQ)*CK
     U=U+EXPF

     DO J=1,3
        UX=-EXPF*ZQ(J)
        UPV(J,IIND)=UPV(J,IIND)+UX
        UPV(J,JIND)=UPV(J,JIND)-UX
        DO K=J,3
           UXX=EXPF*(ZQ(J)*ZQ(K)-2.0D0*Z(J,K))
           UPM(IOF+J,IOF+K,MNUM)=UPM(IOF+J,IOF+K,MNUM)+UXX             ! Fill upper triangle of blocks
           UPM(JOF+J,JOF+K,MNUM)=UPM(JOF+J,JOF+K,MNUM)+UXX
           UPM(IOF+J,JOF+K,MNUM)=UPM(IOF+J,JOF+K,MNUM)-UXX             ! Fill off diagonal blocks
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE RHSSP(IIND,JIND,QIJ,U,UPV,G,UPM,MNUM1,MNUM2,TYPEI,TYPEJ)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K,IIND,JIND,IOF,JOF,MNUM1,MNUM2,NG,TYPEI,TYPEJ
  DOUBLE PRECISION :: QIJ(3),UPV(3,N_ATOM),G(3*LRGMSZE,3*LRGMSZE,N_MOL),DETA,DETAG,ZQ(3),U,QZQ
  DOUBLE PRECISION :: UPM(3*LRGMSZE,3*LRGMSZE,N_MOL),A(3,3),M(3,3),AG(3,3),Z(3,3),EXPF,UX,UXX,AK,CK

  IOF=MOLDATA(3,IIND)
  JOF=MOLDATA(3,JIND)

  DO I=1,3
     DO J=I,3
        A(I,J)=G(IOF+I,IOF+J,MNUM1)+G(JOF+I,JOF+J,MNUM2)
     ENDDO
  ENDDO
        
  CALL INVDET(A,M,DETA)

  DETA=1.0D0/DETA
  A=M*DETA

  NG=INT(GDATA(1,TYPEI,TYPEJ,2))

  DO I=1,NG                                     ! Sum over intra-molecular Gaussians
     AG=A
     AK=GDATA(2*I+1,TYPEI,TYPEJ,2)
     CK=GDATA(2*I,TYPEI,TYPEJ,2)
     DO J=1,3
        AG(J,J)=AG(J,J)+AK
     ENDDO

     CALL INVDET(AG,M,DETAG)

     Z=-(AK**2/DETAG)*M

     DO J=1,3
        Z(J,J)=Z(J,J)+AK
     ENDDO

     Z(3,1)=Z(1,3)
     Z(2,1)=Z(1,2)
     Z(3,2)=Z(2,3)

     ZQ=MATMUL(Z,QIJ)
     QZQ=DOT_PRODUCT(QIJ,ZQ)
     ZQ=2.0D0*ZQ
     
     EXPF=SQRT(DETA/DETAG)*EXP(-QZQ)*CK
     U=U+EXPF

     DO J=1,3
        UX=-EXPF*ZQ(J)
        UPV(J,IIND)=UPV(J,IIND)+UX
        UPV(J,JIND)=UPV(J,JIND)-UX
        DO K=J,3
           UXX=EXPF*(ZQ(J)*ZQ(K)-2.0D0*Z(J,K))
           UPM(IOF+J,IOF+K,MNUM1)=UPM(IOF+J,IOF+K,MNUM1)+UXX             ! Fill upper triangle of blocks
           UPM(JOF+J,JOF+K,MNUM2)=UPM(JOF+J,JOF+K,MNUM2)+UXX
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE FILLBLOCKS(UPM)
  USE VGWM
  IMPLICIT NONE
  INTEGER :: I,J,K,MSZE
  DOUBLE PRECISION :: UPM(3*LRGMSZE,3*LRGMSZE,N_MOL)

  MSZE=3*LRGMSZE

  DO I=1,N_MOL
     DO J=1,LRGMSZE
        DO K=J+1,LRGMSZE
           UPM(3*(J-1)+2,3*(K-1)+1,I)=UPM(3*(J-1)+1,3*(K-1)+2,I)           ! Fill off diagonal blocks
           UPM(3*(J-1)+3,3*(K-1)+1,I)=UPM(3*(J-1)+1,3*(K-1)+3,I)
           UPM(3*(J-1)+3,3*(K-1)+2,I)=UPM(3*(J-1)+2,3*(K-1)+3,I)
        ENDDO
     ENDDO

     DO J=1,MSZE                  ! Fill lower triangle of molecular blocks       
        DO K=J+1,MSZE
           UPM(K,J,I)=UPM(J,K,I)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE GAUSSENERGY(Q,UGAUSS)
  USE VGWM
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(3,N_ATOM),UGAUSS,RSQ,AK,CK,QIJ(3)
  INTEGER :: I,J,K,PCNT,NG,TYPEI,TYPEJ

  UGAUSS=0.0D0

  DO PCNT=1,NPAIRS
     I=PARRAY(1,PCNT)
     J=PARRAY(2,PCNT)
     RSQ=0.0D0
     DO K=1,3
        QIJ(K)=Q(K,I)-Q(K,J)
        IF(PBC) THEN                                            ! Periodic boundary conditions
           IF(QIJ(K).GT.BL2) QIJ(K)=QIJ(K)-BL
           IF(QIJ(K).LT.-BL2) QIJ(K)=QIJ(K)+BL
        ENDIF
        RSQ=RSQ+QIJ(K)**2
     ENDDO
     TYPEI=MOLDATA(2,I)
     TYPEJ=MOLDATA(2,J)
      IF(MOLDATA(1,I).EQ.MOLDATA(1,J)) THEN
         NG=INT(GDATA(1,TYPEI,TYPEJ,1))
         DO K=1,NG
            AK=GDATA(2*K+1,TYPEI,TYPEJ,1)
            CK=GDATA(2*K,TYPEI,TYPEJ,1)
            UGAUSS=UGAUSS+CK*EXP(-AK*RSQ)
         ENDDO
      ELSE
         NG=INT(GDATA(1,TYPEI,TYPEJ,2))
         DO K=1,NG
            AK=GDATA(2*K+1,TYPEI,TYPEJ,2)
            CK=GDATA(2*K,TYPEI,TYPEJ,2)
            UGAUSS=UGAUSS+CK*EXP(-AK*RSQ)
         ENDDO
      ENDIF
   ENDDO
   
   UGAUSS=UGAUSS+E_ZERO

 END SUBROUTINE GAUSSENERGY

