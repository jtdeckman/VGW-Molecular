MODULE VGWMC
  IMPLICIT NONE
  INTEGER, SAVE :: NUPDATES,UFREQ,NEQUIL,NSTEP,NTOTAL,NACCPT,SWPINT,BLCK,VSFREQ,WFREQ,NVOUT
  INTEGER, SAVE :: NTRANS,NTACCP,NACC,NVSTEP,NVACC,NVACCT,MCUF,CSTEP(2),NSWPATT(2),NSWAPS(2)
  DOUBLE PRECISION, SAVE :: STPLEN,VSLEN,RCON,BLMAX
  LOGICAL, SAVE :: CHNGV
END MODULE VGWMC

SUBROUTINE INITIALIZE_VGWMC(NSTP,UFR,WRF,MCF,SI,SLEN,VLEN,RCN,Q,CONT,N,VSF,BLMX,NVO)
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: NSTP,NE,UFR,WRF,MCF,SI,CONT,N,VSF,NVO
  DOUBLE PRECISION :: SLEN,VLEN,RCN,BLMX,Q(3,N)

  NSTEP=NSTP                    ! Number of MC steps to do
  WFREQ=WRF                     ! Observable write frequency (every n updates)
  UFREQ=UFR                     ! Observable update frequency
  MCUF=MCF                      ! Scale step lengths every MCUF MC steps
  SWPINT=SI                     ! PT swap interval
  RCON=RCN                      ! Constraining radius (for non-PBC)
  STPLEN=SLEN                   ! Step length for translational steps
  VSLEN=VLEN                    ! Step length for volume changes
  VSFREQ=VSF                    ! Attempt volume step frequency
  NVOUT=NVO                     ! Volume output frequency

  NACCPT=0                      ! Number of accepted MC steps
  NTOTAL=0                      ! Number of total steps
  NTRANS=0                      ! Number of translational steps
  NVSTEP=0                      ! Number of volume change steps
  NSWPATT=0                     ! Number of PT swap attempts
  NSWAPS=0                      ! Number of swaps made (1 = left node, 2=right node)
  CSTEP=0                       ! Current step (used for PT swapping; swap every CSTEP=SWPINT) 
  BLCK=0                        ! Current block of observable data
  NACC=0                        ! Number of accepted translational steps
  NVACC=0                       ! Number of accepted volume change steps
  NVACCT=0                      ! Total accepted volume steps
  CHNGV=.FALSE.                 ! Change volume (do volume steps)
  BLMAX=BLMX                    ! Maximum box length

  CALL VGWFSETUP(Q,CONT)

END SUBROUTINE INITIALIZE_VGWMC

SUBROUTINE MCSTEP(Q,VWRF,TAG)
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I,J,STYPE,VWRF,REQVOL
  DOUBLE PRECISION :: Q(3,N_ATOM),QVEC(3)
  LOGICAL :: INSIDE,ATEST,TFLAG,UD
  CHARACTER*2 :: TAG
  
  TFLAG=.FALSE.
 
  CALL RANDOM_SEED   
  CALL VGWMQUENCH(Q)
  CALL VGWFINIT   
  
  DO WHILE(.NOT.TFLAG) 
     
     CALL TSTEP(Q)
     
     IF((MOD(NTRANS,VSFREQ).EQ.0).AND.CHNGV) CALL VSTEP(Q)
    
     IF((MOD(NTOTAL,VWRF).EQ.0).AND.(VWRF.GT.0)) CALL VMDCOORDS(Q,TAG)
     
     CALL ATTSWAP(Q,TFLAG)
  ENDDO

END SUBROUTINE MCSTEP

SUBROUTINE TSTEP(Q)                                               ! Translational step
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I,J,IND
  DOUBLE PRECISION :: Q(3,N_ATOM),QT(3,N_ATOM),RNUM,ACCR
  LOGICAL :: INSIDE

  INSIDE=.TRUE.

  NTRANS=NTRANS+1
  IF(mynode.GT.0) CSTEP(1)=CSTEP(1)+1                             ! Parallel
  IF(mynode.LT.(NREPS-1)) CSTEP(2)=CSTEP(2)+1

  QT=Q
  
  CALL RANDOM_NUMBER(RNUM)

  IND=INT(N_ATOM*RNUM)+1                                          ! Pick atom randomly
  IF(IND.GT.N_ATOM) IND=N_ATOM

  DO J=1,3
     CALL RANDOM_NUMBER(RNUM)
     QT(J,IND)=QT(J,IND)+STPLEN*(RNUM-0.5D0)
     IF(PBC) THEN
        IF(QT(J,IND).GE.(0.5D0)) QT(J,IND)=QT(J,IND)-1.0D0        ! Periodic boundary conditions
        IF(QT(J,IND).LE.(-0.5D0)) QT(J,IND)=QT(J,IND)+1.0D0
        IF(ABS(QT(J,IND)).GT.(0.5D0)) INSIDE=.FALSE.
     ENDIF
  ENDDO
  
  IF(.NOT.PBC) CALL RCENTER(QT,N_ATOM,RCON,INSIDE)
  
  IF(INSIDE) THEN
     CALL VGWMQUENCH(QT)                                          ! Get density
     CALL RANDOM_NUMBER(RNUM)
     IF((RNUM < EXP(ROENEW-ROECURR))) THEN                        ! Metropolis
        NACC=NACC+1                                               ! Increment number of accepted steps in current interval
        NTACCP=NTACCP+1                                           ! Increment number of accepted translational steps
        ROECURR=ROENEW                                            ! Set current density to new density (roe)
        ZCURR=ZNEW                                                ! Set current Z to new Z
        Q=QT                                                      ! Set current coordinates to temporary coordinates
     ENDIF
  ENDIF
  
  CALL VGWFUPDATE(Q)

  IF(MOD(NTRANS,MCUF).EQ.0) THEN
     ACCR=FLOAT(NACC)/FLOAT(MCUF)
     NACC=0
      IF(STPLEN.GE.(0.35D0)) THEN
        STPLEN=STPLEN/1.10779652734191D0
     ELSE 
        IF(ACCR.GT.0.42D0) STPLEN=STPLEN*1.12799165273419D0
        IF(ACCR.LT.0.38D0) STPLEN=STPLEN/1.10779652734191D0
     ENDIF
  ENDIF

END SUBROUTINE TSTEP

SUBROUTINE VSTEP(Q)
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I,J,VEC
  DOUBLE PRECISION :: BLTMP,Q(3,N_ATOM),RNUM,ACCR,VTMP

  BLTMP=BL
  NVSTEP=NVSTEP+1

  CALL RANDOM_NUMBER(RNUM)

  BL=BL+VSLEN*(RNUM-0.5D0)
    
  IF(BL.LE.BLMAX) THEN
     BL2=BL/2.0D0
  
     VTMP=VCURR                            ! Store current volume
     VCURR=BL**3                           ! New volume
  
     CALL VGWMQUENCH(Q)
  
     CALL RANDOM_NUMBER(RNUM)

     IF(RNUM < EXP(ROENEW-ROECURR+N_ATOM*LOG(VCURR/VTMP)+2.0D0*LOG(BL/BLTMP))) THEN
        NVACC=NVACC+1                                                                  ! Increment # volume steps
        NVACCT=NVACCT+1                                                                ! Increment # accpeted volume steps
        ROECURR=ROENEW                                                                 ! Set current density to new density 
        ZCURR=ZNEW                                                                     ! Set current Z to new Z
     ELSE
        BL=BLTMP                                                                       ! If not accepted, reset box length and P*V
        VCURR=VTMP                                                                     ! Set current volume to new volume 
        BL2=BL/2.0D0
     ENDIF

  ELSE
     BL=BLTMP
  ENDIF
  
  CALL VGWFUPDATE(Q)
  
  IF(MOD(NVSTEP,MCUF).EQ.0) THEN
     ACCR=FLOAT(NVACC)/FLOAT(MCUF)
     NVACC=0
     IF(VSLEN.GT.(BL/3.0D0)) THEN
        VSLEN=VSLEN/1.10779652734191D0
     ELSE       
        IF(ACCR.GT.0.42D0) VSLEN=VSLEN*1.12799165273419D0
        IF(ACCR.LT.0.38D0) VSLEN=VSLEN/1.10779652734191D0
     ENDIF
  ENDIF
  
END SUBROUTINE VSTEP     

SUBROUTINE RCENTER(Q,N,RC,INSIDE)
  IMPLICIT NONE
  INTEGER :: I,J,N,RN
  DOUBLE PRECISION :: RC,Q(3,N),QC(3),RMAX,RCURR
  LOGICAL :: INSIDE

   DO I=1,N                                  ! Determine center of mass vector
      DO J=1,3
         QC(J)=QC(J)+Q(J,I)                    
      ENDDO
   ENDDO
    
   QC=QC/N
  
   DO I=1,N
      DO J=1,3
         Q(J,I)=Q(J,I)-QC(J)                 ! Zero center of mass
      ENDDO
   ENDDO

   RMAX=-1.0D0

   DO I=1,N
      RCURR=0.0D0
      DO J=1,3
         RCURR=RCURR+Q(J,I)**2
      ENDDO
      IF(RCURR.GT.RMAX) RMAX=RCURR
   ENDDO
   
   IF(RMAX.GT.RC**2) THEN
      INSIDE=.FALSE.
   ELSE
      INSIDE=.TRUE.
   ENDIF

 END SUBROUTINE RCENTER

 SUBROUTINE PBCHECK(Q,N,BL,TVAL)
   IMPLICIT NONE
   INTEGER :: I,J,N
   DOUBLE PRECISION :: Q(3,N),BL,BL2
   LOGICAL :: TVAL
   
   TVAL=.TRUE.

   BL2=BL/2.0D0
 
   DO I=1,N
      DO J=1,3
         IF(Q(J,I).GE.BL2) Q(J,I) = Q(J,I) - BL
         IF(Q(J,I).LE.-BL2) Q(J,I) = Q(J,I) + BL
         IF(ABS(Q(J,I)).GT.BL2) TVAL=.FALSE.                           ! If its still outside the box then reject config
      ENDDO
   ENDDO
   
 END SUBROUTINE PBCHECK
