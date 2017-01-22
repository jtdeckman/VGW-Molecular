SUBROUTINE VGWFUPDATE(Q)
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I
  LOGICAL :: ATST,UD
  DOUBLE PRECISION :: Q(3,N_ATOM),ROEMAX,EFACT

  NTOTAL=NTOTAL+1

  IF(MOD(NTOTAL,UFREQ).EQ.0) THEN
     ROEMAX=ZCURR(ASIZE)
     DO I=1, ASIZE
        EFACT=EXP(ZCURR(I)-ROEMAX)
        ZTOT(I)=ZTOT(I)+EFACT
        VDATA(I)=VDATA(I)+VCURR*EFACT
     ENDDO
     NUPDATES=NUPDATES+1
  ENDIF
  
  IF((NUPDATES.GE.1).AND.(NUPDATES.EQ.WFREQ)) CALL WRTZ(Q)

END SUBROUTINE VGWFUPDATE

SUBROUTINE WRTZ(Q)
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I
  DOUBLE PRECISION :: Q(3,N_ATOM)
   
  BLCK=BLCK+1

  WRITE(FILEIND,"(A8,I6,I8)") "<BLOCK> ",BLCK,NUPDATES
  REWIND(QFILEIND)
  IF(PBC) THEN
     DO I=1,N_ATOM
        WRITE(QFILEIND,"(F9.5,A1,F9.5,A1,F9.5)") BL*Q(1,I)," ",BL*Q(2,I)," ",BL*Q(3,I)
     ENDDO
  ELSE
     DO I=1,N_ATOM
        WRITE(QFILEIND,"(F9.5,A1,F9.5,A1,F9.5)") Q(1,I)," ",Q(2,I)," ",Q(3,I)
     ENDDO
  ENDIF

  CALL FLUSH(QFILEIND)

  WRITE(FILEIND,"(A9,I12)") "<NTOTAL> ", NTOTAL
  WRITE(FILEIND,"(A9,I12)") "<NTACCP> ", NTACCP
  WRITE(FILEIND,"(A9,I12)") "<NTRANS> ", NTRANS
  WRITE(FILEIND,"(A14,F9.4)") "<ACCPTRATIO>  ", FLOAT(NTACCP)/FLOAT(NTRANS)
  WRITE(FILEIND,"(A14,F7.4)") "<STEPLENGTH>  ",STPLEN
  WRITE(FILEIND,"(A10,F12.5)") "<BOXLXYZ> ", BL
  WRITE(FILEIND,"(A8,F12.5,2I12,F9.5)") "<VSTEP> ", VSLEN,NVACCT,NVSTEP,FLOAT(NVACCT)/FLOAT(NVSTEP)
  WRITE(FILEIND,"(A7,F12.5)") "<ENOT> ", ENOT

  IF(mynode.EQ.0) THEN
     WRITE(FILEIND,"(A14,A4,F10.4)") "<PTACCPRATIO> "," -- ", FLOAT(NSWAPS(2))/FLOAT(NSWPATT(2))
  ELSE 
     IF(mynode.EQ.(NREPS-1)) THEN
        WRITE(FILEIND,"(A14,F10.4,A4)") "<PTACCPRATIO> ",FLOAT(NSWAPS(1))/FLOAT(NSWPATT(1)),"  --"
     ELSE
        WRITE(FILEIND,"(A14,F10.4,A1,F10.4)") "<PTACCPRATIO> ",FLOAT(NSWAPS(1))/FLOAT(NSWPATT(1))," ", FLOAT(NSWAPS(2))/FLOAT(NSWPATT(2))
     ENDIF
  ENDIF
  WRITE(FILEIND,"(A8,I8,A1,I8)") "<PTACC> ",NSWAPS(1)," ", NSWAPS(2)
  WRITE(FILEIND,"(A8,I8,A1,I8)") "<PTATT> ",NSWPATT(1)," ", NSWPATT(2)
  
  DO I=1,ASIZE
     WRITE(FILEIND,"(F21.16,A2,ES27.17E4,A2,ES27.17E4)") BETA(I),"  ",ZTOT(I),"  ",VDATA(I)
  ENDDO
  
  CALL FLUSH(FILEIND)

  NUPDATES=0
  ZTOT(:)=0
  VDATA(:)=0

END SUBROUTINE WRTZ

SUBROUTINE VGWFINIT
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  
  NUPDATES=0
  ZCURR=ZNEW
  ROECURR=ROENEW
   
END SUBROUTINE

SUBROUTINE VMDCOORDS(Q,TAG)
  USE VGWM
  USE VGWMC
  IMPLICIT NONE
  INTEGER :: I
  DOUBLE PRECISION :: Q(3,N_ATOM)
  CHARACTER*2 :: TAG

  WRITE(VFILEIND,"(I4)") N_ATOM
  WRITE(VFILEIND,"(A2,I9)") "# ",NTOTAL
 
  IF(PBC) THEN
     DO I=1,N_ATOM
        WRITE(VFILEIND,"(A2,A1,F12.5,A1,F12.5,A1,F12.5)") TAG," ",BL*Q(1,I)," ",BL*Q(2,I)," ",BL*Q(3,I)
     ENDDO
  ELSE
     DO I=1,N_ATOM
        WRITE(VFILEIND,"(A2,A1,F12.5,A1,F12.5,A1,F12.5)") TAG," ",Q(1,I)," ",Q(2,I)," ",Q(3,I)
     ENDDO
  ENDIF

  CALL FLUSH(VFILEIND)

END SUBROUTINE VMDCOORDS