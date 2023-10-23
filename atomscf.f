      SUBROUTINE ATOMSCF(IFNCT)
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (MAXRAD=500)
       PARAMETER (NCOUL=5)
       PARAMETER (MATSZ=MAXLSQF)
       COMMON/TMP1/OVLTAB(MAX_CON,MAX_CON,LDIM),OVLDIA(MAX_CON,LDIM)
     &  ,HKNTAB(MAX_CON,MAX_CON,LDIM),HNLTAB(MAX_CON,MAX_CON,LDIM)
     &  ,BASFU(3,MAX_CON),TABEXP(MAX_BARE)
     &  ,RRAD(MAXRAD),WRAD(MAXRAD),VLOC(MAXRAD),VPOT(2,MAXRAD)
     &  ,VXCSUB(2),VXC(MAXRAD),RHO(5,MAXRAD),RHCL(NCOUL,3),RHOGRD(10)
     &  ,HMAT(MAX_CON,MAX_CON,LDIM),OMAT(MAX_CON,MAX_CON,LDIM)
     &  ,EVALA(MAX_CON,LDIM),SCRV(MAX_CON),RPW(5),ORCNT(3),EXCVEC(4)
     &  ,EVTB(MAX_CON*LDIM),OCCU(MAX_CON*LDIM)
     &  ,AMAT(MATSZ,MATSZ),AVEC(MATSZ),FCOF(MATSZ)
     &  ,NDEG(MAX_CON*LDIM),IVEC(MATSZ),NIGGA(2)
       DATA SCFTATM /1.0D-5/
       DATA SMALL   /1.0D-5/
C
C SETUP
C
       PI= 4*ATAN(1.0D0)
       IF (MAX_PTS .LT. 1) THEN
        PRINT *,'ATOMSCF: MAX_PTS MUST BE > 0'
        CALL STOPIT
       END IF
C
C GET SMALLEST AND LARGEST EXPONENTS
C
       ALMIN=  1.0D30
       ALMAX= -1.0D30
       DO I= 1,N_BARE(IFNCT)
        ALMIN= MIN(ALMIN,BFALP(I,IFNCT))
        ALMAX= MAX(ALMAX,BFALP(I,IFNCT))
       END DO
       ALMIN= 2*ALMIN
       ALMAX= 2*ALMAX
C
C IF ZELC IS ZERO, DEFINE STANDARD POTENTIAL AND DENSITY
C
       IF (ZELC(IFNCT) .LT. SMALL) THEN
        PRINT '(A,I3,A)','ATOM TYPE ',IFNCT,' HAS NO ELECTRONS'
        NRPFIT(IFNCT)= 0
        LDIVR(IFNCT)= 0
        RETURN
       END IF
C
C SAVE VARIABLES THAT HAVE TO BE CHANGED TO MAKE GETVLXC WORK
C
       NIGGA(1)= IGGA(1)
       NIGGA(2)= IGGA(2)
       IGGA(1)= 0
       IGGA(2)= 0
       NNCNT= NCNT
       NCNT= 1
       NIFUCNT= IFUCNT(1)
       IFUCNT(1)= IFNCT
       ORCNT(1)= RCNT(1,1)
       ORCNT(2)= RCNT(2,1)
       ORCNT(3)= RCNT(3,1)
       RCNT(1,1)= 0.0D0
       RCNT(2,1)= 0.0D0
       RCNT(3,1)= 0.0D0
       NNGRP= NGRP
       NGRP= 1
       NNSPN= NSPN
       NSPN= 1
C
C GET OVERLAP, KINETIC ENERGY, AND NONLOCAL POTENTIAL MATRIX ELEMENTS
C
       CALL OVLATM(IFNCT,OVLTAB)
       CALL HKNATM(IFNCT,HKNTAB)
       DO L= 0,LSYMMAX(IFNCT)
        L1= L+1
        DO ICON= 1,N_CON(L1,IFNCT)
         OVLDIA(ICON,L1)= OVLTAB(ICON,ICON,L1)
        END DO
       END DO
C
C GET NONLOCAL POTENTIAL MATRIX ELEMENTS
C
       DO L= 0,LSYMMAX(IFNCT)
        L1= L+1
        DO ICON= 1,N_CON(L1,IFNCT)
         DO JCON= ICON,N_CON(L1,IFNCT)
          HNLTAB(JCON,ICON,L1)= 0.0D0
         END DO
        END DO
       END DO
C
       IF (PSPSYM(IFNCT)(1:3) .NE. 'ALL') THEN
        LPSPMX= MIN(LSYMMAX(IFNCT),LMAXNLO(IFNCT))
        DO IRAD= 1,NRPSP(IFNCT)
         RR= RPSNLO(IRAD,IFNCT) 
         R2= RR*RR
         WFAC= 4*PI*WPSNLO(IRAD,IFNCT)
         DO IBARE= 1,N_BARE(IFNCT)
          TABEXP(IBARE)= EXP(-BFALP(IBARE,IFNCT)*R2)
         END DO
         RPW(1)= 1.0D0
         DO L= 0,LPSPMX
          L1= L+1
          DO ICON= 1,N_CON(L1,IFNCT)
           SUM= 0.0D0
           DO IBARE= 1,N_BARE(IFNCT)
            SUM= SUM+BFCON(IBARE,ICON,L1,IFNCT)*TABEXP(IBARE)
           END DO
           BASFU(1,ICON)= RPW(1)*SUM
          END DO
          DO ICON= 1,N_CON(L1,IFNCT)
           DO JCON= ICON,N_CON(L1,IFNCT)
            HNLTAB(JCON,ICON,L1)= HNLTAB(JCON,ICON,L1)
     &      +WFAC*BASFU(1,ICON)*BASFU(1,JCON)*VPSNLO(L1,IRAD,IFNCT)
           END DO
          END DO
          RPW(1)= RPW(1)*RR
         END DO
        END DO
       END IF

C
C GET ACCURATE RADIAL MESH 
C DEFINE LOCAL AND INITIAL POTENTIAL 
C
       ALONG= 1.0D30
       DO I= 1,N_BARE(IFNCT)
        ALONG= MIN(ALONG,BFALP(I,IFNCT))
       END DO
       ERROR= 1.0D-8
       RINFTY= RCUTOFF(2*LSYMMAX(IFNCT)+2,2*ALONG,ERROR)
       CALL OPTRMSH(MAXRAD,IFNCT,0.0D0,RINFTY,OVLDIA,ERROR,
     &              NRAD,RRAD,WRAD)
       DO IRAD= 1,NRAD
        CALL VLOCAL(1,1,IFNCT,RRAD(IRAD),VLOC(IRAD))
        VPOT(1,IRAD)= VLOC(IRAD)*EXP(-RRAD(IRAD))
       END DO
C
C ITERATIONS
C
       ITER= 0
       TROLD= 1.0D30
  100  CONTINUE
        ITER= ITER+1
C
C FOR EACH ANGULAR MOMENTUM, BUILD HAMILTONIAN MATRIX AND DIAGONALIZE
C
        NSTATES= 0
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          DO JCON= ICON,N_CON(L1,IFNCT)
           HMAT(JCON,ICON,L1)= 0.0D0
          END DO
         END DO
         DO IRAD= 1,NRAD
          RR= RRAD(IRAD)
          R2= RR*RR
          RPW(1)= 1.0D0
          DO I= 1,L
           RPW(1)= RPW(1)*RR
          END DO
          WFAC= 4*PI*WRAD(IRAD)
          DO IBARE= 1,N_BARE(IFNCT)
           TABEXP(IBARE)= EXP(-BFALP(IBARE,IFNCT)*R2)
          END DO
          DO ICON= 1,N_CON(L1,IFNCT)
           SUM= 0.0D0
           DO IBARE= 1,N_BARE(IFNCT)
            SUM= SUM+BFCON(IBARE,ICON,L1,IFNCT)*TABEXP(IBARE)
           END DO
           BASFU(1,ICON)= RPW(1)*SUM
          END DO
          DO ICON= 1,N_CON(L1,IFNCT)
           DO JCON= ICON,N_CON(L1,IFNCT)
            HMAT(JCON,ICON,L1)= HMAT(JCON,ICON,L1)
     &      +WFAC*BASFU(1,ICON)*BASFU(1,JCON)*VPOT(1,IRAD)
           END DO
          END DO
         END DO
         DO ICON= 1,N_CON(L1,IFNCT)
          DO JCON= ICON,N_CON(L1,IFNCT)
           HMAT(JCON,ICON,L1)= HMAT(JCON,ICON,L1)
     &     +HNLTAB(JCON,ICON,L1)+HKNTAB(JCON,ICON,L1)
           HMAT(ICON,JCON,L1)= HMAT(JCON,ICON,L1)
           OMAT(JCON,ICON,L1)= OVLTAB(JCON,ICON,L1)
           OMAT(ICON,JCON,L1)= OVLTAB(JCON,ICON,L1)
          END DO
         END DO
         CALL DIAGGE(MAX_CON,N_CON(L1,IFNCT),HMAT(1,1,L1),OMAT(1,1,L1),
     &               EVALA(1,L1),SCRV,1)
         DO ICON= 1,N_CON(L1,IFNCT)
          NSTATES= NSTATES+1
          EVTB(NSTATES)= EVALA(ICON,L1)
          NDEG(NSTATES)= 4*L+2
         END DO
        END DO
C
C DEFINE OCCUPATION NUMBERS
C
        ELEC= ZELC(IFNCT)
        TEMP= 0.01D0
        CALL FERMILV(NSTATES,ELEC,EFERM,TEMP,EVTB,OCCU,NDEG)
        DO I= 1,NSTATES
         OCCU(I)= OCCU(I)*NDEG(I)
        END DO
C
C GET EIGENVALUE TRACE 
C
        EVTR= 0.0D0
        DO I= 1,NSTATES
         EVTR= EVTR+EVTB(I)*OCCU(I)
        END DO
C
C GET DENSITY AND ITS DERIVATIVES AND INTEGRALS ON MESH
C
        DO IRAD= 1,NRAD
         RHO(1,IRAD)= 0.0D0
         RHO(2,IRAD)= 0.0D0
         RHO(3,IRAD)= 0.0D0
         RR= RRAD(IRAD)
         R2= RR*RR
         DO IBARE= 1,N_BARE(IFNCT)
          TABEXP(IBARE)= EXP(-BFALP(IBARE,IFNCT)*R2)
         END DO
         RPW(1)= 0.0D0
         RPW(2)= 0.0D0
         RPW(3)= 1.0D0
         RPW(4)= RR
         RPW(5)= R2
         IOFS= 0
         DO L= 0,LSYMMAX(IFNCT)
          L1= L+1
          DO ICON= 1,N_CON(L1,IFNCT)
           BASFU(1,ICON)= 0.0D0
           BASFU(2,ICON)= 0.0D0
           BASFU(3,ICON)= 0.0D0
           DO IBARE= 1,N_BARE(IFNCT)
            ALP= BFALP(IBARE,IFNCT)
            FAC= BFCON(IBARE,ICON,L1,IFNCT)*TABEXP(IBARE)
            BASFU(1,ICON)= BASFU(1,ICON)+FAC*RPW(3)
            BASFU(2,ICON)= BASFU(2,ICON)+FAC*(L*RPW(2)-2*ALP*RPW(4))
            BASFU(3,ICON)= BASFU(3,ICON)
     &      +FAC*(L*(L-1)*RPW(1)-2*ALP*((2*L+1)*RPW(3)-2*ALP*RPW(5)))
           END DO
          END DO
          DO I= 1,4
           RPW(I)= RPW(I+1)
          END DO
          RPW(5)= RPW(5)*RR
          DO ISTA= 1,N_CON(L1,IFNCT)
           WF0= 0.0D0
           WF1= 0.0D0
           WF2= 0.0D0
           DO ICON= 1,N_CON(L1,IFNCT)
            WF0= WF0+HMAT(ICON,ISTA,L1)*BASFU(1,ICON)
            WF1= WF1+HMAT(ICON,ISTA,L1)*BASFU(2,ICON)
            WF2= WF2+HMAT(ICON,ISTA,L1)*BASFU(3,ICON)
           END DO
           RHO(1,IRAD)= RHO(1,IRAD)+OCCU(IOFS+ISTA)*WF0*WF0
           RHO(2,IRAD)= RHO(2,IRAD)+OCCU(IOFS+ISTA)*WF0*WF1*2
           RHO(3,IRAD)= RHO(3,IRAD)+OCCU(IOFS+ISTA)*(WF1*WF1+WF0*WF2)*2
          END DO
          IOFS= IOFS+N_CON(L1,IFNCT)
         END DO
C
C INTEGRALS (FOR COULOMB-POT)
C
         RLW= 0.0D0
         IF (IRAD .GT. 1) RLW= RRAD(IRAD-1)
         RUP= RRAD(IRAD)
         CALL GAUSSP(RLW,RUP,NCOUL,RHCL(1,1),RHCL(1,2))
         DO JRAD= 1,NCOUL
          RHCL(JRAD,3)= 0.0D0
          RR= RHCL(JRAD,1)
          R2= RR*RR
          DO IBARE= 1,N_BARE(IFNCT)
           TABEXP(IBARE)= EXP(-BFALP(IBARE,IFNCT)*R2)
          END DO
          RPW(1)= 1.0D0
          IOFS= 0
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           DO ICON= 1,N_CON(L1,IFNCT)
            BASFU(1,ICON)= 0.0D0
            DO IBARE= 1,N_BARE(IFNCT)
             FAC= BFCON(IBARE,ICON,L1,IFNCT)*TABEXP(IBARE)
             BASFU(1,ICON)= BASFU(1,ICON)+FAC*RPW(1)
            END DO
           END DO
           RPW(1)= RPW(1)*RR
           DO ISTA= 1,N_CON(L1,IFNCT)
            WF0= 0.0D0
            DO ICON= 1,N_CON(L1,IFNCT)
             WF0= WF0+HMAT(ICON,ISTA,L1)*BASFU(1,ICON)
            END DO
            RHCL(JRAD,3)= RHCL(JRAD,3)+OCCU(IOFS+ISTA)*WF0*WF0
           END DO
           IOFS= IOFS+N_CON(L1,IFNCT)
          END DO
         END DO
         RHO(4,IRAD)= 0.0D0
         RHO(5,IRAD)= 0.0D0
         DO JRAD= 1,NCOUL
          RHO(4,IRAD)= RHO(4,IRAD)
     &                +RHCL(JRAD,2)*RHCL(JRAD,3)*RHCL(JRAD,1)
          RHO(5,IRAD)= RHO(5,IRAD)
     &                +RHCL(JRAD,2)*RHCL(JRAD,3)*RHCL(JRAD,1)**2
         END DO
        END DO
C
C GET EXCHANGE-CORRELATION POTENTIAL: USE SUBVLXC WITH MODE 1
C WHICH DOESN'T CALCULATE THE LOCAL POTENTIAL
C
        NMSH= 1
        RMSH(1,1)= 0.0D0 
        RMSH(2,1)= 0.0D0 
        WMSH(1)= 0.0D0
        RHOGRD( 2)= 0.0D0
        RHOGRD( 3)= 0.0D0
        RHOGRD( 8)= 0.0D0
        RHOGRD( 9)= 0.0D0
        RHOGRD(10)= 0.0D0
        DO IRAD= 1,NRAD
         RRC= 1.0D0/RRAD(IRAD)
         RMSH(3,1)= RRAD(IRAD)
         RHOGRD( 1)= RHO(1,IRAD)
         RHOGRD( 4)= RHO(2,IRAD)
         RHOGRD( 5)= RHO(2,IRAD)*RRC  
         RHOGRD( 6)= RHO(2,IRAD)*RRC  
         RHOGRD( 7)= RHO(3,IRAD)  
         CALL SUBVLXC(1,0,1,RHOGRD,VXCSUB,VLOC,EXCVEC)
         VXC(IRAD)= VXCSUB(1)
        END DO
C
C COULOMB POTENTIAL
C        
        DO IRAD= 2,NRAD
         RHO(4,NRAD-IRAD+1)= RHO(4,NRAD-IRAD+1)+RHO(4,NRAD-IRAD+2)
         RHO(5,IRAD)= RHO(5,IRAD)+RHO(5,IRAD-1)
        END DO
        DO IRAD= 1,NRAD-1
         VPOT(2,IRAD)= 4*PI*(RHO(4,IRAD+1)+RHO(5,IRAD)/RRAD(IRAD))
        END DO
        VPOT(2,NRAD)= 4*PI*RHO(5,NRAD)/RRAD(NRAD)
C
C PUTTING IT ALL TOGETHER
C
        DO IRAD= 1,NRAD
         VPOT(2,IRAD)= VPOT(2,IRAD)+VXC(IRAD)+VLOC(IRAD)
        END DO
C
C CHECK FOR CONVERGENCE
C
        IF (ABS(EVTR-TROLD) .LT. SCFTATM) GOTO 200
        TROLD= EVTR
C
C SLIGHTLY BETTER THAN SIMPLE MIXING
C
        IF (ITER .EQ. 1) THEN
         DO IRAD= 1,NRAD
          VPOT(1,IRAD)= VPOT(2,IRAD)
         END DO
        ELSE 
         AMIX= 0.10D0
         DO IRAD= 1,NRAD
          VPOT(1,IRAD)= (1.0D0-AMIX)*VPOT(1,IRAD)+AMIX*VPOT(2,IRAD)
         END DO
        END IF
        GOTO 100
  200  CONTINUE
       PRINT '(A,I3,2A)','SCF FOR ATOM TYPE ',IFNCT,' SUCCESSFUL, ',
     &                   'EIGENVALUES: '
       DO L= 0,LSYMMAX(IFNCT)
        L1= L+1
        PRINT 1010,(EVALA(I,L1), I= 1,N_CON(L1,IFNCT))
 1010   FORMAT(5(1X,F12.5))
       END DO
C
C UPDATE ALMIN (TAKE RHO**(1/3) BEHAVIOR OF EXC INTO ACCOUNT)
C
       ALMIN= 0.335D0*ALMIN
       RPFALP(IFNCT)= ALMIN
       I= 0
       X= ALMIN
  210  CONTINUE
        I= I+1
        X= 2*X
        IF (X .LT. ALMAX) GOTO 210
       CONTINUE
       I= MAX(3,I)
       IF (I .GT. MAXLSQF) THEN
        PRINT *,'ATOMSCF: MAXLSQF MUST BE AT LEAST: ',I
        CALL STOPIT
       END IF
       NRPFIT(IFNCT)= I
C
C CHECK FOR SIGN CHANGE IN RHO/POT
C
       DO IRAD= 2,NRAD
        IF ((RHO (1,IRAD)*RHO (1,IRAD-1) .LT. 0.0D0) .OR.
     &      (VPOT(2,IRAD)*VPOT(2,IRAD-1) .LT. 0.0D0)) THEN
         PRINT *,'ATOMSCF: SIGN CHANGE IN RHO/POT'
         PRINT *,'MOST LIKELY YOUR LOCAL POTENTIAL CHANGES ITS SIGN'
         PRINT *,'WHICH IS A BAD, BAD THING'
         CALL STOPIT
        END IF
       END DO
C
C DENSITY FIT (USE SPACE IN VLOC AND IGNORE FIRST TWO EXPONENTS)
C RENORMALIZE DENSITY TO GET CORRECT NORM
C
       DO IRAD= 1,NRAD
        VLOC(IRAD)= RHO(1,IRAD)
       END DO
       CALL ATMFIT(4*ALMIN,NRPFIT(IFNCT)-2,NRAD,RRAD,WRAD,VLOC,
     &             MATSZ,FCOF,AMAT,AVEC,IVEC)
       RPFCMX(IFNCT)= 0.0D0
       RPFCOF(1,1,IFNCT)= 0.0D0
       RPFCOF(1,2,IFNCT)= 0.0D0
       X= 4*ALMIN 
       SUM= 0.0D0
       DO I= 1,NRPFIT(IFNCT)-2
        SUM= SUM+FCOF(I)*(SQRT(PI/X))**3
        X= 2*X
       END DO
       SUM= ZELC(IFNCT)/SUM
       DO I= 1,NRPFIT(IFNCT)-2
        FCOF(I)= FCOF(I)*SUM
        RPFCOF(1,I+2,IFNCT)= FCOF(I)
        RPFCMX(IFNCT)= MAX(RPFCMX(IFNCT),ABS(FCOF(I)))
       END DO
       PRINT '(A,F12.8,A)','CHARGE ERROR IN DENSITY FIT: ',
     &                     ZELC(IFNCT)*(SUM-1.0D0),' (CORRECTED)'
C
C POTENTIAL FIT (USE SPACE IN VLOC)
C IF THERE IS NO 1/R SINGULARITY, FIT -POT, OTHERWISE FIT -POT*R
C
       ZZ= (VPOT(2,1)-VPOT(2,2))/(1.0D0/RRAD(1)-1.0D0/RRAD(2))
       IF (ABS(ZZ) .LT. 0.1D0) THEN
        LDIVR(IFNCT)= 0
        DO IRAD= 1,NRAD
         VLOC(IRAD)= -VPOT(2,IRAD)
        END DO
       ELSE
        LDIVR(IFNCT)= 1
        DO IRAD= 1,NRAD
         VLOC(IRAD)= -VPOT(2,IRAD)*RRAD(IRAD)
         WRAD(IRAD)= WRAD(IRAD)/RRAD(IRAD)**2
        END DO
       END IF
       CALL ATMFIT(ALMIN,NRPFIT(IFNCT),NRAD,RRAD,WRAD,VLOC,
     &             MATSZ,FCOF,AMAT,AVEC,IVEC)
       DO I= 1,NRPFIT(IFNCT)
        RPFCOF(2,I,IFNCT)= FCOF(I)
        RPFCMX(IFNCT)= MAX(RPFCMX(IFNCT),ABS(FCOF(I)))
       END DO
C
C CHECK IF FITTED DENSITY / POTENTIAL HAS SIGN CHANGE
C
       IFAIL= 0
       IF (RPFCOF(1,3,IFNCT) .LT. 0.0D0) IFAIL= IFAIL+1 
       IF (RPFCOF(2,1,IFNCT) .LT. 0.0D0) IFAIL= IFAIL+1 
       IF (IFAIL .EQ. 0) THEN
        RLIMT= 0.0D0
        DO IR= 1,2
         IOFS= 3
         IF (IR .EQ. 2) IOFS= 1
         DO I= IOFS+1,NRPFIT(IFNCT)
          IF (RPFCOF(IR,I,IFNCT) .LT. 0.0D0) THEN
           RLIMT= MAX(RLIMT,
     &           (LOG(-RPFCOF(IR,I,IFNCT))-LOG(RPFCOF(IR,IOFS,IFNCT)))
     &           /(ALMIN*(2**(I-1)-2**(IOFS-1))))
          END IF
         END DO
        END DO
        RLIMT= SQRT(RLIMT)
        DLTX= 0.01D0
        NSTEP= INT(RLIMT/DLTX)
        DO ISTEP= 0,NSTEP
         CALL RPFIT(IFNCT,ISTEP*DLTX,1.0D0,RHOF,POTF)
         IF ((RHOF .LT. 0.0D0) .OR. (POTF .GT. 0.0D0)) IFAIL= IFAIL+1
        END DO
       END IF
       IF (IFAIL .NE. 0) THEN
        PRINT '(A,I8)','WARNING: DENSITY/POTENTIAL FIT CHANGES SIGN'
       ELSE
        PRINT '(A)','POTENTIAL AND DENSITY FITS ARE OK'
       END IF
C
C RESTORE ALTERED VARIABLES AND LEAVE
C
       IGGA(1)= NIGGA(1)
       IGGA(2)= NIGGA(2)
       NCNT= NNCNT
       IFUCNT(1)= NIFUCNT
       RCNT(1,1)= ORCNT(1)
       RCNT(2,1)= ORCNT(2)
       RCNT(3,1)= ORCNT(3)
       NGRP= NNGRP
       NSPN= NNSPN
       RETURN
      END
C
C ********************************************************************
C
       SUBROUTINE OVLATM(IFNCT,OVLTAB)
C
C CREATE TABLE OF ATOMIC OVERLAP INTEGRALS
C
        INCLUDE 'PARAMS'
        INCLUDE 'commons.inc'
        DIMENSION OVLTAB(MAX_CON,MAX_CON,LDIM)
C 
        PI= 4*ATAN(1.0D0) 
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          DO JCON= 1,N_CON(L1,IFNCT)
           OVLTAB(JCON,ICON,L1)= 0.0D0
          END DO
         END DO
        END DO
C
        DO IBARE= 1,N_BARE(IFNCT)
         DO JBARE= 1,N_BARE(IFNCT)
          ALP= BFALP(IBARE,IFNCT)+BFALP(JBARE,IFNCT)
          ALRC= 1.0D0/ALP
          FACO= 2*PI*SQRT(PI*ALRC)
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           FACO= 0.5D0*ALRC*(2*L+1)*FACO
           DO ICON= 1,N_CON(L1,IFNCT)
            DO JCON= 1,N_CON(L1,IFNCT)
             P= BFCON(IBARE,ICON,L1,IFNCT)*BFCON(JBARE,JCON,L1,IFNCT)
             OVLTAB(JCON,ICON,L1)= OVLTAB(JCON,ICON,L1)+FACO*P
            END DO
           END DO
          END DO
         END DO
        END DO
        RETURN
       END
C
C ********************************************************************
C
       SUBROUTINE HKNATM(IFNCT,HKNTAB)
C
C CREATE TABLE OF ATOMIC KINETIC ENERGY INTEGRALS
C
        INCLUDE 'PARAMS'
        INCLUDE 'commons.inc'
        DIMENSION HKNTAB(MAX_CON,MAX_CON,LDIM)
C 
        PI= 4*ATAN(1.0D0) 
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          DO JCON= ICON,N_CON(L1,IFNCT)
           HKNTAB(JCON,ICON,L1)= 0.0D0
          END DO
         END DO
        END DO
C
        DO IBARE= 1,N_BARE(IFNCT)
         DO JBARE= 1,N_BARE(IFNCT)
          ABET= BFALP(JBARE,IFNCT)
          ASUM= BFALP(IBARE,IFNCT)+ABET
          ASRC= 1.0D0/ASUM
          FAC1= 2*PI*SQRT(PI*ASRC)
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           FAC1= 0.5D0*ASRC*(2*L+1)*FAC1
           FAC2= 0.5D0*ASRC*(2*L+3)*FAC1
           DO ICON= 1,N_CON(L1,IFNCT)
            DO JCON= 1,N_CON(L1,IFNCT)
             P= BFCON(IBARE,ICON,L1,IFNCT)*BFCON(JBARE,JCON,L1,IFNCT)
             HKNTAB(JCON,ICON,L1)= HKNTAB(JCON,ICON,L1)
     &                            +ABET*P*((2*L+3)*FAC1-2*ABET*FAC2)
            END DO
           END DO
          END DO
         END DO
        END DO
        RETURN
       END
C
C ********************************************************************
C
       SUBROUTINE RPFIT(IFNCT,R,RRC,RHO,POT)
        INCLUDE 'PARAMS'
        INCLUDE 'commons.inc'
        DATA ACCU /1.0D-20/
C
        FAC= EXP(-RPFALP(IFNCT)*R*R)
        RHO= 0.0D0
        POT= 0.0D0
        DO I= 1,NRPFIT(IFNCT)
         RHO= RHO+RPFCOF(1,I,IFNCT)*FAC
         POT= POT+RPFCOF(2,I,IFNCT)*FAC
         FAC= FAC*FAC
         IF (FAC*RPFCMX(IFNCT) .LT. ACCU*MIN(RHO,POT)) GOTO 100
        END DO
  100   CONTINUE
        IF (LDIVR(IFNCT) .EQ. 1) THEN
         POT= -POT*RRC
        ELSE
         POT= -POT
        END IF
        RETURN
       END
C
C ********************************************************************
C
       SUBROUTINE ATMFIT(ALMIN,NALP,NPTS,RPTS,WPTS,FPTS,
     &                   MATSZ,FCOF,AMAT,AVEC,IVEC)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RPTS(NPTS),WPTS(NPTS),FPTS(NPTS)
        DIMENSION FCOF(MATSZ),AMAT(MATSZ,MATSZ),AVEC(MATSZ),IVEC(MATSZ)
        SAVE
C
C LOW-ACCURACY LEAST-SQUARES FIT OF FUNCTION F(R) GIVEN ON A MESH 
C OF POINTS FPTS(RPTS) WITH WEIGHT WPTS TO THE FUNCTIONAL FORM:
C
C       NALP                        2                   I-1
C F(R)= SUM  FCOF(I) * EXP(-ALP(I)*R )   WHERE ALP(I)= 2   * ALMIN
C       I=1
C
        IF (NALP .LT. 1) RETURN
        IF (NALP .GT. MATSZ) THEN
         PRINT *,'ATMFIT: MATSZ MUST BE AT LEAST: ',NALP
         CALL STOPIT
        END IF
        DO I= 1,NALP
         FCOF(I)= 0.0D0
         DO J= I,NALP
          AMAT(J,I)= 0.0D0
         END DO
        END DO
C
C SETUP AMAT, AVEC
C
        DO IPTS= 1,NPTS
         RR= RPTS(IPTS)
         WR= WPTS(IPTS)
         FR= FPTS(IPTS)
         R2= RR*RR
         REX1= EXP(-ALMIN*R2)
         DO I= 1,NALP
          FAC1= WR*REX1
          FCOF(I)= FCOF(I)+FAC1*FR
          REX2= REX1
          DO J= I,NALP
           AMAT(J,I)= AMAT(J,I)+FAC1*REX2
           REX2= REX2*REX2
          END DO
          REX1= REX1*REX1
         END DO
        END DO
        DO I= 1,NALP
         DO J= I,NALP
          AMAT(I,J)= AMAT(J,I)
         END DO
        END DO
C
C GET SOLUTION
C
        CALL GAUSSPIV(MATSZ,1,NALP,1,AMAT,FCOF,AVEC,IVEC,IER)
        IF (IER .NE. 0) THEN
         PRINT *,'ATMFIT: ERROR IN GAUSSPIV'
         CALL STOPIT
        END IF
        RETURN
       END 
c       
c ***********************************************************
c
        SUBROUTINE GAUSSPIV(NA,NB,NDIM,NCOL,A,B,SCR,ISCR,IER)
c
c ***********************************************************
c  Copyright by Dirk V. Porezag (porezag@dave.nrl.navy.mil) *
c ***********************************************************/
c 
c Gausspiv performs a Gaussian elimination to solve the matrix 
c equation a*x = b. Complete line and row pivoting is applied.
c The result x is stored in b. Gausspiv modifies a and b.
c
c Dimension of the matrices: a(na rows, na columns) 
c                            b(na rows, nb columns)
c
c Input:   na,nb   see above matrix dimensions
c          ndim    actual dimension of the problem (must be <= na)
c          ncol    actual number of vectors in b (must be <= nb)
c          a       matrix a
c          b       matrix b
c          scr     scratch vector of size ndim
c          iscr    scratch vector of size ndim
c Output:  b       matrix x with a*x = b
c          ier     error code
c
c Error codes:      0: no error
c                  -1: na, nb, ncol, and ndim are incompatible
c                  -2: should not happen, can only be caused by
c                      a bug or broken compiler.
c                   n: (n > 0) means that the matrix a has a zero
c                      determinant and that the gauss algorithm
c                      failed in the nth line
c
c *************************************************************** */
c
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER(ZERO=1.0D-30)
         DIMENSION A(NA,NA),B(NA,NB),SCR(NDIM),ISCR(NDIM)
         SAVE
c
c special cases and setup
c
         IER= 0
         IF ((NDIM .LT. 1) .OR. (NCOL .LT. 1)) RETURN
         IER= -1
         IF ((NDIM .GT. NA) .OR. (NCOL .GT. NB)) RETURN
         DO I= 1,NDIM
          ISCR(I)= I
         END DO
c
c Gauss elimination: construction of triangular form
c
         DO I= 1,NDIM
c
c look for largest element in submatrix
c
          DMAX= ABS(A(I,I))
          INDJ= I
          INDK= I
          DO J= I,NDIM
           DO K= I,NDIM
            SAV= ABS(A(K,J))
            IF (SAV .GT. DMAX) THEN
             DMAX= SAV
             INDJ= J
             INDK= K
            END IF
           END DO
          END DO
          XREF= ZERO*ABS(A(1,1))
          IF (I .EQ. 1) XREF= 0.0D0
          IF (DMAX .LE. XREF) THEN
           IER= I
           RETURN
          END IF
c
c exchange rows
c
          IF (I .NE. INDK) THEN
           DO J= I,NDIM
            SAV= A(I,J)
            A(I,J)= A(INDK,J)
            A(INDK,J)= SAV
           END DO
           DO J= 1,NCOL
            SAV= B(I,J)
            B(I,J)= B(INDK,J)
            B(INDK,J)= SAV
           END DO
          END IF
c
c exchange columns
c
          IF (I .NE. INDJ) THEN
           DO J= 1,NDIM
            SAV= A(J,I)
            A(J,I)= A(J,INDJ)
            A(J,INDJ)= SAV
           END DO
           ISAV= ISCR(I)
           ISCR(I)= ISCR(INDJ)
           ISCR(INDJ)= ISAV 
          END IF
c
c one step towards triangular form, invert diagonal element
c
          AREC= 1.0D0/A(I,I)
          A(I,I)= AREC
          DO J= I+1,NDIM
           SCR(J)= AREC*A(J,I)
          END DO
          DO K= I+1,NDIM
           FAC= A(I,K)
           DO J= I+1,NDIM
            A(J,K)= A(J,K)-FAC*SCR(J)
           END DO
          END DO
          DO K= 1,NCOL
           FAC= B(I,K)
           DO J= I+1,NDIM
            B(J,K)= B(J,K)-FAC*SCR(J)
           END DO
          END DO
         END DO
c
c Gauss elimination: backward substitution 
c
         DO I= NDIM,1,-1
          AREC= A(I,I)
          DO K= I+1,NDIM
           SCR(K)= A(I,K)
          END DO
          DO J= 1,NCOL
           DO K= I+1,NDIM
            B(I,J)= B(I,J)-SCR(K)*B(K,J)
           END DO
           B(I,J)= B(I,J)*AREC
          END DO
         END DO
c
c restore old ordering
c 
         IER= -2
         DO I= 1,NDIM
          DO INDI= I,NDIM
           IF (ISCR(INDI) .EQ. I) GOTO 10
          END DO
          RETURN
   10     ISCR(INDI)= ISCR(I)
          DO J= 1,NCOL
           SAV= B(I,J)
           B(I,J)= B(INDI,J)
           B(INDI,J)= SAV
          END DO
         END DO
         IER= 0
         RETURN
        END
