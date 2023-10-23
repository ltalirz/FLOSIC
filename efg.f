       SUBROUTINE EFG
C
C ORIGINALLY WRITTEN BY ALAN JACKSON FALL/1999
C MINOR CHANGES JK 01/2000
C
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (MAXSPH=200) 
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (NRHOG=6*MXSPN-2)
       PARAMETER (EFAK=9.717366384D0)
C                  Ha/bohr^2 -> V/m^2
C
       LOGICAL LMKFIL,EXIST 
       CHARACTER*80 LINE       
       DIMENSION RREL(3),RVECA(3,MX_GRP)
       DIMENSION EFGRAD(3,3,MAXSPH)
       DIMENSION RHOFERMI(MAXSPH)
       DIMENSION EFGVAL(3,MAXSPH),SC(3)
       DIMENSION CENTER(3,MAXSPH)
C
C
C RETURN IF INPUT FILE DOES NOT EXIST
C
        PRINT '(A)','CALCULATING ELECTRIC FIELD GRADIENTS'
        INQUIRE(FILE='EFG',EXIST=EXIST)
        IF (.NOT.EXIST) THEN
         PRINT '(2A)','EFG: FILE EFG DOES NOT EXIST ',
     &                '--> NOTHING TO DO'
         RETURN
        END IF         
C
C CREATE A STANDARD INPUT FILE IF THE CURRENT INPUT FILE IS EMPTY
C
        CALL GTTIME(TIME1)
        LMKFIL=.TRUE.
        OPEN(74,FILE='EFG',FORM='FORMATTED',STATUS='OLD')
        REWIND(74)
        READ(74,*,END=5,ERR=5) ISWITCH
        IF (ISWITCH .NE. 0) LMKFIL=.FALSE.
    5   CLOSE(74)
C                            
        IF (LMKFIL) THEN
         OPEN(74,FILE='EFG',FORM='FORMATTED',STATUS='OLD')
         REWIND(74)
         WRITE(74,*) '0  auto=0, otherwise user-defined'
         WRITE(74,1010) NIDENT
 1010    FORMAT(' ',I5,' NUMBER OF CENTERS')
         DO IATOM=1,NIDENT 
          WRITE(74,1022)(RIDT(J,IATOM),J=1,3),
     &        (SYMATM(L,IATOM),L=1,10)
 1022     FORMAT(3(1X,F10.5),5X,10A)
         END DO
         CLOSE(74)
        END IF               
C
C READ INPUT FILE
C CENTER CONTAINS THE COORDINATES (X,Y,Z) 
C
        OPEN(74,FILE='EFG',FORM='FORMATTED',STATUS='OLD')
        REWIND(74)
        READ(74,*,END=10) ISWITCH
        READ(74,*,END=10) NSPHERE
        IF (NSPHERE.GT.MAXSPH) THEN
         PRINT *,'EFG: MAXSPH MUST BE AT LEAST: ',NSPHERE
         GOTO 20
        END IF                         
        DO I=1,NSPHERE
         READ(74,33) LINE
         READ(LINE,*)(CENTER(J,I),J=1,3)
         DO J=80,1,-1
          IF (LINE(J:J).NE.' ') THEN
           DO L=1,10
            IH=J-10+L
            SYMATM(L,I)=LINE(IH:IH)
           ENDDO
           GOTO 13
          ENDIF
         ENDDO
   13   CONTINUE
        END DO
        GOTO 30
   10   PRINT *,'EFG: INPUT FILE IS INVALID'
   20   CLOSE(74)
        GOTO 900
   30   CONTINUE
   33   FORMAT(A80)
C                              
C
        PI= 4*ATAN(1.0D0)
C         CALL MOSS1(RHOFERMI,EFGRAD,CENTER,NSPHERE)
        RTOL = 0.0001D0
        DO IAT = 1,NSPHERE
         PRINT1022,(CENTER(J,IAT),J=1,3),(SYMATM(L,IAT),L=1,10)
C
C LOOP OVER ALL OTHER SITES
C
        DO IID=1,NIDENT
          ZCHRG = ZNUC(IFUIDT(IID))
          CALL GASITES(1,RIDT(1,IID),MTOT,RVECA,MSITES)
          DO JATOM=1,MTOT
             RREL(1) = CENTER(1,IAT) - RVECA(1,JATOM)
             RREL(2) = CENTER(2,IAT) - RVECA(2,JATOM)
             RREL(3) = CENTER(3,IAT) - RVECA(3,JATOM)
             R2 = RREL(1)**2 + RREL(2)**2 + RREL(3)**2 
             R5 = R2**2.5
             IF(R2.LT.RTOL) GO TO 100
             DO I=1,3
             DO J= 1,3
               EFGRAD(I,J,IAT) = EFGRAD(I,J,IAT) 
     &                - ZCHRG*3.0D0*RREL(I)*RREL(J)/R5
               IF(I.EQ.J) THEN
                  EFGRAD(I,I,IAT) = EFGRAD(I,I,IAT) + ZCHRG*R2/R5
               END IF
             END DO
             END DO
100          CONTINUE
          END DO
         END DO
       END DO
C
C
C PRINT RESULTS
C
        WRITE(74,*)
        WRITE(74,*)
        IEV=1
        DO J=1,NSPHERE
         WRITE(74,*) J,'  ',(SYMATM(L,J),L=1,10),('=',K=1,50)
         WRITE(74,*) RHOFERMI(J)
         DO K=1,3
           WRITE(74,90) (EFGRAD(I,K,J),I=1,3)
         END DO
         WRITE(74,*) ' '
         CALL DIAGSP(3,3,EFGRAD(1,1,J),EFGVAL(1,J),SC,IEV)
         WRITE(74,*) 'AFTER DIAG : '
         WRITE(74,*) (EFGVAL(I,J),I=1,3)
         WRITE(74,*) ' '
         DO K=1,3
           WRITE(74,90) (EFGRAD(I,K,J),I=1,3)
         END DO             
         WRITE(74,*) ' '
C EFG AND ETA
         EFGMAX=-100.d0
         EFGMIN=999999.9d0
         DO K=1,3
          IF(ABS(EFGVAL(K,J)).GT.EFGMAX) THEN
             KMAX=K
             EFGMAX=ABS(EFGVAL(K,J))
          ENDIF
          IF(ABS(EFGVAL(K,J)).LT.EFGMIN) THEN
             KMIN=K
             EFGMIN=ABS(EFGVAL(K,J))
          ENDIF
         ENDDO 
         KY=6-KMAX-KMIN
         VZZ=EFGVAL(KMAX,J)
         ETA=0.0D0
         IF(ABS(VZZ).GT.RTOL) ETA=(EFGVAL(KY,J)-EFGVAL(KMIN,J))/VZZ
         WRITE(74,91) (SYMATM(L,J),L=1,10),VZZ,VZZ*EFAK
         WRITE(74,92) (SYMATM(L,J),L=1,10),ETA
         WRITE(74,*)' '
        END DO
        CLOSE(74)
  90    FORMAT(3F12.6)
  91    FORMAT(1X,'EFG : ',10A,2X,'= ',F12.6,3X,'= ',
     &        F12.6,' *10**21 V/m**2')
  92    FORMAT(1X,'ETA : ',10A,2X,'= ',F12.6)
C
  900   CONTINUE
        CALL GTTIME(TIME2)
        CALL TIMOUT('ELECTRIC FIELD GRADIENTS:          ',TIME2-TIME1)    
        RETURN
       END
C
       SUBROUTINE  MOSS1(RHOFERMI,EFGRAD,CENTER,NSPHERE)
C
C COMPUTE EXPECTATION VALUES OF X**2/R**5, X*Y/R**5, ETC. FOR USE
C IN FINDING ELECTRIC FIELD GRADIENTS AT THE NUCLEAR SITES THIS IS USED
C FOR CALCULATING MOSSBAUER PARAMETERS.
C
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (LMXX=20)
       PARAMETER (LSIZ=(LMXX+1)**2)
       PARAMETER (MAXANG=((2*LMXX+1)*(LMXX+1)))
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MXCHG=56)
       CHARACTER*80 LINE
C
       LOGICAL LMKFIL,EXIST
       LOGICAL IUPDAT,ICOUNT
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP1/PTS(NSPEED,3),PSIL(MAX_OCC,LSIZ,2)
     &  ,RANG(3,MAXANG),QL(LMXX+2),QTOT(LMXX+2),DOS(LMXX+2)
     &  ,QLDS(LMXX+2,MAX_OCC),YLM(MAXANG,LSIZ)
     &  ,EMN(MAXSPH),EMX(MAXSPH),XRAD(MAXRAD),WTRAD(MAXRAD)
     &  ,CCCCCC(6,MAXSPH),ANGLE(3,MAXANG),DOMEGA(MAXANG)
     &  ,RVECA(3,MX_GRP),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3)
       COMMON/TMP2/PSIG(NMAX,MAX_OCC),RHOG(MAX_PTS,NVGRAD,MXSPN)
C
       DIMENSION  CENTER(3,NSPHERE)
C
       DIMENSION ISIZE(3),MSITES(1),SPN(2)
       DIMENSION XMAT(3,3,MAX_IDENT), PT(3), PTREL(3)
       DIMENSION EFGRAD(3,3,MAX_IDENT),RHOFERMI(MAX_IDENT)
       DATA ISIZE/1,3,6/
       DATA AU2ANG/0.529177D0/

C
C CALCULATE CHARGE DENSITIES AT EACH POINT IN SPACE:
C
        IF (DEBUG) THEN
           PRINT *, 'IN NEW MOSS1'
           PRINT *, 'NMSH NSPN', NMSH, NSPN
        END IF
           NXTRA = NSPHERE
           DO J = 1, NSPHERE
              RMSH(1,NMSH+J) = CENTER(1,J)
              RMSH(2,NMSH+J) = CENTER(2,J)
              RMSH(3,NMSH+J) = CENTER(3,J)
           END DO
           NMSH = NMSH+NXTRA
           DO ISPN=1,NSPN
            SPN(ISPN) = 0.0D0
            DO IPTS=1,NMSH
             RHOG(IPTS,ISPN,1)=0.0D0
            END DO
           END DO
C
C POINTS LOOP
C
         LPTS=0
   40    CONTINUE
          MPTS=MIN(NMAX,NMSH-LPTS)
          LPBEG=LPTS
C
C INITIALIZE PSIG 
C
          DO IWF=1,NWF
           DO IPTS=1,MPTS
            PSIG(IPTS,IWF)=0.0D0
           END DO  
          END DO  
          ISHELLA=0
          DO 86 IFNCT=1,NFNCT
           LMX1=LSYMMAX(IFNCT)+1
           DO 84 I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
            DO 82 J_POS=1,M_NUC
             CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                    RVECA,L_NUC,1)
             IF(L_NUC.NE.M_NUC)THEN
              PRINT *,'DECOMP: PROBLEM IN UNRAVEL'
              CALL STOPIT 
             END IF
             LPTS=LPBEG
             DO 80 JPTS=1,MPTS,NSPEED
              NPV=MIN(NSPEED,MPTS-JPTS+1)
              DO LPV=1,NPV
               PTS(LPV,1)=RMSH(1,LPTS+LPV)-RVECA(1,J_POS)
               PTS(LPV,2)=RMSH(2,LPTS+LPV)-RVECA(2,J_POS)
               PTS(LPV,3)=RMSH(3,LPTS+LPV)-RVECA(3,J_POS)
              END DO
C
C GET ORBITS AND DERIVATIVES
C
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
C
C UPDATING ARRAY PSIG
C
              IF (IUPDAT) THEN
               IPTS=JPTS-1
               ILOC=0
               DO 78 LI=1,LMX1
                DO MU=1,ISIZE(LI)
                 DO ICON=1,N_CON(LI,IFNCT)
                  ILOC=ILOC+1
                  IF (ICOUNT(ICON,LI)) THEN
                   DO IWF=1,NWF
                    FACTOR=PSI(ILOC,IWF,1)
                    DO LPV=1,NPV
                     PSIG(IPTS+LPV,IWF)=PSIG(IPTS+LPV,IWF)
     &               +FACTOR*GRAD(LPV,1,MU,ICON,LI)
                    END DO
                   END DO  
                  END IF
                 END DO  
                END DO  
   78          CONTINUE
              END IF
              LPTS=LPTS+NPV
   80        CONTINUE
   82       CONTINUE
   84      CONTINUE
   86     CONTINUE
C
C UPDATE CHARGE DENSITY:
C
           IWF=0
           DO ISPN=1,NSPN
           DO JWF=1,NWFS(ISPN)
           IWF=IWF+1
            DO IPTS=1,MPTS
             RHOG(IPTS+LPBEG,ISPN,1)=
     &       RHOG(IPTS+LPBEG,ISPN,1)+PSIG(IPTS,IWF)**2
            END DO
           END DO
           END DO
C
C CHECK IF ALL POINTS DONE
C
          LPTS=LPBEG+MPTS
          IF(LPTS.GT.MAX_PTS)THEN
           PRINT *,'DECOMP: ERROR: LPTS >',MAX_PTS
           CALL STOPIT
          ELSE
           IF(LPTS.LT.NMSH) GOTO 40
          END IF
  500    CONTINUE
C
C  REMOVE EXTRA POINTS FROM MESH.  FILL UP RHOFERMI ARRAY.  
C
          NMSH = NMSH - NXTRA
          DO J=1,NSPHERE
            IF(NSPN.EQ.1) THEN
              RHOFERMI(J) = 2.0D0*RHOG(NMSH+J,1,1) 
            ELSE
              RHOFERMI(J) = RHOG(NMSH+J,1,1) + RHOG(NMSH+J,2,1)
            END IF
          END DO
C 
C COMPUTE EXPECTATION VALUES
C 
          DO IDT=1,NSPHERE
          DO I=1,3
          DO J=1,3
            EFGRAD(I,J,IDT)=0.0D0
          END DO
          END DO
          END DO
          DO ISPN=1,NSPN
          DO IPTS=1,NMSH
            DO IROT = 1,NGRP
              DO I=1,3
                PT(I)=0.0
                DO J=1,3
                   PT(I) = PT(I) + RMSH(J,IPTS)*RMAT(J,I,IROT)
                END DO
              END DO
C
C  EVALUATE MATRIX ELEMENTS FOR ALL IDENTITY MEMBERS
C
              DO IDT = 1,NSPHERE
                PTREL(1) = PT(1) - CENTER(1,IDT)
                PTREL(2) = PT(2) - CENTER(2,IDT)
                PTREL(3) = PT(3) - CENTER(3,IDT)
                R2 = PTREL(1)**2 + PTREL(2)**2 + PTREL(3)**2
                R1 = SQRT(R2)
                R5 = R2**2.5D0
                DO I=1,3
                DO J=1,3
                  EFGRAD(I,J,IDT) = EFGRAD(I,J,IDT) + 
     &              3.0D0*PTREL(I)*PTREL(J)/R5*2.0D0/FLOAT(NSPN)*
     &              RHOG(IPTS,ISPN,1)*WMSH(IPTS)/FLOAT(NGRP)
                END DO
                END DO
                IF(IDT.EQ.1) VTEST = VTEST + RHOG(IPTS,ISPN,1)/R1
     &                *WMSH(IPTS)/FLOAT(NGRP)*2.0D0/FLOAT(NSPN)
              END DO
              SPN(ISPN)=SPN(ISPN)+RHOG(IPTS,ISPN,1)*WMSH(IPTS)
     &             /FLOAT(NGRP)
            END DO
          END DO
          END DO
         IF (DEBUG) THEN
          PRINT *, 'VTEST', VTEST
          DNS=SPN(1)+SPN(NSPN)
          RMN=SPN(NSPN)-SPN(1) 
          SPN(1)=DNS
          SPN(2)=RMN
         PRINT *, 'CHARGES'
         PRINT *, 'NSPN', NSPN
         PRINT *, DNS, RMN
         END IF 
         DO IDT = 1,NSPHERE
           R2MAT = EFGRAD(1,1,IDT) + EFGRAD(2,2,IDT) + EFGRAD(3,3,IDT)
           R2MAT = R2MAT/3.0D0
           DO I=1,3
              EFGRAD(I,I,IDT) = EFGRAD(I,I,IDT) - R2MAT
           END DO
          END DO
        RETURN
       END


