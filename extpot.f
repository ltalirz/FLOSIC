C
C this defines an external potential of the form
C
C Vext= Sum_N Vext_i
C Vext_i= C(i)/r**M(i) * (x-Ax(i))**Nx(i)*exp(-Gx(i)*(x-Ax(i))**2) *
C                      * (y-Ay(i))**Ny(i)*exp(-Gy(i)*(y-Ay(i))**2) *
C                      * (z-Az(i))**Nz(i)*exp(-Gz(i)*(z-Az(i))**2)  
C
C===========EXTPOT==================================================
C NTERMS                # NUMBER OF TERMS'
C C M                   # PREFACTOR C(i), POWER M(i) of 1/r**M(i)')
C Nx Ny Nz              # POWER OF POLYNOMIALS NX(i), NY(i),NZ(i)')
C Ax Ay Az              # CENTER OF EXT POT Ax(i),Ay(i),Az(i)')
C gx gy gz              # EXP OF ENVELOPE Gx(i),Gy(i),Gz(i)')
C the last 4 lines are repeated for every NTERM
C==================================================================
C
      SUBROUTINE READEXT
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXT=10)
      COMMON/EXPOT/C(MAXT),M(MAXT),N(3,MAXT),A(3,MAXT),G(3,MAXT),NTERMS
      LOGICAL EXIST
      CHARACTER*1 AXIS(3)
      DATA AXIS/'X','Y','Z'/
      INQUIRE(FILE='EXTPOT',EXIST=EXIST)
      IF(.NOT.EXIST)RETURN
      OPEN(90,FILE='EXTPOT',FORM='FORMATTED')
      READ(90,*)NTERMS
      IF(NTERMS.GT.MAXT) THEN
       PRINT*,'MAXT IN EXTPOT TOO SMALL : ',MAXT,NTERMS
       CALL STOPIT
      ENDIF
      PRINT*,'EXTERNAL POTENTIAL FROM EXTPOT USED'
      DO I=1,NTERMS
       READ(90,*)C(I),M(I)
       READ(90,*)(N(J,I),J=1,3)
       READ(90,*)(A(J,I),J=1,3)
       READ(90,*)(G(J,I),J=1,3)
      END DO
      REWIND(90)
      WRITE(90,'(I4,23X,A)')NTERMS,'# NUMBER OF TERMS'
      DO I=1,NTERMS
       WRITE(90,10)C(I),M(I)
       WRITE(90,11)(N(J,I),J=1,3)
       WRITE(90,12)(A(J,I),J=1,3)
       WRITE(90,13)(G(J,I),J=1,3)
      END DO
      WRITE(90,*)
      WRITE(90,*)
      DO I=1,NTERMS
       WRITE(90,20)C(I),M(I),(AXIS(J),A(J,I),N(J,I),G(J,I),
     &  AXIS(J),A(J,I),J=1,3)
      END DO
 10   FORMAT(F8.4,I4,15X,'# PREFACTOR C(i), POWER M(i) of 1/r**M(i)')
 11   FORMAT(3I4,15X,'# POWER OF POLYNOMIALS NX(i), NY(i),NZ(i)')
 12   FORMAT(3F8.4,3X,'# CENTER OF EXT POT Ax(i),Ay(i),Az(i)')
 13   FORMAT(3F8.4,3X,'# EXP OF ENVELOPE Gx(i),Gy(i),Gz(i)')
 20   FORMAT(' +',F6.2,'/R**',I1,3(/,'*(',A1,'-',F6.2,')**',I1,
     &   '*EXP(-',F6.2,'*(',A1,'-',F6.2,')**2)'))
      CLOSE(90)
      RETURN
      END                                      
C
      SUBROUTINE EXTPOT(R,POT,DERIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXT=10)
      COMMON/EXPOT/C(MAXT),M(MAXT),N(3,MAXT),A(3,MAXT),G(3,MAXT),NTERMS
      LOGICAL FIRST
      DIMENSION R(3),DERIV(3)
      DATA FIRST/.TRUE./
      DATA EPS/ 1.0D-7/
      IF(FIRST)THEN
        NTERMS=0
        CALL READEXT
        FIRST=.FALSE.
      END IF
      POT=0.0D0
      DO J=1,3
      DERIV(J)=0.0D0
      END DO
      DO I=1,NTERMS
          XX=R(1)-A(1,I)
          YY=R(2)-A(2,I)
          ZZ=R(3)-A(3,I)
          RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
              IF(RR.GT.1.0D-6)THEN
                RP=1.0D0
                DO J=1,M(I)
                RP=RP/RR
                END DO
                ELSE
                RP=0.0D0
               ENDIF
          POLX=1.0D0
          POLY=1.0D0
          POLZ=1.0D0
          DERX=1.0D0
          DERY=1.0D0
          DERZ=1.0D0
          DO J=1,N(1,I)
          POLX=POLX*XX
             IF(J.GT.1)DERX=DERX*XX
          END DO
          DO J=1,N(2,I)
          POLY=POLY*YY
             IF(J.GT.1)DERY=DERY*YY
          END DO
          DO J=1,N(3,I)
          POLZ=POLZ*ZZ
             IF(J.GT.1)DERZ=DERZ*ZZ
          END DO
          DERX=N(1,I)*DERX
          DERY=N(2,I)*DERY
          DERZ=N(3,I)*DERZ
          X2=XX*XX
          Y2=YY*YY
          Z2=ZZ*ZZ
          ENV=EXP(-G(1,I)*X2-G(2,I)*Y2-G(3,I)*Z2)
          POT=POT+C(I)*POLX*POLY*POLZ*ENV*RP
          DERIV(1)=DERIV(1)+C(I)*DERX*POLY*POLZ*ENV*RP
          DERIV(2)=DERIV(2)+C(I)*POLX*DERY*POLZ*ENV*RP
          DERIV(3)=DERIV(3)+C(I)*POLX*POLY*DERZ*ENV*RP
          DERIV(1)=DERIV(1)-C(I)*POLX*POLY*POLZ*ENV*2.0*G(1,I)*XX*RP
          DERIV(2)=DERIV(2)-C(I)*POLX*POLY*POLZ*ENV*2.0*G(2,I)*YY*RP
          DERIV(3)=DERIV(3)-C(I)*POLX*POLY*POLZ*ENV*2.0*G(3,I)*ZZ*RP
          IF(RR.GT.EPS)THEN
          DERIV(1)=DERIV(1)-C(I)*POLX*POLY*POLZ*ENV*M(I)*XX*RP/RR/RR
          DERIV(2)=DERIV(2)-C(I)*POLX*POLY*POLZ*ENV*M(I)*YY*RP/RR/RR
          DERIV(3)=DERIV(3)-C(I)*POLX*POLY*POLZ*ENV*M(I)*ZZ*RP/RR/RR
          END IF
      END DO
      END 
C
C 
      SUBROUTINE AFPOT(NSPN,RHOG,POT,POTIN)
C AFPOT adds a spin dependent potential of the form
C V(s_up)=V(s_up)+c(s_up)*exp(-alpha*(r-a)**2)
C V(s_dn)=V(s_dn)+c(s_dn)*exp(-alpha(r-a)**2)
C at sites given in AFPOTSYM/AFPOT
C
      INCLUDE 'PARAMS'
      PARAMETER(MAXT=10)
      COMMON/MESH/WMSH(MAX_PTS),R(3,MAX_PTS),NMSH
      DIMENSION C(2,MAXT),A(3,MAXT),G(MAXT)
      DIMENSION POT(*),B(3,MX_GRP),MSITES(1)
      DIMENSION POTIN(*),RHOG(MAX_PTS,NVGRAD,MXSPN)
      DIMENSION AHELP(2)
      LOGICAL FIRST,EXIST
      DATA FIRST/.TRUE./
      DATA EPS/ 1.0D-7/
      SAVE
      IF(FIRST)THEN
        NTERMS=0
        INQUIRE(FILE='AFPOT',EXIST=EXIST)
       IF(EXIST)THEN
        OPEN(43,FILE='AFPOT')
        READ(43,*)NTERMS
         DO I=1,NTERMS
         READ(43,*)(A(J,I),J=1,3)
         READ(43,*)G(I),(C(ISPN,I),ISPN=1,NSPN)
         END DO
        CLOSE(43)
       END IF                     
C        FIRST=.FALSE.
      END IF
C
      DO I=1,NTERMS
      CALL GASITES(1,A(1,I),NATOMS,B,MSITES)
        DO J=1,NATOMS
        PRINT*,' AFPOT: ',I,(B(K,J),K=1,3),C(1,I),C(2,I),NSPN
          DO IPTS=1,NMSH
          XX=R(1,IPTS)-B(1,J)
          YY=R(2,IPTS)-B(2,J)
          ZZ=R(3,IPTS)-B(3,J)
          RR=XX*XX+YY*YY+ZZ*ZZ
          DO ISPN=1,NSPN
           AHELP(ISPN)=EXP(-G(I)*RR)*C(ISPN,I)
           POT(IPTS+(ISPN-1)*NMSH)=
     &     POT(IPTS+(ISPN-1)*NMSH)+AHELP(ISPN)
          END DO
          IF(NSPN.EQ.1) THEN
           AHELP(2)=AHELP(1)
          ENDIF
          AFTMP=AHELP(1)*(RHOG(IPTS,1,1)-RHOG(IPTS,1,NSPN))
          AFTMP=AFTMP+AHELP(2)*RHOG(IPTS,1,NSPN)
          IF(RHOG(IPTS,1,1).GT.1D-6) THEN
          POTIN(IPTS)=POTIN(IPTS)+AFTMP/RHOG(IPTS,1,1)
          ENDIF
          END DO
        END DO                    
      END DO                    
      RETURN
      END      
