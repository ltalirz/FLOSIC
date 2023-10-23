C
C *********************************************************************
C
       SUBROUTINE GTRHOCR(ISGGA,NPTS,RPTS,RHOC,XTMP,DTMP,RTMP)
C
C DIRK POREZAG, FEBRUARY 1998
C GTRHOCR CALCULATES THE CORE DENSITY NEEDED FOR PSEUDOPOTENTIAL
C CALCULATIONS USING NONLINEAR CORE CORRECTIONS
C
        INCLUDE 'PARAMS'
        INCLUDE 'commons.inc'
        LOGICAL   ISGGA
        DIMENSION RPTS(3,*),RHOC(10,*)
        DIMENSION XTMP(3,*),DTMP(3,*),RTMP(*)
C
C SETUP
C
        NGRAD=1
        NDRV=1
        IF (ISGGA) THEN
         NGRAD=10
         NDRV=3
        END IF
        DO IPTS=1,NPTS
         DO I=1,NGRAD
          RHOC(I,IPTS)= 0.0D0
         END DO
        END DO
        IF (ISNLCC .NE. 1) RETURN
C
C LOOP OVER ALL ATOMS
C
        DO 100 ICNT=1,NCNT
         IFU= IFUCNT(ICNT)
         IF (NLCC(IFU) .EQ. 1) THEN
          DO IPTS=1,NPTS 
           XTMP(1,IPTS)= RPTS(1,IPTS)-RCNT(1,ICNT)
           XTMP(2,IPTS)= RPTS(2,IPTS)-RCNT(2,ICNT)
           XTMP(3,IPTS)= RPTS(3,IPTS)-RCNT(3,ICNT)
           RTMP(IPTS)= SQRT(XTMP(1,IPTS)**2
     &                     +XTMP(2,IPTS)**2
     &                     +XTMP(3,IPTS)**2)
          END DO
          CALL FINTPOL(8,NPTS,RTMP,0.1D0,NRADTAB(IFU),3,NDRV,
     &                 RRADTAB(1,IFU),RHOCOR(1,1,IFU),DTMP)
          DO IPTS=1,NPTS
           RHOC(1,IPTS)= RHOC(1,IPTS)+DTMP(1,IPTS)
          END DO
          IF (ISGGA) THEN
           DO IPTS=1,NPTS
            D1= DTMP(2,IPTS)
            D2= DTMP(3,IPTS)
            RRC= 1.0D0/RTMP(IPTS)
            DO I=1,3
             FAC= XTMP(I,IPTS)*RRC
             RHOC(I+1,IPTS)= RHOC(I+1,IPTS)+D1*FAC
             RHOC(I+4,IPTS)= RHOC(I+4,IPTS)+(D2-D1*RRC)*FAC*FAC+D1*RRC
             DO J=I+1,3
              RHOC(I+J+5,IPTS)= RHOC(I+J+5,IPTS)
     &                         +(D2-D1*RRC)*FAC*XTMP(J,IPTS)*RRC
             END DO
            END DO
           END DO
          END IF
         END IF
  100   CONTINUE
        RETURN
        END
C
C *********************************************************************
C
       SUBROUTINE RHCDRV(IID,NPTS,RTMP,DTMP,RHOCDR)
C
C DIRK POREZAG, MARCH 1998
C RHCDRV CALCULATES THE RADIAL DERIVATIVES OF THE CORE DENSITY 
C NEEDED FOR NONLINEAR CORE CORRECTIONS
C RETURN: RHOCDR= dRHO/dR * 1/R WHERE R= RTMP(IPTS)
C
        INCLUDE 'PARAMS'
        INCLUDE 'commons.inc'
        LOGICAL     SKIP 
        DIMENSION   RTMP(*),DTMP(3,*),RHOCDR(*)
C
C SETUP
C
        IFU= IFUIDT(IID)
        SKIP= (ISNLCC .NE. 1) 
        IF (.NOT. SKIP) THEN
         SKIP= (NLCC(IFU) .NE. 1)
        END IF
        IF (SKIP) THEN
         DO IPTS=1,NPTS
          RHOCDR(IPTS)= 0.0D0
         END DO
         RETURN
        END IF
C
C INTERPOLATION
C
        CALL FINTPOL(8,NPTS,RTMP,0.1D0,NRADTAB(IFU),3,2,
     &               RRADTAB(1,IFU),RHOCOR(1,1,IFU),DTMP)
        DO IPTS=1,NPTS
         RHOCDR(IPTS)= DTMP(2,IPTS)/RTMP(IPTS)
        END DO
        RETURN
        END
               
