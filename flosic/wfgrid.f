C
C **************************************************************
C
C WFGRID DIRK POREZAG JUNE 1998 
C CALCULATION OF WAVEFUNCTION DENSITY ON A GRID OF POINTS
C
C MODIFIED TO WRITE THE OUTPUT IN THE GAUSSIAN CUBE FORMAT
C THE HOMO DENSITY IS PRINTED OUT BY DEFAULT  -- TB 04/03
C
       SUBROUTINE WFGRID
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (NTOTWF=MAX_OCC)
       LOGICAL EXIST,IUPDAT,ICOUNT,WFF
       CHARACTER*20 FORMSTR
       CHARACTER*20 BAND, TRASH
       CHARACTER*6 FILE
       CHARACTER*9 FILENAME(20)
       COMMON/TMP2/PSIG(MPBLOCK,MAX_OCC)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,RVECA(3,MX_GRP),RBAS(3,4),RGRID(3,MPBLOCK)
     &  ,NGRID(3),INFOWF(4,MAX_OCC),ICOUNT(MAX_CON,3)
       DIMENSION ISIZE(3)
       DIMENSION INDWRTHOMO(NTOTWF),ISP(NTOTWF),IRP(NTOTWF),IBS(NTOTWF)
       DATA ISIZE/1,3,6/
       DATA FACT/0.148203857088/
C
C READ FILE WFGRID WHICH CONTAINS: 
C * MODE (1=UNFORMATTED, 2=FORMATTED)
C * NUMBER OF POINTS FOR THE THREE BASIS VECTORS
C * ORIGIN
C * BASIS VECTORS
C * NUMBER OF STATES
c * COR/VAL STATES; VAL REQUIRES LEVELS ORDERED WRT HOMO: E.G. +1, -1
C * SPIN, REPRESENTATION, INDEX FOR EACH STATE (ORDERED) FOR OPTION COR
C
       INQUIRE(FILE='WFGRID',EXIST=EXIST)
       FORMSTR= ' '
       IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
       PRINT '(2A)','WAVEFUNCTION GRID',FORMSTR
       IF (.NOT. EXIST) RETURN

       NWAVF=1
       NSTORE=NMSH
       OPEN(80,FILE='WFGRID',FORM='FORMATTED',STATUS='OLD')
       REWIND(80)
       READ(80,*,END=110) ITYPE,IFORM
       GO TO 120
***********************TB******************************
 110   CONTINUE
       REWIND(80) 
       
C      CREATE A DEFAULT INPUT FILE   - TB
       ITYPE=1
       IFORM=2
       DO I=1,3
          RBAS(I,1)=1.0D30
          RBAS(I,2)=-1.0D30
       END DO
       DO ICNT=1,NCNT
          DO I=1,3
            RBAS(I,1)=MIN(RBAS(I,1),RCNT(I,ICNT))
            RBAS(I,2)=MAX(RBAS(I,2),RCNT(I,ICNT))
          END DO
       END DO
       DO I=1,3
          RBAS(I,1)=RBAS(I,1)-5.0D0
          RBAS(I,2)=RBAS(I,2)+5.0D0
          NGRID(I)=(RBAS(I,2)-RBAS(I,1))/0.5D0+2
          RBAS(I,2)=(RBAS(I,2)-RBAS(I,1))/(NGRID(I)-1)
       END DO
       WRITE(80,1010) 1,2, 'GRID MODE, FORMATTED FILE IN CUBE FORMAT'
       WRITE(80,1020) (NGRID(I),I=1,3),   'NUMBER OF GRID POINTS'
       WRITE(80,1030) (RBAS(I,1),I=1,3),   'ORIGIN'
       WRITE(80,1030) RBAS(1,2), 0.0D0, 0.0D0,  'BASIS VECTOR 1'
       WRITE(80,1030) 0.0D0, RBAS(2,2), 0.0D0,  'BASIS VECTOR 2'
       WRITE(80,1030) 0.0D0, 0.0D0, RBAS(1,2),  'BASIS VECTOR 1'
 1010  FORMAT( 2(I6,1X), 8X,A)
 1020  FORMAT( 3(I6,1X), 1X,A)
 1030  FORMAT( 3(F7.3,1X), 1X,A)
       REWIND(80)
       READ(80,*) ITYPE,IFORM
 120   CONTINUE
***********************TB******************************

       IF (ITYPE .GT. 2) ITYPE= 2
       IF (ITYPE .LT. 1) ITYPE= 1
       IF (IFORM .GT. 2) IFORM= 2
       IF (IFORM .LT. 1) IFORM= 1
       IF (ITYPE .EQ. 2) IFORM= 2
       FORMSTR='UNFORMATTED'
       IF (IFORM .EQ. 2) FORMSTR='FORMATTED'
C
C GET NUMBER OF GRID POINTS, ORIGIN, AND BASIS VECTORS
C
       READ(80,*,END=110)(NGRID(J), J=1,3)
       IF (ITYPE .EQ. 1) THEN
        DO I=1,4
         READ(80,*,END=110)(RBAS(J,I), J=1,3)
        END DO
       END IF
       DO I=1,3
        IF (NGRID(I) .LT. 1) THEN
         PRINT *,'WFGRID: NUMBER OF GRID POINTS MUST BE >= 1'
         GOTO 920
        END IF
       END DO
C
C READ DATA FOR WAVEFUNCTIONS
C PRINT OUT THE HOMO DENSITY BY DEFAULT
C
       READ(80,*,END=210) NWAVF
       IF (NWAVF.GT.NTOTWF) THEN
        WRITE(6,*)'WFGRID: CHANGE NTOTWF TO ',NWAVF
        WRITE(6,*)'AND RECOMPILE AND RERUN'
        WRITE(6,*) 'TILL THEN, BYE BYE FROM WGRID'
        RETURN
       END IF
       GO TO 220
 210   WRITE(80,'(I3,A)') NWAVF, '    NUMBER OF ORBITALS '
 220   READ(80,240,END=310) BAND
 240   FORMAT(A20)
       GO TO 320
 310   BAND='VAL' 
       INDWRTHOMO(NWAVF)=0
       WRITE(80,'(A)') 'VAL    IF VALENCE GIVE 
     &   INDEX WITH RESPECT TO HOMO'
       WRITE(80,'(I3,2X,A)') INDWRTHOMO(NWAVF), 
     &   '   FOR CORE STATES GIVE SPIN, REPRESENATION, BASIS INDEX'
       GO TO 230
 320   CONTINUE

       IF(BAND(1:3).NE.'VAL') THEN
       NUNSYM=0
       ISPO=0
       IRPO=0
       IBSO=0
       DO I=1,NWAVF
        READ(80,*,END=910) ISPN,IREP,IBAS
        IF ((ISPN .LT. 1) .OR. (ISPN .GT. NSPN)) THEN
         PRINT *,'WFGRID: SPIN FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        IF ((IREP .LT. 1) .OR. (IREP .GT. N_REP)) THEN
         PRINT *,'WFGRID: REPRESENTATION FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        IF ((IBAS .LT. 1) .OR. (IBAS .GT. NS_TOT(IREP))) THEN
         PRINT *,'WFGRID: INDEX FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        NDIM=NDMREP(IREP)
        NUNSYM=NUNSYM+NDIM
        IF (NUNSYM .GT. MAX_OCC) THEN
         PRINT *,'WFGRID: NUMBER OF STATES IS TOO LARGE'
         PRINT *,'        INCREASE MAX_OCC TO AT LEAST: ',NUNSYM
         GOTO 920
        END IF


C
C CHECK FOR ORDERING
C
        IF (ISPN .LT. ISPO) GOTO 10
        IF (ISPN .GT. ISPO) THEN
         ISPO=ISPN
         IRPO=0
        END IF
        IF (IREP .LT. IRPO) GOTO 10
        IF (IREP .GT. IRPO) THEN
         IRPO=IREP
         IBSO=0
        END IF
        IF (IBAS .LT. IBSO) GOTO 10
        IBSO=IBAS
        INFOWF(1,I)=ISPN
        INFOWF(2,I)=IREP
        INFOWF(3,I)=IBAS
        INFOWF(4,I)=NDIM
       END DO
       GOTO 30
C
C WRONG ORDERING
C
   10  PRINT *,'WFGRID: STATES MUST BE ORDERED ACCORDING TO: '
       PRINT *,'      * SPIN           (HIGHEST PRIORITY)'
       PRINT *,'      * REPRESENTATION (SECOND HIGHEST PRIORITY)'
       PRINT *,'      * INDEX          (LOWEST PRIORITY)'
       GOTO 920

       ELSE
         NUNSYM=0
         READ(80,*,END=910) (INDWRTHOMO(I),I=1,NWAVF)
 230     CALL FINDSTATE(INDWRTHOMO,NWAVF,ISP,IRP,IBS,IERR)
         DO I=1,NWAVF
           WRITE(6,*)'STATE :', ISP(I),IRP(I),IBS(I)
           NDIM=NDMREP(IRP(I))
          IF ((ISP(I) .LT. 1) .OR. (ISP(I) .GT. NSPN)) THEN
           PRINT *,'WFGRID: SPIN FOR STATE ',I,' IS INVALID'
           GOTO 920
          END IF
          IF ((IRP(I).LT. 1) .OR. (IRP(I).GT. N_REP)) THEN
           PRINT *,'WFGRID: REPRESENTATION FOR STATE ',I,' IS INVALID'
           GOTO 920
          END IF
          IF ((IBS(I) .LT. 1) .OR. (IBS(I) .GT. NS_TOT(IRP(I)))) THEN
           PRINT *,'WFGRID: INDEX FOR STATE ',I,' IS INVALID'
           GOTO 920
          END IF
          NDIM=NDMREP(IRP(I))
          NUNSYM=NUNSYM+NDIM
          IF (NUNSYM .GT. MAX_OCC) THEN
           PRINT *,'WFGRID: NUMBER OF STATES IS TOO LARGE'
           PRINT *,'        INCREASE MAX_OCC TO AT LEAST: ',NUNSYM
           GOTO 920
          END IF
           INFOWF(1,I)=ISP(I)
           INFOWF(2,I)=IRP(I)
           INFOWF(3,I)=IBS(I)
           INFOWF(4,I)=NDIM
         END DO
       END IF
         
************************************************************************
C
C DIAGONALIZE FOR EACH SPIN AND REPRESENTATION 
C STORE RELEVANT STATES IN PSI_COEF 
C
C  THIS PART IS SKIPPED, THE STORED WF ARE USED -- TB
   30  PRINT '(A)','SETUP OF WAVEFUNCTIONS ... '
c       GO TO 200
c
c      IWOFS=0
c      NWF=0
c      DO 200 ISPN=1,NSPN
c        NWFS(ISPN)=0
c        IF (IWOFS .GE. NWAVF) GOTO 190
c        IF (INFOWF(1,IWOFS+1) .GT. ISPN) GOTO 190
c        CALL OVERLAP(1)
c        CALL OVERLAP(2)
C
C LOOP OVER REPRESENTATIONS
C GET MATRIX ELEMENTS AND DIAGONALIZE
C
c        KBAS=0
c        DO 100 IREP=1,N_REP
c         NBAS=NS_TOT(IREP)
c         IF ((N_OCC(IREP,ISPN) .LT. 1) .OR. (NBAS .LT. 1)) THEN
c          KBAS=KBAS+(NBAS*(NBAS+1))/2
c          GOTO 90
c         END IF
c         IF (NBAS .GT. NDH) THEN
c          PRINT *,'WFGRID: NDH MUST BE AT LEAST: ',NBAS
c          CALL STOPIT
c         END IF
c         DO IBAS=1,NBAS
c          DO JBAS=IBAS,NBAS
c           KBAS=KBAS+1
c           OVER(JBAS,IBAS)=HSTOR(KBAS,1)
c           HAM (JBAS,IBAS)=HSTOR(KBAS,2)
c          END DO
c         END DO
c         CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
c
C
C SAVE CORRECT EIGENSTATES
C
c         DO I=1,N_OCC(IREP,ISPN)
c          INDX=INFOWF(3,IWOFS+I)
c          EVLOCC(IWOFS+I)=EVAL(INDX)
c          OCCUPANCY(IWOFS+I)= 1.0D0
c          DO IB=1,NBAS
c           PSI_COEF(IB,I,IREP,ISPN)=HAM(IB,INDX)
c          END DO
c          PRINT '(A,I1,A,I2,A,I5,A,F15.6)',
c     &          'SPIN ',ISPN,', REPRESENTATION ',IREP,', STATE ',
c     &           INDX,', EIGENVALUE: ',EVAL(INDX)
c         END DO
c         IWOFS=IWOFS+N_OCC(IREP,ISPN)
c         NWFS(ISPN)=NWFS(ISPN)+N_OCC(IREP,ISPN)*NDMREP(IREP)
c   90    CONTINUE
cc  100   CONTINUE
c  190   CONTINUE
c        NWF=NWF+NWFS(ISPN)
c
c  200  CONTINUE

        NWF=0
        DO ISPN=1,NSPN
         NWF=NWF+NWFS(ISPN)
        END DO
C
C OPEN OUTPUT FILE
C
       IF  (BAND(1:3).EQ.'VAL') THEN
          FILE='WFHOMO'
       ELSE
          FILE='WFCORE'
       END IF

       DO I=1,NWAVF
        IF(BAND(1:3).EQ.'VAL') THEN
           II=INDWRTHOMO(I)
        ELSE
           II=I
        END IF
        IF (II.LT.0) THEN
           WRITE(FILENAME(I),'(A,I3.2)') FILE, II
        ELSE
           WRITE(FILENAME(I),'(A,I2.2)') FILE, II
        END IF

c       OPEN(90,FILE='WFGROUT',FORM=FORMSTR,STATUS='UNKNOWN')
c       REWIND(90)
C
C WRITE FILE HEADER
C
       IUNIT=90+I
       OPEN(IUNIT,FILE=FILENAME(I),FORM=FORMSTR,STATUS='UNKNOWN')
       REWIND(IUNIT)
       IF (IFORM .EQ. 1) THEN
        WRITE(IUNIT) ITYPE,NSPN
        WRITE(IUNIT)(NGRID(J), J=1,3),MPBLOCK
        IF (ITYPE .EQ. 1) THEN
         WRITE(IUNIT)((RBAS(J,K), J=1,3), K=1,4)
        END IF
        WRITE(IUNIT) NWAVF
        WRITE(IUNIT)(INFOWF(J,I), J=1,4),EVLOCC(I)
       ELSE
        JSP=INFOWF(1,I)
        JRP=INFOWF(2,I)
        JST=INFOWF(3,I)
        WRITE(IUNIT,*) 'ORBITAL DENSITY FOR WF'
        WRITE(IUNIT,'(A,I1,A,I2,A,I5)')'SPIN ',JSP,
     &     ' REPRESENTATION ',JRP, ' STATE ',JST
         OPEN(77,FILE='XMOL.DAT')
         REWIND(77)
         READ(77,*) NATOM
         READ(77,*)
         WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
         DO K=1,3
         WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(K),(RBAS(J,K+1),J=1,3)
         ENDDO
         DO K=1,NATOM
           READ(77,*)IZ, X, Y, Z
           CHR=REAL(IZ)
           WRITE(IUNIT,2002)IZ, CHR, X, Y, Z
         END DO
         CLOSE(77)
 2002    FORMAT(I6,4F16.10)
       END IF
       END DO
C
C LOOP OVER ALL POINTS IN BLOCKS
C
       ITOTWF=0
       DO IW=1,NWAVF
        K_REP=INFOWF(2,IW)
        NDM=NDMREP(K_REP)
        ITOTWF=ITOTWF+NDM
       END DO

       DO 800 IX=1,NGRID(1)
        DO 790 IY=1,NGRID(2)
         NMSH=NGRID(3)
         DO 780 IOFS=0,NMSH-1,MPBLOCK
          MPTS=MIN(MPBLOCK,NMSH-IOFS)
C
C SETUP GRID POINTS AND INITIALIZE PSIG
C
          IF (ITYPE .EQ. 1) THEN
           DO IPTS=1,MPTS
            IZ=IOFS+IPTS
            DO I=1,3
             RGRID(I,IPTS)=RBAS(I,1)+(IX-1)*RBAS(I,2)
     &                    +(IY-1)*RBAS(I,3)+(IZ-1)*RBAS(I,4)
            END DO
           END DO
          ELSE
           DO IPTS=1,MPTS
            READ(80,*,END=900)(RGRID(I,IPTS), I=1,3)
           END DO
          END IF
          DO IWF=1,MAX_OCC
           DO IPTS=1,MPTS
            PSIG(IPTS,IWF)=0.0D0
           END DO
          END DO
C
C LOOP OVER ALL FUNCTION SETS, THEIR POSITIONS, EQUIVALENT SITES
C
          ISHELLA=0
          DO IFNCT=1,NFNCT
           LMAX1=LSYMMAX(IFNCT)+1
           DO I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
            DO J_POS=1,M_NUC
             CALL WFRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &           RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
             IF(L_NUC .NE. M_NUC)THEN
              PRINT *,'WFGRID: PROBLEM IN UNRAVEL'
              CALL STOPIT
             END IF
C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
             DO JOFS=0,MPTS-1,NSPEED
              NPV=MIN(NSPEED,MPTS-JOFS)
              DO IPV=1,NPV
               PTS(IPV,1)=RGRID(1,JOFS+IPV)-RVECA(1,J_POS)
               PTS(IPV,2)=RGRID(2,JOFS+IPV)-RVECA(2,J_POS)
               PTS(IPV,3)=RGRID(3,JOFS+IPV)-RVECA(3,J_POS)
              END DO
C
C GET VALUE OF BASIS FUNCTIONS
C
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
C
C UPDATE PSIG
C

              IF (IUPDAT) THEN
              DO IWF=1,ITOTWF
               ILOC=0
               DO LI=1,LMAX1
                DO MU=1,ISIZE(LI)
                 DO ICON=1,N_CON(LI,IFNCT)
                  ILOC=ILOC+1
c                  IF (ICOUNT(ICON,LI)) THEN
                    FACTOR=PSI(ILOC,IWF,1)
                    DO IPV=1,NPV
                     PSIG(JOFS+IPV,IWF)=PSIG(JOFS+IPV,IWF)
     &                                 +FACTOR*GRAD(IPV,1,MU,ICON,LI)
                    END DO
c                  END IF
                   END DO
                 END DO
                END DO
               END DO

              END IF
             END DO
            END DO
           END DO
          END DO
C
C GET ORBITAL DENSITY FROM WAVEFUNCTION 
C

          IWF=0
          DO IW=1,NWAVF
            ISPN=INFOWF(1,IW)
            IREP=INFOWF(2,IW)
            IBAS=INFOWF(3,IW)
            NDIM=NDMREP(IREP)
            IWF=IWF+1
             DO IPTS=1,MPTS
              PSIG(IPTS,IW)=PSIG(IPTS,IWF)**2
             END DO
             DO IDIM=2,NDIM
              IWF=IWF+1
              DO IPTS=1,MPTS
               PSIG(IPTS,IW)=PSIG(IPTS,IW)
     &                              +PSIG(IPTS,IWF)**2
              END DO
             END DO
           END DO
           IF(IWF.GT.ITOTWF) THEN
            WRITE(6,*) 'WFGRID: THE TOTAL NUMBER OF STATES DO NOT MATCH'
            WRITE(6,*)  IWF, ' VS ', NWF
            CALL STOPIT
           END IF
C
C WRITE OUTPUT
C
                DO IWF=1,NWAVF
                 IUNIT=90+IWF

          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(IUNIT)(PSIG(IPTS,IWF)*FACT, IPTS=1,MPTS)
           ELSE
             WRITE(IUNIT,'(3(1X,E20.12))')(PSIG(IPTS,IWF)*FACT,
     &        IPTS=1,MPTS)
           END IF
          ELSE
           DO IPTS=1,MPTS
            WRITE(IUNIT,'(3(1X,F20.12))')(RGRID(I,IPTS), I=1,3)
            WRITE(IUNIT,'(3(1X,E20.12))')(PSIG(IPTS,IWF))
           END DO
          END IF

c                   GO TO 770
c                  END IF
c                END DO
c             END DO
c  770        CONTINUE

           END DO
  780    CONTINUE
  790   CONTINUE
  800  CONTINUE
       CLOSE(80)
       DO I=1,NWAVF
         IUNIT=90+I
         CLOSE(IUNIT)
       END DO
       NMSH=NSTORE
       RETURN
C
C ERROR HANDLING
C
  900  CLOSE(80)
  910  PRINT *,'WFGRID: ERROR READING FILE WFGRID' 
  920  CLOSE(90)
       NMSH=NSTORE
       RETURN
      END


********************************************************************
C
C      THIS ROUTINE FINDS SPIN, REPRESENTATION AND THE INDEX OF THE BASIS
C      HOMO OF THE MOLECULE AND THE HOMO+INDEX STATES.
C                                          TB 04/03
C
      SUBROUTINE FINDSTATE_old(INDEX,N,IS,IR,IB)
      INCLUDE 'PARAMS' 
      INCLUDE 'commons.inc' 
      DIMENSION INDEX(N),IS(N),IR(N),IB(N)
      DIMENSION ISP(MAX_OCC),IRP(MAX_OCC),IBS(MAX_OCC),EVL(MAX_OCC)
      
        WRITE(6,*)'FINDSTATE:' ,NSPN,N_REP
        IWF=0
        DO ISPIN=1,NSPN
          DO IREPR=1,N_REP
            JVIRT=0
            WRITE(6,*)'N_OCC', IREPR,ISPIN,N_OCC(IREPR,ISPIN)
            DO IVIRT=1,N_OCC(IREPR,ISPIN)
               DO IDEG=1,NDMREP(IREPR)
                IWF=IWF+1
                JVIRT=JVIRT+1
                EVL(IWF)=EVLOCC(IWF)
                ISP(IWF)=ISPIN
                IRP(IWF)=IREPR
                IBS(IWF)=JVIRT
               END DO
            END DO
          END DO
        END DO
        WRITE(6,*) IWF
        DO I=1,IWF
          DO J=1,IWF
           IF(EVL(J).GT.EVL(I)) THEN  
             CALL SWAP(EVL(I),EVL(J))
             CALL ISWAP(ISP(I),ISP(J))
             CALL ISWAP(IRP(I),IRP(J))
             CALL ISWAP(IBS(I),IBS(J))
           END IF
          END DO
        END DO
        CALL FINDHOMO(IWF,EVL,IWFHOMO) 
        DO I=1,N
           KWF=IWFHOMO+INDEX(I)
           IS(I)=ISP(KWF)
           IR(I)=IRP(KWF)
           IB(I)=IBS(KWF)
        END DO
      RETURN
      END
     
********************************************************************

      SUBROUTINE FINDHOMO(IWF,EVL,IWFHOMO)
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      DIMENSION EVL(IWF) 
      EF=EFERMI(1)
      IF(NSPN.EQ.2) EF=MAX(EF,EFERMI(NSPN))
      DIFF=1.0D30
       DO I=1,IWF
         DIFFNEW=EF-EVL(I)
         IF(DIFFNEW.GE.-1.0D-4) THEN
           IF(DIFFNEW.LT.DIFF) THEN
            DIFF=DIFFNEW
            IWFHOMO=I
           END IF 
         END IF 
       END DO 
      RETURN
      END

********************************************************************
C
C      THIS ROUTINE FINDS SPIN, REPRESENTATION AND THE INDEX OF THE BASIS
C      HOMO OF THE MOLECULE AND THE HOMO+INDEX STATES.
C                                          TB 04/03
C
      SUBROUTINE FINDSTATE(INDEX,N,IS,IR,IB,IERR)
      INCLUDE 'PARAMS' 
      INCLUDE 'commons.inc' 
      LOGICAL EXIST
      CHARACTER*20 TRASH
      DIMENSION INDEX(N),IS(N),IR(N),IB(N),EV(MXSPN*MAX_OCC)
      DIMENSION ISP(MXSPN*MAX_OCC),IRP(MXSPN*MAX_OCC),
     &           IBS(MXSPN*MAX_OCC),EVL(MXSPN*MAX_OCC)
 
        IERR=0
        INQUIRE(FILE='EVALUES',EXIST=EXIST) 
        IF(EXIST) THEN
          OPEN(90,FILE='EVALUES',STATUS='OLD')
          REWIND(90)
        ELSE
          WRITE(6,*)'EVALUES DOES NOT EXIST, CALLING'
          WRITE(6,*)'OLD ROUTINE TO FIND HOMO AND OCCUPIED STATES'
          WRITE(6,*)'NEED EVALUES TO FIND LUMO'
          IERR=1
          RETURN
        END IF
        IWF=0

        DO 300 ISPIN=1,NSPN
          READ(90,240) TRASH
          DO 200 IREPR=1,N_REP
            IF(NS_TOT(IREPR).LE.0) GO TO 200
            READ(90,240) TRASH
            READ(90,*) NDIM, NBASE
             JVIRT=0
             READ(90,*) (EV(IVIRT),IVIRT=1,NBASE)
              DO IV=1,NBASE
                IWF=IWF+1
                JVIRT=JVIRT+1
                EVL(IWF)=EV(IV)
                ISP(IWF)=ISPIN
                IRP(IWF)=IREPR
                IBS(IWF)=JVIRT
              END DO
 200      CONTINUE
 300    CONTINUE
 240    FORMAT(A20)
        DO I=1,IWF
          DO J=1,IWF
           IF(EVL(J)-EVL(I).GT.1.0D-6) THEN  
             CALL SWAP(EVL(I),EVL(J))
             CALL ISWAP(ISP(I),ISP(J))
             CALL ISWAP(IRP(I),IRP(J))
             CALL ISWAP(IBS(I),IBS(J))
           END IF
          END DO
        END DO
        CLOSE(90)
        CALL FINDHOMO(IWF,EVL,IWFHOMO) 
        DO I=1,N
            KWF=IWFHOMO+INDEX(I)
            IS(I)=ISP(KWF)
            IR(I)=IRP(KWF)
            IB(I)=IBS(KWF)
        END DO
      RETURN
      END

C
C******************************************************************
C
       SUBROUTINE WFRAVEL(IFNCT,ISHELLA,I_SITE,RVEC,RVECI,N_NUC,ILOC,
     &            INFO,NW,ITOTW)
C ORIGINALLY WRITTEN BY MARK R PEDERSON (1985)
C SLIGHTLY MODIFIED TO INCLUDE UNOCCUPIED STATES FOR THE PURPOSE OF CALCULATION 
C OF ORBITAL DENSITY               - TB 04/03
C
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       LOGICAL WF
       DIMENSION NDEG(3),IND_SALC(ISMAX,MAX_CON,3)
       DIMENSION RVECI(3,MX_GRP),RVEC(3)
       DIMENSION ISHELLV(2)
       DIMENSION INFO(4,NW)
       DATA NDEG/1,3,6/
       DATA ICALL1,ICALL2,ISHELL/0,0,0/
C
       IF(I_SITE.EQ.1)THEN
        ICALL1=ICALL1+1
        CALL GTTIME(TIMER1)
        CALL OBINFO(1,RVEC,RVECI,N_NUC,ISHELLV(ILOC))
        CALL GSMAT(ISHELLV(ILOC),ILOC)
        CALL GTTIME(TIMER2)
        T1UNRV=T1UNRV+TIMER2-TIMER1
        IF (DEBUG.AND.(1000*(ICALL1/1000).EQ.ICALL1)) THEN
         PRINT*,'WASTED1=',T1UNRV,ISHELL
         PRINT*,'ICALL2,AVERAGE:',ICALL1,T1UNRV/ICALL2
        END IF
       END IF
       CALL GTTIME(TIMER1)
       ISHELL=ISHELLV(ILOC)
C
C UNSYMMETRIZE THE WAVEFUNCTIONS....
C
       IWF=0
       DO IW=1,NW
       DO 1020 ISPN=1,NSPN
        KSALC=0
        DO 1010 K_REP=1,N_REP
C
C CALCULATE ARRAY LOCATIONS:
C
         DO 5 K_ROW=1,NDMREP(K_REP)
          KSALC=KSALC+1
    5    CONTINUE
         INDEX=INDBEG(ISHELLA,K_REP)
         DO 20 LI =0,LSYMMAX(IFNCT)
          DO 15 IBASE=1,N_CON(LI+1,IFNCT)
           DO 10 IQ=1,N_SALC(KSALC,LI+1,ISHELL)
            INDEX=INDEX+1
            IND_SALC(IQ,IBASE,LI+1)=INDEX
   10      CONTINUE
   15     CONTINUE
   20    CONTINUE
C
C END CALCULATION OF SALC INDICES FOR REPRESENTATION K_REP
C
C         DO 1000 IOCC=1,N_OCC(K_REP,ISPN)

         DO 1000 IOCC=1,NS_TOT(K_REP)
          WF=.FALSE.
              IF(INFO(1,IW).EQ.ISPN) THEN
               IF (INFO(2,IW).EQ.K_REP) THEN
                IF(INFO(3,IW).EQ.IOCC) THEN
                  WF=.TRUE.
                  GO TO 112
                END IF
               END IF
              END IF
 112      IF(.NOT.WF) GO TO 1000
    
          I_SALC=KSALC-NDMREP(K_REP)
          DO 950 IROW=1,NDMREP(K_REP)
           I_SALC=I_SALC+1
           IWF=IWF+1
           I_LOCAL=0
           DO 900 LI=0,LSYMMAX(IFNCT)
            DO 890 MU=1,NDEG(LI+1)
             IMS=MU+NDEG(LI+1)*(I_SITE-1)
             DO 880 IBASE=1,N_CON(LI+1,IFNCT)
              I_LOCAL=I_LOCAL+1
              PSI(I_LOCAL,IWF,ILOC)=0.0D0
              IQ_BEG=IND_SALC(1,IBASE,LI+1)-1
              DO 800 IQ=1,N_SALC(KSALC,LI+1,ISHELL)
               PSI(I_LOCAL,IWF,ILOC)=PSI(I_LOCAL,IWF,ILOC)+
     &         PSI_COEF(IQ+IQ_BEG,IOCC,K_REP,ISPN)*
     &         U_MAT(IMS,IQ,I_SALC,LI+1,ILOC)
  800         CONTINUE
  880        CONTINUE
  890       CONTINUE
  900      CONTINUE
           IF(I_LOCAL.GT.MAXUNSYM)THEN
            PRINT*,'UNRAVEL: MAXUNSYM MUST BE AT LEAST:',I_LOCAL
            CALL STOPIT
           END IF
            FACTOR=SQRT(1.0D0/REAL(NDMREP(K_REP)))
           DO 25 J_LOCAL=1,I_LOCAL
            PSI(J_LOCAL,IWF,ILOC)=FACTOR*PSI(J_LOCAL,IWF,ILOC)
   25      CONTINUE
  950     CONTINUE
 1000    CONTINUE
 1010   CONTINUE
 1020  CONTINUE
           END DO
C
       IF (IWF.GT.ITOTW) THEN
        PRINT *,'UNRAVEL: BIG BUG: NUMBER OF STATES IS INCORRECT'
        PRINT *,'IWF,ITOTWF:',IWF,ITOTW
        CALL STOPIT
       END IF
       CALL GTTIME(TIMER2)
       T2UNRV=T2UNRV+(TIMER2-TIMER1)
       ICALL2=ICALL2+1
       IF (DEBUG.AND.(1000*(ICALL2/1000).EQ.ICALL2)) THEN
        PRINT*,'WASTED2:',T2UNRV
        PRINT*,'ICALL2,AVG:',ICALL2,T2UNRV/ICALL2
       END IF
       RETURN
       END
C
C *****************************************************************
       SUBROUTINE POTRHOSER 
C
C POTRHOGRID VERSION DIRK POREZAG JUNE 1998.
C CALCULATION OF CHARGES WITHIN A SPHERE moved into atomsph JK 3/99
C * DENSITY AND POTENTIAL ON A GRID OF POINTS
C
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MAXANG=200)
C
       LOGICAL ICOUNT,EXIST
       CHARACTER*20 FNAME(2),FORMSTR
       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
     &  ,RDIS(NSPEED),VLOC(NSPEED)
     &  ,XRAD(MAXRAD),WTRAD(MAXRAD),CENTER(5,MAXSPH)
     &  ,ANGLE(3,MAXANG),DOMEGA(MAXANG)
     &  ,RVECA(3,MX_GRP),PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3)
       DIMENSION NGRID(3),DERIV(3)
       DIMENSION RBAS(3,4)
C
C MODE:  1 ... CALCULATE DENSITY ON A MESH OF POINTS
C        2 ... CALCULATE POTENTIAL ON A MESH OF POINTS
C
C FOR DENSITY EVALUATIONS, USE OLD SCHEME FOR DENSITY CALCULATION
C OTHERWISE, GET DENSITY FROM COUPOT
C
       OPEN(65,FILE='SYMBOL')
       READ(65,65)FORMSTR
       CLOSE(65)
 65    FORMAT(A)
             IRET=0
             DO I=1,17
             IF(FORMSTR(I:I+2).EQ.'SCF')IRET=1
             END DO
             IF(IRET.EQ.0)THEN
             PRINT*,'RETURNING FROM POTRHOGRD/POTRHOSER'
             PRINT*,'MUST USE SCF-ONLY MODE IN SYMBOL'
             RETURN
             END IF
       PRINT '(A)',' '
       PRINT '(A)','POTENTIAL/DENSITY GRIDS'
       CALL GTTIME(TIME1)
       DO 900 MODE=1,2
C
C READ IN NECESSARY INPUT DATA
C
        IF (MODE.EQ.1) THEN
         FNAME(1)='RHOGRID'
        ELSE
         FNAME(1)='POTGRID'
        END IF
        INQUIRE(FILE=FNAME(1),EXIST=EXIST)
C
C DETERMINE IF THERE IS ANYTHING TO DO
C
        FORMSTR= ' '
        IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
        IF (MODE .EQ. 1) PRINT '(2A)','DENSITY GRID  ',FORMSTR
        IF (MODE .EQ. 2) PRINT '(2A)','POTENTIAL GRID',FORMSTR
        IF (.NOT.EXIST) GOTO 900
C
C READ INPUT DATA
C
        OPEN(72,FILE=FNAME(1),FORM='FORMATTED',STATUS='OLD')
        REWIND(72)
C
C CHECK IF THIS FILE IS EMPTY
C IF YES, CREATE A DEFAULT ONE
C
        I=1
        READ(72,*,END=10) ITYPE,IFORM
        I=0
        REWIND(72)
   10   IF (I .EQ. 1) THEN
         DO I=1,3
          RBAS(I,1)=  1.0D30
          RBAS(I,2)= -1.0D30
         END DO
         DO ICNT=1,NCNT
          DO I=1,3
           RBAS(I,1)= MIN(RBAS(I,1),RCNT(I,ICNT))
           RBAS(I,2)= MAX(RBAS(I,2),RCNT(I,ICNT))
          END DO
         END DO
         DO I=1,3
          RBAS(I,1)= RBAS(I,1)-5.0D0
          RBAS(I,2)= RBAS(I,2)+5.0D0
          NGRID(I)= (RBAS(I,2)-RBAS(I,1))/0.5D0+2
          RBAS(I,2)= (RBAS(I,2)-RBAS(I,1))/(NGRID(I)-1)
         END DO
         REWIND(72)
         WRITE(72,1010) 1,2,'Grid mode, formatted file'
         WRITE(72,1020) (NGRID(I), I=1,3),    'Number of grid points'
         WRITE(72,1030) (RBAS(I,1), I=1,3),   'Origin'
         WRITE(72,1030) RBAS(1,2),0.0D0,0.0D0,'Basis vector 1'
         WRITE(72,1030) 0.0D0,RBAS(2,2),0.0D0,'Basis vector 1'
         WRITE(72,1030) 0.0D0,0.0D0,RBAS(3,2),'Basis vector 1'
 1010    FORMAT(2(I6,1X),8X,A)
 1020    FORMAT(3(I6,1X),1X,A)
 1030    FORMAT(3(F7.3,1X),1X,A)
         CLOSE(72)
         OPEN(72,FILE=FNAME(1),FORM='FORMATTED',STATUS='OLD')
         REWIND(72)
        END IF
C
C GRID INPUT
C
        READ(72,*,END=880) ITYPE,IFORM
        IF (ITYPE .GT. 2) ITYPE= 2
        IF (ITYPE .LT. 1) ITYPE= 1
        IF (IFORM .GT. 2) IFORM= 2
        IF (IFORM .LT. 1) IFORM= 1
        IF (ITYPE .EQ. 2) IFORM= 2
        FORMSTR='UNFORMATTED'
        IF (IFORM .EQ. 2) FORMSTR='FORMATTED'
C
C GET NUMBER OF GRID POINTS, ORIGIN, AND BASIS VECTORS
C
        READ(72,*,END=880)(NGRID(J), J=1,3)
        write(6,*)(NGRID(J), J=1,3)
        IF (ITYPE .EQ. 1) THEN
         DO I=1,4
          READ(72,*,END=880)(RBAS(J,I), J=1,3)
         END DO
        END IF
        DO I=1,3
         IF (NGRID(I) .LT. 1) THEN
          PRINT *,'POTRHOGRID: NUMBER OF GRID POINTS MUST BE >= 1'
          GOTO 890
         END IF
        END DO
        IF (NGRID(3) .GT. MAX_PTS) THEN
         PRINT *,'POTRHOGRID: MAX_PTS MUST BEAT LEAST: ',NGRID(3)
         PRINT *,'SKIPPING GRID EVALUATION FOR MODE ',MODE
         GOTO 890
        END IF
        NLOOP=NGRID(1)*NGRID(2)
        IF(NSPN.EQ.1) THEN
        FNAME(1)=FNAME(1)(1:3)//'GROUT'
        ELSE
         IF (MODE.EQ.1) THEN
           FNAME(1)=FNAME(1)(1:3)//'TOT'
           FNAME(2)=FNAME(1)(1:3)//'SPN'
         ELSE
           FNAME(1)=FNAME(1)(1:3)//'MAJ'
           FNAME(2)=FNAME(1)(1:3)//'MIN'
         ENDIF
        ENDIF
C
C OPEN OUTPUT FILES, WRITE HEADER
C
        DO IS=1,NSPN
        IUNIT=73+IS
   20   OPEN(IUNIT,FILE=FNAME(IS),FORM=FORMSTR,STATUS='UNKNOWN')
        REWIND(IUNIT)
        IF (IFORM .EQ. 1) THEN
         WRITE(IUNIT) ITYPE,NSPN
         WRITE(IUNIT)(NGRID(J), J=1,3),NGRID(3)
         IF (ITYPE .EQ. 1) THEN
          WRITE(IUNIT)((RBAS(J,I), J=1,3), I=1,4)
         END IF
         WRITE(IUNIT) NSPN
         DO ISPN=1,NSPN
          X=ISPN
          WRITE(IUNIT) ISPN,ISPN,ISPN,ISPN,X
         END DO
        ELSE
         WRITE(IUNIT,*)'CLUSTER OUTPUT'
         IF(MODE.EQ.1) THEN
             IF(NSPN.EQ.1) WRITE(IUNIT,*)'SCF TOTAL DENSITY (ANG)'
             IF(NSPN.EQ.2) WRITE(IUNIT,*)'SCF DENSITY (ANG)'
         ELSE
             WRITE(IUNIT,2001)'SCF POTENTIAL'
         END IF
         OPEN(77,FILE='XMOL.DAT')
         READ(77,*) NATOM
         READ(77,*)
         WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
         DO I=1,3
         WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(I),(RBAS(J,I+1),J=1,3)
         ENDDO
         DO I=1,NATOM
           READ(77,*)IZ, X, Y, Z
           CHR=REAL(IZ)
           WRITE(IUNIT,2002)IZ, CHR, X, Y, Z
         END DO
         CLOSE(77)
c         WRITE(IUNIT,'(4(1X,I10))')(NGRID(J), J=1,3),NGRID(3)
c         IF (ITYPE .EQ. 1) THEN
c          WRITE(IUNIT,'(3(1X,F20.12))')((RBAS(J,I), J=1,3), I=1,4)
c         END IF
        ENDIF
        ENDDO

2001    FORMAT(A25)
2002    FORMAT(I6,4F16.10)
C
C LOOP FOR EACH PILE
C
        DO 850 ILOOP=1,NLOOP
C
C SETUP POINTS
C
c         NMSH=NGRID(3)
         NPTS=NGRID(3)
         IF (ITYPE .EQ. 1) THEN
          IY=MOD(ILOOP-1,NGRID(2))
          IX=(ILOOP-1)/NGRID(2)
          DO IZ=1,NPTS
           I=IZ-1
           RMSH(1,IZ)=RBAS(1,1)+IX*RBAS(1,2)+IY*RBAS(1,3)+I*RBAS(1,4)
           RMSH(2,IZ)=RBAS(2,1)+IX*RBAS(2,2)+IY*RBAS(2,3)+I*RBAS(2,4)
           RMSH(3,IZ)=RBAS(3,1)+IX*RBAS(3,2)+IY*RBAS(3,3)+I*RBAS(3,4)
           WMSH(IZ)=0.0D0
          END DO
         ELSE
          DO IZ=1,NPTS
           READ(72,*,END=870)(RMSH(I,IZ), I=1,3)
           WMSH(IZ)=0.0D0
          END DO
         END IF
C
C NOW: CALCULATE ELECTRONIC COULOMB POTENTIAL AND/OR DENSITY
C DENSITY WILL BE STORED IN RHOG
C
         NGRAD=1
         MODDEN=2
         IF (MODE .EQ. 2) THEN
          MODDEN=1
          IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
          DO IAT=1,NIDENT
           GAUSS_CUT(IAT)=1.0D30
          END DO
c
c
          CALL COUPOT1
         ELSE
          I1=IGGA(1)
          I2=IGGA(2)
          IGGA(1)=0
          IGGA(2)=0
          CALL DENSOLD(VOL)
          IGGA(1)=I1
          IGGA(2)=I2
         END IF
C
C UPDATE DATA IN RHOG
C
         DO IGRAD=1,NGRAD
          DO IPTS=1,NPTS
           RHOG(IPTS,IGRAD,1)=RHOG(IPTS,IGRAD,1)+RHOG(IPTS,IGRAD,NSPN)
          END DO
         END DO
C
C DENSITY GRID
C
C        DENSITY IN  ANG^-3 FOR THE GAUSSIAN CUBE FORMAT
C
        FACTOR=0.5292D0**3
        IF (MODE .EQ. 1) THEN

         IF (NSPN.EQ.1) THEN
           IF (ITYPE .EQ. 1) THEN

            IF (IFORM .EQ. 1) THEN
             WRITE(74)(RHOG(IPTS,1,1), IPTS=1,NPTS)
            ELSE
             WRITE(74,9010)(RHOG(IPTS,1,1)*FACTOR, IPTS=1,NPTS)
            END IF
           ELSE

           DO IPTS=1,NPTS
             WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(74,9010) RHOG(IPTS,1,1)
            END DO
          END IF
         ELSE     !IF SPIN-POLARISED
          DO IS=1,NSPN
           IUNIT=73+IS
           IF (ITYPE .EQ. 1) THEN
            IF (IFORM .EQ. 1) THEN
             WRITE(IUNIT)(RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN), IPTS=1,NPTS)
            ELSE
            IF(IS.EQ.1) WRITE(IUNIT,9010) (RHOG(IPTS,1,1)
     &               *FACTOR,IPTS=1,NPTS)
            IF(IS.EQ.2) WRITE(IUNIT,9010) ((RHOG(IPTS,1,1)-2*
     &               RHOG(IPTS,1,NSPN)) *FACTOR,IPTS=1,NPTS)
            END IF
           ELSE
            DO IPTS=1,NPTS
             WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
             WRITE(74,9010) RHOG(IPTS,1,1),
     &                 RHOG(IPTS,1,1)-2*RHOG(IPTS,1,NSPN)
            END DO
            END IF
          END DO
         END IF
         GOTO 850
         END IF
C
C THE FOLLOWING PART IS ONLY DONE IF (MODE .EQ. 2)
C CALCULATING KOHN-SHAM POTENTIAL POT
C
         CALL GETVLXC(MAX_PTS,RHOG,POT,POTIN)
C
C ADD EFIELD POTENTIAL TO LOCAL POTENTIAL
C
         DO IPTS=1,NPTS
          CALL EXTPOT(RMSH(1,IPTS),EXHERE,DERIV)
          EFHERE=EFIELD(1)*RMSH(1,IPTS)+EFIELD(2)*RMSH(2,IPTS)
     &          +EFIELD(3)*RMSH(3,IPTS)
          POTIN(IPTS)=POTIN(IPTS)+EFHERE+EXHERE
         END DO
C
C WRITE OUTPUT: POTIN   CONTAINS LOCAL
C               COULOMB CONTAINS COULOMB
C               POT     CONTAINS EXCHANGE-CORRELATION
C
         IF (NSPN.EQ.1) THEN
          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(74)(POTIN(IPTS),COULOMB(IPTS),POT(IPTS),
     &                IPTS=1,NPTS)
           ELSE
            WRITE(74,9010)(POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS), IPTS=1,NPTS)
           END IF
          ELSE
           DO IPTS=1,NPTS
            WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
            WRITE(74,9010) POTIN(IPTS),COULOMB(IPTS),POT(IPTS)
           END DO
          END IF
         ELSE
          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(74)(POTIN(IPTS),COULOMB(IPTS),
     &                POT(IPTS),POT(IPTS+NPTS), IPTS=1,NPTS)
           ELSE
            WRITE(74,9010)(POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS),POT(IPTS+NPTS), IPTS=1,NPTS)
           END IF
          ELSE
           DO IPTS=1,NPTS
            WRITE(74,9020)(RMSH(I,IPTS), I=1,3)
            WRITE(74,9010) POTIN(IPTS),COULOMB(IPTS),
     &                     POT(IPTS),POT(IPTS+NPTS)
           END DO
          END IF
         END IF
  850   CONTINUE
        CLOSE(74)
        GOTO 890
C
C ERROR HANDLING
C
  870   CLOSE(74)
  880   PRINT *,'ERROR IN INPUT FILE, SKIPPING MODE ',MODE
  890   CLOSE(72)
  900  CONTINUE
       RETURN
 9010  FORMAT(3(1X,E20.12))
 9020  FORMAT(3(1X,F20.12))
       END
