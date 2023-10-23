C
C MINIMUM FUNCTIONALITY PREPROCESSOR FOR CONDITIONAL COMPILATION
C READS FROM STDIN AND WRITES TO STDOUT
C DIRK POREZAG, MAY 1998
C 
C CALL: condcomp -Ddef1 -Ddef2 ... 
C
C *****************************************************************
C
C SUPPLEMENTARY FUNCTION TO FIND NONBLANK LENGTH OF STRING
C
      INTEGER FUNCTION LENSTR(S)
       CHARACTER*300 S
       LENSTR= 0
       DO I= 300, 1, -1
        IF (S(I:I) .NE. ' ') THEN
         LENSTR= I
         RETURN
        END IF
       END DO
       RETURN
      END
C
C *****************************************************************
C
       PROGRAM CONDCOMP
        PARAMETER (MAXDEF=100)
        PARAMETER (MXNEST=100)
        LOGICAL       LERR,LPRT(MXNEST),LELS(MXNEST)
        CHARACTER*300 ARGUM,LINE,DEFLIST(MAXDEF)
        character*80 sInput,sOutputA
c       OPEN(50,FILE='fortran.f')
        NARG=IARGC()
        NDEF=0
        DO IARG=1,NARG
         CALL GETARG(IARG,ARGUM)
         ILEN= LENSTR(ARGUM)
         IF ((ARGUM(1:2) .EQ. '-D') .AND. (ILEN .GE. 3)) THEN
          NDEF= NDEF+1
          IF (NDEF .GT. MAXDEF) THEN
           PRINT*, 'COMPERR: MAXDEF MUST BE LARGER'
           STOP
          END IF
          DEFLIST(NDEF)=ARGUM(3:ILEN)
         ELSE
          PRINT*, 'COMPERR: OPTIONS MUST HAVE FORM -Dxxx'
          STOP
         END IF
        END DO
C
C PROCESSING
C
        NEST=0
        NOPR=0
        NLINE=0
  100   CONTINUE
         MODE= 0
         NLINE= NLINE+1
         READ(*,'(A300)',END=200,ERR=200) LINE
         IF (LINE(1:6) .EQ. '%ifdef')  MODE= 7
         IF (LINE(1:7) .EQ. '%ifndef') MODE= 8
         IF (LINE(1:7) .EQ. '%else')   MODE= 10
         IF (LINE(1:7) .EQ. '%endif')  MODE= 11
C
C PRINT
C        
         IF ((MODE .EQ. 0) .AND. (NOPR .EQ. 0)) THEN
          WRITE(6,'(A)') LINE(1:LENSTR(LINE))
         END IF
C
C DEAL WITH %ifdef AND %ifndef  
C
         IF ((MODE .EQ. 7) .OR. (MODE .EQ. 8)) THEN
          NEST= NEST+1
          IF (NEST .GT. MXNEST) THEN
           PRINT*, 'COMPERR: MXNEST MUST BE LARGER'
           STOP
          END IF
          I=MODE
  110     CONTINUE
           IF ((I .EQ. 300) .OR. (LINE(I:I) .NE. ' ')) GOTO 120
           I= I+1
           GOTO 110
  120     CONTINUE
          INDX=I
  130     CONTINUE
           IF ((I .EQ. 300) .OR. (LINE(I:I) .EQ. ' ')) GOTO 140
           I= I+1
           GOTO 130
  140     CONTINUE
          I= I-1
          ARGUM= LINE(INDX:I)
          IGOT=0
          ILEN= LENSTR(ARGUM)
          DO IDEF=1,NDEF
           JLEN= LENSTR(DEFLIST(IDEF))
           IF (ARGUM(1:ILEN) .EQ. DEFLIST(IDEF)(1:JLEN)) IGOT=IGOT+1
          END DO
          LELS(NEST)= .FALSE.
          IF (((MODE .EQ. 7) .AND. (IGOT .NE. 0)) .OR.
     &        ((MODE .EQ. 8) .AND. (IGOT .EQ. 0))) THEN
           LPRT(NEST)= .TRUE.
          ELSE
           LPRT(NEST)= .FALSE.
           NOPR=NOPR+1
          END IF
         END IF
C
C DEAL WITH %else
C
         IF (MODE .EQ. 10) THEN
          LERR= (NEST .LT. 1)
          IF (.NOT. LERR) LERR= LELS(NEST)
          IF (LERR) THEN
           PRINT*, 'COMPERR: EXTRA %else IN LINE ',NLINE
           STOP
          ELSE
           LELS(NEST)= .TRUE.
          END IF
          IF (LPRT(NEST)) THEN
           LPRT(NEST)= .FALSE.
           NOPR=NOPR+1
          ELSE
           LPRT(NEST)= .TRUE.
           NOPR= NOPR-1
          END IF
         END IF
C
C DEAL WITH %endif
C
         IF (MODE .EQ. 11) THEN
          IF (NEST .LT. 1) THEN
           PRINT*, 'COMPERR: EXTRA %endif IN LINE ',NLINE
           STOP
          END IF
          IF (.NOT. LPRT(NEST)) NOPR= NOPR-1
          NEST= NEST-1
         END IF
C
C CHECK IF NOPR IS >= 0
C
         IF (NOPR .LT. 0) THEN
          IF (MODE .EQ. 10) THEN
           PRINT*, 'COMPERR: EXTRA %else IN LINE ',NLINE
          ELSE
           PRINT*, 'COMPERR: EXTRA %endif IN LINE ',NLINE
          END IF
          STOP
         END IF
         GOTO 100
C
  200    CONTINUE
        close(6)
        END

