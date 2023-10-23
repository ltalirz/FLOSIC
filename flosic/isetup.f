        subroutine perfect_cluster
        implicit real*8 (a-h,o-z)
        dimension r(3),a(3),b(3),grp(3,3,121)
        character*100 xchgcp
        character*50 symm
c       r(1)=1.97800933760767    
c       r(2)=1.97800953225088    
c       r(3)=1.97800929452247   
       open(90,file='CLUSTER')
       read(90,'(a20)')xchgcp 
                 if(xchgcp(1:1).eq.'@')then
                 close(90)
                 return
                 end if
       read(90,'(a)')symm
       open(95,file='tmp')
       write(95,'(a20)')xchgcp
       write(95,'(a)')symm
       read(90,*)natms
       write(95,*)natms, 'No. of inequivalent atoms'
       do i=1,natms
         read(90,'(a100)')xchgcp
c         read(90,*)(r(j),j=1,3),nz
         do j=98,1,-1
         if(xchgcp(j:j+2).eq.'ALL')lb=j
         end do 
         do j=1,100
         if(xchgcp(j:j).ne.' ')le=j
         end do
         read(xchgcp,*)(r(j),j=1,3),nz
         call repair_cluster(r)
         write(95,98)(r(j),j=1,3),nz,xchgcp(lb:le)
98        format(3f20.14,i3,' ',A)
       end do
       read(90,*)chrg,spin
       write(95,99) chrg,spin
99     format(2f6.2)
       close(90)
       close(95)
       call system('cp CLUSTER CLUSTER_ORI')
       call system('mv tmp CLUSTER')
       end 
       subroutine repair_cluster(r)
        implicit real*8 (a-h,o-z)
        dimension r(3),a(3),b(3),g(3,3,121)
        open(91,file='GRPMAT')
        read(91,*)nn
        a=0.0d0
        mm=0
          do igp=1,nn
          read(91,*)((g(j,i,igp),j=1,3),i=1,3)
          end do
        close(91)
        do igp=1,nn
          do k=1,3
          b(k)=0.0d0
           do l=1,3
           b(k)=b(k)+g(k,l,igp)*r(l)
           end do
          end do
          err=abs(b(1)-r(1))
     &       +abs(b(2)-r(2))
     &       +abs(b(3)-r(3))
          if (err.lt.0.001)then
c             print '(i5,7f12.4)',mm,err,b,r
              mm=mm+1
              do j=1,3
              a(j)=a(j)+b(j)
              end do
          end if
        end do
            do j=1,3
            r(j)=a(j)/mm
            end do
        return
        end 

       SUBROUTINE GETGRP(NAME)
C ORIGINAL VERSION BY MARK R PEDERSON 1995-1996
       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*3 NAME
       DIMENSION R(3,3,121),PXYZ(3,3,6),RXYZ(3,3,3)
       SAVE
C
       IF (NAME .EQ. 'GRP') RETURN
C
C PXYZ(1): XYZ=>YXZ
C PXYZ(2): XYZ=>ZYX
C PXYZ(3): XYZ=>XZY
C PXYZ(4): XYZ=>YZX
C PXYZ(5): XYZ=>ZXY
C PXYZ(6): XYZ=>XYZ
C
       DO IP=1,3
        DO I=1,3
         DO J=1,3
          PXYZ(J,I,IP  )=0.0D0
          PXYZ(J,I,IP+3)=0.0D0
          RXYZ(J,I,IP  )=0.0D0
         END DO
         RXYZ(I,I,IP  )=1.0D0
        END DO
        RXYZ(IP,IP,IP)=-1.0D0
        IF (IP.EQ.1) THEN
         IX=1
         IY=2
         IZ=3
        ELSE IF (IP.EQ.2) THEN
         IX=3
         IY=1
         IZ=2
        ELSE IF (IP.EQ.3) THEN
         IX=2
         IY=3
         IZ=1 
        END IF
        PXYZ(IY,IX,IP)=1.0D0
        PXYZ(IX,IY,IP)=1.0D0
        PXYZ(IZ,IZ,IP)=1.0D0
       END DO
C
       PXYZ(1,2,4)=1.0D0
       PXYZ(2,3,4)=1.0D0
       PXYZ(3,1,4)=1.0D0
       PXYZ(1,3,5)=1.0D0
       PXYZ(2,1,5)=1.0D0
       PXYZ(3,2,5)=1.0D0
       PXYZ(1,1,6)=1.0D0
       PXYZ(2,2,6)=1.0D0
       PXYZ(3,3,6)=1.0D0
       MGP=1
       DO I=1,3
        DO J=1,3
         R(J,I,1)=PXYZ(J,I,6)
        END DO
       END DO
C
C START WITH GROUPS
C
       IF(NAME.EQ.'C3V')THEN
        MGP=6
        DO IGP=1,6
         DO I=1,3
          DO J=1,3
           R(J,I,IGP)=PXYZ(J,I,MGP)
          END DO
         END DO
         MGP=MGP-1
        END DO
        MGP=6
       END IF
C
       IF (NAME.EQ.'IH') THEN
        DO IGP=1,3
         MGP=MGP+1
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=RXYZ(J,I,IGP)
          END DO
         END DO
        END DO
        CALL CLSGRP(MGP,R) 
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=PXYZ(J,I,4)
         END DO
        END DO
        CALL CLSGRP(MGP,R) 
        A= 0.3090169943749474D0     
        B= 0.5000000000000000D0
        C= 0.8090169943749474D0  
        MGP=MGP+1
        R(1,1,MGP)= A
        R(1,2,MGP)= B
        R(1,3,MGP)=-C
        R(2,1,MGP)=-B
        R(2,2,MGP)= C
        R(2,3,MGP)= A
        R(3,1,MGP)= C
        R(3,2,MGP)= A
        R(3,3,MGP)= B 
       END IF
C
       IF ((NAME.EQ.'TD').OR.(NAME.EQ.'OH')) THEN
        DO I=1,3
         DO J=1,3
          R(J,I,2)=PXYZ(J,I,1)
         END DO
        END DO
        DO I=1,3
         DO J=1,3
          R(J,I,3)=PXYZ(J,I,4)
         END DO
        END DO
        DO I=1,3
         DO J=1,3
          R(J,I,4)=0.0D0       
         END DO
         R(I,I,4)=-1.0D0
        END DO
        R(1,1,4)=1.0D0
        MGP=4
C
        IF (NAME.EQ.'OH') THEN
         DO I=1,3
          DO J=1,3
           R(J,I,5)=0.0D0
          END DO
          R(I,I,5)=-1.0D0
         END DO
         MGP=5
        END IF
       END IF 
C
C CHECK FOR REFLECTIONS:
C
       IF ((NAME.EQ.'X').OR.(NAME.EQ.'Y').OR.(NAME.EQ.'Z')) THEN
        IF(NAME.EQ.'X')IC=1
        IF(NAME.EQ.'Y')IC=2
        IF(NAME.EQ.'Z')IC=3
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=0.0D0
         END DO
         R(I,I,MGP)=1.0D0
        END DO
        R(IC,IC,MGP)=-1.0D0
       END IF
C
       IF ((NAME(1:2).EQ.'XY').OR.(NAME(1:2).EQ.'YX').OR.
     &     (NAME(1:2).EQ.'XZ').OR.(NAME(1:2).EQ.'ZX').OR.
     &     (NAME(1:2).EQ.'YZ').OR.(NAME(1:2).EQ.'ZY')) THEN
        IF(NAME(1:1).EQ.'X')I1=1
        IF(NAME(1:1).EQ.'Y')I1=2
        IF(NAME(1:1).EQ.'Z')I1=3
        IF(NAME(2:2).EQ.'X')I2=1
        IF(NAME(2:2).EQ.'Y')I2=2
        IF(NAME(2:2).EQ.'Z')I2=3
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=0.0D0
         END DO
         R(I,I,MGP)=1.0D0
        END DO
        R(I1,I1,MGP)=-1.0D0
        MGP=MGP+1
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=0.0D0
         END DO
         R(I,I,MGP)=1.0D0
        END DO
        R(I2,I2,MGP)=-1.0D0
        IF ((NAME(3:3).EQ.'X').OR.(NAME(3:3).EQ.'Y').OR.
     &      (NAME(3:3).EQ.'Z')) THEN
         IF(NAME.EQ.'X')IC=1
         IF(NAME.EQ.'Y')IC=2
         IF(NAME.EQ.'Z')IC=3
         MGP=MGP+1
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=0.0D0
          END DO
          R(I,I,MGP)=1.0D0
         END DO
         R(IC,IC,MGP)=-1.0D0
        END IF
       END IF
C 
       IF ((NAME.EQ.'D4H').OR.(NAME.EQ.'C4V')) THEN
        DO MGP=1,3
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=0.0D0
          END DO
          R(I,I,MGP)=1.0D0
         END DO
         IF (MGP.GT.1) R(MGP-1,MGP-1,MGP)= -1.0D0
        END DO
        MGP=4
        DO I=1,3
         DO J=1,3
          R(J,I,MGP)=PXYZ(J,I,1)
         END DO
        END DO
        IF (NAME.EQ.'D4H') THEN
         MGP=MGP+1
         DO I=1,3
          DO J=1,3
           R(J,I,MGP)=0.0D0
          END DO
          R(I,I,MGP)=1.0D0
         END DO
         R(3,3,MGP)= -1.0D0
        END IF
       END IF
C
       CALL CLSGRP(MGP,R)
       OPEN(60,FILE='GRPMAT',FORM='FORMATTED',STATUS='UNKNOWN')
       WRITE(60,*) MGP,' ',NAME
       DO IGP=1,MGP
        DO I=1,3
         WRITE(60,60)(R(J,I,IGP),J=1,3)
        END DO
        WRITE(60,*)' '
       END DO
 60    FORMAT(' ',3G25.16)
       CLOSE(60)
       RETURN
       END
C
C ***************************************************************
C
        SUBROUTINE CLSGRP(MGP,R)
C ORIGINAL VERSION BY MARK R PEDERSON (1996)
C
C CLSGROUP FINDS ALL THE GROUP OPERATIONS GENERATED BY THE MGP
C INPUT GROUP ELEMENTS
C
         IMPLICIT REAL*8  (A-H,O-Z)
         DIMENSION R(3,3,121)
         SAVE
C
C GENERATE GROUP MATRIX:
C
         NGP=MGP 
   10    CONTINUE
          MGP=NGP
          IF (NGP.GT.120) THEN
           PRINT *,'CLSGRP: NUMBER OF GROUP ELEMENTS > 120'
           CALL STOPIT
          END IF
          DO IGP=1,MGP
           DO JGP=1,MGP
            NGP=NGP+1
                    IF(NGP.GT.121)THEN
                    PRINT*,'OUT OF SPACE IN CLSGRP'
                    PRINT*,'PROBABLE ERROR'
                    CALL STOPIT
                    END IF
            DO I=1,3
             DO J=1,3
              R(J,I,NGP)=0.0D0
              DO K=1,3
               R(J,I,NGP)=R(J,I,NGP)+R(J,K,IGP)*R(K,I,JGP)
              END DO
             END DO
            END DO
            ITOT=0
            DO KGP=1,NGP-1
             ERROR=0.0D0
             DO I=1,3
              DO J=1,3
               ERROR=ERROR+ABS(R(J,I,NGP)-R(J,I,KGP))
              END DO
             END DO
             IF(ERROR .LT. 1.0D-4) ITOT=ITOT+1
            END DO
            IF (ITOT.GT.0) NGP=NGP-1
           END DO
          END DO
          IF (NGP.GT.MGP) GOTO 10
         CONTINUE
         RETURN
        END
C
C *******************************************************************
       SUBROUTINE NONSYM(GENGRP)
       IMPLICIT REAL*8 (A-H,O-Z)
C GENGRP = 2PI009 MEANS ADD ROTATION by 2PI/9
       DIMENSION RMAT(3,3,121)
       CHARACTER*8  GENGRP
       READ(GENGRP(4:6),*)NRT
C         NAXIS=3
C         DO I=1,8
C         IF(GENGRP(I:I).EQ.'X')NAXIS=1
C         IF(GENGRP(I:I).EQ.'Y')NAXIS=2
C         IF(GENGRP(I:I).EQ.'Z')NAXIS=3
C         END DO
       RMAT=0.0D0
       DO I=1,3
       RMAT(I,I,1)=1.0D0
       END DO
            DO J=1,3
            RMAT(J,J,2)=1.0D0
            END DO
C      IF(NAXIS.EQ.3)THEN
       IX=1
       IY=2
       IZ=3
C      ELSE IF(NAXIS.EQ.2)THEN
C      IX=1
C      IY=3
C      IZ=2
C      ELSE IF(NAXIS.EQ.1)THEN
C      IX=3
C      IY=2
C      IZ=1
C      END IF
       RMAT(IX,IX,2)= COS(8.D0*ATAN(1.D0)/NRT)
       RMAT(IX,IY,2)= SIN(8.D0*ATAN(1.D0)/NRT)
       RMAT(IY,IX,2)=-SIN(8.D0*ATAN(1.D0)/NRT)
       RMAT(IY,IY,2)= COS(8.D0*ATAN(1.D0)/NRT)
       MGP=2
       IF(MGP.GT.121)THEN
         PRINT*,'121   MUST BE>',MGP
         CALL STOPIT
       END IF
C      CALL CLSGRP(MGP,RMAT)
       OPEN(60,FILE='GRPMAT',FORM='FORMATTED',STATUS='UNKNOWN')
       WRITE(60,*) MGP,' ',GENGRP
       DO IGP=1,MGP
        DO I=1,3
         PRINT    60,(RMAT(J,I,IGP),J=1,3)
         WRITE(60,60)(RMAT(J,I,IGP),J=1,3)
 60      FORMAT(3F20.15)
        END DO
        WRITE(60,*)' '
       END DO
       CLOSE(60)
       RETURN
       END  
       SUBROUTINE ISETUP  
C ORIGINAL VERSION BY MARK R PEDERSON (1996)
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER(MXCHG=56)
       PARAMETER(MXATM=1000)
       PARAMETER(MXSET=100)
       LOGICAL   EXIST,PSEUDO
       CHARACTER*60 NAMFUNCT
       CHARACTER*3  NAMGRP,TABLE(0:MXCHG)
       CHARACTER*3  NAMPSP(MXSET),SPNDIR(MXSET)
       CHARACTER*6  CHARINF
       CHARACTER*8  GENGRP
       DIMENSION IZNUC(MXSET),IZELC(MXSET)
       DIMENSION R(3,MXATM),IDXSET(MXATM)
       DIMENSION ALPSET(MAX_BARE),CONSET(MAX_BARE,MAX_CON,3),NBASF(2,3)
       DIMENSION RNUC(3,MX_GRP),MSITES(1)
       DIMENSION V(3),ELECTRONS(2)
       DATA TABLE/
     & 'PSE',
     & 'HYD','HEL',
     & 'LIT','BER','BOR','CAR','NIT','OXY','FLU','NEO',
     & 'SOD','MAG','ALU','SIL','PHO','SUL','CHL','ARG',
     & 'POT','CAL','SCA','TIT','VAN','CHR','MAN','IRO','COB','NIC',
     & 'COP','ZIN','GAL','GER','ARS','SEL','BRO','KRY',
     & 'RUB','STR','YTR','ZIR','NIO','MOL','TEC','RHU','RHO','PAL',
     & 'SLV','CAD','IND','TIN','ANT','TEL','IOD','XEN',
     & 'CES','BAR'/
       DATA TOLER/1.0D-5/
       INQUIRE(FILE='GRPMAT',EXIST=EXIST)
C
       PSEUDO= .FALSE.
       INQUIRE(FILE='CLUSTER',EXIST=EXIST)
       IF (.NOT. EXIST) GOTO 800
       INQUIRE(FILE='GRPMAT',EXIST=EXIST)
       if(EXIST)call perfect_cluster
       OPEN(90,FILE='CLUSTER',FORM='FORMATTED',STATUS='OLD')
       REWIND(90)
       READ(90,'(A)', END=200) NAMFUNCT
            IF(NAMFUNCT(1:1).EQ.'@')THEN
                CLOSE(90)
                CALL EZSTART
                CALL STOPIT
            END IF
C      READ(90,'(A3)',END=200) NAMGRP
       READ(90,'(A6)',END=200) GENGRP
       NAMGRP(1:3)=GENGRP(1:3)
   30  FORMAT(A3)
C
C GET GROUP MATRIX IF DESIRED
C
       IF(GENGRP(1:3).EQ.'2PI')THEN
       CALL NONSYM(GENGRP)
       ELSE
       CALL GETGRP(NAMGRP)
       END IF
C
C READ GROUP INFORMATION, CHECK GROUP, AND STORE IT IN COMMON /GROUP/
C
       CALL FGMAT
C
C DEAL WITH ATOMS
C
       READ(90,*,END=200) MMATOMS
       IF (MMATOMS .LE. 0) THEN
        PRINT *,'ISETUP: NUMBER OF ATOMS IS ZERO'
        GOTO 900
       END IF
       NSETS=0
       NATOMS=0
       DO IATRD=1,MMATOMS
        NATOMS=NATOMS+1
        IDXSET(NATOMS)=0
        NSETS=NSETS+1  
        IF (NATOMS .GT. MXATM) THEN
         PRINT *,'ISETUP: MXATM MUST BE AT LEAST: ',NATOMS
         GOTO 900
        END IF
        IF (NSETS .GT. MXSET) THEN
         PRINT *,'ISETUP: MXSET MUST BE AT LEAST: ',NSETS
         GOTO 900
        END IF
C
C READ INFO FOR ONE ATOM, CHECK IF CHARGE IN BOUNDS
C
        CALL IGETATM(90,R(1,NATOMS),IZNUC(NSETS),CHARINF,IREAD)
                NAMPSP(NSETS)(1:3)=CHARINF(1:3)
                SPNDIR(NATOMS)(1:3)=CHARINF(4:6)
                print*,CHARINF(1:3),CHARINF(4:6)
        IF (IREAD .NE. 0) GOTO 200
        IF (NAMPSP(NSETS) .NE. 'ALL') PSEUDO= .TRUE.
        IF (IZNUC(NSETS) .GT. MXCHG) THEN
         PRINT *,'ISETUP: ATOMS WITH Z > ',MXCHG,' ARE NOT SUPPORTED'
         GOTO 900
        END IF
        IF (IZNUC(NSETS) .LT. 0) THEN
         PRINT *,'ISETUP: ATOMS WITH Z < 0 ARE NOT SUPPORTED'
         GOTO 900
        END IF
C
C CHECK FOR EQUIVALENT ATOMS
C
        ITOT=0
        DO IATOMS=1,NATOMS-1
         JTOT=0
         DO IGP=1,NGRP
          DO I=1,3
           V(I)=0.0D0
           DO J=1,3
            V(I)=V(I)+RMAT(I,J,IGP)*R(J,IATOMS)
           END DO
          END DO
          ERROR=ABS(V(1)-R(1,NATOMS)) 
     &         +ABS(V(2)-R(2,NATOMS)) 
     &         +ABS(V(3)-R(3,NATOMS)) 
          IF (ERROR .LT. TOLER) JTOT=JTOT+1
         END DO
         IF ((JTOT .NE. 0) .AND. 
     &       ((IZNUC (IDXSET(IATOMS)) .NE. IZNUC (NSETS)) .OR.
     &        (NAMPSP(IDXSET(IATOMS)) .NE. NAMPSP(NSETS)))) THEN
          PRINT *,'ISETUP: FOUND ATOMS OF DIFFERENT TYPE AT THE'
          PRINT *,'SAME COORDINATES - ASK A NUCLEAR PHYSICIST'
          GOTO 900
         END IF
         ITOT=ITOT+JTOT
        END DO
C
C FORGET THIS ATOM IF IT CAN BE CONSTRUCTED BY SYMMETRY
C
        IF (ITOT .NE. 0) THEN 
         NATOMS=NATOMS-1
         NSETS=NSETS-1
        ELSE
C
C CHECK IF NEW ATOM TYPE
C
         ITOT=NSETS
         DO ISET=1,NSETS-1 
          IF ((IZNUC (ISET) .EQ. IZNUC (NSETS)) .AND.
     &        (NAMPSP(ISET) .EQ. NAMPSP(NSETS))) ITOT=ISET
         END DO
         IF (ITOT .NE. NSETS) NSETS=NSETS-1
         IDXSET(NATOMS)=ITOT
        END IF
       END DO
C
C DEFINE NET CHARGE AND SPIN
C IPSCHG RETURNS THE ACTUAL NUMBER OF ELECTRONS ON THE ATOM
C
       READ(90,*,END=200) CHNET,SPNNET
       ZTOT=0.0D0
       DO IATOMS=1,NATOMS
        ISET=IDXSET(IATOMS)
        IZELC(ISET)=IPSCHG(NAMPSP(ISET),IZNUC(ISET))  
        CALL GASITES(1,R(1,IATOMS),NNUC,RNUC,MSITES)
        ZTOT=ZTOT+NNUC*IZELC(ISET) 
       END DO
       ELECTRONS(1)= -(CHNET-SPNNET-ZTOT)/2.0D0
       ELECTRONS(2)= -(CHNET+SPNNET-ZTOT)/2.0D0
       CLOSE(90)
       GOTO 210
C
C ERROR
C
  200  PRINT *,'ISETUP: FILE CLUSTER IS INVALID'
       CLOSE(90)
       CALL STOPIT
C
C GENERATE INPUT FILES
C FIRST, WRITE HEADER OF ISYMGEN, SYMBOL
C
  210  PRINT '(A)','GENERATING INPUT FILES FROM DATA IN FILE CLUSTER'
       OPEN(90,FILE='ISYMGEN',FORM='FORMATTED',STATUS='UNKNOWN')
       OPEN(92,FILE='SYMBOL' ,FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(90)
       REWIND(92)
       IF (PSEUDO) THEN
        OPEN(94,FILE='PSPINP' ,FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(94)
       END IF
C      
C CHANGE DEFAULT TO LBFGS....
       IF(NATOMS.EQ.1) THEN
         WRITE(92,'(A)') 'CONJUGATE-GRADIENT'
       ELSE
         WRITE(92,'(A)') 'LBFGS'                 
       ENDIF
       WRITE(92,'(A)') NAMFUNCT
       WRITE(92,'(A)') 'OLDMESH'
       WRITE(92,'(A)') '1    NUMBER OF SYMBOLIC FILES'
       WRITE(92,'(A)') 'ISYMGEN = INPUT'
       WRITE(92,'(I3,A)') NATOMS+2,'  NUMBER OF SYMBOLS IN LIST'
       WRITE(92,'(I3,A)') NATOMS,  '  NUMBER OF NUCLEI'
       DO IATOMS=1,NATOMS 
        WRITE(92,'(A)') '1.0  1 1 1'
       END DO
       WRITE(92,'(A)') '  1  CALCULATION SETUP BY ISETUP'
C
C UPDATE ISYMGEN, SYMBOL, PSPINP
C
       WRITE(90,'(1X,I3,10X,A)') NSETS,'TOTAL NUMBER OF ATOM TYPES'
       DO ISET=1,NSETS  
        IZN=IZNUC(ISET)
        IZE=IZELC(ISET)
        PSPSYM(1)(1:3)= NAMPSP(ISET)
        PSPSYM(1)(4:4)= '-'
        PSPSYM(1)(5:7)= TABLE(IZN)
        CALL SETPSP(94,PSPSYM(1),IZN)
        CALL SETBAS(IZN,NAMPSP(ISET),ALPSET,CONSET,NALP,NBASF)
        NFIND=0
        DO IATOMS=1,NATOMS
         IF (IDXSET(IATOMS) .EQ. ISET) NFIND=NFIND+1
        END DO 
        PRINT 1020,ISET,PSPSYM(1),IZN,NFIND
 1020   FORMAT('ATOM TYPE ',I3,' (',A,'): NUCLEAR CHARGE= ',I3,', ',
     &         I3,' ATOM(S)')
        WRITE(90,'(2(1X,I3),6X,A)') IZE,IZN,
     &           'ELECTRONIC AND NUCLEAR CHARGE'
        IF (NAMPSP(ISET) .EQ. 'ALL') THEN
         WRITE(90,'(A3,11X,A)') 'ALL','ALL-ELECTRON ATOM TYPE'
        ELSE
         WRITE(90,'(A7,7X,A)') PSPSYM(1),'PSP-SYMBOL'
        END IF
        WRITE(90,'(1X,I3,10X,A,A3)') NFIND,
     &           'NUMBER OF ATOMS OF TYPE ',TABLE(IZN)
C
C ADD TO ISYMGEN, SYMBOL
C
        NFIND=0
        DO IATOMS=1,NATOMS
         IF (IDXSET(IATOMS) .EQ. ISET) THEN
          NFIND=NFIND+1
          WRITE(90,1030) PSPSYM(1),NFIND
 1030     FORMAT(A7,I3.3)
                 IOK=0
                 IF(SPNDIR(IATOMS).EQ.'SUP')IOK=1
                 IF(SPNDIR(IATOMS).EQ.'SDN')IOK=1
                 IF(SPNDIR(IATOMS).EQ.'UPO')IOK=1
                 IF(IOK.EQ.0)SPNDIR(IATOMS)='UPO'
          WRITE(92,1040) PSPSYM(1),NFIND,(R(J,IATOMS),J=1,3),
     &      SPNDIR(IATOMS) 
 1040     FORMAT(A7,I3.3,' =',3(1X,F20.10),' ',A3)
         END IF
        END DO
C
C ADD BASIS SET TO ISYMGEN
C 
        WRITE(90,'(2A)') 'EXTRABASIS    CONTROLS USAGE OF ',
     &                   'SUPPLEMENTARY BASIS FUNCTIONS'
        WRITE(90,1050) NALP,'NUMBER OF BARE GAUSSIANS'
 1050   FORMAT(1X,I3,10X,A) 
        WRITE(90,1060)(NBASF(1,L), L=1,3),'NUMBER OF S,P,D FUNCTIONS'
        WRITE(90,1060)(NBASF(2,L), L=1,3),
     &                'SUPPLEMENTARY S,P,D FUNCTIONS'
 1060   FORMAT(3(1X,I3),2X,A) 
        WRITE(90,1070)(ALPSET(J), J=1,NALP)
        WRITE(90,*)
        DO L=1,3
         DO IL=1,NBASF(1,L)+NBASF(2,L)
          WRITE(90,1070)(CONSET(J,IL,L), J=1,NALP)
          WRITE(90,*)
         END DO
        END DO 
       END DO
 1070  FORMAT(3(1X,D20.8))
       WRITE(90,'(A)') 'ELECTRONS'
       WRITE(90,'(A)') 'WFOUT'
C
C ADD SOME STUFF TO SYMBOL
C
       WRITE(92,1080) ELECTRONS
       WRITE(92,1090)
 1080  FORMAT('ELECTRONS  =',2(1X,F15.6))
 1090  FORMAT('EXTRABASIS = 0')
       CLOSE(90)
       CLOSE(92)
       IF (PSEUDO) CLOSE(94)
       PRINT '(A)',' '
       RETURN
C
C SETUP DEFAULT CLUSTER FILE
C
  800  OPEN(90,FILE='CLUSTER',FORM='FORMATTED',STATUS='NEW')
       REWIND(90)
       WRITE(90,'(2A)') 'GGA-PBE*GGA-PBE          ',
     &                  '(DF TYPE EXCHANGE*CORRELATION)'
       WRITE(90,'(2A)') 'TD                       ',
     &                  '(TD, OH, IH, X, Y, XY, ... OR GRP)'
       WRITE(90,'(2A)') '2                        ',
     &                  '(NUMBER OF INEQUIV. ATOMS IN CH4)'
       WRITE(90,'(2A)') '0.00  0.00  0.00  6  ALL ',
     &                  '(R, Z, PSEUDOPOT. TYPE FOR CARBON)'
       WRITE(90,'(2A)') '1.20  1.20  1.20  1  ALL ',
     &                  '(R, Z, PSEUDOPOT. TYPE FOR HYDROGEN)'
       WRITE(90,'(2A)') '0.0 0.0                  ',
     &                  '(NET CHARGE AND NET SPIN)' 
       WRITE(90,*)'--------------OR-------------------'
       WRITE(90,'(1A)')'@XMOL.DAT'
       WRITE(90,*)'IF YOU WISH TO START FROM AN XYZ XMOL FILE'
       CLOSE(90)
       PRINT '(A)','LOOK AT AND EDIT FILE NAMED CLUSTER, THEN RERUN'
       RETURN
C
C ERROR
C
  900  CLOSE(90)
       CALL STOPIT
       END
C
C ******************************************************************
C
        SUBROUTINE IGETATM(IUNIT,R,IZ,NAME,IREAD)
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION    R(3)
         SAVE
C
         CHARACTER*6  NAME
         CHARACTER*1  CCH
         CHARACTER*80 LINE 
         IREAD=1
         READ(IUNIT,'(A80)',END=100) LINE
                 K=81
                 DO I=80,1,-1
                 IF(LINE(I:I).EQ.'!')K=I
                 END DO
                 DO I=K,80
                 LINE(I:I)=' '
                 END DO
         READ(LINE,*,END=100) R,IZ
         IPOS=0
         IBLANK=1
         DO IJUMP=1,5
   10     CONTINUE
           IPOS=IPOS+1
           IF (IPOS .GT. 80) GOTO 100
           CCH=LINE(IPOS:IPOS)
           IF (IBLANK .EQ. 1) THEN
            IF ((CCH .NE. ' ') .AND. (CCH .NE. '\t') 
     &                         .AND. (CCH .NE. ',')) THEN
             IF (IJUMP .EQ. 5) GOTO 20
             IBLANK=0
            END IF
            GOTO 10
           ELSE
            IF ((CCH .NE. ' ') .AND. (CCH .NE. '\t') 
     &                         .AND. (CCH .NE. ',')) THEN
             GOTO 10
            ELSE
             IBLANK=1
            END IF
           END IF
          CONTINUE
         END DO
   20    CONTINUE 
               K=0
               DO I=IPOS,80
                IF(LINE(I:I).NE.' ')THEN
                K=K+1
                LINE(K:K)=LINE(I:I)
                END IF
               END DO  
            NAME=LINE(1:6)
         IREAD=0
  100    RETURN
        END
