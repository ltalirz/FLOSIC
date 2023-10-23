       SUBROUTINE PART_HOLE(NDH,ISP,IRP,NBS,NCH,HAM,OVR,EV,SC,ND)
       IMPLICIT REAL*8 (A-H,O-Z)
       LOGICAL EXIST,IFRET
       DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),EV(NDH),SC(2*NDH)
       DIMENSION OO(1000),HH(1000)
       INQUIRE (FILE='PART_HOLE',EXIST=EXIST)
       OPEN(20,FILE='PART_HOLE')               
       IF(.NOT.EXIST)THEN
       WRITE(20,*).FALSE.,' T or F to do DSCF'
       WRITE(20,*)'1 1 1 !ISPN, IREP, IEIG'
       WRITE(20,*)0,' NCH=0 For "NEUTRAL" NCH=1 for "Ionized"'
       CLOSE(20)
       END IF
       REWIND(20)
       READ(20,*)IFRET
       READ(20,*)JSP,JRP,JBS
       READ(20,*)NCH 
         IERR=ABS(JSP-ISP)+ABS(JRP-IRP)
         IF(ND.EQ.2.AND.NCH.EQ.0.AND.IERR.EQ.0.AND.NCH.EQ.0)THEN    
         WRITE(20,*)EV(JBS)
         WRITE(20,'(6F20.12)')(HAM(J,JBS),J=1,NBS)
         END IF
       CLOSE(20)
       OPEN(20,FILE='PART_HOLE')               
       READ(20,*)IFRET
       READ(20,*)JSP,JRP,JBS
       READ(20,*)NCH
         IF(JSP.NE.ISP.OR.JRP.NE.IRP)THEN
            CLOSE(20)
            RETURN
         END IF
       IF(ND.EQ.1.AND.NCH.EQ.1)THEN    
       READ (20,*)EV_NP
       READ (20,*)(SC(J),J=1,NBS)
C CONSTRUCT (1-|SC><SC|)H(1-|SC><SC|)=H - |SC><SC|H - H|SC><SC|+
C                                              <SC|H|SC>|SC><SC|
           DO M=1,NBS
           OO(M)=0.0D0
           HH(M)=0.0D0
             DO N=1,NBS
             OO(M)=OO(M)+SC(N)*OVR(N,M)
             END DO
             DO N=1,NBS
             HH(M)=HH(M)+SC(N)*HAM(N,M)
             END DO
           END DO
           FK=0.0D0
           DO N=1,NBS
           DO M=1,NBS
           FK=FK+SC(N)*SC(M)*HAM(M,N)
           END DO
           END DO
           FK=FK-5000.
           DO M=1,NBS
           DO N=1,NBS
           HAM(N,M)=HAM(N,M)-OO(N)*HH(M)-OO(M)*HH(N)+FK*OO(N)*OO(M)
           END DO
           END DO
       ELSE IF(ND.EQ.2.AND.NCH.EQ.1)THEN
       READ (20,*)EV_NP
       READ (20,*)(SC(J),J=1,NBS)
           DO I=1,JBS-1
                DO J=1,NBS
                HAM(J,I)=HAM(J,I+1)
                END DO
                EV(I)=EV(I+1)
           END DO
                DO J=1,NBS
                HAM(J,JBS)=SC(J)
                END DO
                EV(JBS)=EV_NP        
       END IF
       CLOSE(20)
       RETURN      
       END   
