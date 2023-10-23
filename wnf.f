       SUBROUTINE WANNIER(MODE)
       END
       
C
C***********************************************************************
C
         SUBROUTINE MATMUL(NDIM,NWV,PSI,HAM,OVER)
         include 'PARAMS'
         DIMENSION PSI(NDH,*)
         DIMENSION OVER(NDH,NDH),HAM(NDH,NDH)
         DIMENSION TMP(NDH,NDH)

          DO I=1,NDH 
            DO J=1,NDH 
              HAM(I,J)=0.0D0
              TMP(I,J)=0.0D0
            END DO
          END DO

          DO IW=1,NWV
            DO K=1,NDIM
              DO L=1,NDIM
               TMP(K,IW)= TMP(K,IW)+PSI(L,IW)*OVER(K,L)
              END DO
            END DO
          END DO

          DO IW=1,NWV
            DO JW=1,NWV
              DO K=1,NDIM
               HAM(JW,IW)= HAM(JW,IW)+PSI(K,JW)*TMP(K,IW)
              END DO
            END DO
          END DO
5034      FORMAT (5F13.5)
         RETURN
         END 

C
C***********************************************************************
C

           SUBROUTINE LOWDEN(NDH,M,OVER,HAM,EVAL,SC1)
C   Mark Pederson 14 September 2001
C TO USE
C PLACE OVERLAP MATRIX IN OVER.
C NEW WAVEFUNCTIONS ARE RETURNED IN HAM...
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION HAM(NDH,NDH),OVER(NDH,NDH),EVAL(NDH),SC1(*)
           LOGICAL CHECKIT
           DATA CHECKIT/.FALSE./
           IF(CHECKIT)THEN
           PRINT*," NOTE THAT IF CHECKIT=.TRUE. THERE IS N**4 SCALING."
           PRINT*," NORMAL CALLS SHOULD BE WITH CHECKIT=.FALSE."
           END IF
              OPEN(32,FILE='SCRLOW',FORM='UNFORMATTED')
              WRITE(32)((OVER(J,I),J=1,M),I=1,M) 
                  DO I=1,M
                  DO J=1,M
                  HAM(J,I)=0.0D0
                  END DO  
                  HAM(I,I)=1.0D0
                  END DO
           CALL DIAGGE(NDH,M,OVER,HAM,EVAL,SC1,1)
           PRINT*,' LOWDEN OVERLAP EIGENVALUES:'
           PRINT 20,(EVAL(I),I=1,M)
                   DO I=1,M
                   EVAL(I)=1.0D0/SQRT(EVAL(I))
                   END DO
C CHECK FOR ORTHOGONALITY:
          IF(CHECKIT)THEN
               DO I=1,M
               DO J=1,M
                 SC1(J)=0.0D0
                 DO K=1,M
                 SC1(J)=SC1(J)+OVER(K,I)*OVER(K,J)
                 END DO 
               END DO
               PRINT 20,(SC1(J),J=1,M)
               END DO
           END IF
C PERFORM BACK TRANSFORMATION:
               DO J=1,M
                  DO I=1,M
                  HAM(I,J)=0.0D0
                   DO N=1,M
                   HAM(I,J)=HAM(I,J)+OVER(J,N)*OVER(I,N)*EVAL(N)
                   END DO
                  END DO
               END DO
C RELOAD OVERLAP MATRIX:
              REWIND(32)
              READ(32)((OVER(J,I),J=1,M),I=1,M) 
              CLOSE(32,STATUS='DELETE')
C CHECK AGAINST 3-ORDER LOWDEN:
              IF(CHECKIT)THEN
               PRINT*,'COMPARISON TO 3RD ORDER LOWDEN...'
               DO I=1,M
                  DO J=1,M
                  SC1(J)=-OVER(J,I)/2.0D0
                  END DO
                  SC1(I)=1.0D0
                  DO J=1,M
                  DO K=1,M
                  IF(K.NE.J.AND.K.NE.I)THEN
                  SC1(J)=SC1(J)+OVER(K,J)*OVER(K,I)*(3./8.)
                  END IF
                  END DO 
                  END DO
               PRINT 20,(HAM(J,I),J=1,M)
               PRINT 20,(SC1(J  ),J=1,M)
               PRINT*,' '
               END DO
C CHECK FOR ORTHOGONALITY:
C
               PRINT*,'TRANSFORMED OVERLAP MATRIX:'
               DO I=1,M
               DO J=1,M
                   SC1(J)=0.0D0
                   DO K=1,M
                   DO L=1,M
                   SC1(J)=SC1(J)+HAM(L,I)*HAM(K,J)*OVER(L,K)
                   END DO
                   END DO
               END DO
               PRINT 20,(SC1(J),J=1,M)
               END DO
              END IF
 20            FORMAT(' ',5g15.6)
         END

        subroutine window(nspn,evl,occ)
        character*80 line
        dimension evl(*),occ(*)
        dimension part(500)
        open(43,file='WANWIN')
        rewind(43)
        emax=-1.0d30
        open(90,file='EVALUES')
        open(83,file='SCRATCH')
        do 200 ispn=1,nspn
        iwf=0
        rewind(90)
        rewind(83)
  100   continue
        read(90,90)line
 90     format(a80)
        if(line(1:5).ne.' SUMM')go to 100
        print*,line
        do 115 iii=1,100000
        rewind(83)
        read(90,90)line 
c       print*,line
              do i=1,72
              if(line(i:i+3).eq.'REP:')line(i:i+3)='    '
              if(line(i:i+3).eq.'DEG:')line(i:i+3)='    '
              if(line(i:i+3).eq.'OCC:')line(i:i+3)='    '
              if(line(i:i+4).eq.'SPIN:')line(i:i+4)='     '
              if(line(i:i+6).eq.'ENERGY:')line(i:i+6)='       '
              end do
c       print*,line
        write(83,*)line
        rewind(83)
        read(83,*)n,m,l,is,e,o
        emax=max(e,emax)
           if(is.eq.ispn.and.o.gt.0.01)then
           iwf=iwf+1
           evl(iwf)=e
           occ(iwf)=o
           end if
 10        format(' ',i5,3f12.6)
        if(o.lt.0.01)go to 120
 115    continue
 120    continue
           print 10,iwf,evl(iwf),occ(iwf),emax
           print*,'ispn:',ispn
           npt=0
           npt=npt+1  
           part(npt)=(emax+evl(iwf))/2.
           do jwf=iwf-1,1,-1
           if(abs(evl(jwf+1)).lt.0.8*abs(evl(jwf)) )then
           npt=npt+1
           part(npt)=(evl(jwf+1)+evl(jwf))/2.
           end if
           end do
           npt=npt+1
           part(npt)=-10000000.
           print *,npt
           write(43,*)npt-1
           do ipt=npt,2,-1
           print 20,part(ipt),part(ipt-1)
           write(43,20)part(ipt),part(ipt-1)
           end do
 20     format(' ',2f20.6)
        close(83,status='delete')
 200    continue
        close(90)
        close(43)
        end
