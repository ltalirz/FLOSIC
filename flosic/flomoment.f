       SUBROUTINE CONDENSE(NSPN)   
       INCLUDE 'PARAMS'
       COMMON/MESH/WMSH(MAX_PTS),RMSH(3,MAX_PTS),NMSH
       COMMON/SPARSE/CNTFLO(3),ALPFLO(3),NANSWER
       COMMON/TMP1/CMB    (MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
C      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
       DIMENSION IPTR(MAX_PTS)
C      EQUIVALENCE (IPTR(1),POTIN(1))
       LOGICAL SKIPIT
       DATA SKIPIT/.FALSE./   !!!!
         IF(SKIPIT)THEN
         PRINT*,'SKIPPING CONDENSE'
      !!!   NMSH=0
         CALL COUPOT1
         RETURN
         END IF
       OPEN(22,FILE='DIRMSH',FORM='UNFORMATTED')
       WRITE(22)NMSH
       WRITE(22) (WMSH(IMSH)         ,IMSH=1,NMSH)
       WRITE(22)((RMSH(J,IMSH),J=1,3),IMSH=1,NMSH) 
       CLOSE(22)
           MMSH=0
           
           DO IMSH=1,NMSH
           ARG=ALPFLO(1)*SQRT(((RMSH(1,IMSH)-CNTFLO(1))**2
     &                   +(RMSH(2,IMSH)-CNTFLO(2))**2
     &                   +(RMSH(3,IMSH)-CNTFLO(3))**2))
           IF(ARG.LT.20.0D0)then
           MMSH=MMSH+1
           IPTR(MMSH)=IMSH
C IPTR(MMSH)>MMSH
           RMSH(1,MMSH)=RMSH(1,IMSH)
           RMSH(2,MMSH)=RMSH(2,IMSH)
           RMSH(3,MMSH)=RMSH(3,IMSH)
           WMSH(  MMSH)=WMSH(  IMSH)
           END IF
           END DO
       PRINT 60,(CNTFLO(J),J=1,3),(ALPFLO(J),J=1,3),
     &FLOAT(MMSH)/FLOAT(NMSH)
 60    FORMAT('CONDENSE:',3F12.6,'  ',3F8.3,'  ',F12.4)
       NMSH=MMSH
       CMB =0.0D0
       RHOG=0.0D0
       CALL COUPOT1
       OPEN(22,FILE='DIRMSH',FORM='UNFORMATTED')
       READ (22)NMSH
       READ (22) (WMSH(IMSH)         ,IMSH=1,NMSH) 
       READ (22)((RMSH(J,IMSH),J=1,3),IMSH=1,NMSH)
       CLOSE(22,STATUS='DELETE')
           DO JMSH=MMSH,1,-1       
           DO JSPN=1,NSPN
           DO IG=1,NVGRAD !POTENTIAL BUG
           RHOG(IPTR(JMSH),IG,JSPN)=RHOG(JMSH,IG,JSPN)
           IF(JMSH.NE.IPTR(JMSH))THEN
           RHOG(     JMSH,IG,JSPN)=0.0D0
           END IF
           END DO
           END DO
           END DO
           DO JMSH=MMSH,1,-1   
           CMB(IPTR(JMSH))=CMB(JMSH)
           IF(JMSH.NE.IPTR(JMSH))THEN
                           CMB(JMSH)=0.0D0
           END IF
           END DO
       RETURN
       END    
       SUBROUTINE FLOMOMENT
       INCLUDE 'PARAMS'
       include 'commons.inc'
       COMMON/TMP4/PHI    (MAX_PTS)! FLO amplitude on grid
       COMMON/FORRx/NFLO,KSPX!!,TMAT(NDH,NDH,2)
       COMMON/SPARSE/CNTFLO(3),ALPFLO(3),NANSWER
       dimension cnt(3)
c
c  COMPUTING moments for PHI = abs(NFLO)
c
       print *, 'in FLOMOMENT, NFLO = ', NFLO
       print *, 'nmsh', nmsh
       call flonase(time)
       chg = 0.0d0
       cnt = 0.0d0
       do i = 1, nmsh
        chg = chg + phi(i)**2*wmsh(i)
        do j = 1,3
           cnt(j) = cnt(j) + rmsh(j,i)*phi(i)**2*wmsh(i)
        end do
       end do
       print *, 'FLOMOMENT CHG: ', CHG
       print *, 'FLOMOMENT centroid ', (cnt(j),j=1,3)
       CNTFLO=CNT
       r2 = 0.0d0
       r4 = 0.0d0
       r6 = 0.0d0
       do i = 1, nmsh
        dist = 0.0d0
        do j = 1,3
           dist = dist + (rmsh(j,i) - cnt(j))**2
        end do
        r2 = r2 + dist   *phi(i)**2*wmsh(i)
        r4 = r4 + dist**2*phi(i)**2*wmsh(i)
        r6 = r6 + dist**3*phi(i)**2*wmsh(i)
       end do
       print *, 'FLOMOMENT r2, r4, r6', r2, r4, r6
c
c  KAJ model FLO as a Slater and use moments to determine a cutoff range
c    for Slater of exponent a: <r^2> = 3/a^2; <r^4> = 45/(2*a^4); <r^6> = 9*7*5/a^6
c    -  e^-23 = 1E-10
c    -  cut off:  a*r_c = 23  -->  r_c = 23/a
c 
       a2 = sqrt(3.0d0/r2)
       a4 = (45.0d0/(2.0d0*r4))**(1.0d0/4.0d0)
       a6 = (315.0d0/r6)**(1.0d0/6.0d0)
       ALPFLO(1)=A2
       ALPFLO(2)=A4
       ALPFLO(3)=A6
       print *, 'Predicted Slater exponents a2, a4, a6 ', a2, a4, a6
       print *, 'predicted cutoffs ', 23.0/(2.0*a2), 23.0/(2.0*a4), 
     &   23.0/(2.0*a6)
       return
       end
        

