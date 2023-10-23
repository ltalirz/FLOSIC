       SUBROUTINE EXCITE(MODE)
C ORIGINAL VERSION BY MARK R PEDERSON (August 27 2004)  
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       PARAMETER (MAXOSQ=(MAX_OCC*(MAX_OCC+1))/2)
       PARAMETER (MXPR=MXPOISS)
       PARAMETER (MXLG=3)
       COMMON/COUPDATA/
     &  AIV(3,MXPR,MX_CNT+2,MXLG),AJV(3,MXPR,MX_CNT+2,MXLG)
     &  ,DMTV(10,10,MXPR,MX_CNT+2,MXLG),ALPIV(MXPR,MX_CNT+2,MXLG)
     &  ,ALPJV(MXPR,MX_CNT+2,MXLG),CENTER(3,MX_CNT+2)
     &  ,XDD(MAXUNSYM,MAXUNSYM,2),RVECI(3,MX_GRP),RVECJ(3,MX_GRP)
     &  ,NPAIRS(MX_CNT+2),IP(MX_CNT+2,MXLG)
       CHARACTER*80 LINE
       LOGICAL FIRST_ALPJ,CALLEDJ,NWRD,MWRD
       COMMON/TMP1/H(MAXOSQ,4),VEC(MAXUNSYM),HOL4(MAXUNSYM,MAXUNSYM,4)
       DIMENSION MAP(MAX_OCC),DOS(500,2),ENG(500),EXT(500),PHZ(500)
       DIMENSION NAP(MAX_OCC),LINE(2)
       DIMENSION IBEG(3),IEND(3),NDEG(3)
       DIMENSION AI(3),AJ(3),DIPTOT(3)
       DIMENSION SS(10,10),W4(10,10,4),C(3)
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA ZED/1.0D-30/
       DATA NDEG/1,3,6/
                     OPEN(80,FILE='DIPOLE')
                            DO K=1,3
                            DIPTOT(K)=0.0D0
                            END DO
                     READ(80,*,END=80)(DIPTOT(K),K=1,3)
 80                  CONTINUE
       IF(MODE.EQ.2)THEN   
C DEFIND ENERGY WINDOW
       call system('grep -i FERMI EVALUES >fermi.out')
       open(55,file='fermi.out')
       rewind(55)
         do jspn=1,1     
         read(55,55)line(jspn) 
         end do
       rewind(55)
         do jspn=1,1       
         write(55,55)line(jspn)(14:80)
         end do
       rewind(55)
       do jspn=1,1           
       read(55,*)efermi(jspn)
       print*,jspn,efermi(jspn),' read fermi'
       end do
       close(55)
 55    format(a)
                         EFERMI(NSPN)=EFERMI(1)
       EFMIN=MIN(EFERMI(1),EFERMI(NSPN))
       EFMAX=MAX(EFERMI(1),EFERMI(NSPN))
       EMIN=EFMIN-0.5D0 
       EMAX=EFMAX+0.5D0 
       EMIN=-0.25
       EMAX= 0.25
       CALL WFWIND(EMIN,EMAX,.TRUE.,.TRUE.,IFAIL)
       END IF
             KWF=(NWFS(1)+NWFS(NSPN))/(2/NSPN)
                    DO IWF=1,KWF
                    PRINT 22,IWF,EVLOCC(IWF)
                    END DO
 22                 FORMAT(' EXCITE:',I5,G15.6)
             PRINT*,'EXCITE: KWF,NWFS:',KWF,NWFS
             MAXOS=(KWF*(KWF+1))/2
                       DO I=1,4
                       DO J=1,MAXOS  
                       H(J,I)=0.0D0
                       END DO
                       END DO
       CALL GTTIME(TIME1)
       CHARGE=0.0D0
       ISITE=0
       ISHELLA=0
       DO 120 IFNCT=1,NFNCT
       DO 120 I_POS=1,N_POS(IFNCT)
        JSITE_BEG=ISITE
        ISHELLA=ISHELLA+1
        CALL OBINFO(1,RIDT(1,ISHELLA),RVECI,IST,ISHDUM)
C
C  TRANSLATE SHELLS TO ATOMS; GET GAUSS CUTOFF FOR ATOM I
C
        IATOM=ISHELLA
        GAUSS_CUTI=GAUSS_CUT(IATOM)
        DO 118 K_SITEI=1,IST
         ISITE=ISITE+1
         CALL UNRAVEL(IFNCT,ISHELLA,K_SITEI,RIDT(1,ISHELLA),
     &                RVECI,LST,1)
         IF(LST.NE.IST)THEN
          PRINT *,'COUPOT: PROBLEM IN UNRAVEL'
          CALL STOPIT
         END IF
         I_LOC_BEG=0
         AI(1)=RVECI(1,K_SITEI)
         AI(2)=RVECI(2,K_SITEI)
         AI(3)=RVECI(3,K_SITEI)
         JSITE=JSITE_BEG
         JSHELLA=0
         DO 116 JFNCT=1,NFNCT
         DO 116 J_POS=1,N_POS(JFNCT)
          JSHELLA=JSHELLA+1
          IF(JSHELLA.LT.ISHELLA) GO TO 116
          CALL OBINFO(1,RIDT(1,JSHELLA),RVECJ,JST,JSHDUM)
C
C  TRANSLATE SHELLS TO ATOMS;  GET GAUSS CUT OFF  FOR ATOM J
C
          JATOM=JSHELLA
          GAUSS_CUTJ=GAUSS_CUT(JATOM)
          JJCALL=0 
          DO 114 K_SITEJ=1,JST
           CALLEDJ=.FALSE.
           JSITE=JSITE+1
               CALL UNRAVEL(JFNCT,JSHELLA,K_SITEJ,RIDT(1,JSHELLA),
     &                      RVECJ,MST,2)
                                    IF(JST.NE.MST)THEN   
                                    PRINT*,'PROBLEM:IST.NE.JST'
                                    CALL STOPIT
                                    END IF
           IF(JSITE.LT.ISITE) GO TO 114
           J_LOC_BEG=0
           IF(JSITE.EQ.ISITE)J_LOC_BEG=I_LOC_BEG
           AJ(1)=RVECJ(1,K_SITEJ)
           AJ(2)=RVECJ(2,K_SITEJ)
           AJ(3)=RVECJ(3,K_SITEJ)
C ***************EXTRACTED FROM OVERLAP*****************
C
C calculate local indices:
c
           MAXIND=0
           DO LI=0,LSYMMAX(IFNCT)
            DO IBASE=1,N_CON(LI+1,IFNCT)
             DO MUI=1,NDEG(LI+1)
              MAXIND=MAXIND+1
C             LNDX(MUI,IBASE,LI+1,1)=MAXIND
             END DO
            END DO
           END DO
           MAXJND=0
           DO LJ=0,LSYMMAX(JFNCT)
            DO JBASE=1,N_CON(LJ+1,JFNCT)
             DO MUJ=1,NDEG(LJ+1)
              MAXJND=MAXJND+1
C             LNDX(MUJ,JBASE,LJ+1,2)=MAXJND
             END DO
            END DO
           END DO
                         DO K=1,4
                         DO INDEX=1,MAXIND
                         DO JNDEX=1,MAXJND
                         HOL4(JNDEX,INDEX,K)=0.0D0
                         END DO
                         END DO
                         END DO
           DO 230 IALP=1,N_BARE(IFNCT)
            ALPHAI=BFALP(IALP,IFNCT)
            DO 220 JALP=1,N_BARE(JFNCT)
             ALPHAJ=BFALP(JALP,JFNCT)
             ARG=(ALPHAI*ALPHAJ/(ALPHAI+ALPHAJ))
     &          *((AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2)
             IF (ARG .GT. CUTEXP) GOTO 220
              CALL GINTED(ALPHAI,ALPHAJ,AI,AJ,W4)
C             CALL OVMXSF(ALPHAI,ALPHAJ,AI,AJ,SS)
C                   DO K=1,4
C                   ERR=0.0D0
C                     DO I=1,10
C                     DO J=1,10 
C                     ERR=ERR+ABS(SS(J,I)-W4(J,I,K))
C                     END DO
C                     END DO
C                   PRINT*,'K, ERR:',K,ERR
C                   END DO
             DO 215 K=1,4
C K=1 <psi_i|x|psi_j>
C K=2 <psi_i|y|psi_j>
C K=3 <psi_i|z|psi_j>
C K=4 <psi_i|1|psi_j>
             INDEX=0
             DO 210 LI =0,LSYMMAX(IFNCT)
             DO 210 MUI=IBEG(LI+1),IEND(LI+1)
             DO 210 IBASE=1,N_CON(LI+1,IFNCT)
                INDEX=INDEX+1
                JNDEX=0
             DO 210 LJ =0,LSYMMAX(JFNCT)
             DO 210 MUJ=IBEG(LJ+1),IEND(LJ+1)
             DO 210 JBASE=1,N_CON(LJ+1,JFNCT)
              PROD=BFCON(IALP,IBASE,LI+1,IFNCT)
     &            *BFCON(JALP,JBASE,LJ+1,JFNCT)
              JNDEX=JNDEX+1
              HOL4(JNDEX,INDEX,K)=HOL4(JNDEX,INDEX,K)+PROD*W4(MUI,MUJ,K)
  210       CONTINUE
  215       CONTINUE
  220       CONTINUE
  230      CONTINUE
C          PRINT *,'AI:',AI
C          PRINT *,'AI:',AJ
C          PRINT*,'DIAGONAL OVERLAP MATRIX IN EXCITE:' 
C          PRINT 50,(ABS(HOL4(I,I)),I=1,MAX(MAXJND,MAXIND))
C50        FORMAT(' ',7F10.3)
****************END OF OVERLAP EXTRACTION******************
C    
               DO K=1,4
               IPTS=0 
               DO IWF=1  ,KWF                   
                  DO J_LOC=1,MAXJND
                   VEC(J_LOC)=0.0D0
                   DO I_LOC=1,MAXIND
                   VEC(J_LOC)=VEC(J_LOC)+
     &                 HOL4(J_LOC,I_LOC,K)*PSI(I_LOC,IWF,1)
                   END DO
                  END DO 
               DO JWF=IWF,KWF
                IPTS=IPTS+1
                DO J_LOC=1,MAXJND   
                H(IPTS,K)=H(IPTS,K)+VEC(J_LOC)*PSI(J_LOC,JWF,2)
                END DO
               END DO
               END DO
              IF(JSITE.GT.ISITE)THEN
               IPTS=0 
               DO IWF=1  ,KWF             
                  DO I_LOC=1,MAXIND
                   VEC(I_LOC)=0.0D0 
                   DO J_LOC=1,MAXJND
                   VEC(I_LOC)=VEC(I_LOC)+ 
     &                 HOL4(J_LOC,I_LOC,K)*PSI(J_LOC,IWF,2)
                   END DO
                  END DO
               DO JWF=IWF,KWF
                IPTS=IPTS+1
                DO I_LOC=1,MAXIND    
                H(IPTS,K)=H(IPTS,K)+VEC(I_LOC)*PSI(I_LOC,JWF,1)
                END DO
               END DO
               END DO
              END IF
              END DO
C
  105       CONTINUE
  110      CONTINUE
  114     CONTINUE
  116    CONTINUE
  118   CONTINUE
  120  CONTINUE
       CALL GTTIME(TIME2)
C CREATE TABLE FOR JOINT DOS:
       
C
C EMPTY TABLE OF PAIRS
C WRITE OUT OVERLAPS:
                     CHARGE=0.0D0       
                     IPTS=0
                     DO IWF=1  ,KWF
                     DO JWF=IWF,KWF
                     IPTS=IPTS+1
                     CHARGE=CHARGE+H(IPTS,4)
                     END DO
                     END DO
       ERRTOT=0.0D0
               PRINT*,'CHARGE IN EXCITE:',CHARGE,KWF
               IPTS=0
               DO IWF=1  ,KWF
               ERR=0.0D0
               DO JWF=IWF,KWF
               IPTS=IPTS+1
C                 PRINT*,'EXCITE:',IWF,JWF,H(IPTS,1)
                   IF(IWF.EQ.JWF)THEN
                   ERR=ERR+ABS(1.0D0-H(IPTS,4))
                   ELSE
                   ERR=ERR+ABS(      H(IPTS,4))
                   END IF
               END DO
C              PRINT*,'EXCITE IWF:',IWF,' ERROR:',ERR
               ERRTOT=ERRTOT+ERR
               END DO
C ORDER THE WAVEFUNCTIONS:
                DO IWF=1,KWF
                 MAP(IWF)=IWF
                 DO JWF=1,IWF-1   
                 IF(EVLOCC(MAP(JWF)).GT.EVLOCC(MAP(IWF)))THEN
                                 MM=MAP(IWF)
                                 MAP(IWF)=MAP(JWF)
                                 MAP(JWF)=MM
                 END IF
                 END DO
                END DO
                DO IWF=1,KWF
                PRINT*,MAP(IWF),EVLOCC(MAP(IWF)),' 4-27-05'
                END DO 
       PRINT*,'TIME IN EXCITE:',TIME2-TIME1,' ERR,CHG:',ERRTOT,CHARGE
       OPEN(40,FILE='DIPINF')
       REWIND(40)
       WRITE(40,*)'1 DIPOLES FOR EACH STATE'
            EFRM=MAX(EFERMI(1),EFERMI(NSPN))
            TEMP=0.0001
            DO IWF=1  ,KWF
            DO JWF=IWF,KWF
            JPTS=JPTS+1
            IF(IWF.EQ.JWF)THEN
                    IVF=MAP(IWF)
                    JVF=MAP(JWF)
                 MXK=MAX(JVF,IVF)
                 MNK=MIN(JVF,IVF)
                     MBG=1
                     DO MWF=1,MNK-1
                     MBG=MBG+KWF-(MWF-1)
                     END DO
                 IPTS=(MXK-MNK)+MBG             
            OCC=FERMFN(EVLOCC(IVF),EFRM,TEMP)
            WRITE(40,400)IVF,EVLOCC(IVF),OCC,(H(IPTS,K),K=1,3)
            END IF
            END DO
            END DO
 400   FORMAT(' ',I8,5F12.6,' E,Q,P')
       IPTS=0
       DO IWF=1  ,KWF
       DO JWF=IWF,KWF
       IPTS=IPTS+1
               H(IPTS,4)=0.0D0                        
       END DO
       END DO
       FF=2.5415803**2 !(AU TO DEBYE)
       DO K=1,3
       IPTS=0
       DO IWF=1  ,KWF
       DO JWF=IWF,KWF
       IPTS=IPTS+1
               H(IPTS,4)=H(IPTS,4)+FF*H(IPTS,K)**2
       END DO
       END DO
       END DO
       WRITE(40,*)'1 DIPOLE  MATRIX ELEMENTS:'   
       WRITE(40,*)'E1, E2, Q1, Q2, TM, DX, DY, DZ'
       DO IWF=1  ,KWF
       DO JWF=IWF,KWF 
                    IVF=MAP(IWF)
                    JVF=MAP(JWF)
                 MXK=MAX(JVF,IVF)
                 MNK=MIN(JVF,IVF)
                     MBG=1
                     DO MWF=1,MNK-1
                     MBG=MBG+KWF-(MWF-1)
                     END DO
                 IPTS=(MXK-MNK)+MBG             
        OCI=FERMFN(EVLOCC(IVF),EFRM,TEMP)
        OCJ=FERMFN(EVLOCC(JVF),EFRM,TEMP)
       PRINT     23,EVLOCC(IVF),EVLOCC(JVF),OCI,OCJ,
     &      H(IPTS,4),(H(IPTS,J),J=1,3)
       WRITE(40,23)EVLOCC(IVF),EVLOCC(JVF),OCI,OCJ,
     &      H(IPTS,4),(H(IPTS,J),J=1,3)
       END DO 
       END DO 
 23    FORMAT(' ',2F12.3,3F8.3,3G15.6)
       CLOSE(40)
       DO 550 IPLT=1,2
       IF(IPLT.EQ.1)THEN
        OPEN(40,FILE='QJNTOUT')
        DENG=0.01 !0.27 eV broadening for joint dos 
       ELSE
        OPEN(40,FILE='EXTCOEF')
             DNU=0.05 ! (1/20) PetaHz Broadening
             PHZ2AU=6.626/1.602/27.2118
        DENG=PHZ2AU*DNU
       END IF
       REWIND(40)
                EEE=0.0D0    
                NNN=400
                DLTE=1.0D0/(NNN-1)
                DO I=1,NNN   
                EEE=EEE+DLTE 
                ENG(I)=EEE
                DOS(I,1)=0.0D0
                DOS(I,2)=0.0D0
                END DO
        BPOL=0.0D0
        IPTS=0
        DO IWF=1  ,KWF
        PRINT*,'IWF,EVLOCC(IWF):',IWF,EVLOCC(IWF)
        DO JWF=IWF,KWF
        IPTS=IPTS+1
                 IF(IWF.EQ.JWF)NAP(IWF)=IPTS
        OCI=FERMFN(EVLOCC(IWF),EFRM,TEMP)
        OCJ=FERMFN(EVLOCC(JWF),EFRM,TEMP)
                                ISN=2
                                JSN=2
              IF(IWF.LE.NWFS(1))ISN=1
              IF(JWF.LE.NWFS(1))JSN=1
               FCT=OCI*(1.0D0-OCJ)
               IF(ISN.NE.JSN)FCT=0.0D0
               DIP=FCT*H(IPTS,4)*(2/NSPN)
C UNITS ARE liters/mole/cm for extinction coefficient
               DFF=EVLOCC(IWF)-EVLOCC(JWF)  
               IF(FCT.GT.0.00001)THEN
               BPOL=BPOL+DIP/DFF
               PRINT 2313,IWF,JWF,DIP,DFF,BPOL,OCI,OCJ
               PRINT*,'NNN:',NNN
 2313   format(2i5,5g15.6)
                  DO I=1,NNN  
                  ARG=((DFF+ENG(I))/DENG)**2
                  ARG=EXP(-ARG)
                       IF(IPLT.EQ.2)THEN
                            ARG=(ARG/FF)*327.65*ENG(I)/DENG
                       END IF
                  DOS(I,1)=DOS(I,1)+DIP*ARG
                  DOS(I,2)=DOS(I,2)+FCT*ARG
                  END DO
               END IF
        END DO
        END DO
        IF(IPLT.EQ.1)THEN
        PRINT*,'BARE POLARIZABILITY, FERMI-LEVEL:',BPOL,EFRM,NSPN
C       BPOL=BPOL*(2/NSPN)
        PRINT*,'BARE POLARIZABILITY, FERMI-LEVEL:',BPOL,EFRM,NSPN
        END IF
        DO I=1,NNN  
                  IF(IPLT.EQ.1)XENG=ENG(I)*27.2118
                  IF(IPLT.EQ.2)XENG=ENG(I)/PHZ2AU 
        WRITE(40,500)XENG,DOS(I,1)!,DOS(I,2)
        END DO
        CLOSE(40)
 550    CONTINUE
        OPEN(40,FILE='DIPDFF')
                  DO I=1,NNN  
                  DOS(I,1)=0.0D0                  
                  DOS(I,2)=0.0D0                 
                  END DO
        DO IWF=1  ,KWF
        IVF=MAP(IWF)
        IPI=NAP(IVF)
        DO JWF=IWF,KWF
        JVF=MAP(JWF)
        IPJ=NAP(JVF)
        OCI=FERMFN(EVLOCC(IVF),EFRM,TEMP)
        OCJ=FERMFN(EVLOCC(JVF),EFRM,TEMP)
            FCT=OCI*(1.0D0-OCJ)
            DO J=1,3
            DIPOLE(J)=H(IPI,J)-H(IPJ,J)
            END DO
            DPT=DIPOLE(1)**2+ 
     &          DIPOLE(2)**2+
     &          DIPOLE(3)**2 
            DPT=SQRT(DPT)
        WRITE(40,510)EVLOCC(IVF),EVLOCC(JVF),FCT,DPT,(DIPOLE(J),J=1,3)
               DFF=EVLOCC(IVF)-EVLOCC(JVF)
               IF(FCT.GT.0.00001)THEN
               DPT=DPT*FCT
                  DO I=1,NNN  
                  ARG=((DFF+ENG(I))/0.01)**2
                  DOS(I,1)=DOS(I,1)+DPT*EXP(-ARG)
                  END DO
               END IF
        END DO
        END DO
        CLOSE(40)
        OPEN(40,FILE='DIPPLT')
        DO I=1,NNN  
        WRITE(40,500)27.2110*ENG(I),DOS(I,1)!,DOS(I,2)
        END DO
        CLOSE(40)
        call dipstr
 500    FORMAT(' ',3F12.4) 
 510    FORMAT(' ',7F12.4) 
       RETURN       
       END
       FUNCTION FERMFN(EI,EF,TEMP)
       IMPLICIT REAL*8 (A-H,O-Z)
            ARG=MIN(ABS((EI-EF)/TEMP),30.0D0)
            IF(EI.LT.EF)ARG=-ARG
            FERMFN=1.0D0/(1.0D0+EXP(ARG))             
       RETURN
       END 
       subroutine dipstr      
       parameter (mx=300)
       implicit real*8 (A-H,O-Z)
           dimension d(mx),e(mx)
           character*80 line
           open(80,file='DIPINF')
           rewind(80)
  7        continue 
           read(80,5)line
  5        format(a)
                ifnd=0
                do i=1,79
                if(line(i:i+1).eq.'DX')ifnd=1
                end do
           if(ifnd.eq.0)go to 7
           k=1
  10       continue       
           read(80,*,end=20)e1,e2,q1,q2,dum,dx,dy,dz 
           factor=q1*(1.-q2)
           if(abs(factor).gt.0.1.and.e2-e1.lt.1.0)then
           e(k)=(e2-e1)*27.2118
           d(k)=(dx*dx+dy*dy+dz*dz)*2.5415803**2
               k=k+1
           end if
           if(k.eq.mx)then 
                      do i=1,k
                       do j=i+1,k
                                if(e(j).lt.e(i))then
                                    es=e(j)
                                    e(j)=e(i)
                                    e(i)=es
                                    ds=d(j)
                                    d(j)=d(i)
                                    d(i)=ds
                                end if
                       end do
                      end do
                  k=0.6*k
                  go to 10
           else
                 go to 10
           end if
 20        continue
           open(60,file='DIPSTR')
           rewind(60)
           print*,'found:',k,' excitations'
                       do i=1,k
                       do j=i+1,k
                            if(e(j).lt.e(i))then
                                ds=d(j)
                                d(j)=d(i)
                                d(i)=ds
                                ds=e(j)
                                e(j)=e(i)
                                e(i)=ds
                            end if
                       end do
                       ff=1.602/6.626
                       write(60,30)e(i),e(i)*ff,d(i)
                       end do
           close(80)
           close(60)
 30        format(' ',3g15.6,' Eng (eV), Eng (Hz), D(debye^2)')
           end
c  -18.806  -18.806    1.000    1.000    8.307   -7.030    2.057    0.000
c  -18.806  -18.799    1.000    1.000    0.000    0.000    0.000    0.000
c  -18.806  -18.790    1.000    1.000    0.000    0.000    0.000    0.000
c  -18.806  -18.740    1.000    1.000    0.000    0.000    0.000    0.000
c  -18.806  -14.040    1.000    1.000    0.000    0.000    0.000    0.000
