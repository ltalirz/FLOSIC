      SUBROUTINE UPDFOV(MXFOD,NSPN,NFOD,FOD,GAM,ALPHAI,ALPHAJ,AI,AJ,FO)
C WRITTEN BY MARK R PEDERSON 24-FEB 2022
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (MPFOD=1000)
        DIMENSION FO(MXFOD,MXFOD,2)
        DIMENSION FOD(3,MXFOD,2)
        DIMENSION AI(3),AJ(3),P1(10,MPFOD,2),P2(10,MPFOD,2),GAM(10,10,2)
        DIMENSION ARG1(MPFOD,2),ARG2(MPFOD,2)
        DIMENSION NFOD(2),P3(10,MPFOD)
C
C NOTE THE MU AND NU CORRESPOND TO INDICES OF FODS --- NORMALLY WRITTEN AS {\bf a_i} and {\bf a_j}.
C
C FO SHOULD BE EQUAL TO ZERO ON FIRST CALL
C
C MXFOD      -- PHYSICAL DIMENSION OF FOD
C NFOD(ISPN) -- NUMBER OF FODS FOR SPIN ISPN
C NSPN       -- 1 or 2 for unpolarized or polarized
C FOD(J,IFOD,ISPN) -- Coordinates of FODs
C GAM -- density matrix for two bare gaussians (alphi,ai) and (alphaj, aj)
C FO -- density matrix for each spin at positions alphai and alphaj
C
C
                      IF(NFOD(1).GT.MPFOD.OR.NFOD(NSPN).GT.MPFOD)THEN
                      PRINT*,'FIX PARAMETER',MPFOD,' IN UPDATE_DMT'
                      PRINT*,'MUST BE GREATER THAN:',NFOD(1),NFOD(NSPN) 
                      END IF
        !PRINT*,'FOD' ,FOD
C        PRINT*,' ALPHAI ALPHAJ',  ALPHAI, ALPHAJ
        DO ISPN=1,NSPN
        DO MU=1,NFOD(ISPN)
            P1(1,MU,ISPN)=1.0D0                       
            P2(1,MU,ISPN)=1.0D0
            DO J=2,4
            P1(J,MU,ISPN)=FOD(J-1,MU,ISPN)-AI(J-1)  !px, py, pz
            P2(J,MU,ISPN)=FOD(J-1,MU,ISPN)-AJ(J-1)  !px, py, pz
            END DO
            DO J=5,7
            P1(J,MU,ISPN)=P1(J-3,MU,ISPN)*P1(J-3,MU,ISPN) !dxx,dyy,dzz
            P2(J,MU,ISPN)=P2(J-3,MU,ISPN)*P2(J-3,MU,ISPN) !dxx,dyy,dzz
            END DO
           P1(8 ,MU,ISPN)=P1(2,MU,ISPN)*P1(3,MU,ISPN)     !dxy
           P1(9 ,MU,ISPN)=P1(2,MU,ISPN)*P1(4,MU,ISPN)     !dxz
           P1(10,MU,ISPN)=P1(3,MU,ISPN)*P1(4,MU,ISPN)     !dyz
           P2(8 ,MU,ISPN)=P2(2,MU,ISPN)*P2(3,MU,ISPN)     !dxy
           P2(9 ,MU,ISPN)=P2(2,MU,ISPN)*P2(4,MU,ISPN)     !dxz
           P2(10,MU,ISPN)=P2(3,MU,ISPN)*P2(4,MU,ISPN)     !dyz
           ARG1(MU,ISPN)=ALPHAI*(P1(5,MU,ISPN)
     &                          +P1(6,MU,ISPN)
     &                          +P1(7,MU,ISPN)) !r^2
           ARG2(MU,ISPN)=ALPHAJ*(P2(5,MU,ISPN)
     &                          +P2(6,MU,ISPN)
     &                          +P2(7,MU,ISPN)) !r^2
           ARG1(MU,ISPN)=EXP(-ARG1(MU,ISPN))      !gaussian
           ARG2(MU,ISPN)=EXP(-ARG2(MU,ISPN))      !gaussian
        END DO
        END DO
        DO ISPN=1,NSPN
             P3=0.0D0
             DO I=1,10
             DO NU=1,NFOD(ISPN)
             DO J=1,10
               P3(I,NU)=P3(I,NU)+P2(J,NU,ISPN)*GAM(I,J,ISPN)
             END DO
             END DO
             END DO
        DO MU=1,NFOD(ISPN)
        DO NU=1,NFOD(ISPN)
           ENV=ARG1(MU,ISPN)*ARG2(NU,ISPN)
           TT=0.0D0
           DO I=1,10
C           DO J=1,10
C           TT=TT+P1(I,MU,ISPN)*P2(J,NU,ISPN)*GAM(I,J,ISPN)           
            TT=TT+P1(I,MU,ISPN)*P3(I,NU)
           END DO
C           END DO
           FO(MU,NU,ISPN)=FO(MU,NU,ISPN)+TT*ENV 
        END DO
        END DO
        END DO


       RETURN
       END   

C DENSOLD BASED ON OLD VERSION OF APOTNL BY M. PEDERSON AND D. POREZAG
C ATTENTION: FIRST TWO ARRAYS OF COMMON BLOCK TMP1 MUST BE IDENTICAL IN 
C DENSOLD AND APOTNL SINCE THEY ARE USED TO PASS DENSITY AND COULOMB POT
C KW: THis gives Fermi orbital transfodmation 04/15/2022
       SUBROUTINE GETTMAT(TOTQNUM)!(FOD,NFOD,MXFOD)       
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       COMMON/FORRx/NFLO,KSPX!!,TMAT(NDH,NDH,2)
       COMMON/FORRX2/MEQV
       COMMON/TSTCHG/TOTCHG,COUSIC
       DIMENSION FMAT(NDH,NDH)
C      DIMENSION OVER2(NDH,NDH)
C      REAL*8 FOD(3,MXFOD,2)!,TMAT(MXFOD,MXFOD,2)
C      INTEGER NFOD(2),MAP(2,2*MXFOD)
       INTEGER         MAP(2,MAX_OCC)
       PARAMETER (NMAX=MPBLOCK)
C
C RETURN:
C RHOG(IPTS,1, 1)= rho_up   
C RHOG(IPTS,2, 1)= d rho_up/dx
C RHOG(IPTS,3, 1)= d rho_up/dy
C RHOG(IPTS,4, 1)= d rho_up/dz
C RHOG(IPTS,5, 1)= d^2 rho_up/dx^2
C RHOG(IPTS,6, 1)= d^2 rho_up/dy^2
C RHOG(IPTS,7, 1)= d^2 rho_up/dz^2
C RHOG(IPTS,8, 1)= d^2 rho_up/dxdy
C RHOG(IPTS,9, 1)= d^2 rho_up/dxdz
C RHOG(IPTS,10,1)= d^2 rho_up/dydz
C RHOG(IPTS,1, 2)= rho_dn   
C RHOG(IPTS,2, 2)= d rho_dn/dx
C RHOG(IPTS,3, 2)= d rho_dn/dy
C RHOG(IPTS,4, 2)= d rho_dn/dz
C RHOG(IPTS,5, 2)= d^2 rho_dn/dx^2
C RHOG(IPTS,6, 2)= d^2 rho_dn/dy^2
C RHOG(IPTS,7, 2)= d^2 rho_dn/dz^2
C RHOG(IPTS,8, 2)= d^2 rho_dn/dxdy
C RHOG(IPTS,9, 2)= d^2 rho_dn/dxdz
C RHOG(IPTS,10,2)= d^2 rho_dn/dydz
C
       LOGICAL ICOUNT,EXIST
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,NVGRAD,MXSPN)
       COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
        COMMON/FLOINFO/FOD(3,MAX_OCC,2),NFOD(2),MFOD(2)
        COMMON/WTSFOD/WTFOD
        COMMON/SIC_ENERGY/SICENERGY,TOT_SIC
        COMMON/FLOMESH/RABCD(3,1000)
        LOGICAL FEXIST,FMESH
        DATA FMESH/.TRUE./
        DIMENSION NT(MAX_OCC,2),IFND(MX_GRP,MAX_OCC,2)

C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       LOGICAL LGGA,IUPDAT
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
        
       call symfrm(1)
                 NFLO=0
        INQUIRE(FILE='FRMORB',EXIST=EXIST)
        IF(EXIST)THEN
        OPEN(80,FILE='FRMORB')
        READ(80,*)(MFOD(I),NFOD(I),I=1,2)
        PRINT*,'NFOD:',NFOD
        DO I=1,NSPN
        DO J=1,MFOD(I)
        READ(80,*)(FOD(K,J,I),K=1,3),NT(J,I),(IFND(K,J,I),K=1,NT(J,I))
        PRINT 40,(FOD(K,J,I),K=1,3),NT(J,I),(IFND(K,J,I),K=1,NT(J,I))
        END DO
C       END DO
C       DO I=1,NSPN
        KOUNT=MFOD(I)
        DO J=1,MFOD(I)
        DO L=2,NT(J,I)
        KOUNT=KOUNT+1
        READ(80,*)(FOD(K,IFND(L,J,I),I),K=1,3)
        PRINT 40, (FOD(K,IFND(L,J,I),I),K=1,3)
        END DO
        END DO
            IF(KOUNT.NE.NFOD(I))THEN 
            PRINT*,'FIX FRMORB FILE'
            CALL STOPIT
            END IF
        END DO
        CLOSE(80)
        END IF
 40    FORMAT(3F12.6,I5,60I4)
C
       IF(FMESH) THEN
       RABCD=0.0D0
       INQUIRE(FILE='ABCD',EXIST=FEXIST) 
       IF(FEXIST) THEN
        OPEN(80,FILE='ABCD')
        DO KSPN=1,NSPN
        DO I=1,MFOD(KSPN)
         READ(80,*)(RABCD(J,I+MFOD(1)*(KSPN-1)),J=1,2)
        END DO
        END DO
        CLOSE(80)
       ELSE
         RABCD(1,:)=1700.0D0
         RABCD(2,:)= 150.0D0
       END IF 
       END IF
       FMESH=.FALSE.



       TIMEGORB=0.0D0
       TMAT=0.0D0
       CALL GTTIME(APT1)
       LGGA= .FALSE.
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) THEN
        LGGA= .TRUE.
        NGRAD=10
       END IF
C
C LOOP OVER ALL POINTS
C
       NMSH_save=NMSH
       NMSH=0
       PRINT*,'NWF NWFS',NWF,NWFS 
C KW: Read FODs as RMSH        
       DO KSPN=1,NSPN
       DO KFOD=1,NFOD(KSPN)
        NMSH=NMSH+1
        MAP(1,NMSH)=KFOD
        MAP(2,NMSH)=KSPN
        RMSH(1,NMSH)=FOD(1,KFOD,KSPN)
        RMSH(2,NMSH)=FOD(2,KFOD,KSPN)
        RMSH(3,NMSH)=FOD(3,KFOD,KSPN)
       END DO
       END DO
        



       LPTS=0
 10    CONTINUE
        IF(LPTS+NMAX.LT.NMSH)THEN
         MPTS=NMAX
        ELSE
         MPTS=NMSH-LPTS
        END IF

C
C INITIALIZE PSIG AND RHOB
C
         
C        DO IWF=1,NWF
C         DO IGR=1,NGRAD
C          DO IPTS=1,MPTS
C           PSIG(IPTS,IGR,IWF)=0.0D0
C          END DO
C         END DO  
C        END DO  
C        DO ISPN=1,NSPN
C         DO IGR=1,NGRAD
C          DO IPTS=1,MPTS
C           RHOG(LPTS+IPTS,IGR,ISPN)=0.0D0
C          END DO
C         END DO  
C        END DO  
        ISHELLA=0
C
        PSIG=0.0D0
        RHOG=0.0D0
C FOR ALL CENTER TYPES
C
        DO 86 IFNCT=1,NFNCT
         LMAX1=LSYMMAX(IFNCT)+1
C
C FOR ALL POSITIONS OF THIS CENTER
C
         DO 84 I_POS=1,N_POS(IFNCT)
          ISHELLA=ISHELLA+1
C
C GET SYMMETRY INFO
C
          CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
          IF(NWF.GT.MAX_OCC)THEN
           PRINT *,'APTSLV: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
          DO 82 J_POS=1,M_NUC
C
C UNSYMMETRIZE 
C
           CALL RxRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                  RVECA,L_NUC,1)
C          CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
C    &                  RVECA,L_NUC,1)
           IF(L_NUC.NE.M_NUC)THEN
            PRINT *,'APTSLV: PROBLEM IN UNRAVEL'
            CALL STOPIT
           END IF
C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
           KPTS=0
           DO 80 JPTS=1,MPTS,NSPEED
            NPV=MIN(NSPEED,MPTS-JPTS+1)
            DO LPV=1,NPV
             KPTS=KPTS+1
             PTS(LPV,1)=RMSH(1,LPTS+KPTS)-RVECA(1,J_POS)
             PTS(LPV,2)=RMSH(2,LPTS+KPTS)-RVECA(2,J_POS)
             PTS(LPV,3)=RMSH(3,LPTS+KPTS)-RVECA(3,J_POS)
            END DO
C
C GET ORBITS AND DERIVATIVES
C
            NDERV=0
            IF (LGGA) NDERV=2
            CALL GTTIME(TIME3)
            CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            CALL GTTIME(TIME4)
            TIMEGORB=TIMEGORB+TIME4-TIME3
C
C UPDATING ARRAY PSIG
C
            IF (IUPDAT) THEN
             IPTS=JPTS-1
             ILOC=0
             DO 78 LI=1,LMAX1
              DO MU=1,ISIZE(LI)
               DO ICON=1,N_CON(LI,IFNCT)
                ILOC=ILOC+1
                IF (ICOUNT(ICON,LI)) THEN
                 DO IWF=1,NWF
                  FACTOR=PSI(ILOC,IWF,1)
                  DO IGR=1,NGRAD
                   DO LPV=1,NPV
                    PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF)
     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
                   END DO
                  END DO  
                 END DO  
                END IF
               END DO  
              END DO  
   78        CONTINUE
            END IF
   80      CONTINUE
   82     CONTINUE
   84    CONTINUE
   86   CONTINUE
C
C UPDATING RHOG, START WITH DENSITY 
C
        DO ISPN=1,NSPN
         JBEG= (ISPN-1)*NWFS(1) 
         DO JWF=1,NWFS(ISPN)
          JLOC=JWF+JBEG
          DO IPTS=1,MPTS
           RHOG(LPTS+IPTS,1,ISPN)=RHOG(LPTS+IPTS,1,ISPN)
     &     +PSIG(IPTS,1,JLOC)**2
          END DO
         END DO
         PRINT*,'RHOG',(RHOG(LPTS+IPTS,1,ISPN),IPTS=1,MPTS)
        END DO
C
C KW : TMAT (Fermi orbital transformation)
          DO IPTS=1,MPTS
           KSPN=MAP(2,IPTS+LPTS)
           KFOD=MAP(1,IPTS+LPTS)
           KFOD_MAX=MAX(KFOD,KFOD_MAX)
           DO KWF=1,NWFS(KSPN)
            KLOC=KWF+(KSPN-1)*NWFS(1)
            KLOC_MAX=MAX(KLOC,KLOC_MAX)
C            TMAT(KFOD,KLOC,KSPN)=SQRT(RHOG(LPTS+IPTS,1,KSPN))
            TMAT(KFOD,KWF,KSPN)=PSIG(IPTS,1,KLOC)/
     &                    SQRT(RHOG(LPTS+IPTS,1,KSPN))
           END DO
          END DO
C
        LPTS=LPTS+MPTS
        IF (LPTS .LT. NMSH) GOTO 10
       CONTINUE
       NMSH=NMSH_save
       PRINT*,'NWF NWFS NFOD ',NWF, NWFS,NFOD
       PRINT*,'KFOD_MAX KLOC_MAX',KFOD_MAX,KLOC_MAX
       DO ISPN=1,NSPN
       OVER=0.0D0
       PRINT*,'NFOD:',NFOD(ISPN)
          DO MU=1,NFOD(ISPN)
          DO NU=1,NFOD(ISPN)
           DO IW=1,NWFS(ISPN)
           OVER(MU,NU)=OVER(MU,NU)
     &    +TMAT(MU,IW,ISPN)*TMAT(NU,IW,ISPN)
           END DO
           END DO
         PRINT 20,(OVER(MU,NU),NU=1,5)
         END DO
       CALL LOWSIC(ISPN,NDH,NFOD(ISPN),OVER,HAM,FMAT,EVAL,SC1)
C FLO_NU = SUM(MU) HAM(MU,NU)*FO_MU
C FLO_NU = SUM(MU)SUM(IW) HAM(MU,NU)*TMAT(MU,IW,ISPN)*|KS_IW>
             PRINT *, 'IN FODMAT, ISPN, NWFS(ISPN)', ISPN,NWFS(ISPN)
             OVER=0.0D0
             DO NU=1,NFOD(ISPN)
             DO IW=1,NWFS(ISPN)
                   DO MU=1,NFOD(ISPN)
                   OVER(IW,NU)=OVER(IW,NU)+HAM(MU,NU)*TMAT(MU,IW,ISPN) 
                   END DO
             END DO
             PRINT 20,(OVER(IW,NU),IW=1,NWFS(ISPN))
             END DO
             DO NU=1,NFOD(ISPN)
             DO IW=1,NWFS(ISPN)
             TMAT(IW,NU,ISPN)=OVER(IW,NU)
             END DO
             END DO
C FOR NORMAL CASES THIS IS A UNITARY MATRIX:
             HAM=0.0D0
             DO NU=1,NFOD(ISPN)
             DO MU=1,NFOD(ISPN)
                DO IW=1,NWFS(ISPN)
                HAM(NU,MU)=HAM(NU,MU)+TMAT(IW,NU,ISPN)*TMAT(IW,MU,ISPN)
                END DO
             END DO
             PRINT 20,(HAM(NU,MU),MU=1,5)
             END DO
       END DO
C REPAIR DAMAGE DONE TO VMESH
          PRINT*,'NMSH:',NMSH
          OPEN(99,FILE='VMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
          REWIND(99)
          READ(99) NMSH,JCALC
          READ(99)((RMSH(J,I),J=1,3),I=1,NMSH)
          CLOSE(99)
          CALL APOTNL(TOTQNUM)
C Save XC DFT energies, because GETVLXC is called in SIC, copy them back
C at the end 
          E1=ERGXL
          E2=ERGXN
          E3=ERGCL
          E4=ERGCN
         IDEBUG=0  !!!
       IF(IDEBUG.EQ.1)THEN
         PRINT*,'CHECKING SIC CALC...SET IDEBUG=0 To SAVE TIME'
         SIC_ENG=0.0D0
         TOT_CHG=0.0D0
         DO LSPX=1,NSPN
         KSPX=LSPX
         DO IFLO=1,MFOD(LSPX)
         ORB_CHG=0.0D0
         ORB_COU=0.0D0
         ORB_EXC=0.0D0
c        print*,'NT',NT(IFLO,LSPX)
         DO JFLO=1,NT(IFLO,LSPX)
              NFLO=-IFND(JFLO,IFLO,LSPX)
                COULOMB=0.0D0
                RHOG   =0.0D0
                CALL COUPOT1
              PI=4.0D0*ATAN(1.0D0)
              DO IPTS=1,NMSH
                        RV=RHOG(IPTS,1,LSPX)*WMSH(IPTS)
         ORB_CHG=ORB_CHG+RV                         
         ORB_COU=ORB_COU+RV*0.25*COULOMB(IPTS)
         ORB_EXC=ORB_EXC+RV*0.75*((6.0D0/PI)*RHOG(IPTS,1,LSPX))**(1./3.)
              XEKS=0.75D0*((6.0d0/PI)*RHOG(IPTS,1,LSPX))**(1./3.)
              SIC_ENG=SIC_ENG-(0.25*COULOMB(IPTS)-0.0*XEKS)*RV
 
              TOT_CHG=TOT_CHG+RV!           RHOG(IPTS,1,LSPX)*WMSH(IPTS)
              END DO
         END DO
                 ORB_COU=ORB_COU/NT(IFLO,LSPX)
                 ORB_EXC=ORB_EXC/NT(IFLO,LSPX)
                 ORB_CHG=ORB_CHG/NT(IFLO,LSPX)
         PRINT 44,LSPX,KSPX,ABS(NFLO),NMSH,
     &  SIC_ENG,TOT_CHG,ORB_CHG,ORB_COU,ORB_EXC
         END DO
 44      FORMAT(3I3,I8,' DEBUG (OK BY PUBLIC FLOSIC 6/22/2022)',6F12.6)
         END DO
C      END DO
       END IF
       OPEN(50,FILE='SICSCI',FORM='UNFORMATTED',STATUS='UNKNOWN')
       DO KSPX=1,NSPN
          TOTCHG=0.0D0
          COUSIC=0.0D0
             DO INDX=1,NDH_TOT
             HSTOR(INDX,2)=0.0D0
C
C kaj 6-24
c  initialize HSTOR(INDX,1) to compute fod forces
C
             HSTOR(INDX,1)=0.0d0
             END DO
       FOUND=0.0D0
       DO KFLO=1,MFOD(KSPX)  !!!MFOD(KSPX)  !!fix this later
          IFLO=KFLO
          NFLO=-IFLO
          PRINT*,'IFLO:',IFLO, "MADE IT HERE"
                       COULOMB=0.0D0
                        RHOG   =0.0D0
                        IF(IFLO.EQ.1.AND.KSPX.EQ.1) SICENERGY=0.0D0
                        MEQV=NT(IFLO,KSPX)
                        call coupot1
                        SICENERGY=SICENERGY+TOT_SIC*MEQV
       END DO !kflo
C
                print*,'CONVERGED',CONVERGENCE,ITSCF,MAXSCF
c  kaj 7-12-23
c
c uncomment to get fod forces  ***************************************************
                IF(CONVERGENCE.OR.ITSCF.EQ.MAXSCF) THEN
                 CALL GET_EPS(kspx)  ! ,MFOD(kspnx))
                END IF
c end change ***************************************************
                 DO INDX=1,NDH_TOT 
                   HSTOR(INDX,2)=(HSTOR(INDX,2)+HSTOR(INDX,1))*0.5D0
                 END DO
                WRITE(50)(HSTOR(INDX,2),INDX=1,NDH_TOT)
                call overlap(1)
c
c  kaj 7-12-23   get fod forces
c
c  uncomment  to get forces **********************************
       IF(CONVERGENCE.OR.ITSCF.EQ.MAXSCF) THEN
       call frmorb2(kspx)
       end if  
       END DO   !kspx
       print *,'finished getting forces'
c      call stopit
                call overlap(1)
c end change ***************************************************
c KAJ  remove when working 
C Adjusting parameters for the ellipsoid 
       OPEN(51,FILE='ABCD')
       DO KSPN=1,NSPN
       DO KFLO=1,MFOD(KSPN)
       IF((1.0D0-RABCD(3,KFLO+MFOD(1)*(KSPN-1))).GT.0.0000009D0) THEN
       DO J=1,2
       RABCD(J,KFLO+MFOD(1)*(KSPN-1))=RABCD(J,KFLO+MFOD(1)*(KSPN-1))
     &                              *1.1D0
       END DO
       END IF
       IF((1.0D0-RABCD(3,KFLO+MFOD(1)*(KSPN-1))).LT.0.0000001D0) THEN
       DO J=1,2
       RABCD(J,KFLO+MFOD(1)*(KSPN-1))=RABCD(J,KFLO+MFOD(1)*(KSPN-1))
     &                              *0.95D0
       END DO
       END IF
       WRITE(51,*)(RABCD(J,KFLO+MFOD(1)*(KSPN-1)),J=1,3)
       END DO
       END DO
       CLOSE(51)


       JNDX=0
       DO KREP=1,N_REP
       NBAS=NS_TOT(KREP)
          HADD=0.0D0
          OADD=0.0D0
          FKEV=0.01D0
        DO IBAS=1   ,NBAS        
        DO JBAS=IBAS,NBAS            
        JNDX=JNDX+1
        HAM (IBAS,JBAS)=HSTOR(JNDX,2)
        HAM (JBAS,IBAS)=HSTOR(JNDX,2)
        OVER(IBAS,JBAS)=HSTOR(JNDX,1)
        OVER(JBAS,IBAS)=HSTOR(JNDX,1)
        HAM(IBAS,JBAS)=HAM(IBAS,JBAS)+FKEV*OVER(IBAS,JBAS)
        HAM(JBAS,IBAS)=HAM(IBAS,JBAS)
        HADD=HADD+ABS(HAM (IBAS,JBAS))
        OADD=OADD+ABS(OVER(IBAS,JBAS))
        END DO
C       PRINT 556,(HAM(IBAS,JBAS),JBAS=1,NBAS)
        END DO
        IF(OADD.NE.OADD) THEN 
         PRINT*,'OADD IS NAN'
         CALL STOPIT 
        END IF
        IF(HADD.NE.HADD) THEN 
         PRINT*,'HADD IS NAN'
         CALL STOPIT 
        END IF
CCCC           CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
           DO IBAS=1,NBAS
           EVAL(IBAS)=EVAL(IBAS)-FKEV
           FOUND=FOUND+EVAL(IBAS)*NDMREP(KREP)/LDMREP(KREP)
           END DO
       PRINT*,'EIGENVALUES:',IREP,NBAS,JNDX,HADD,OADD
       PRINT 556,(EVAL(IBAS),IBAS=1,NBAS)
       END DO
       CLOSE(50)
          PRINT*,'TOTCHG:',TOTCHG,FOUND
          call flush(6)
 20    FORMAT(20F9.3)
 556   FORMAT(10F12.6)
       CALL FLUSH(6)
CC      STOP   !!!!!! remove this later
       ERGXL=E1
       ERGXN=E2
       ERGCL=E3
       ERGCN=E4
       RETURN
       END
C
       SUBROUTINE LOWSIC(ISPN,NDH,M,OVER,HAM,FMAT,EVAL,SC1)
C   Mark Pederson 14 September 2001
C TO USE
C PLACE OVERLAP MATRIX IN OVER.
C NEW WAVEFUNCTIONS ARE RETURNED IN HAM...
CKW       USE MPIDAT1,only : IRANK
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION HAM(NDH,NDH),OVER(NDH,NDH),FMAT(NDH,NDH)
     &,EVAL(NDH),SC1(NDH)
!      DIMENSION OVER2(NDH,NDH)
       REAL*8,allocatable :: OVER2(:,:)
       LOGICAL CHECKIT
       LOGICAL LOWEVAL
       LOGICAL AUTOCORR
       DATA CHECKIT/.FALSE./
       DATA LOWEVAL/.FALSE./
       DATA AUTOCORR/.FALSE./
       DATA ISPN_SAV/0/
       SAVE
       IF(CHECKIT)THEN
         PRINT*," NOTE IF CHECKIT=.TRUE. THERE IS N**4 SCALING."
         PRINT*," NORMAL CALLS SHOULD BE WITH CHECKIT=.FALSE."
       END IF
C THa: Check if we want to have very small LOWDEN 
C      eigenvalues autocorrected
C      this may be usefull during startup and bad
C      FO positions 
C      CAUTION: *may* lead to undefined results !
       LOWEVAL  = .FALSE.
       AUTOCORR = .FALSE.
       INQUIRE(FILE='LOWAUTO',EXIST=AUTOCORR)
           
       ALLOCATE(OVER2(NDH,NDH))
       FORALL(I=1:M, J=1:M)
         HAM(J,I) = 0.0D0
         OVER2(J,I)=OVER(J,I)
       END FORALL
       DO I=1,M
         HAM(I,I)=1.0D0
       END DO

       WRITE(6+IRANK,*)'CALLING DIAGGE:',NDH,M
       CALL DIAGGE(NDH,M,OVER,HAM,EVAL,SC1,1)
       FORALL(I=1:M, J=1:M)
         FMAT(J,I)=OVER(J,I)
       END FORALL
           
       WRITE(6+IRANK,*)' LOWDEN OVERLAP EIGENVALUES:'
       WRITE(6+IRANK,20)(EVAL(I),I=1,M)

       DO I=1,M
C        check for too small eigenvalues
C        if there are some, stop after printing them all out
         IF(EVAL(I).LT.1.0D-9) THEN
           LOWEVAL=.TRUE.
C          PRINT *, "> SET <"
         ENDIF
       END DO
       IF(LOWEVAL .eqv. .TRUE.) THEN
         WRITE(6+IRANK,*)
     &    "-----------------------------------------------"
         WRITE(6+IRANK,*)
     &    "> LOWSIC: LOWDEN EIGENVALUES TOO SMALL (< 1e-8)"
         DO I=1,M
C          PRINT *, '>> ', EVAL(I)
           IF(EVAL(I).LT.1.0D-8) THEN
             WRITE(6+IRANK,*) ">  EV: ", I , "=", EVAL(I)
             IF (AUTOCORR .eqv. .TRUE.) THEN
               EVAL(I) = 1.0D-8
               WRITE(6+IRANK,*) ">    corrected to", EVAL(I)
             END IF
           END IF
         END DO
         WRITE(6+IRANK,*)
     &    ">         BAD POSITIONS IN FRMORB ??           "
         WRITE(6+IRANK,*)
     &    "-----------------------------------------------"
         IF(.NOT.AUTOCORR) CALL STOPIT
       END IF
C      normalize
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
C          PRINT 20,(SC1(J),J=1,M)
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
       DO J=1,M
         DO I=1,M
           OVER(J,I)=OVER2(J,I)
         END DO
       END DO
!YY
       DEALLOCATE(OVER2)
C CHECK AGAINST 3-ORDER LOWDEN:
       IF(CHECKIT)THEN
         WRITE(6+IRANK,*)'COMPARISON TO 3RD ORDER LOWDEN...'
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
C          PRINT 20,(HAM(J,I),J=1,M)
C          PRINT 20,(SC1(J  ),J=1,M)
           PRINT*,' '
         END DO
C CHECK FOR ORTHOGONALITY:
C
C        PRINT*,'TRANSFORMED OVERLAP MATRIX:'
         DO I=1,M
           DO J=1,M
             SC1(J)=0.0D0
             DO K=1,M
               DO L=1,M
                 SC1(J)=SC1(J)+HAM(L,I)*HAM(K,J)*OVER(L,K)
               END DO
             END DO
           END DO
C          PRINT 20,(SC1(J),J=1,M)
         END DO
       END IF
 20    FORMAT(' ',5g15.6)
       END
       SUBROUTINE RxRAVEL
     &   (IFNCT,ISHELLA,I_SITE,RVEC,RVECI,N_NUC,ILOC)
C ORIGINALLY WRITTEN BY MARK R PEDERSON (1985)
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       COMMON/FORRx/NFLO,KSPX!!,TMAT(NDH,NDH,2)
       COMMON/FLOINFO/ FOD(3,MAX_OCC,MXSPN),NFOD(MXSPN),MFOD(MXSPN)
       DIMENSION NDEG(3),IND_SALC(ISMAX,MAX_CON,3)
       DIMENSION RVECI(3,MX_GRP),RVEC(3)
       DIMENSION ISHELLV(2)
       DIMENSION FLO(MAX_OCC)
       DATA NDEG/1,3,6/
       DATA ICALL1,ICALL2,ISHELL/0,0,0/
       logical exist
C
c      print *, 'RxRAVEL TMAT'
c      do i = 1,5
c        print 101, (TMAT(j,i,1),j=1,5)
c      end do
c       PSI = 0.0d0
 101    format(5(f10.4))
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
       ITOT=0
       IWF=0
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
         DO 1000 IOCC=1,N_OCC(K_REP,ISPN)
          ITOT=ITOT+1
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
           FACTOR=SQRT(OCCUPANCY(ITOT))
           DO 25 J_LOCAL=1,I_LOCAL
            PSI(J_LOCAL,IWF,ILOC)=FACTOR*PSI(J_LOCAL,IWF,ILOC)
   25      CONTINUE
  950     CONTINUE
 1000    CONTINUE
 1010   CONTINUE
C END OF WAVEFUNCTION LOOP
c           print *, 'IN RXRAVEL, NWFS', NWFS
       IF(NFLO.NE.0)THEN
         IDID=0
               IOFS=(ISPN-1)*NWFS(1)
            !  MFOD=NFOD(ISPN)
          DO JL=1,I_LOCAL
              DO MU=ABS(NFLO),ABS(NFLO)!1,NFOD(ISPN)     
              FLO(MU)=0.0D0
                  DO IW=1,NWFS(ISPN)
                  FLO(MU)=FLO(MU)+TMAT(IW,MU,ISPN)*PSI(JL,IW+IOFS,ILOC)
                  END DO
              END DO
C                 DO IW=1,NWFS(ISPN)
C                 PSI(JL,IW+IOFS  ,ILOC)=0.0D0
C                 END DO
              IF(NFLO.EQ.1)THEN
              IDID=1
                 DO MU=ABS(NFLO),ABS(NFLO)!NFOD(ISPN)     
                 PSI(JL,MU+IOFS  ,ILOC)=FLO(MU)
                 END DO
              ELSE IF(NFLO.LT.0.AND.ISPN.EQ.KSPX)THEN
                 IDID=1       
                 PSI(JL,ABS(NFLO)+IOFS,ILOC)=FLO(ABS(NFLO))
              END IF
          END DO
C         IF(IDID.EQ.0.AND.NFLO.LT.0)STOP'IDID'
       END IF
 1020  CONTINUE
C
       IF (IWF.NE.NWF) THEN
        PRINT *,'UNRAVEL: BIG BUG: NUMBER OF STATES IS INCORRECT'
        PRINT *,'IWF,NWF:',IWF,NWF
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
C
