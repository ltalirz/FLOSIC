C
C define maximum number of symbols
C
       PARAMETER (MXSYMBS=MAX_IDENT+20)
C
C DEBUG:  determines if debug info will be printed
C MPIDAT: used by the MPI version to keep track of CPU usage
C
       LOGICAL DEBUG
       COMMON/DEBUG/DEBUG
       COMMON/MPIDAT/INUSE(MX_PROC),NPROC,NCALLED,IRANK
C
C ALPCOR: largest exponent of NLCC core density
C PSPANG: lmax of different radial zones used for nonlocal PSPs
C PSPINF: general PSP info 
C PSPNLO: tabulated nonlocal PSPs
C BHSPSP: data for local potential of type BHS
C TABPSP: data for local potential of type TAB (including NLCC)
C
       CHARACTER*7 PSPSYM(MAX_FUSET)
       COMMON/ALPCOR/ALPCOR(MAX_FUSET)
       COMMON/PSPANG/PSRZONE(4),LMXPSRZ(5)
       COMMON/PSPINF/ISITPSP,ISNLCC,PSPSYM
       COMMON/PSPNLO/RPSNLO(MXRPSP,MAX_FUSET),WPSNLO(MXRPSP,MAX_FUSET)
     &  ,VPSNLO(MXLPSP+1,MXRPSP,MAX_FUSET)
     &  ,LMAXNLO(MAX_FUSET),NRPSP(MAX_FUSET)
       COMMON/BHSPSP/BHSALP(2,MAX_FUSET),BHSCOF(2,MAX_FUSET)
       COMMON/TABPSP/RRADTAB(MXPTAB,MAX_FUSET)
     &  ,VLRTAB(2,MXPTAB,MAX_FUSET),RHOCOR(3,MXPTAB,MAX_FUSET)
     &  ,NRADTAB(MAX_FUSET),NLCC(MAX_FUSET)
C
C NUCLEI:   information about atom locations and types
C BASET:    basis set / nuclear charge data for atom types
C FERMIONS: number of electrons, wavefunction output file
C DFTYP:    information about density functional
C SPIN:     spin info
C DIPOLE:   dipole moment and external homegeneous electric field
C ENERG:    stores a lot of different energies
C FORCES:   force information
C NUMFRC:   used for calculation Pulay corrections to the forces
C OPTGEOM:  data needed in geometry optimizations
C SYMBNAME: contains the symbols used for the identity members
C
       CHARACTER*40 WFFILE
       CHARACTER*80 SYMBS(2,MXSYMBS)
       CHARACTER SYMATM(10,MXSYMBS)
       COMMON/NUCLEI/RIDT(3,MAX_IDENT),RCNT(3,MX_CNT)
     &  ,IFUIDT(MAX_IDENT),IFUCNT(MX_CNT),NIDENT,NCNT
       COMMON/BASET/ZELC(MAX_FUSET),ZNUC(MAX_FUSET)
     &  ,BFCON(MAX_BARE,MAX_CON,LDIM,MAX_FUSET)
     &  ,BFALP(MAX_BARE,MAX_FUSET)
     &  ,N_BARE(MAX_FUSET),N_CON(LDIM,MAX_FUSET)
     &  ,LSYMMAX(MAX_FUSET),N_POS(MAX_FUSET),NFNCT
       COMMON/FERMIONS/E_UP,E_DN,WFFILE
       COMMON/DFTYP/IGGA(2),IDFTYP(2)
       COMMON/SPIN/ISPN,NSPN
       COMMON/DIPOLE/DIPOLE(3),EFIELD(3)
       COMMON/ENERG/ETOTAL,ENNUC,ELOCAL,ECOUL,ERGFLD,EKINONL,ENONLO
     &  ,ERGXL,ERGXN,ERGCL,ERGCN
       COMMON/FORCES/FHELLF(3,MAX_IDENT),FNONL(3,MAX_IDENT)
     &  ,FTOT(3,MAX_IDENT)
       COMMON/NUMFRC/FRC1(3,MAX_IDENT),FRC2(3,MAX_IDENT)
       COMMON/OPTGEOM/ATMSYM(MAX_IDENT),RIDSYM(3,MAX_IDENT),GTOL
     &  ,MOVING(3,MAX_IDENT),NSCHAR(2,MXSYMBS),NIDSYM,NSYMBS,SYMBS
       COMMON/SYMBNAME/SYMATM
C
C GROUP: point group operations
C
       COMMON/GROUP/RMAT(3,3,MX_GRP),NGRP,MULTAB(MX_GRP,MX_GRP)
C
C STPOT:   data for spin dependence of start potential
C RHOPFIT: data for fitted atomic potentials and density
C
       COMMON/STPOT/ISPIDT(MAX_IDENT)
       COMMON/RHOPFIT/RPFALP(MAX_FUSET),RPFCMX(MAX_FUSET)
     &  ,RPFCOF(2,MAXLSQF,MAX_FUSET),NRPFIT(MAX_FUSET),LDIVR(MAX_FUSET)
C
C SCFDAT:   SCF and startup data
C CBLK4:    contains wavefunctions and their occupation numbers
C UNSYM:    contains unsymmetrized wavefunctions
C CBLK7:    used to store the hamiltonian and overlap matrices
C CBLK30:   hamiltonian submatrix and temporary array needed during H setup
C EFRMI:    Fermi energy for each spin system
C ELEVELS:  contains occupied eigenvalues
C FOR_DIAG: used during matrix diagonalization
C
       LOGICAL CONVERGENCE,HAVEHAM
       COMMON/SCFDAT/ISTSCF,IHIPOL,ITSCF,MAXSCF,CONVERGENCE,HAVEHAM
       COMMON/CBLK4/PSI_COEF(NDH,MAX_VIRT_PER_SYM,MAX_REP,MXSPN)
     &  ,OCCUPANCY(MAX_VIRT_PER_SYM*MAX_REP*MXSPN),N_OCC(MAX_REP,MXSPN)
       COMMON/UNSYM/PSI(MAXUNSYM,MAX_OCC,2),NWF,NWFS(MXSPN)
C      COMMON/TMP1/HSTOR(NDH_TOT,2)
        COMMON/CBLK7/HSTOR(NDH_TOT,2)
       COMMON/CBLK30/HOLD(MAXUNSYM,MAXUNSYM),HTEMP(MTEMP)
       COMMON/EFRMI/EFERMI(2)
       COMMON/ELEVELS/EVLOCC(MAX_OCC)
       COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),EVAL(NDH)
     &  ,SC1(NDH),SC2(NDH)
C
C MESH:   mesh points and weights
C SPH:    data for spheres that contains subsets of mesh points
C
       COMMON/MESH/WMSH(MAX_PTS),RMSH(3,MAX_PTS),NMSH
       COMMON/MSPH/TSPH(4,MX_SPH),LIMSPH(2,MX_SPH),NSPHERES
C
C FORDEN:   determines if the density is created in COUPOT or DENSOLD
C CUTOFF:   defines cutoff for exponents in COUPOT
C WASTED:   time wasted in UNRAVEL
C
       COMMON/FORDEN/MODDEN
       COMMON/CUTOFF/GAUSS_CUT(MAX_IDENT)
       COMMON/WASTED/T1UNRV,T2UNRV
C
C symmetry stuff 
C
       COMMON/BASE_REP/S_REP(MX_GRP),P_REP(3,3,MX_GRP)
     &  ,D_REP(6,6,MX_GRP)
       COMMON/CBLK11/REP(5,5,MX_GRP,MAX_REP),N_REP,NDMREP(MAX_REP)
       COMMON/CBLK12/U_MAT(LOCMAX,ISMAX,MAXSYMSALC,3,2)
     &  ,N_SALC(MAXSYMSALC,3,MAX_IDENT)
       COMMON/CBLK18/IGEN(MX_GRP)
       COMMON/INFOR/RDENT(3,MAX_IDENT),NUMSITES(MAX_IDENT)
     &  ,IGGEN(MX_GRP,MAX_IDENT),N_IDNT
       COMMON/INTS_INDEX/INDBEG(MAX_IDENT,MAX_REP),NS_TOT(MAX_REP)
       COMMON/REDREP/LDMREP(MAX_REP)

C
C SIC STUFF (TMAT FERMI-LOWDIN TRNASFORMATION)
       COMMON/FTMAT/TMAT(NDH,NDH,2)
       SAVE
