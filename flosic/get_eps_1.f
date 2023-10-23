      SUBROUTINE GET_EPS(KSPNX)
C
C
      INCLUDE 'PARAMS'
      INCLUDE 'commons.inc'
      character*3 spin
      character*7 line
c     COMMON /FTMAT/TMAT(NDH,NDH,2)  !in commons.inc
      COMMON/FLOINFO/FOD(3,MAX_OCC,2),NFOD(2),MFOD(2)

      common /sicmat/SIC(max_occ,max_occ,2)
      dimension ham1(200,200,10),ham3(200,200),ham2(200,200)
      dimension knt(200,4),npsi(10),nstart(10)
c     dimension ham1(100,100),ham2(100,100),ham3(100,100)
c      DATA JBEG,JEND/1,2,5,1,4,10/
c      DATA ISIZE/1,3,6/
C
c  get matrix elements of SIC matrices -- needed for FOD forces
C
c
c  convert eigenvectors to FLO representation
c

c
c  FIX THIS -- GET INFO RE # OF FODs INTO the subroutine
c
c       nfod = 8
        nbasis = 0 
 
        do irep = 1, n_rep
        call flush(6)
        do irow = 1, ndmrep(irep)
        do ibas = 1, ns_tot(irep)
          nbasis = nbasis+1
        end do
        end do  !krow
        end do  !irep
c       print *, 'nbasis', nbasis
c       print *, 'tmat'
c       do i = 1, nfod(kspnx)
c          test = 0.00
c          do jbas = 1,nbasis   !nfod(kspnx)
c            test = test + tmat(jbas,i,kspnx)**2
c          end do
c          print 100, (tmat(jbas,i,kspnx),jbas=1,nfod(kspnx)+10)
c          print *, 'test', test
c       end do        
 100  format(8(f10.4))
        call flush(6)
c
c  print HAM
c
      indx = 0
      sic = 0.0d0
      ham1 = 0.0d0
      do i = 1, n_rep
c
c  ham3 stores hamiltonian in fixed basis rep
c
        ham3 = 0.0d0
         print *, 'ham in fixed basis'
         do j = 1,ns_tot(i)
           do k = j,ns_tot(i)
             indx = indx+1
             if(j.eq.k) then
                    diff = hstor(indx,2) - hstor(indx,1)
                    if(diff.gt.0.000001) then
                             print *, 'get_eps: hstor diag not same'
                             print *, hstor(indx,2), hstor(indx,1)
                    end if
             end if
             ham3(k,j) = hstor(indx,2)
             ham3(j,k) = hstor(indx,1)
           end do
          end do
          do j = 1,8 
           print 100, (ham3(k,j), k=1,8)
          end do
c
c  convert HAM to psi basis
c      HAM2 holds intermediate array -- first half of the transformation
c
         ham2 = 0.0d0
        do l = 1, n_occ(i,kspnx)
          do n = 1, ns_tot(i)
          do k = 1, ns_tot(i)
             ham2(l,n) = ham2(l,n) + 
     &                         ham3(k,n)*psi_coef(k,l,i,kspnx)
           end do   !k
           end do   !n
c        print *, 'second half'
         do j = 1, n_occ(i,kspnx)
           do k = 1, ns_tot(i)
              ham1(l,j,i) = ham1(l,j,i) + 
     &                               ham2(l,k)*psi_coef(k,j,i,kspnx)
           end do
         end do
        end do
      call flush(6)
      print *, 'ham in psi rep'
      do l = 1, n_occ(i,kspnx)
         print 100, (ham1(j,l,i),j = 1, n_occ(i,kspnx))
      end do
      call flush(6)

      end do       !irep
c     stop
c
c  READ IN INFO ABOUT PSI's FROM EVALUES
c
      open(unit=15,file='EVALUES')
      do j = 1,1000
      read(15,*,end=300) line
c       print *, line
        if(line.eq.'PROJECT') go to 201
      end do
 300  continue
      print *, 'cannot find evalues list'
      stop
 201  continue
c
c     create array to find KSO indices for TMAT
c
      ntot = 0
      knt = 0
      do irep = 1,n_rep
        nstart(irep) = ntot
        ntot = ntot + n_occ(irep,kspnx)*ndmrep(irep)
      end do
      print *, 'nstart', (nstart(i),i=1,n_rep)
      npsi = 0
      ipsi = 0
      do ii = 1,100
      read(15,10,end=101) i ,irep, ideg, spin
      if(ii.ne.i) then
          print *, 'ii ne i', ii, i
          go to 101
      end if
      ispn = 0
      if(spin.eq.' u') ispn = 1
      if(spin.eq.' d') ispn = 2
          if(ispn.eq.0) then
              print *, 'spin stop', spin
              stop
          end if
       if(ispn.eq.kspnx) then    
         npsi(irep) = npsi(irep) + 1
         do jrow = 1,ideg
          ipsi = ipsi + 1
          knt(ipsi,1) = irep
          knt(ipsi,2) = jrow
          knt(ipsi,3) = npsi(irep)
          knt(ipsi,4) = nstart(irep) + (npsi(irep)-1)*ideg + jrow
         end do
       end if
c   1 Rep:  1 Deg:  1 Spin: u Energy     -11.152799 Occ:       1.000000 <==                         NO-FLIP        1  1  1
 10    format(i5,5x,i3,5x,i3,6x,a3) 
c     print *, 'i, irep, ideg, spin', i, irep, ideg, spin
      end  do
 101  continue
      close(15)
      print *, 'finished reading evalues list'
c     print *, 'knt array, ipsi', ipsi
c     do iwf = 1, ipsi
c        print *, (knt(iwf,j),j=1,4)
c     end do
c
c  convert HAM to phi basis; use HAM2 for intermediate step
c  knt array links psi index to irrep, row, etc.
c
c
       SIC = 0.0d0
       HAM2 = 0.0d0
       do iflo = 1,NFOD(kspnx)
       do jflo = 1,NFOD(kspnx)
           do l1 = 1,nfod(kspnx)     !psi on left
             lrep = knt(l1,1)
             lrow = knt(l1,2)
             lwf  = knt(l1,3)
             lpsi = knt(l1,4)
           do k1 = 1,nfod(kspnx)     !psi on right
             krep = knt(k1,1)
             krow = knt(k1,2)
             kwf  = knt(k1,3)
             kpsi = knt(k1,4)
             if(lrep.eq.krep.and.lrow.eq.krow) then
c                    ham2(l1,iflo) = ham2(l1,iflo) +
c    &                  ham1(lwf,kwf,lrep)*tmat(k1,iflo,kspnx)
             SIC(iflo,jflo,kspnx) = SIC(iflo,jflo,kspnx) + 
     & tmat(kpsi,iflo,kspnx)*ham1(kwf,lwf,lrep)*tmat(lpsi,jflo,kspnx)
             end if
           end do
           end do
        end do
c        do iflo = 1,nfod(kspnx)
c          do jflo = 1,NFOD(kspnx)
c          do k1 = 1,nfod(kspnx)
             
c            SIC(iflo,jflo,kspnx) = sic(iflo,jflo,kspnx) - 
c    &                         tmat(k1,jflo,kspnx)*ham2(k1,iflo)
c          end do
c          end do
       end do
         print *, 'SIC  nfod', nfod
         do if = 1,nfod(kspnx)
          print 100, (SIC(jf,if,kspnx),jf=1,nfod(kspnx))
         end do

      RETURN
      END
