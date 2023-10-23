           subroutine ezstart
           implicit real*8 (a-h,o-z)
           character*30,name,rname
           dimension r(5,3000) 
           dimension g(3,3,121),ind(121)
           data au2ang/0.529177d0/
            open(80,file='CLUSTER')
            rewind(80)
                       do i=1,30
                       name(i:i)=' '
                       end do
            read(80,10)name
 10         format(a30) 
            if(name(1:1).ne.'@')then
               print*,'improper call to eztart'
               call stopit
            endif
                 if(name(30:30).ne.' ')then
                     print*,'xmol file must be'
                     print*,'less than 30 characters'
                  end if
                     do i=2,30
                     name(i-1:i-1)=name(i:i)
                     end do
            close(80)
            do i=1,48
            do j=1,3
            do k=1,3
            g(k,j,i)=0.0d0
            end do
            end do
            end do
             g(1,1,1)= 1.0d0
             g(2,2,1)= 1.0d0
             g(3,3,1)= 1.0d0
             g(1,1,4)=-1.0d0
             g(2,2,4)= 1.0d0
             g(3,3,4)= 1.0d0
             g(1,2,2)= 1.0d0
             g(2,3,2)= 1.0d0
             g(3,1,2)= 1.0d0
             g(1,2,3)= 1.0d0
             g(2,1,3)= 1.0d0
             g(3,3,3)= 1.0d0
             mgp=4
                call clsgrp(mgp,g)
                ngp=mgp
                   print*,'ngp:',ngp
                open(90,file=name)
                    read(90,*)nat
                    read(90,*)
                           if(nat.gt.300)then
                                 print*,'nat too big in ezstart'
                            end if
                    do iat=1,nat
                    read(90,*)nuc,(r(j,iat),j=1,3)
                    print 80,(r(j,iat),j=1,3),nuc
                          r(5,iat)=dfloat(nuc)
                    end do
                    call prnaxes(nat,r)
                    ngp=0
                    do kgp=1,mgp
                       call chkrot(nat,r,g(1,1,kgp),noyes)
                       if(noyes.eq.1)then
                       ngp=ngp+1
                       ind(ngp)=kgp
                       end if
                    end do
                    do igp=1,ngp
                               do i=1,3
                               do j=1,3
                               g(i,j,igp)=g(i,j,ind(igp))
                               end do
                               end do
                    end do
                    print*,'Oh Elements:',ngp
                                          noh=ngp
c check for additional odd rotations about principle axes...
                    do is=1,3
                       if(is.eq.1)then
                       s1=1.0d0
                       s2=1.0d0
                       rname='-FOLD ROTATION'
                       else if (is.eq.2)then
c probably redundant...
                       s1=-1.0d0
                       s2=-1.0d0
                       rname='-FOLD ROTATION'
                       else if (is.eq.2)then
                       rname='-FOLD PSEUDO ROTATION'
                       s1=-1.0d0
                       s2= 1.0d0
                       end if
                    do ix=1,3
                       iy=ix+1
                       iz=iy+1
                       if(iy.gt.3)iy=iy-3
                       if(iz.gt.3)iz=iz-3
                        do n=2,7   
                        ngp=ngp+1
                        ang=s2*8.0d0*atan(1.0d0)/n
                              do j=1,3
                              do k=1,3
                              g(j,k,ngp)=0.0d0
                              end do
                              end do
                              g(ix,ix,ngp)=    1.0d0*s1
                              g(iy,iy,ngp)= cos(ang)*s1
                              g(iy,iz,ngp)= sin(ang)*s1
                              g(iz,iy,ngp)=-sin(ang)*s1
                              g(iz,iz,ngp)= cos(ang)*s1
c check to see if this is a new group operation...
                              new=1
                              do igp=1,ngp-1
                                diff=0.0d0
                                do i=1,3
                                do j=1,3
                                  diff=diff+abs(g(j,i,igp)-g(j,i,ngp))
                                end do
                                end do
                                if(diff.lt.0.1)new=0
                              end do
                       if(new.eq.1)then
                       call chkrot(nat,r,g(1,1,ngp),noyes)
                       else
                              noyes=0
                       end if
                          if(noyes.eq.0)then
                          ngp=ngp-1
                          else 
                          print*,'FOUND',n,rname,ix,iy,iz
                          call clsgrp(ngp,g)
                          end if
                        end do
                    end do
                    end do
                    call rotall(noh,ngp,nat,r,g)
                    open(80,file='GRPMAT')
                    write(80,*)ngp,' Number of group elements'
                    do igp=1,ngp
                    write(80,75)((g(i,j,igp),j=1,3),i=1,3)
                    write(80,*)' '
                    end do
                    close(80)
 75                 format(' ',3g21.14)
                mat=nat+1
                do iat=1,nat
                   r(4,iat)=1.0d0
                   errmin=1.0d0
                   do igp=1,ngp
                     do j=1,3
                     r(j,mat)=0.0d0
                       do k=1,3
                       r(j,mat)=r(j,mat)+g(j,k,igp)*r(k,iat)
                       end do
                     end do
                     do jat=1,iat-1   
                     dd=distance(r(1,mat),r(1,jat))
     &                 +abs(r(5,jat)-r(5,iat))
                     errmin=min(dd,errmin)
                     end do
                   end do
                     if(errmin.lt.0.001)r(4,iat)=-1.0d0
                end do
                mat=0
                do iat=1,nat
                   if(r(4,iat).gt.0.0d0)mat=mat+1
                end do
                open(80,file='CLUSTER')
                write(80,78)
                write(80,79)
                write(80,*)mat,' Number of inequivalent atoms'
                do iat=1,nat
                   if(r(4,iat).gt.0.0d0)then
                   write(80,80)(r(j,iat)/au2ang,j=1,3),nint(r(5,iat))
                   end if
                end do
                write(80,81)chrg,spin
                close(80)
                call reorient_molecule
                call perfect_cluster
 30             format(a20)
 78             format('GGA-PBE*GGA-PBE')
 79             format('GRP')
 80             format(' ',3f20.14,I5,' ALL UPO ')
 81             format(' ',2f12.4,' Net Charge and Moment')
                end
                function distance(x,y)
                implicit real*8 (a-h,o-z)
                dimension x(*),y(*)
                distance=(x(1)-y(1))**2+(x(2)-y(2))**2+(x(3)-y(3))**2
                distance=sqrt(distance)
                return
                end
                subroutine prnaxes(nat,r)
                include 'PARAMS'            
                include 'commons.inc'
                dimension r(5,*),p(3)
c center molecule:
                          do i=1,3
                          ham(i,1)=0.0d0
                              do j=1,nat
                              ham(i,1)=ham(i,1)+r(i,j)
                              end do
                          ham(i,1)=ham(i,1)/nat
                          end do
                print*,'center of charge:',(ham(j,1),j=1,3)
                              do i=1,3
                              do j=1,nat
                              r(i,j)=r(i,j)-ham(i,1)
                              end do
                              end do
                do i=1,3
                do j=1,3
                    ham (i,j)=0.0d0
                    over(i,j)=0.0d0
                    do k=1,nat
                    ham(i,j)=ham(i,j)-r(i,k)*r(j,k)*r(5,k)
                    end do
                end do
                    over(i,i)=1.0d0
                end do
                   do i=1,3
                   print*,(ham(i,j),j=1,3)
                   end do
                call diagge(ndh,3,ham,over,eval,sc1,1)
                     do i=1,3  
                     print*,i,eval(i),(ham(j,i),j=1,3)
                     end do
                     do i=1,nat
                     print 10,(r(j,i),j=1,5)
                       do j=1,3
                       p(j)=0.0d0
                         do k=1,3
                         p(j)=p(j)+ham(k,j)*r(k,i)
                         end do
                       end do
                       do j=1,3
                       r(j,i)=p(j)
                       end do
                     print 10,(r(j,i),j=1,5)
                     end do
 10             format(' ',5f12.4)
                end
                subroutine rotall(noh,ngp,nat,r,g)
                implicit real*8 (a-h,o-z)
                dimension r(5,*),g(3,3,121),rot(3,3),grt(3,3)
                dimension bet(3,3),p(3)
                dimension ind(3)
c this subroutine tries to rotate molecule to find the maximum
c number of Oh matrices => better effiency for mesh construction 
c and vib modes. 
c
c primarily for problems with 3-fold and 5-fold rotations...
c 5-fold rotation fixer still needs to be completed.
               nrot=2
               nbst=noh
               do irot=2,2   
               do ip=1,6
                   if(ip.eq.1)then
                      ind(1)=1
                      ind(2)=2
                      ind(3)=3
                   else if(ip.eq.2)then
                      ind(1)=3
                      ind(2)=1
                      ind(3)=2
                   else if(ip.eq.3)then
                      ind(1)=2
                      ind(2)=3
                      ind(3)=1
                   else if(ip.eq.4)then
                      ind(1)=2
                      ind(2)=1
                      ind(3)=3
                   else if(ip.eq.5)then
                      ind(1)=3
                      ind(2)=2
                      ind(3)=1
                   else if(ip.eq.6)then
                      ind(1)=1
                      ind(2)=3
                      ind(3)=2
                   end if
                   if(irot.eq.1)then
                       do i=1,3
                       do j=1,3
                       rot(ind(j),ind(i))=0.0d0 
                       end do
                       rot(ind(i),ind(i))=1.0d0 
                       end do
                   else if(irot.eq.2)then
                    rot(1,ind(1))= 1.0d0/sqrt(3.0d0)
                    rot(2,ind(1))= 1.0d0/sqrt(3.0d0)
                    rot(3,ind(1))= 1.0d0/sqrt(3.0d0)
                    rot(1,ind(2))= 1.0d0/sqrt(2.0d0)
                    rot(2,ind(2))=-1.0d0/sqrt(2.0d0)
                    rot(3,ind(2))= 0.0d0           
                    rot(1,ind(3))= 1.0d0/sqrt(6.0d0)
                    rot(2,ind(3))= 1.0d0/sqrt(6.0d0)
                    rot(3,ind(3))=-2.0d0/sqrt(6.0d0)
 10                 format(' ',3f12.4)
                   end if
                       ngd=0
                   do igp=1,ngp
                       add=0.0d0
                       do i=1,3
                       do j=1,3
                        grt(i,j)=0.0d0
                        do k=1,3
                        do l=1,3
                        grt(i,j)=grt(i,j)+g(k,l,igp)*rot(i,k)*rot(j,l)
                        end do
                        end do
                        add=add+abs(grt(i,j))
                       end do
                       end do
                       if(abs(add-3.0d0).le.1.0D-6)ngd=ngd+1
                   end do
                      if(ngd.gt.nbst)then
                      nbst=ngd
                      print*,irot,ngd,nbst,' irot, ngd...'
                      do i=1,3
                      do j=1,3
                      bet(j,i)=rot(j,i)
                      end do
                      end do
                    do i=1,3
                    print 10,(rot(i,j),j=1,3)
                    end do
                      end if
               end do
               end do
                  if(noh.ge.ngd)return
                      do i=1,3
                      do j=1,3
                      rot(j,i)=bet(j,i)
                      end do
                      end do
                   do iat=1,nat
                    do j=1,3
                    p(j)=0.0d0
                    do k=1,3
                    p(j)=p(j)+rot(j,k)*r(k,iat)
                    end do
                    end do
                    do j=1,3
                    r(j,iat)=p(j)
                    end do
                   end do
                   do igp=1,ngp
                    add=0.0d0
                       do i=1,3
                       do j=1,3
                        grt(i,j)=0.0d0
                        do k=1,3
                        do l=1,3
                        grt(i,j)=grt(i,j)+g(k,l,igp)*rot(i,k)*rot(j,l)
                        end do
                        end do
                        add=add+abs(grt(i,j))
                       end do
                       end do
                    print*,igp,add
                        do k=1,3
                        do l=1,3
                        g(k,l,igp)=grt(k,l)
                        end do
                        end do
                   end do
                   do igp=1,ngp
                    do jgp=igp+1,ngp
                      add1=0.0d0
                      add2=0.0d0
                         do i=1,3
                         do j=1,3
                         add1=add1+abs(g(j,i,igp))
                         add2=add2+abs(g(j,i,jgp))
                         end do
                         end do
                         if(add2.lt.add1)then
                              do i=1,3
                              do j=1,3
                              gs=g(j,i,igp)
                              g(j,i,igp)=g(j,i,jgp)
                              g(j,i,jgp)=gs          
                              end do
                              end do
                         end if
                    end do
                   end do
c some asthetics...
                  do igp=1,ngp
                   do i=1,3
                   do j=1,3
                       if(abs(g(j,i,igp)).le.1.0d-14)then
                              g(j,i,igp)=0.0d0
                       end if
                   end do
                   end do
                  end do
               return
               end  
                subroutine chkrot(nat,r,g,noyes)
                implicit real*8 (a-h,o-z)
                dimension r(5,*),g(3,3)
                       kat=nat
                       do iat=1,nat
                       kat=kat+1
                        do j=1,3
                         r(j,kat)=0.0d0
                         do k=1,3
                         r(j,kat)=r(j,kat)+g(j,k)*r(k,iat)
                         end do
                        end do
                        errmin=1.0d0
                        do jat=1,nat
                        dd=distance(r(1,kat),r(1,jat))
     &                    +abs(r(5,iat)-r(5,jat))
                        errmin=min(dd,errmin)
                        end do
                        if(errmin.lt.0.001)kat=kat-1
                       end do
                       if(kat.eq.nat)then
                         noyes=1
                       else
                         noyes=0
                       end if
                return
                end
C Sample Input File:
C
C   4
C              
C 25     1.71819    1.71819    1.71819
C 25    -1.71819   -1.71819    1.71819
C 25     1.71819   -1.71819   -1.71819
C 25    -1.71819    1.71819   -1.71819

C FOR_DIAG: used during matrix diagonalization
       SUBROUTINE REORIENT_MOLECULE
       IMPLICIT REAL*8 (A-H,O-Z)
       LOGICAL DOIT 
       CHARACTER*80 line
       DIMENSION O(3,3), U(3,3), V(3,3), G(3,3), RI(3), RO(3)
              print*,'in REORIENT_MOLECULE'
C written by Karma Dema- April 2020
C this subroutine looks at ...      
       U(1,1)=-1.0d0/SQRT(2.0D0)
       U(1,2)= 1.0d0/SQRT(2.0D0)
       U(1,3)= 0.d0
       U(2,1)=-1.0d0/SQRT(6.0D0)
       U(2,2)=-1.0d0/SQRT(6.0D0)
       U(2,3)= 2.0d0/SQRT(6.0D0)
       U(3,1)= 1.0d0/SQRT(3.0D0)
       U(3,2)= 1.0d0/SQRT(3.0D0)
       U(3,3)= 1.0d0/SQRT(3.0D0)

       V = TRANSPOSE(U)

      O(1,1)= -0.500000000000000
      O(1,2)=0.86602540378444
      O(1,3)=0.0000000000000
      O(2,1)=-0.86602540378444
      O(2,2)=-0.50000000000000
      O(2,3)=0.0000000000000
      O(3,1)=0.0000000000000
      O(3,2)=0.0000000000000
      O(3,3)=1.0000000000000


        open(2, file='GRPMAT',status='unknown')
        doit=.false.
        read(2,*)ngp
        do igp=1,ngp
        read(2,*)O
        do i=1,3
        do j=1,3
        if(abs(o(j,i)-0.5d0).lt.1.0d-6)doit=.true.
        if(abs(o(j,i)+0.5d0).lt.1.0d-6)doit=.true.
        end do
        end do
        end do 
                  print*,'REORIENT:',doit        
        if (doit)then
        rewind(2)
        else
        close(2)
        return
        end if 

        open(1, file='data1.dat',status='unknown')

        do i=1,3
                !print*, (      U(i,j), j=1,3)
        end do


        do i=1,3
                !print*, (V(i,j), j=1,3)
        end do


        do i=1,3
                !print*,  (      O(i,j), j=1,3)
        end do
        read(2,*)ngp
        write(1,*)ngp
        do 10 it=1,2
        rewind(2)
        read(2,*)ngp         
        do igp=1,ngp
        read(2,*)O



        G= MATMUL(O,U)
        G= MATMUL(V,G)
        
        nt=1        
        do i=1,3
        do j=1,3
        if(abs(G(i,j)).le.0.99999.and.abs(G(i,j)).gt.0.000001)nt=2 
        end do
        end do
        
        if(nt.eq.it)then
        if(nt.eq.1)then
        do i=1,nt
        do j=1,nt
        g(i,j)=dfloat(nint(g(i,j)))
        end do
        end do
        end if

        do i=1,3
        write(1,100)  (G(i,j), j=1,3)
 100    format(3F20.13)
        end do
        write(1,*)' '
        end if
                end do
 10     continue
        close(1)
        close(2)

        open(1,file='CLUSTER')
        open(2,file='CLUSALT')

        read(1,200)line
        write(2,200)line
        read(1,200)line
        write(2,200)line
 200    format(A80)       
        read(1,*)nat
        write(2,*)nat

        do iat=1,nat
        read(1,200)line
        do j=1,77
        if(line(j:j+2).eq.'ALL')l=j-3
        end do
        print*,line(l:80)
        read(line,*)RI
                
        do j=1,3
        RO(j)=0.0d0
        do k=1,3
        RO(j)=RO(J)+V(J,K)*RI(K)
        end do
        end do
                print*, (RO(J),J=1,3)
                write(2,400) (RO(J),J=1,3),line(l:80)
 400    format(3f20.12,1x,A)        
        end do
        read(1,200)line
        write(2,200)line
        write(2,*)'Molecule reoriented by Karma Dema'
        print   *,'Molecule reoriented by Karma Dema'
        close(1)
        close(2)
       
       CALL SYSTEM ('cp GRPMAT GRPMAT_EZSTART')
       CALL SYSTEM ('mv data1.dat GRPMAT') 
       CALL SYSTEM ('cp CLUSTER CLUSTER_EZSTART') 
       CALL SYSTEM ('mv CLUSALT CLUSTER') 
       CALL SYSTEM ('cat GRPMAT >>CLUSTER') 
       
       END
