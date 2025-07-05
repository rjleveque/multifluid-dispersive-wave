c
c---------------------------------------------------------------
c     
      subroutine out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &           dx,dy,q,t,iframe,aux,maux)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension qloc(20),qp(20)
      character*10 fname1,fname2,fname3
      logical outaux
c
c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>

      outaux = .false.
c
c     # first create the file name and open file
c
      fname1 = 'fort.qxxxx'
      fname2 = 'fort.txxxx'
      fname3 = 'fort.axxxx'
      nstp = iframe
      do ipos=10,7,-1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         fname2(ipos:ipos) = char(ichar('0') + idigit)
         fname3(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp/10
         enddo    
c
      open(unit=50,file=fname1,status='unknown',form='formatted')
      open(unit=60,file=fname2,status='unknown',form='formatted')
c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:
      ngrids = 1
      mptr   = 1
      level  = 1
c
      write(50,1001) mptr,level,mx,my
      write(50,1002) xlower,ylower,dx,dy
c
      do j=1,my
         do i=1,mx
            do m=1,meqn
               qloc(m) = q(i,j,m)
               enddo
c
            call ctomfd(qloc,qp(1),qp(2),qp(3),qp(4),rhoh0,
     &           qp(5),grue0,zeta0,qp(6),zfb0,c20)
c
            do m=1,6   
               if (dabs(qp(m)) .lt. 1.d-40) qp(m) = 0.d0
               enddo
c
            write(50,1005) (qp(m),m=1,6)
            enddo   
         enddo  
c
      if (outaux) then 
c         # also output the aux arrays:
          open(unit=70,file=fname3,status='unknown',form='formatted')
          write(70,1001) mptr,level,mx,my
          write(70,1002) xlower,ylower,dx,dy
          do j=1,my
             do i=1,mx
                do m=1,maux
c                  # exponents with more than 2 digits cause problems reading
c                  # into matlab... reset tiny values to zero:
                   if (dabs(aux(i,j,m)) .lt. 1d-99) aux(i,j,m) = 0.d0
                   enddo
c
                write(70,1006) (aux(i,j,m),m=1,maux)
                enddo    
             enddo     
          close(unit=70)
      endif
c
      meqout = 6
      write(60,1000) t,meqout,ngrids,maux
c
      close(unit=50)
      close(unit=60)
c
      call out1(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,t,iframe,aux,maux)
c
      call outgrd(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &     dx,dy,iframe,aux,maux)
      return
c
 1000 format(f26.16,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 maux'/,/)
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')
 1002 format(f26.16,'    xlow', /,
     &       f26.16,'    ylow', /,
     &       f26.16,'    dx', /,
     &       f26.16,'    dy',/)
 1005 format(6e18.10)
 1006 format(10f26.16)
      end
