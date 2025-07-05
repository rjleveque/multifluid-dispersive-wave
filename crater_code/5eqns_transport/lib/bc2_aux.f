c
c ----------------------------------------------------------
c
      subroutine bc2_aux(maxmx,maxmy,ma1,ma2,mbc,mx,my,
     &           xlower,ylower,dx,dy,maux,aux,t,dt,mthbc)
      implicit double precision (a-h,o-z)
      integer mthbc(4)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension rot(4),uv(2)
c
c     # left boundary
      go to (100,110,120,110) mthbc(1)+1
      goto 199
c
 100  continue
      go to 199
c
 110  continue
c     # zero-order extrapolation
      do ma=ma1,ma2
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               aux(1-ibc,j,ma) = aux(1,j,ma)
               enddo
            enddo
         enddo
      go to 199

 120  continue
c     # periodic
      do ma=ma1,ma2
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               aux(1-ibc,j,ma) = aux(mx+1-ibc,j,ma)
               enddo
            enddo
         enddo
      go to 199
c
 199  continue
c     # right boundary
      go to (200,210,220,210) mthbc(2)+1
      goto 299
c
 200  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

 210  continue
c     # zero-order extrapolation
      do ma=ma1,ma2
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               aux(mx+ibc,j,ma) = aux(mx,j,ma)
               enddo
            enddo
         enddo
      go to 299

 220  continue
c     # periodic
      do ma=ma1,ma2
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               aux(mx+ibc,j,ma) = aux(ibc,j,ma)
               enddo
            enddo
         enddo
      go to 299
c
 299  continue
c
c     # bottom boundary
      go to (300,310,320,310) mthbc(3)+1
      goto 399
c
 300  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
 310  continue
c     # zero-order extrapolation
      do ma=ma1,ma2
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               aux(i,1-jbc,ma) = aux(i,1,ma)
               enddo
            enddo
         enddo
      go to 399
c
 320  continue
c     # periodic
      do ma=ma1,ma2
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               aux(i,1-jbc,ma) = aux(i,my+1-jbc,ma)
               enddo
            enddo
         enddo
      go to 399
c
 399  continue
c     # top boundary
      go to (400,410,420,410) mthbc(4)+1
      goto 499
c
 400  continue
      go to 499

 410  continue
c     # zero-order extrapolation
      do ma=ma1,ma2
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               aux(i,my+jbc,ma) = aux(i,my,ma)
               enddo
            enddo
         enddo
      go to 499
c
 420  continue
c     # periodic
      do ma=ma1,ma2
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               aux(i,my+jbc,ma) = aux(i,jbc,ma)
               enddo
            enddo
         enddo
      go to 499
c
 499  continue
      return
c
      end
