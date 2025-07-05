c
c ----------------------------------------------------------
c
      subroutine copyq2(maxmx,maxmy,meqn,mbc,mx,my,q1,q2)
      implicit double precision (a-h,o-z)
      dimension q1(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension q2(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
c
c     # copy the contents of q1 into q2
c
      do m=1,meqn
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               q2(i,j,m) = q1(i,j,m)
               enddo
            enddo
         enddo
      return
c
      end
