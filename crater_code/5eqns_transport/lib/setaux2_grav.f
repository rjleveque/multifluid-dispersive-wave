c
c -----------------------------------------------------------
c
      subroutine setaux_misc(maxmx,maxmy,mbc,mx,my,maux,
     &           xlower,ylower,dx,dy,aux,work,mwork,mthbc)
      implicit double precision (a-h,o-z)
      common /cgrav/ grav
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension work(mwork)
      dimension mthbc(4)
c
c     # set auxiliary state
c     # gravitational potential
c     # psi = g y
c
      do j=1-mbc,my+mbc
         yc = ylower+dble(j-1)*dy+0.5d0*dy
         do i=1-mbc,mx+mbc
            aux(i,j,2) = grav*yc             
            enddo
         enddo
      return
c
      end
