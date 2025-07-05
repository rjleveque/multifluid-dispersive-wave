c
c -----------------------------------------------------------
c
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &           dx,dy,maux,aux)
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
c
c     # radial distance
      do i=1-mbc,mx+mbc
         xc = xlower+dble(i-1)*dx+0.5d0*dx
         do j=1-mbc,my+mbc
            aux(i,j,1) = xc
            enddo
         enddo
      return
c
      end
