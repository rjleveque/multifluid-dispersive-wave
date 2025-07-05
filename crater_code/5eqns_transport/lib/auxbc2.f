c
c ----------------------------------------------------------
c
      subroutine auxbc2(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &           dx,dy,maux,aux,t,dt,mthbc)
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      integer mthbc(4)
c
c     # dummy routine
c
      return
c
      end
