c
c ----------------------------------------------------------
c
      subroutine auxbc2_lr(ixy,maxm,mbc,mx,xlower,dx,
     &           maux,auxl,auxr,t,dt,mthbc)
      implicit double precision (a-h,o-z)
      common /auxbc_info/ ma1,ma2
      dimension auxl(1-mbc:maxm+mbc,maux)
      dimension auxr(1-mbc:maxm+mbc,maux)
      dimension rot(4),uv(2)
      dimension mthbc(4)
c
c     # set cell-edge boundary condition
c     # auxiliary states
c
c     # dummy routine
c
      return
c
      end
