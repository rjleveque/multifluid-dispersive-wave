c
c -------------------------------------------------------------
c
      subroutine dg2(maxmx,maxmy,meqn,mbc,mx,my,
     &           xlower,ylower,dx,dy,time,q,qold,
     &           aux,auxold,maux,nstep,work,mwork)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension work(mwork)
c
c     # data diagnosis 
c
      return
c
      end
