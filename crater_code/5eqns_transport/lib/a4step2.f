c
c -----------------------------------------------------------
c
      subroutine a4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &           xlower,ylower,dx,dy,t,dt,q,qold,aux,auxold,
     &           dq,work,mwork,mthlim,mwaves,mthbc)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension dq(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension work(mwork)
      dimension mthlim(mwaves),mthbc(4)
c
      return
c
      end
