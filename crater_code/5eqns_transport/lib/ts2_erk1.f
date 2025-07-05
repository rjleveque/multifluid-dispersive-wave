c
c --------------------------------------------------------------------
c
      subroutine ts2_erk1(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &           mx,my,xlower,ylower,q,qold,aux,auxold,
     &           dx,dy,told,dt,method,mthlim,cfl,
     &           mthbc,qrk1,dq1d,
     &           q1d,dtdx1d,dtdy1d,aux1,aux2,aux3,
     &           work,mwork,cflv2)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension qrk1(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension q1d(1-mbc:maxm+mbc,meqn)
      dimension dq1d(1-mbc:maxm+mbc,meqn)
      dimension aux1(1-mbc:maxm+mbc,maux)
      dimension aux2(1-mbc:maxm+mbc,maux)
      dimension aux3(1-mbc:maxm+mbc,maux)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension method(10)
      dimension mthlim(mwaves)
      dimension mthbc(4)
      dimension work(mwork)
c
c     # first order forward Euler
c
      call b4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,told,dt,
     &     q,qold,aux,auxold,
     &     work,mwork,mthlim,mwaves,mthbc,cflv2)
c
      call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &     mx,my,q,qold,aux,auxold,qrk1,
     &     dx,dy,dt,method,mthlim,cfl,
     &     dq1d,q1d,dtdx1d,dtdy1d,aux1,aux2,aux3,
     &     work,mwork,mthbc)
c
c     # q = q^n+dt*L(q^n)
      do m=1,meqn
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               q(i,j,m)= q(i,j,m)+qrk1(i,j,m)
               enddo
            enddo
         enddo
c
c     # after step2
      call a4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,told,dt,q,qold,aux,auxold,
     &     qrk1,work,mwork,mthlim,mwaves,mthbc)
      return
c
      end
