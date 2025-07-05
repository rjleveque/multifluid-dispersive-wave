c
c --------------------------------------------------------------------
c
      subroutine ts2_erk2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
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
c     # SSPRK22 time stepping (2-stage,2nd order)
c     # 2 memory registers
c
c     # constants for RK time stepping
      c2 = 1.d0
c
c      write(66,*) 'stage 1'
c
c     # stage 1
c     # before step2
      call b4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,told,dt,
     &     q,qold,aux,auxold,
     &     work,mwork,mthlim,mwaves,mthbc,cflv2)
c
c     # store (dt)*L(q^n) = (dt)*d/dt(q^n(t^n) in qrk1
      call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &     mx,my,q,qold,aux,auxold,qrk1,
     &     dx,dy,dt,method,mthlim,cfl,
     &     dq1d,q1d,dtdx1d,dtdy1d,aux1,aux2,aux3,
     &     work,mwork,mthbc)
c
c     # q^(1) = q^n+dt*L(q^n) => q
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
c
c     # stage 2
      that = told+c2*dt
c
c     # extend data from grid to bordering boundary cells:
      call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux,that,dt,mthbc)
c
c     # before step2
      call b4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,that,dt,
     &     q,qold,aux,auxold,
     &     work,mwork,mthlim,mwaves,mthbc,cflv2)
c
c      write(66,*) 'stage 2'
c
c     # store (dt)*L(q^1) = (dt)*d/dt(q^(1)(t^(1)) in qrk1
      call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &     mx,my,q,qold,aux,auxold,qrk1,
     &     dx,dy,dt,method,mthlim,cfl,
     &     dq1d,q1d,dtdx1d,dtdy1d,aux1,aux2,aux3,
     &     work,mwork,mthbc)
c
c     # q = 1/2*u^n+1/2*u^(1)+1/2*dt*L(q^1)
      do m=1,meqn
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               q(i,j,m) = 0.5d0*(qold(i,j,m)+q(i,j,m)+qrk1(i,j,m))
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
