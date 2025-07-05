c
c ---------------------------------------------------------------------
c
      subroutine ts2_erk3(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
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
c     # SSPRK33 time stepping (3-stage, 3rd order)
c     # 2 memory registers
c
c     Constants for RK time stepping
      third = 1.d0/3.d0
      c2    = 1.d0
      c3    = 0.5d0
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
c     # extend data from grid to bordering boundary cells:
      call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux,that,dt,mthbc)
c
c     # before step2
      call b4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,told,dt,
     &     q,qold,aux,auxold,
     &     work,mwork,mthlim,mwaves,mthbc,cflv2)
c
c     # store (dt)*L(q^1) = (dt)*d/dt(q^(1)(t^(1)) in qrk1
      call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &     mx,my,q,qold,aux,auxold,qrk1,
     &     dx,dy,dt,method,mthlim,cfl,
     &     dq1d,q1d,dtdx1d,dtdy1d,aux1,aux2,aux3,
     &     work,mwork,mthbc)
c
c     # q^(2) = 3/4*q^n+1/4*(q^(1)+dt*L(q^1)) => q
      do m=1,meqn
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               q(i,j,m) = 0.75d0*qold(i,j,m)+
     &                    0.25d0*(q(i,j,m)+qrk1(i,j,m))
               enddo
            enddo
         enddo
c
c     # after step2
      call a4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,told,dt,q,qold,aux,auxold,
     &     qrk1,work,mwork,mthlim,mwaves,mthbc)
c
c     # stage 3
      that = told+c3*dt
c     # extend data from grid to bordering boundary cells:
      call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux,that,dt,mthbc)
c
c     # before step2
      call b4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &     xlower,ylower,dx,dy,told,dt,
     &     q,qold,aux,auxold,
     &     work,mwork,mthlim,mwaves,mthbc,cflv2)
c
c     # store (dt)*L(q^2) = (dt)*d/dt(q^(2)(t^(2)) in qrk1
      call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &     mx,my,q,qold,aux,auxold,qrk1,
     &     dx,dy,dt,method,mthlim,cfl,
     &     dq1d,q1d,dtdx1d,dtdy1d,aux1,aux2,aux3,
     &     work,mwork,mthbc)
c
c     # q^(n+1) = 1/3*q^n+2/3*(q^(2)+dt*L(q^2)) => q
      do m=1,meqn
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               q(i,j,m) = third*(qold(i,j,m)+
     &                    2.d0*(q(i,j,m)+qrk1(i,j,m)))
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
