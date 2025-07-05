c
c----------------------------------------------------------------
c
      subroutine src2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &           xlower,ylower,dx,dy,q,qold,aux,auxold,
     &           dq,work,mwork,t,dt)
      implicit double precision(a-h,o-z)
      common /cgrav/ grav
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension dq(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension work(mwork)
      dimension qloc(20)
      dimension psi(20)
c
c     # source terms for cylindrical symmetry in 2d Euler equations
c     # about y-axis (so x=radius)
c
c     # aux(i,j,1) stores x coordinate of cell center
c
      ndim = 2

      do i=1,mx
         do j=1,my
            do m=1,meqn
               qloc(m) = q(i,j,m)
               enddo
c
            call ctomfd(qloc,rhoa0,rhob0,vx0,vy0,
     &           rhoh0,p0,grue0,zeta0,zfa0,zfb0,c20)
c
            gaxis  = -dble(ndim-1)/aux(i,j,1)
            psi(1) = dt*gaxis*qloc(1)*vx0
            psi(2) = dt*gaxis*qloc(2)*vx0
            psi(3) = dt*gaxis*(qloc(1)+qloc(2))*vx0**2
            psi(4) = dt*gaxis*(qloc(1)+qloc(2))*vx0*vy0
            psi(5) = dt*gaxis*rhoh0*vx0
            psi(6) = 0.d0
c
            do m=1,meqn
               dq(i,j,m) = dq(i,j,m)+psi(m)
               enddo   
            enddo   
         enddo   
c
c     # source terms for gravitational force
      do i=1,mx
         do j=1,my
c           # vertical momentum
            psi(4) = -dt*grav*(q(i,j,1)+q(i,j,2))
c
c           # total energy
            psi(5) = -dt*grav*q(i,j,4)
c     
            do m=4,5
               dq(i,j,m) = dq(i,j,m)+psi(m)
               enddo
            enddo
         enddo
      return
c
      end
