c
c ----------------------------------------------------------
c
      subroutine bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &           dx,dy,q,maux,aux,t,dt,mthbc)
      implicit double precision (a-h,o-z)
      parameter(rhoeps = 1.0d-8)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      common /bcs/ qin(20,4),qpin(20,4),qref(20)
      common /eos_local/ geos0,bn0
      common /grav_steady_info/ q0,ubc,H0,gS0,phibc
      integer mthbc(4)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension qbc(20)
      external grav_steady_sg,dgrav_steady_sg
c
c     # left boundary:
      go to (100,110,120,130) mthbc(1)+1
      goto 199
c
 100  continue
c     # supersonic inflow
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = qin(m,1)
               enddo
            enddo
         enddo
      go to 199
c
 110  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = q(1,j,m)
               enddo
            enddo
         enddo
      go to 199

 120  continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = q(mx+1-ibc,j,m)
               enddo
            enddo
         enddo
      go to 199
c
 130  continue
c     # solid wall
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
      do ibc=1,mbc
         do j=1-mbc,my+mbc
            q(1-ibc,j,3) = -q(1-ibc,j,3)
            enddo
         enddo
      go to 199
c
 199  continue
c     # right boundary:
      go to (200,210,220,230) mthbc(2)+1
      goto 299
c
 200  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

 210  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
               enddo
            enddo
         enddo
      go to 299

 220  continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(mx+ibc,j,m) = q(ibc,j,m)
               enddo
            enddo
         enddo
      go to 299

 230  continue
c     # solid wall
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
      do ibc=1,mbc
         do j=1-mbc,my+mbc
            q(mx+ibc,j,3) = -q(mx+ibc,j,3)
            enddo
         enddo
      go to 299

 299  continue
c     # bottom boundary:
      go to (300,310,320,330) mthbc(3)+1
      goto 399
c
 300  continue
c     # gravitational source terms
c     # steady-state equilibrium conditions
c
      do jbc=1,mbc
         do i=1-mbc,mx+mbc
            do m=1,meqn
               qbc(m) = q(i,jbc,m)
               enddo
c
            call ctomfd(qbc,rhoa0,rhob0,u0,v0,rhoh0,p0,
     &           grue0,dum2,zfa0,zfb0,c20)
c
c           # total density
            rho0 = rhoa0+rhob0
c
c           # mass fraction
            Y10  = rhoa0/rho0
            Y20  = rhob0/rho0
c
c           # constant momentum
            q0 = rho0*v0
            w0 = rho0*u0*v0
c
c           # constant entropy
            S0 = (p0+bn0)/rho0**geos0
c
c           # constant modified enthalpy
            phi0 = aux(i,jbc,2)
            H0   = (qbc(5)+p0)/rho0+phi0   
c
c           # transverse velocity
            ubc = u0
c
            gS0   = geos0/(geos0-1.0d0)*S0
            phibc = aux(i,jbc-1,2)
c
c           # solve for density
            rho00 = rho0
            rhobc = zero1d_aux(rho00,grav_steady_sg,
     &                dgrav_steady_sg,rhoeps,iflag)
c
            if (iflag .ne. 1) then
                write(66,*) 'error in bottom bc: iflag=',iflag
                write(66,*) 'qbc=',(qbc(m),m=1,meqn)
                write(66,*) rho0,u0,v0,p0
                stop
            endif
c
            vbc = q0/rhobc
            pbc = S0*rhobc**geos0-bn0
c
            rho1_bc = Y10*rhobc/zfa0
            rho2_bc = Y20*rhobc/zfb0
c
            call prmtoc(qbc,rho1_bc,rho2_bc,ubc,vbc,
     &           pbc,pbc,zfa0,zfb0)
c
            do m=1,meqn
               q(i,1-jbc,m) = qbc(m)
               enddo
            enddo
         enddo
      go to 399
c
 310  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
               enddo
            enddo
         enddo
      go to 399

 320  continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,my+1-jbc,m)
               enddo
            enddo
         enddo
      go to 399

 330  continue
c     # solid wall
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
      do jbc=1,mbc
         do i=1-mbc,mx+mbc
            q(i,1-jbc,4) = -q(i,1-jbc,4)
            enddo
         enddo
      go to 399

 399  continue
c     # top boundary:
      go to (400,410,420,430) mthbc(4)+1
      goto 499
c
 400  continue
c     # gravitational source terms
c     # steady-state equilibrium conditions
c
      do jbc=1,mbc
         do i=1-mbc,mx+mbc
            do m=1,meqn
               qbc(m) = q(i,my+1-jbc,m)
               enddo
c
            call ctomfd(qbc,rhoa0,rhob0,u0,v0,rhoh0,p0,
     &           grue0,zeta0,zfa0,zfb0,c20)
c
c           # total density
            rho0 = rhoa0+rhob0
c
c           # mass fraction
            Y10  = rhoa0/rho0
            Y20  = rhob0/rho0
c
c           # constant momentum
            q0 = rho0*v0
            w0 = rho0*u0*v0
c
c           # constant entropy
            S0 = (p0+bn0)/rho0**geos0
c
c           # constant modified enthalpy
            phi0 = aux(i,my+1-jbc,2)
            H0   = (qbc(5)+p0)/rho0+phi0
c
c           # transverse velocity
            ubc = u0
c
            gS0   = geos0/(geos0-1.0d0)*S0
            phibc = aux(i,my+1-jbc+1,2)
c
c           # solve for density
            rho00 = rho0
            rhobc = zero1d_aux(rho00,grav_steady_sg,
     &                dgrav_steady_sg,rhoeps,iflag)
c
            if (iflag .ne. 1) then
c                write(66,*) 'error in top bc: time=',t,
c     &                       ', iflag=',iflag
c                write(66,*) 'qbc=',(qbc(m),m=1,meqn)
c                write(66,*) 'rho0=',rho0,u0,v0,p0
c
c                stop
            else 
                vbc = q0/rhobc
                pbc = S0*rhobc**geos0-bn0
c
                rho1_bc = Y10*rhobc/zfa0
                rho2_bc = Y20*rhobc/zfb0
c
                call prmtoc(qbc,rho1_bc,rho2_bc,ubc,vbc,
     &               pbc,pbc,zfa0,zfb0)
            endif
c
            do m=1,meqn
               q(i,my+jbc,m) = qbc(m)
               enddo
            enddo
         enddo
      go to 499
c
 410  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
               enddo
            enddo
         enddo
      go to 499
c
 420  continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
               enddo
            enddo
         enddo
      go to 499
c
 430  continue
c     # solid wall
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
      do jbc=1,mbc
         do i=1-mbc,mx+mbc
            q(i,my+jbc,4) = -q(i,my+jbc,4)
            enddo
         enddo
      go to 499
c
 499  continue
      return
c
      end
c
c ------------------------------------------------------------
c
      double precision function grav_steady_sg(rho)
      implicit double precision (a-h,o-z)
      common /eos_local/ geos0,bn0 
      common /grav_steady_info/ q0,u0,H0,gS0,phi
c
c     # gravitational source terms
c     # steady equilibrium condition
c     # stiffened gas EOS case
c
      rhoK = 0.5d0*(q0**2/rho+rho*u0**2)
      grav_steady_sg = rho*(H0-phi)-
     &             (gS0*rho**geos0+rhoK)
c
c      write(66,*) 'rho=',rho,
c     &      ', grav_steady_sg=',grav_steady_sg
      return
c
      end
c
c ------------------------------------------------------------
c
      double precision function dgrav_steady_sg(rho)
      implicit double precision (a-h,o-z)
      common /eos_local/ geos0,bn0 
      common /grav_steady_info/ q0,u0,H0,gS0,phi
c
c     # gravitational source terms
c     # steady equilibrium condition
c     # stiffened gas EOS case
c
      drhoK = 0.5d0*(-q0**2/rho**2+u0**2)
      dgrav_steady_sg = (H0-phi)-
     &       (geos0*gS0*rho**(geos0-1.0d0)+drhoK)
c
c      write(66,*) 'rho=',rho,
c     &      ', dgrav_steady_sg=',dgrav_steady_sg
      return
      end
