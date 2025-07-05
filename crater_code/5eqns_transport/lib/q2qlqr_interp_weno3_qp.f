c
c--------------------------------------------------------------------
c     
      subroutine q2qlqr_interp(ixy,maxm,meqn,mwaves,mbc,maux,mx,
     &           dt,dx,q,aux,ql,qr,auxl,auxr,s,wave,qp,dq,uu,hh,
     &           evl,evr,mthlim)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxm+mbc,meqn)
      dimension aux(1-mbc:maxm+mbc,maux)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension auxl(1-mbc:maxm+mbc,maux)
      dimension auxr(1-mbc:maxm+mbc,maux)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension qp(1-mbc:maxm+mbc,meqn)
      dimension dq(1-mbc:maxm+mbc,meqn)
      dimension uu(1-mbc:maxm+mbc,meqn,2)
      dimension hh(1-mbc:maxm+mbc,-2:2)
      dimension evl(1-mbc:maxm+mbc,meqn,mwaves)
      dimension evr(1-mbc:maxm+mbc,meqn,mwaves)
      dimension mthlim(mwaves)
c
c     # third order WENO reconstruction
c     # primitive variables case
c
      epweno   = 1.0d-16
      onethird = 1.d0/3.d0
      twothird = 2.d0/3.d0
c
c     # initialization: piecewise constant
      do m=1,meqn
         do i=1-mbc,mx+mbc
            ql(i,m) = q(i,m)
            qr(i,m) = q(i,m)
            enddo
         enddo
c
      do ma=1,maux
         do i=2-mbc,mx+mbc
            auxr(i-1,ma) = aux(i-1,ma)
            auxl(i,ma)   = aux(i,ma)
            enddo
         enddo
c
c     # change conservative to primitive
      call c2prm_1d(maxm,meqn,mbc,maux,mx,q,qp,aux)
c
      do m=1,meqn
         do i=2-mbc,mx+mbc
            dq(i,m) = qp(i,m)-qp(i-1,m)
            enddo
         enddo
c
c     # component-by-component reconstruction
      do m=1,meqn
         do i=1-mbc,mx+mbc
c           # smoothness indicator
            beta0 = dq(i+1,m)**2
            beta1 = dq(i,m)**2
c
c           # left-hand side
            d0        = onethird   
            d1        = twothird  
            alpha0    = d0/(epweno+beta0)**2
            alpha1    = d1/(epweno+beta1)**2
            omega0    = alpha0/(alpha0+alpha1)
            omega1    = alpha1/(alpha0+alpha1)
            uu0       = 0.5d0*(3.d0*qp(i,m)-qp(i+1,m))
            uu1       = 0.5d0*(qp(i-1,m)+qp(i,m))
            uu(i,m,1) = omega0*uu0+omega1*uu1
c
c           # right-hand side
            d0      = twothird  
            d1      = onethird   
            alpha0  = d0/(epweno+beta0)**2
            alpha1  = d1/(epweno+beta1)**2
            omega0  = alpha0/(alpha0+alpha1)
            omega1  = alpha1/(alpha0+alpha1)
            uu0     = 0.5d0*(qp(i,m)+qp(i+1,m))
            uu1     = 0.5d0*(-qp(i-1,m)+3.d0*qp(i,m))
            uu(i,m,2) = omega0*uu0+omega1*uu1
            enddo      
         enddo      
c
c     # change primitive to conservative
      call prm2c_edge(maxm,meqn,mbc,maux,mx,
     &     ql,qr,uu,q,aux)
      return
c
      end
