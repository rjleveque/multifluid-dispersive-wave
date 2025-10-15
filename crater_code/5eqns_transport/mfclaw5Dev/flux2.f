c
c ---------------------------------------------------------------
c
      subroutine flux2(ixy,maxm,meqn,maux,mbc,mx,
     &           q1d,dtdx1d,aux1,aux2,aux3,
     &           qadd,fadd,gadd,
     &           cfl1d,wave,s,amdq,apdq,
     &           sigma,
     &           cqxx,bmasdq,bpasdq,rpn2,rpt2)
     
      implicit double precision (a-h,o-z)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      dimension q1d(1-mbc:maxm+mbc,meqn)
      dimension amdq(1-mbc:maxm+mbc,meqn)
      dimension apdq(1-mbc:maxm+mbc,meqn)
      dimension bmasdq(1-mbc:maxm+mbc,meqn)
      dimension bpasdq(1-mbc:maxm+mbc,meqn)
      dimension cqxx(1-mbc:maxm+mbc,meqn)
      dimension qadd(1-mbc:maxm+mbc,meqn)
      dimension fadd(1-mbc:maxm+mbc,meqn)
      dimension gadd(1-mbc:maxm+mbc,meqn,2)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension aux1(1-mbc:maxm+mbc,maux+1)
      dimension aux2(1-mbc:maxm+mbc,maux+1)
      dimension aux3(1-mbc:maxm+mbc,maux+1)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension sigma(1-mbc:maxm+mbc,meqn,mwaves)
      dimension method(10),mthlim(mwaves)
      external rpn2,rpt2

      mwaves = NEED TO SEt
      method  = NEED TO SEt,
      mthlim = NEED TO SET
c
c     # initialize flux increments
      do m=1,meqn
         do i=1-mbc,mx+mbc
            qadd(i,m) = 0.d0
            fadd(i,m) = 0.d0
            gadd(i,m,1) = 0.d0
            gadd(i,m,2) = 0.d0
            enddo    
         enddo    
c
c     # solve Riemann problem at each interface and compute Godunov updates
      call rpn2(ixy,maxm,meqn,mwaves,mbc,maux,mx,q1d,q1d,aux2,aux2,
     &     wave,s,amdq,apdq)
c
c     # Set qadd for the donor-cell upwind method (Godunov)
      do m=1,meqn
         do i=1,mx+1
            qadd(i,m)   = qadd(i,m)-dtdx1d(i)*apdq(i,m)
            qadd(i-1,m) = qadd(i-1,m)-dtdx1d(i-1)*amdq(i,m)
            enddo
         enddo 
c
c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.d0
      do mw=1,mwaves
         do i=1,mx+1
c           # if s>0 use dtdx1d(i) to compute CFL,
c           # if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, dtdx1d(i)*s(i,mw), 
     &                          -dtdx1d(i-1)*s(i,mw))
            enddo
         enddo
c
      if (method(2) .eq. 1) go to 130
c
c     # modify F fluxes for second order q_{xx} correction terms:
c     # apply limiter to waves:
      call limiter(maxm,meqn,mwaves,mbc,mx,sigma,wave,s,mthlim)
c
      do i=2-mbc,mx+mbc
c        # For correction terms below, need average of dtdx in cell
c        # i-1 and i.  Compute these and overwrite dtdx1d:
c
c        # modified in Version 4.3 to use average only in cqxx, not transverse
c         dtdxave = 0.5d0*(dtdx1d(i-1)+dtdx1d(i))
         dtdxave = 2.0d0/(1.d0/dtdx1d(i-1)+1.d0/dtdx1d(i))
c
         do m=1,meqn
            cqxx(i,m) = 0.d0
            do mw=1,mwaves
c
c              # second order corrections:
               cqxx(i,m) = cqxx(i,m)+dabs(s(i,mw))*
     &                     (1.d0-dabs(s(i,mw))*dtdxave)*sigma(i,m,mw)
               enddo
            fadd(i,m) = fadd(i,m)+0.5d0*cqxx(i,m)
            enddo 
         enddo 
c
  130  continue
c
      if (method(3).le.0) go to 999   !# no transverse propagation
c
      if (method(2) .gt. 1 .and. method(3) .eq. 2) then
c         # incorporate cqxx into amdq and apdq so that it is split also.
          do m=1,meqn
             do i=1,mx+1
                amdq(i,m) = amdq(i,m)+cqxx(i,m)
                apdq(i,m) = apdq(i,m)-cqxx(i,m)
                enddo
             enddo
      endif
c
c     # modify G fluxes for transverse propagation
c     # split the left-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,maux,mx,q1d,q1d,aux1,aux2,aux3,
     &     1,wave,s,amdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^- Delta q and  B^+ A^- Delta q:
      do m=1,meqn
         do i=1,mx+1
            gadd(i-1,m,1) = gadd(i-1,m,1)- 
     &                      0.5d0*dtdx1d(i-1)*bmasdq(i,m)
            gadd(i-1,m,2) = gadd(i-1,m,2)-
     &                      0.5d0*dtdx1d(i-1)*bpasdq(i,m)
            enddo
         enddo 
c
c     # split the right-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,maux,mx,q1d,q1d,aux1,aux2,aux3,
     &     2,wave,s,apdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
      do m=1,meqn
         do i=1,mx+1
            gadd(i,m,1) = gadd(i,m,1)- 
     &                    0.5d0*dtdx1d(i)*bmasdq(i,m)
            gadd(i,m,2) = gadd(i,m,2)- 
     &                    0.5d0*dtdx1d(i)*bpasdq(i,m)
            enddo
         enddo   
  999 continue
      return
c
      end
