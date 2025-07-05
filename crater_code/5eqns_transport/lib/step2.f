c
c ---------------------------------------------------------------
c
      subroutine step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &           mx,my,q,qold,aux,auxold,dq,dx,dy,dt,
     &           method,mthlim,cfl,dq1d,q1d,dtdx1d,dtdy1d,
     &           aux1,aux2,aux3,work,mwork,mthbc)
      implicit double precision (a-h,o-z)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension dq(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
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
c     # Evaluate (delta t) *dq/dt
c     # On entry, q should give the initial data for this step
c     # dimension-by-dimension version
c     # (No genuinely multi-d terms)
c    
c     # dq1d is used to return increments to q from flux2
c     # See the flux2 documentation for more information.
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave  = 1
      i0s     = i0wave+(maxm+2*mbc)*meqn*mwaves
      i0amdq  = i0s+(maxm+2*mbc)*mwaves
      i0apdq  = i0amdq+(maxm+2*mbc)*meqn
      i0amdq2 = i0apdq+(maxm+2*mbc)*meqn
      i0apdq2 = i0amdq2+(maxm+2*mbc)*meqn
      i0ql    = i0apdq2+(maxm+2*mbc)*meqn
      i0qr    = i0ql+(maxm+2*mbc)*meqn
      i0auxl  = i0qr+(maxm+2*mbc)*meqn
      i0auxr  = i0auxl+(maxm+2*mbc)*maux
      i0evl   = i0auxr+(maxm+2*mbc)*maux
      i0evr   = i0evl+(maxm+2*mbc)*meqn*mwaves
      i0uu    = i0evr+(maxm+2*mbc)*meqn*mwaves
      i0hh    = i0uu+(maxm+2*mbc)*meqn*2   
      i0q1l   = i0hh+(maxm+2*mbc)*5
      i0q1r   = i0q1l+(maxm+2*mbc)*meqn
      i0q2l   = i0q1r+(maxm+2*mbc)*meqn
      i0q2r   = i0q2l+(maxm+2*mbc)*meqn
      i0next  = i0q2r+(maxm+2*mbc)*meqn
      i0used  = i0next-1
      if (i0used.gt.mwork) then
          write(6,*) '*** not enough work space in step2'
          write(6,*) '*** iused = ', iused, '   mwork =',mwork
          stop 
      endif
c
      mcapa = method(6)
      cfl   = 0.d0
      dtdx  = dt/dx
      dtdy  = dt/dy
c
      if (mcapa .eq. 0) then
c        # no capa array:
         do i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            enddo 
      endif
c
      do m=1,meqn
         do i=1-mbc,mx+mbc
            do j=1-mbc,my+mbc
               dq(i,j,m)=0.d0
               enddo
            enddo
         enddo
c
c     First evaluate source terms
c     Evaluate (delta t) * psi(q,t) and store in dq
      if (method(5) .eq. 1) then
          call src2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &         xlower,ylower,dx,dy,q,qold,aux,auxold,
     &         dq,work,mwork,t,dt)
      endif
c
c     # perform x-sweeps
      do j=0,my+1
c        # copy data along a slice into 1d arrays
         do m=1,meqn
            do i=1-mbc,mx+mbc
               q1d(i,m) = q(i,j,m)
               enddo
            enddo
c
         if (mcapa .gt. 0)  then
            do i=1-mbc,mx+mbc
               dtdx1d(i) = dtdx/aux(i,j,mcapa)
               enddo 
         endif
c
         if (maux .gt. 0)  then
             do ma=1,maux
               do i=1-mbc,mx+mbc
                  aux1(i,ma) = aux(i,j-1,ma)
                  aux2(i,ma) = aux(i,j,ma)
                  aux3(i,ma) = aux(i,j+1,ma)
                  enddo
               enddo 
         endif
c
c        # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c           
c        # compute modification dq1d along this slice:
         call flux2(1,maxm,meqn,mwaves,mbc,maux,mx,dt,dx,
     &        q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,dq1d,cfl1d,
     &        work(i0ql),work(i0qr),
     &        work(i0auxl),work(i0auxr),work(i0wave),work(i0s),
     &        work(i0amdq),work(i0apdq),work(i0amdq2),work(i0apdq2),
     &        work(i0evl),work(i0evr),work(i0uu),work(i0hh),
     &        work(i0q1l),work(i0q1r),work(i0q2l),work(i0q2r),
     &        mthbc)
c
         cfl = dmax1(cfl,cfl1d)
c
         do m=1,meqn
            do i=1,mx
               dq(i,j,m) = dq(i,j,m)+dq1d(i,m)
               enddo
            enddo
         enddo 
c
c     # perform y sweeps
      do i=0,mx+1
c        # copy data along a slice into 1d arrays:
         do m=1,meqn
            do j=1-mbc,my+mbc
               q1d(j,m) = q(i,j,m)
               enddo
            enddo 
c
         if (mcapa .gt. 0)  then
             do j=1-mbc,my+mbc
                dtdy1d(j) = dtdy/aux(i,j,mcapa)
                enddo
         endif
c
         if (maux .gt. 0)  then
             do mau=1,maux
                do j=1-mbc,my+mbc
                   aux1(j,mau) = aux(i-1,j,mau)
                   aux2(j,mau) = aux(i,  j,mau)
                   aux3(j,mau) = aux(i+1,j,mau)
                   enddo
                enddo
         endif
c
c        # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c           
         call flux2(2,maxm,meqn,mwaves,mbc,maux,my,dt,dy,
     &        q1d,dtdy1d,aux1,aux2,aux3,method,mthlim,dq1d,cfl1d,
     &        work(i0ql),work(i0qr),
     &        work(i0auxl),work(i0auxr),work(i0wave),work(i0s),
     &        work(i0amdq),work(i0apdq),work(i0amdq2),work(i0apdq2),
     &        work(i0evl),work(i0evr),work(i0uu),work(i0hh),
     &        work(i0q1l),work(i0q1r),work(i0q2l),work(i0q2r),
     &        mthbc)
c
         cfl = dmax1(cfl,cfl1d)
c
         do m=1,meqn
            do j=1,my
               dq(i,j,m) = dq(i,j,m)+dq1d(j,m)
               enddo
            enddo
         enddo  
      return
c
      end
