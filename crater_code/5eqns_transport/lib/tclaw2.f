c
c ---------------------------------------------------------------
c
      subroutine tclaw2(maxmx,maxmy,meqn,mwaves,mbc,maux,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,rkwork,mrkwork,
     &           work,mwork,info)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension work(mwork)
      dimension rkwork(mrkwork)
      dimension mthlim(mwaves)
      dimension method(10)
      dimension dtv(5),cflv(4),nv(2),mthbc(4)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c     # Solves a hyperbolic system of conservation laws in 
c     # two space dimensions of the general form
c  
c     capa * q_t + A q_x + B q_y = psi
c
c     The "capacity function" capa(x,y) and source term psi are optional 
c     # Beginning of tclaw2 code
c
      maxm = max0(maxmx,maxmy)
      info = 0
      t    = tstart
      maxn = nv(1)
      dt   = dtv(1)   
      cflmax = 0.d0
      dtmin = dt
      dtmax = dt
      nv(2) = 0
c
c     # check for errors in data:
      if (mx .gt. maxmx .or. my .gt. maxmy .or.
     &    mbc .lt. 2) then
          info = 1
          write(6,*) 'tclaw2 ERROR...  check mx,maxmx,my,maxmy,mbc'
          go to 900
      endif
c
      if (method(1) .eq. 0) then
c         # fixed size time steps.  Compute the number of steps:
          if (tend .lt. tstart) then
c             # single step mode 
              maxn = 1
          else
              maxn = (tend - tstart + 1d-10) / dt
              if (dabs(maxn*dt - (tend-tstart)) .gt.
     &                          1d-5*(tend-tstart)) then
c                # dt doesn't divide time interval integer number of times
                 info = 2
                 write(6,*) 
     &               'tclaw2 ERROR... dt does not divide (tend-tstart)'
                 go to 900
                 endif
          endif
      endif
c
      if (method(1) .eq. 1 .and. 
     &    cflv(2) .gt. cflv(1)) then
          info = 3
          write(6,*) 'tclaw2 ERROR...  cflv(2) > cflv(1)'
          go to 900
      endif
c
      if (method(6) .gt. method(7)) then
          info = 5
          write(6,*) 'tclaw2 ERROR...  method(6) > method(7)'
          go to 900
      endif
c
      if (method(5) .gt. 1) then
          write(*,*) 'Strang splitting not allowed'
          stop
      endif
c
      mrkwork0 = (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
      if (mrkwork .lt. mrkwork0) then
          info = 4
          write(6,*) 'tclaw2 ERROR... mrkwork should be increased to ',
     &               mrkwork0
          go to 900
      endif
c
c     # partition work array into pieces needed for local storage in 
c     # step2 routine. Find starting index of each piece:
c
      i0dq1d  = 1
      i0q1d   = i0dq1d+(maxm+2*mbc)*meqn
      i0dtdx  = i0q1d+(maxm+2*mbc)*meqn  
      i0dtdy  = i0dtdx+(maxm+2*mbc)
      i0qwork = i0dtdy+(maxm+2*mbc)
      i0awork = i0qwork+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
      i0aux1  = i0awork+(maxmx+2*mbc)*(maxmy+2*mbc)*maux 
      i0aux2  = i0aux1+(maxm+2*mbc)*maux
      i0aux3  = i0aux2+(maxm+2*mbc)*maux
      i0next  = i0aux3+(maxm+2*mbc)*maux  
      mused   = i0next-1                  
      mwork1  = mwork-mused               
c
c     # partition rkwork array into pieces for passing into ts_*
      i0qrk1 = 1
      if (method(2).eq.1 .or. method(2).eq.3) then
          i0rkend = i0qrk1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn-1
      endif
c
c     # main loop
c
      if (maxn.eq.0) go to 900
      do 100 n=1,maxn
         told = t   
c
c        # adjust dt to hit tend exactly if we're near end of computation
c        #  (unless tend < tstart, which is a flag to take only a single step)
         if (told+dt.gt.tend .and. tstart.lt.tend) dt = tend - told
c
 9       continue
c
c        # store dt and t in the common block comxyt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
         dycom = dy
c
c        # main steps
c
c        # extend data from grid to bordering boundary cells
         call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &        dx,dy,q,maux,aux,told,dt,mthbc)
c
         call auxbc2(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &        dx,dy,maux,aux,told,dt,mthbc)
c
         call copyq2(maxmx,maxmy,meqn,mbc,mx,my,q,work(i0qwork))
c
         call copyaux2(maxmx,maxmy,maux,mbc,mx,my,aux,work(i0awork))
c
c        # take one step on the conservation law
         go to (10,20,30) method(2)
c
 10      continue
c        # explicit first-order Runge-Kutta method 
         call ts2_erk1(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &        mx,my,xlower,ylower,q,work(i0qwork),
     &        aux,work(i0awork),dx,dy,told,dt,method,mthlim,cfl,
     &        mthbc,rkwork(i0qrk1),work(i0dq1d),
     &        work(i0q1d),work(i0dtdx),work(i0dtdy),
     &        work(i0aux1),work(i0aux2),work(i0aux3),
     &        work(i0next),mwork1,cflv(2))
         go to 299
c
 20      continue
c        # explicit second-order Runge-Kutta method 
         call ts2_erk2(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &        mx,my,xlower,ylower,q,work(i0qwork),
     &        aux,work(i0awork),dx,dy,told,dt,method,mthlim,cfl,
     &        mthbc,rkwork(i0qrk1),work(i0dq1d),
     &        work(i0q1d),work(i0dtdx),work(i0dtdy),
     &        work(i0aux1),work(i0aux2),work(i0aux3),
     &        work(i0next),mwork1,cflv(2))
         go to 299
c
 30      continue
c        # explicit third-order Runge-Kutta method 
         call ts2_erk3(maxm,maxmx,maxmy,meqn,mwaves,mbc,maux,
     &        mx,my,xlower,ylower,q,work(i0qwork),
     &        aux,work(i0awork),dx,dy,told,dt,method,mthlim,cfl,
     &        mthbc,rkwork(i0qrk1),work(i0dq1d),
     &        work(i0q1d),work(i0dtdx),work(i0dtdy),
     &        work(i0aux1),work(i0aux2),work(i0aux3),
     &        work(i0next),mwork1,cflv(2))
         go to 299
c
 299     continue
         t = told+dt
c
         call dg2(maxmx,maxmy,meqn,mbc,mx,my,
     &        xlower,ylower,dx,dy,t,
     &        q,work(i0qwork),aux,work(i0awork),
     &        maux,n,work(i0next),mwork1)
c
         if (method(4).eq.1) then
c            # verbose mode
             write(6,601) n,cfl,dt,t
         endif
c
c        # check to see if the Courant number was too large:
         if (cfl .le. cflv(1)) then
c            # accept this step
             cflmax = dmax1(cfl,cflmax)
         else
c            # Reject this step.  Reset q to qwork from previous time:
c            # Note that i0qwork points to work space where previous
c            # solution is stored.
             t = told
             call copyq2(maxmx,maxmy,meqn,mbc,mx,my,work(i0qwork),q)
             call copyaux2(maxmx,maxmy,maux,mbc,mx,my,work(i0awork),aux)
c
             if (method(4) .eq. 1) then
c                # verbose mode
                 write(6,602)
             endif
c
             write(66,*) 'reject this step: cfl=',cfl
             stop
c
             if (method(1) .eq. 1) then
c                # if variable dt, go back and take a smaller step.
                 dt = dmin1(dtv(2),dt*cflv(2)/cfl)
                 go to 9 
             else
c                # if fixed dt, give up and return
                 cflmax = dmax1(cfl,cflmax)
                 go to 900
             endif
         endif
c
         if (method(1) .eq. 1) then
c            # choose new time step if variable time step
             if (cfl.eq.0.d0) then 
                 dt = dtv(2)
             else
                 dt = dmin1(dtv(2),dt*cflv(2)/cfl)
             endif
             dtmin = dmin1(dt,dtmin)
             dtmax = dmax1(dt,dtmax)
         endif
c
c        # see if we are done:
         nv(2) = nv(2)+1
         if (t .ge. tend) go to 900
  100    continue
c
  900  continue
c 
c      # return information
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          write(6,*) 'tclaw2 ERROR...  too many timesteps'
          info = 11
       endif
c
       if (method(1).eq.0 .and. cflmax .gt. cflv(1)) then
c         # Courant number too large with fixed dt
          write(6,*) 'tclaw2 ERROR...  Courant number too large'
          info = 12
       endif
c
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt
       return 
c
 601   format('tclaw2... Step',i6,
     &        '   Courant number =',f6.3,'  dt =',d12.4,
     &        '  t =',d12.4)
 602   format('tclaw2 rejecting step... Courant number too large')
       end
