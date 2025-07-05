c
c ---------------------------------------------------------------------
c
      subroutine tclaw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,
     &           mthlim,q,work,mrkwork,rkwork,aux)
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension work(mwork),rkwork(mrkwork)
      dimension mthlim(mwaves)
      dimension method(10),dtv(5),cflv(4),nv(2),mthbc(4)
      dimension tout(100)
      common /restrt_block/ tinitial,iframe
      common /tinit/ t0
      logical rest
      character*10 fname
c
c     # driver routine 
c     # semi-discretized scheme for time-dependent PDEs
c
      open(55,file='claw2ez.data',status='old',form='formatted')
      open(10,file='fort.info',status='unknown',form='formatted')
c
c     # Read the input in standard form from claw2ez.data:

c     domain variables
      read(55,*) mx
      read(55,*) my

c     i/o variables
      read(55,*) nout
      read(55,*) outstyle
      if (outstyle.eq.1) then
          read(55,*) tfinal
          nstepout = 1
      elseif (outstyle.eq.2) then
          read(55,*) (tout(i), i=1,nout)
          nstepout = 1
      elseif (outstyle.eq.3) then
          read(55,*) nstepout, nstop
          nout = nstop
      endif
c
c     timestepping variables
      read(55,*) dtv(1)
      read(55,*) dtv(2)
      read(55,*) cflv(1)
      read(55,*) cflv(2)
      read(55,*) nv(1)
c
c     # input parameters for clawpack routines
      read(55,*) method(1)
      read(55,*) method(2)
      read(55,*) method(3)
      read(55,*) method(4)
      read(55,*) method(5)
      read(55,*) method(6)
      read(55,*) method(7)
      read(55,*) method(8)
      read(55,*) method(9)
      read(55,*) method(10)
c
      read(55,*) meqn1
      read(55,*) mwaves1
      read(55,*) (mthlim(mw),mw=1,mwaves1)
      read(55,*) t0
      read(55,*) xlower
      read(55,*) xupper
      read(55,*) ylower
      read(55,*) yupper
      read(55,*) mbc1
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)
c
c     # check to see if we are restarting:
      rest = .false.
c     # The next two lines may not exist in old versions of claw2ez.data.
c     # Jump over the second read statement if the 1st finds an EOF:
      read(55,*,end=199,err=199) rest
      read(55,*) iframe   !# restart from data in fort.qN file, N=iframe
 199  continue
c
      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
          write(6,*) '*** ERROR ***  periodic boundary conditions'
          write(6,*) 'require mthbc(1) and mthbc(2) BOTH be set to 2'
          stop
      endif
c
      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
          write(6,*) '*** ERROR ***  periodic boundary conditions'
          write(6,*) 'require mthbc(3) and mthbc(4) BOTH be set to 2'
          stop
      endif
c
c     # These values were passed in, but check for consistency:
c
      if (method(7) .ne. maux) then
          write(6,*) '*** ERROR ***  method(7) should equal maux'
          stop
      endif
      if (meqn1 .ne. meqn) then
          write(6,*) '*** ERROR ***  meqn set wrong in input or driver'
          stop
      endif
      if (mwaves1 .ne. mwaves) then
          write(6,*) 
     &     '*** ERROR ***  mwaves set wrong in input or driver'
          stop
      endif
      if (mbc1 .ne. mbc) then
          write(6,*) 
     &     '*** ERROR ***  mbc set wrong in input or driver'
          stop
      endif
c
c     # check storage allocation                     
      call claw_maloc(mx,my,maxmx,maxmy,meqn,mwaves,mbc,
     &     maux,mwork,mrkwork,method)
c
      write(6,*) 'running...'
      write(6,*) ' '
c
c     # grid spacing
      dx = (xupper-xlower)/dble(mx)
      dy = (yupper-ylower)/dble(my)
c
c     # time increments between outputing solution:
      if (outstyle .eq. 1) then
          dtout = (tfinal-t0)/dble(nout)
      endif
c
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      call setprob()
c
c     # set aux array
      if (maux .gt. 0)  then
c         # mesh
          call setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &         maux,aux)
c
c         # miscellaneous geometrical quantity
          call setaux_misc(maxmx,maxmy,mbc,mx,my,maux,
     &         xlower,ylower,dx,dy,aux,work,mwork,mthbc)
      endif
c
c     # set initial conditions:
      if (rest) then
          call restart(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &         dx,dy,q)
          t0 = tinitial
      else
c         # physical variables
          call qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &          dx,dy,q,maux,aux)
c
c         # boundary conditions
          call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &        dx,dy,q,maux,aux,told,dt,mthbc)
c
          iframe = 0
      endif
c
      if (.not. rest) then
c         # output initial data
          call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &         q,t0,iframe,aux,maux)
          write(6,601) iframe, t0
      endif
c
      if (method(1) .eq. 1) then
          call estdt(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &         dx,dy,dt0,cflv(2),q,aux)
c
          dtv(1) = dmin1(dtv(1),dt0)
      endif
c
c     # data diagnosis
      fname = 'fort.d0010'
      open(unit=11,file=fname,status='unknown',form='formatted')
c
      call dg2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,t0,q,q,aux,aux,maux,0,work,mwork)
c
c     # main loop
      tend = t0
      n0   = iframe*nstepout + 1
      do 100 n=n0,nout
         tstart = tend
         if (outstyle .eq. 1) tend = tstart+dtout
         if (outstyle .eq. 2) tend = tout(n)
         if (outstyle .eq. 3) tend = tstart-1.d0  !# single-step mode
c
         if (method(3) .eq. 0) then 
c            # dim-by-dim reconstruction
             call tclaw2(maxmx,maxmy,meqn,mwaves,mbc,maux,mx,my,
     &            q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &            cflv,nv,method,mthlim,mthbc,rkwork,mrkwork,
     &            work,mwork,info)
         else                       
c            # multi-d reconstruction
             write(66,*) 'error: multi-d reconstruction'
c             call tclaw2md(maxmx,maxmy,meqn,mwaves,mbc,maux,mx,my,
c     &            q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
c     &            cflv,nv,method,mthlim,mthbc,rkwork,mrkwork,
c     &            work,mwork,info,bc2,rpn2,tfluct2,src2,
c     &            a4step2,b4step2)
         endif
c
c        # check to see if an error occured:
         if (info .ne. 0) then
            write(6,*) 'tclaw2ez aborting: Error return from tclaw2',
     &                 ' with info =',info
            go to 999
         endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
c        # if outstyle=1 or 2, then nstepout=1 and we output every time
c        # we reach this point, since claw1 was called for the entire time
c        # increment between outputs.
c
c        # if outstyle=3 then we only output if we have taken nstepout
c        # time steps since the last output.

c        # iframe is the frame number used to form file names in out1
         iframe = n/nstepout
         if (iframe*nstepout .eq. n) then
             call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &            dx,dy,q,tend,iframe,aux,maux)
c
            write(6,601) iframe,tend
            write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &           cflv(3),cflv(4),nv(2)
         endif
c
c        # formats for writing out information about this call to claw:
c
  601    format('tclaw2ez: Frame ',i4,
     &           ' matlab plot files done at time t =',
     &           d12.4,/)
c
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
  100    continue
c
  999 continue
      return
c
      end
