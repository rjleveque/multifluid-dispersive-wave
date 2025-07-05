c
c ---------------------------------------------------------------------
c
      subroutine claw_maloc(mx,my,maxmx,maxmy,meqn,mwaves,mbc,
     &           maux,mwork,mrkwork,method)
      implicit double precision (a-h,o-z)
      dimension method(10)
c
c     # check that enough storage has been allocated:
c
      if (method(5) .ge. 2) then
          write(*,*) 'Strang splitting not allowed with RK stepping!'
          stop
      endif
c
      maxm = max0(maxmx,maxmy)
c
c     # storage for Runge-Kutta stages:
      mrkwork1 = (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
c
      if (method(2) .eq. 4) then
c         #SSPRK95 - 11 registers
          mrkwork1 = 10*mrkwork1
      endif
c 
      mwork1 = 0
      if (method(3) .eq. 0) then 
c         # dim-by-dim reconstruction
c
c         # allocated in tclaw2
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*maux
          mwork1 = mwork1+(maxm+2*mbc)*maux
          mwork1 = mwork1+(maxm+2*mbc)*maux
          mwork1 = mwork1+(maxm+2*mbc)*maux
c
c         # allocated in step2
          mwork1 = mwork1+(maxm+2*mbc)*meqn*mwaves
          mwork1 = mwork1+(maxm+2*mbc)*mwaves
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*maux
          mwork1 = mwork1+(maxm+2*mbc)*maux
          mwork1 = mwork1+(maxm+2*mbc)*meqn*mwaves
          mwork1 = mwork1+(maxm+2*mbc)*meqn*mwaves
          mwork1 = mwork1+(maxm+2*mbc)*meqn*2
          mwork1 = mwork1+(maxm+2*mbc)*5
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
c
          if (method(10) .eq. 1) then
c             # in a4step2 (anti-diffusion interface-sharpening step) 
c
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+
     &                 (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
              mwork1 = mwork1+
     &                 (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
              mwork1 = mwork1+
     &                 (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
              mwork1 = mwork1+(maxm+2*mbc)*meqn
          elseif (method(10) .eq. 2) then
c             # in an2step (diffusive step for drift flux)
c
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
              mwork1 = mwork1+(maxm+2*mbc)*meqn
          elseif (method(10) .eq. 3) then
c             # in an2step (CSF step for surface tension)
c
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*2
          elseif (method(10) .eq. 4) then
c             # in b4step2 (VSF step for surface tension)
c
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          elseif (method(10) .eq. 5) then
c             # in b4step2 (split VSF step for surface tension)
c
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxm+2*mbc)*meqn
          elseif (method(10) .eq. 6) then
c             # in b4step2 (interface direction)
c
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
              mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          endif
      else
          write(6,*) 'error: fully multi-dimensional case'
          stop
      endif
c
      if ((mx .gt. maxmx) .or. (my .gt. maxmy) .or. 
     &    (mwork .lt. mwork1) .or. (mrkwork .lt. mrkwork1)) then
c         # insufficient storage
          maxmx1 = max0(mx,maxmx)
          maxmy1 = max0(my,maxmy)
          maxm1  = max0(maxmx1,maxmy1)
c
          mwork1 = 0
          if (method(3) .eq. 0) then 
c             # do dim-by-dim reconstruction

c             # allocated in tclaw2
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*maux
              mwork1 = mwork1+(maxm1+2*mbc)*maux
              mwork1 = mwork1+(maxm1+2*mbc)*maux
              mwork1 = mwork1+(maxm1+2*mbc)*maux
c
c             # allocated in step2
              mwork1 = mwork1+(maxm1+2*mbc)*meqn*mwaves
              mwork1 = mwork1+(maxm1+2*mbc)*mwaves
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*maux
              mwork1 = mwork1+(maxm1+2*mbc)*maux
              mwork1 = mwork1+(maxm1+2*mbc)*meqn*mwaves
              mwork1 = mwork1+(maxm1+2*mbc)*meqn*mwaves
              mwork1 = mwork1+(maxm1+2*mbc)*meqn*2
              mwork1 = mwork1+(maxm1+2*mbc)*5
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
c
              if (method(10) .eq. 1) then
c                 # in a4step2 (anti-diffusion interface-sharpening step) 
c
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
                  mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
                  mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
                  mwork1 = mwork1+(maxm1+2*mbc)*meqn
              elseif (method(10) .eq. 2) then
c                 # in an2step (diffusive step for drift flux)
c
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
                  mwork1 = mwork1+(maxm1+2*mbc)*meqn
              elseif (method(10) .eq. 3) then
c                 # in an2step (CSF step for surface tension)
c
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*2
              elseif (method(10) .eq. 4) then
c                 # in b4step2 (VSF step for surface tension)
c
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              elseif (method(10) .eq. 5) then
c                 # in b4step2 (split VSF step for surface tension)
c
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxm1+2*mbc)*meqn
              elseif (method(10) .eq. 6) then
c                 # in b4step2 (interface direction)
c
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
                  mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              endif
          else
              write(6,*) 'error: fully multi-dimensional case'
              stop
          endif
c
          mrkwork1 = (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
c
          if (method(2) .eq. 4) then
c             # SSPRK95 - 11 registers
              mrkwork1 = 10*mrkwork1
          endif
c
          write(6,*) ' '
          write(6,*) '*** ERROR *** Insufficient storage allocated'
          write(6,*) 'Recompile after increasing values in driver.f:'
          write(6,611) maxmx1
          write(6,612) maxmy1
          write(6,613) mwork1
          write(6,615) mrkwork1
          stop
      endif
      return
c
 611  format(/,'parameter (maxmx = ',i5,')')
 612  format('parameter (maxmy = ',i5,')')
 613  format('parameter (mwork = ',i9,')',/)
 615  format('parameter (mrkwork = ',i9,')',/)
      end
