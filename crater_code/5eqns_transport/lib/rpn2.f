c
c ---------------------------------------------------------------
c
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,maux,mx,ql,qr,
     &           auxl,auxr,wave,s,amdq,apdq)
      implicit double precision (a-h,o-z)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension amdq(1-mbc:maxm+mbc,meqn)
      dimension apdq(1-mbc:maxm+mbc,meqn)
      dimension auxl(1-mbc:maxm+mbc,maux)
      dimension auxr(1-mbc:maxm+mbc,maux)
      dimension ql_state(20),qr_state(20)
      dimension aux1(30),aux2(30)
      dimension delta(20)
      dimension wave_local(20,20)
      dimension s_local(20)
      dimension speeds(20,2)
c
c     # Riemann solver in the normal direction for the Euler
c     # equations on a cartesian grid
c
      if (ixy .eq. 1) then
          mu = 3
          mv = 4
      else
          mu = 4
          mv = 3
      endif
c
      do i=2-mbc,mx+mbc-1
         do m=1,meqn
            ql_state(m) = qr(i-1,m)
            qr_state(m) = ql(i,m)
            enddo
c
         ql_state(3) = qr(i-1,mu)
         ql_state(4) = qr(i-1,mv)
         qr_state(3) = ql(i,mu)
         qr_state(4) = ql(i,mv)
c
         do ma=1,maux
            aux1(ma) = auxr(i-1,ma)
            aux2(ma) = auxl(i,ma)
            enddo
c
         do m=1,meqn
            delta(m) = qr_state(m)-ql_state(m)
            enddo
c
         call rp_solver(ql_state,qr_state,aux1,aux2,
     &        delta,s_local,wave_local,meqn)
c
         area = 1.0d0               
         do mw=1,mwaves
            wave_3loc = wave_local(3,mw)
            wave_4loc = wave_local(4,mw)
c
            wave_local(mu,mw) = wave_3loc
            wave_local(mv,mw) = wave_4loc
c
            speeds(mw,1) = area*dmin1(s_local(mw),0.d0)
            speeds(mw,2) = area*dmax1(s_local(mw),0.d0)
            s(i,mw) = speeds(mw,1)+speeds(mw,2)
            enddo

         do m=1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw=1,mwaves
               wave(i,m,mw) = wave_local(m,mw)
               amdq(i,m) = amdq(i,m)+speeds(mw,1)*wave(i,m,mw)
               apdq(i,m) = apdq(i,m)+speeds(mw,2)*wave(i,m,mw)
               enddo
            enddo
         enddo
      return
c
      end
