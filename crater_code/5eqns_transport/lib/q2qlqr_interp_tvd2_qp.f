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
c     # Second order TVD reconstruction 
c     # based on limited cell-average   
c
c     # initialization: piecewise constant
      do m=1,meqn
         do i=2-mbc,mx+mbc
            qr(i-1,m) = q(i-1,m)
            ql(i,m)   = q(i,m)
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

c     # compute limited slope
      do m=1,meqn
         do i=1-mbc,mx+mbc
            if (i .eq. 1-mbc) then
c               # one-sided difference
c
                dq(i,m) = qp(i+1,m)-qp(i,m)
            else if (i .eq. mx+mbc) then
c               # one-sided difference
c
                dq(i,m) = qp(i,m)-qp(i-1,m)
            else
c               # central difference
c
                dqL = qp(i,m)-qp(i-1,m)
                dqR = qp(i+1,m)-qp(i,m)
c
                if (dabs(dqL) .lt. 1.0d-16) then
                    dq(i,m) = 0.0d0
                else
                    dq(i,m) = philim(dqL,dqR,mthlim(1))*dqL
                endif
            endif
            enddo
         enddo
c
c     # piecewise linear reconstruction
      do m=1,meqn
         do i=1-mbc,mx+mbc
            uu(i,m,1) = qp(i,m)-0.5d0*dq(i,m)
            uu(i,m,2) = qp(i,m)+0.5d0*dq(i,m)
            enddo    
         enddo    
c
c     # change primitive to conservative
      call prm2c_edge(maxm,meqn,mbc,maux,mx,
     &     ql,qr,uu,q,aux)
      return
c
      end
