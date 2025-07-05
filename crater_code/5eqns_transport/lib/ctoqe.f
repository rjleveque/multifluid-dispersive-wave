c
c ----------------------------------------------------------
c
      subroutine c2prm_1d(maxm,meqn,mbc,maux,mx,q,qp,aux)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxm+mbc,meqn)
      dimension aux(1-mbc:maxm+mbc,maux)
      dimension qp(1-mbc:maxm+mbc,meqn)
c
c     # 5-equation transport model
c     # change conservative to primitive
c     # 2-phase case
c
      do i=1-mbc,mx+mbc
         qp(i,1) = q(i,1)
         qp(i,2) = q(i,2)
         qp(i,3) = q(i,3)/(q(i,1)+q(i,2))
         qp(i,4) = q(i,4)/(q(i,1)+q(i,2))
         qp(i,5) = q(i,5)-0.5d0*(q(i,3)**2+q(i,4)**2)/
     &             (q(i,1)+q(i,2))
         qp(i,6) = q(i,6)
         enddo
      return
c
      end
c
c --------------------------------------------------------------
c
      subroutine prm2c_edge(maxm,meqn,mbc,maux,mx,
     &           ql,qr,uu,q,aux)
      implicit double precision (a-h,o-z)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension q(1-mbc:maxm+mbc,meqn)
      dimension aux(1-mbc:maxm+mbc,maux)
      dimension uu(1-mbc:maxm+mbc,meqn,2)
c
c     # 5-equation transport model
c     # change primitive variables to conservative variables
c     # 2-phase case
c
      do i=0,mx+1
c        # left cell-edge 
c
         ql(i,1) = uu(i,1,1)
         ql(i,2) = uu(i,2,1)
         ql(i,3) = (ql(i,1)+ql(i,2))*uu(i,3,1)
         ql(i,4) = (ql(i,1)+ql(i,2))*uu(i,4,1)
         ql(i,5) = uu(i,5,1)+
     &             0.5d0*(ql(i,3)**2+ql(i,4)**2)/
     &             (ql(i,1)+ql(i,2))
         ql(i,6) = uu(i,6,1)
c
c        # right cell-edge
         qr(i,1) = uu(i,1,2)
         qr(i,2) = uu(i,2,2)
         qr(i,3) = (qr(i,1)+qr(i,2))*uu(i,3,2)
         qr(i,4) = (qr(i,1)+qr(i,2))*uu(i,4,2)
         qr(i,5) = uu(i,5,2)+
     &             0.5d0*(qr(i,3)**2+qr(i,4)**2)/
     &             (qr(i,1)+qr(i,2))
         qr(i,6) = uu(i,6,2)
         enddo
c
c     # positivity-preserving check
      do i=0,mx+1      
         if ((ql(i,1) .lt. 0.d0) .or. 
     &       (ql(i,2) .lt. 0.d0) .or.
     &       (qr(i,1) .lt. 0.d0) .or.
     &       (qr(i,2) .lt. 0.d0)) then
             write(66,*) 'i=',i
             write(66,*) 'left-state'
             write(66,*) (qr(i-1,m),m=1,meqn)
             write(66,*) 'right-state'
             write(66,*) (ql(i,m),m=1,meqn)
             write(66,*) 'orig-state'
             write(66,*) (q(i,m),m=1,meqn)
             stop
         endif
         enddo
      return
c
      end
