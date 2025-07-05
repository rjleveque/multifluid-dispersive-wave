c
c ----------------------------------------------------------------
c
      subroutine q2qlqr(ixy,maxm,meqn,mbc,maux,mx,dt,dx,
     &           q1d,aux1,aux2,aux3,ql,qr,auxl,auxr,qp,iflag)
      implicit double precision (a-h,o-z)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      common /thinc0/ thinc_beta,thinc_eta,gp1,gp2,vofmin,
     &        intf_test1,intf_test2
      dimension q1d(1-mbc:maxm+mbc,meqn)
      dimension aux1(1-mbc:maxm+mbc,maux)
      dimension aux2(1-mbc:maxm+mbc,maux)
      dimension aux3(1-mbc:maxm+mbc,maux)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension auxl(1-mbc:maxm+mbc,maux)
      dimension auxr(1-mbc:maxm+mbc,maux)
      dimension qp(1-mbc:maxm+mbc,meqn)
      dimension qloc(20),dq(20)
c
c     # THINC-M interface-sharpening reconstruction
c     # Ref: Cassidy, Edwards, Tian, JCP 228 (2009) 5628-5649
c
      zeps = 1.0d-12
      beta = thinc_beta
c
      if (iflag .eq.  1) then
c         # default: piecewise constant reconstruction
          do m=1,meqn
             do i=2-mbc,mx+mbc
                qr(i-1,m) = q1d(i-1,m)
                ql(i,m)   = q1d(i,m)
                enddo
             enddo
c
          do ma=1,maux 
             do i=2-mbc,mx+mbc
                auxr(i-1,ma) = aux2(i-1,ma)
                auxl(i,ma)   = aux2(i,ma)
                enddo
             enddo
      else
          do i=2-mbc,mx+mbc
             zfa0 = q1d(i,6)
             zfa1 = q1d(i-1,6)
             zfa2 = q1d(i+1,6)
c
c            # test for interface
             if (dmin1(zfa0,1.0d0-zfa0) .gt. vofmin) then 
                 if ((zfa0-zfa1)*(zfa2-zfa0) .gt. 0.d0) then
c                    # monotone profile
c
c                    # reconstruction: geometrical variable
                     call thinc_LR(zfaL,zfaR,zfa1,zfa0,zfa2,beta,zeps)
c
                     do m=1,meqn
                        qloc(m) = q1d(i,m)
                        enddo
c
                     call ctomfd(qloc,rhoa0,rhob0,vx0,vy0,rhoh0,
     &                    pr0,grue0,zeta0,zfa0,zfb0,c20)
c
                     rho10 = rhoa0/zfa0
                     rho20 = rhob0/zfb0
                     rhoe1 = rhoe_prho_ideal(pr0,geos(1))
                     rhoe2 = rhoe_prho_sg(pr0,geos(2),bn(2))
                     engk  = 0.5d0*(vx0**2+vy0**2)
c
c                    # cell edge on the left
                     dq(6) = zfaL-q1d(i,6)  
                     dq(1) = rho10*dq(6)
                     dq(2) = -rho20*dq(6)
                     dq(3) = vx0*(dq(1)+dq(2))
                     dq(4) = vy0*(dq(1)+dq(2))
                     dengk = engk*(dq(1)+dq(2))
                     drhoe = (rhoe1-rhoe2)*dq(6)
                     dq(5) = drhoe+dengk
c
                     do m=1,meqn
                        ql(i,m) = q1d(i,m)+dq(m)
                        enddo
c
c                    # cell edge on the right
                     dq(6) = zfaR-q1d(i,6) 
                     dq(1) = rho10*dq(6)
                     dq(2) = -rho20*dq(6)
                     dq(3) = vx0*(dq(1)+dq(2))
                     dq(4) = vy0*(dq(1)+dq(2))
                     dengk = engk*(dq(1)+dq(2))
                     drhoe = (rhoe1-rhoe2)*dq(6)
                     dq(5) = drhoe+dengk
c
                     do m=1,meqn
                        qr(i,m) = q1d(i,m)+dq(m)
                        enddo
c
c                    # positivity-preserving check
                     if ((ql(i,1) .lt. 0.d0) .or.
     &                   (ql(i,2) .lt. 0.d0) .or.
     &                   (qr(i,1) .lt. 0.d0) .or.
     &                   (qr(i,2) .lt. 0.d0)) then
                          write(66,*) 'i=',i
                          write(66,*) 'left-state'
                          write(66,*) (ql(i,m),m=1,meqn)
                          write(66,*) 'right-state'
                          write(66,*) (qr(i,m),m=1,meqn)
                          write(66,*) 'orig-state'
                          write(66,*) (q1d(i,m),m=1,meqn)
                          stop
                     endif          
                 endif
             endif
             enddo
      endif
      return
c
      end
