c
c--------------------------------------------------------------------
c     
      subroutine q2qlqr_interp(ixy,maxm,meqn,mwaves,mbc,maux,mx,
     &           dt,dx,q,aux,ql,qr,s,wave,qp,dq,uu,hh,
     &           evl,evr,mthlim)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxm+mbc,meqn)
      dimension aux(1-mbc:maxm+mbc,maux)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
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
c     # fifth order WENO reconstruction
c
      epweno = 1.e-16
c
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
         do m1=1,2
c           # m1=1: construct ql
c           # m1=2: construct qr
c
            im    = (-1)**(m1+1)
            ione  = im
            inone = -im
            intwo = -2*im
c 
            do i=1,mx+1
               t1  = dble(im)*(dq(i+intwo,m)-dq(i+inone,m))
               t2  = dble(im)*(dq(i+inone,m)-dq(i,m))
               t3  = dble(im)*(dq(i,m)-dq(i+ione,m))
               tt1 = 13.d0*t1**2+
     &               3.d0*(dq(i+intwo,m)-3.d0*dq(i+inone,m))**2
               tt2 = 13.d0*t2**2+
     &               3.d0*(dq(i+inone,m)+dq(i,m))**2
               tt3 = 13.d0*t3**2+
     &               3.d0*(3.d0*dq(i,m)-dq(i+ione,m))**2
               tt1 = (epweno+tt1)**2
               tt2 = (epweno+tt2)**2
               tt3 = (epweno+tt3)**2
               tt1 = tt1**2
               tt2 = tt2**2
               tt3 = tt3**2
               s1  = tt2*tt3
               s2  = 6.d0*tt1*tt3
               s3  = 3.d0*tt1*tt2
               t0  = 1.d0/(s1+s2+s3)
               s1  = s1*t0
               s3  = s3*t0
c 
               uu(i,m,m1) = (s1*(t2-t1)+
     &                  (0.5d0*s3-0.25d0)*(t3-t2))/3.d0
               enddo    
            enddo   
         enddo      
c
      do m=1,meqn
         do i=1,mx+1
            qr(i-1,m) = (-qp(i-2,m)+
     &                  7.d0*(qp(i-1,m)+qp(i,m))-qp(i+1,m))/12.d0
            ql(i,m)   = qr(i-1,m)
c
            qr(i-1,m) = qr(i-1,m)+uu(i,m,1)
            ql(i,m)   = ql(i,m)+uu(i,m,2)
            enddo
         enddo
c
      do m=1,meqn
         do i=2-mbc,mx+mbc
            uu(i,m,1) = ql(i,m)   
            uu(i,m,2) = qr(i,m)   
            enddo
         enddo

c     # change primitive to conservative
      call prm2c_edge(maxm,meqn,mbc,maux,mx,
     &     ql,qr,uu,q,aux)
      return
c
      end
