c
c ---------------------------------------------------------------
c
      subroutine estdt(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &           dx,dy,dt,cfl,q,aux)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension qloc(20)
c
c     # set initial variable time step
c
      cmax = 0.d0
c
      do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
            do m=1,meqn
               qloc(m) = q(i,j,m)
               enddo
c            
            call ctomfd(qloc,rhoa0,rhob0,u0,v0,
     &           rhoh0,p0,grue0,zeta0,zfa0,zfb0,c20)
c
            cmax  = dmax1(cmax,dabs(u0)+dsqrt(c20),
     &                         dabs(v0)+dsqrt(c20))
            enddo
         enddo
      go to 999
c      
 999  continue
      if (cmax .ne. 0.d0) dt = cfl*dmin1(dx,dy)/cmax
      return
c
      end
