c
c-----------------------------------------------------------
c      
      subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,
     &           ylower,dx,dy,q,maux,aux)
      implicit double precision (a-h,o-z)
      common /loops/ xloop(10),yloop(10),zloop(10),rloop(10),
     &        axloop(10),bxloop(10),inbloop(2,2,10)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      common /phypar/ den0,velx0,vely0,velz0,pr0
      common /cgrav/ grav
      common /crater_info/ RC,RD,DC,r1,rdecay,linfactor
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension q1(20),q2(20)
      dimension qloc(20)
c
c     # astroid impact problem         
c
      p0 = 1.01325d5
c
c     # phase 1
c     # gas phase
      rho1 = rref(1)
      vx1  = velx0 
      vy1  = vely0 
      p1   = p0     
c
c     # phase 2
c     # liquid phase                    
      rho2 = rref(2)
      vx2  = velx0
      vy2  = vely0
      p2   = p0    
c
      call prmtoc(q1,rho1,rho2,vx1,vy1,p1,p2,
     &     1.0d0-zf(1),zf(1))
      call prmtoc(q2,rho1,rho2,vx2,vy2,p1,p2,
     &     zf(2),1.0d0-zf(2))

      write(66,*) 'q1=',(q1(m),m=1,meqn)
      write(66,*) 'q2=',(q2(m),m=1,meqn)
c
c     # initialize air-water interface
      call fnpar(1)
c
      do i=1-mbc,mx+mbc
         xc = xlower+dble(i-1)*dx+0.5d0*dx
         do j=1-mbc,my+mbc
            yc = ylower+dble(j-1)*dy+0.5d0*dy
c
c            p1 = p0-rho1*grav*(yc-yloop(1))
c
            p1 = p0
            p2 = p0-rho2*grav*(yc-yloop(1))
c
            if (fdisc(xc,yc) .lt. 0.0d0) then
c               # liquid phase
c
                call prmtoc(qloc,rho1,rho2,vx2,vy2,p1,p2,
     &               zf(2),1.0d0-zf(2))
            else
c               # gas phase
c
                do m=1,meqn
                   qloc(m) = q1(m)
                   enddo

c                call prmtoc(qloc,rho1,rho2,vx2,vy2,p1,p2,
c     &               1.0d0-zf(1),zf(1))
            endif
c
            do m=1,meqn
               q(i,j,m) = qloc(m)
               enddo
            enddo
         enddo
c
c     # initialize crater    
      call fnpar(2)
c
      do j=1-mbc,my+mbc
         ylow = ylower+dble(j-1)*dy
         do i=1-mbc,mx+mbc
            xlow = xlower+dble(i-1)*dx
            call cellave(xlow,ylow,dx,dy,wlin)
c
            if (wlin .ne. 0.d0) then
                do m=1,meqn
                   q(i,j,m) = wlin*q1(m)+(1.d0-wlin)*q2(m)
                   enddo
            endif
            enddo
         enddo
c
c     # cavity: hydrostatic pressure readjustment 
      do i=1-mbc,mx+mbc
         xlow = xlower+dble(i-1)*dx
         if (xlow+0.5d0*dx .lt. RC) then
             do j=my+mbc,1-mbc,-1
                ylow = ylower+dble(j-1)*dy
                call cellave(xlow,ylow,dx,dy,wlin)
c
                if (wlin .ne. 1.0d0 .and.
     &              wlin .ne. 0.0d0) then
c                   # find the interface cell

                    jk = j
                    yk = ylower+dble(jk-1)*dy+0.5d0*dy
c
                    p1 = p0 
c
                    do jj=jk-1,1-mbc,-1
                       yc = ylower+dble(jj-1)*dy+0.5d0*dy
c
                       p2 = p0-rho2*grav*(yc-yk)
c
                       call prmtoc(qloc,rho1,rho2,vx2,vy2,p1,p2,
     &                   zf(2),1.0d0-zf(2))
c
                       do m=1,meqn
                          q(i,jj,m) = qloc(m)
                          enddo
                       enddo
                endif
                enddo
         endif
         enddo
      return
c
      end
c
c -------------------------------------------------
c
      double precision function fdisc(r,z)
      implicit double precision (a-h,o-z)
      common /cominit3/ x0,y0,z0,r0,jfront
      common /crater_info/ RC,RD,DC,r1,rdecay,linfactor
c
      if (jfront .eq. 1) then
c         # flat air-water interface
c
          fdisc = z-z0
      else
c         # crater
c
          eta_crater = -dsqrt(RC**2-r**2)
c
c         # no lip if RD==RC
          if (r > RD) then
              eta_crater = 0.0d0
          endif
c
          fdisc = eta_crater-z
      endif      
      return
c                       
      end   
c
c ------------------------------------------------------------
c
      subroutine fnpar(iloops)
      implicit double precision (a-h,o-z)
      common /loops/ xloop(10),yloop(10),zloop(10),rloop(10),
     &        axloop(10),bxloop(10),inbloop(2,2,10)
      common /cominit3/ x0,y0,z0,r0,jfront
c
      x0 = xloop(iloops)
      y0 = yloop(iloops)
      z0 = zloop(iloops)
      r0 = rloop(iloops)

      jfront = iloops
      return
c
      end
