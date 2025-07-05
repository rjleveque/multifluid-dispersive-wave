c
c --------------------------------------------------------------
c
      subroutine copyaux2(maxmx,maxmy,maux,mbc,mx,my,aux1,aux2)
      implicit double precision (a-h,o-z)
      dimension aux1(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension aux2(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
c
c     # copy the contents of aux1 into aux2
c
      do ma=1,maux
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               aux2(i,j,ma) = aux1(i,j,ma)
               enddo
            enddo
         enddo
      return
c
      end
