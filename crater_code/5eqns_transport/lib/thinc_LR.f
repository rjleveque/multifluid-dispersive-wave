c
c ----------------------------------------------------------------
c
      subroutine thinc_LR(zL,zR,zm,z0,zp,beta,zeps)
      implicit double precision (a-h,o-z)
c
c     # THINC-M interface-sharpening reconstruction
c     # Ref: Cassidy, Edwards, Tian, JCP 228 (2009) 5628-5649
c     # return: zL and zR (volume fractions to the left and right at 
c     #                    the cell edge)
c
      zmin = dmin1(zm,zp)
      zmax = dmax1(zm,zp)
      tsgn = dsign(1.d0,zp-zm)
      BB   = dexp(tsgn*beta*(2.d0*(z0-zmin+zeps)/
     &           (zmax-zmin+zeps)-1.d0))
      CC   = (BB/dcosh(beta)-1.d0)/dtanh(beta)
      DD   = (dtanh(beta)+CC)/(1.d0+CC*dtanh(beta))
      zL   = zmin+0.5d0*(zmax-zmin)*(1.d0+tsgn*CC)
      zR   = zmin+0.5d0*(zmax-zmin)*(1.d0+tsgn*DD)
      return
c
      end
