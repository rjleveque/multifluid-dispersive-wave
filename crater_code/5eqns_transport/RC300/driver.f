c
c ----------------------------------------------------------------
c
      program driver
      implicit double precision (a-h,o-z)
      parameter(maxmx = 5000)
      parameter(maxmy = 800)
      parameter(mwork = 33069636)
      parameter(mrkwork = 24209016)
      parameter(mbc = 3)
      parameter(meqn = 6)
      parameter(mwaves = 3)
      parameter(maux = 2)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension mthlim(mwaves)
      dimension work(mwork)
      dimension rkwork(mwork)
      real*4 cpusec,cputm(2)
c
      call cpu_time(cpustart)
c
      call tclaw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,
     &     mthlim,q,work,mrkwork,rkwork,aux)
c
      call cpu_time(cpuend)
      cpusec = cpuend-cpustart
c
      write(23,*) cpusec
c
      write(10,*)
     &    'Program has been running for', cpusec, 'seconds.'
      stop
c
      end
