c@**********************************************************************
c                                                                      *
c     function dpmpar                                                  *
c                                                                      *
c     this function provides double precision machine parameters       *
c     when the appropriate set of data statements is activated (by     *
c     removing the c from column 1) and all other data statements are  *
c     rendered inactive. most of the parameter values were obtained    *
c     from the corresponding bell laboratories port library function.  *
c                                                                      *
c     the function statement is                                        *
c                                                                      *
c       function dpmpar (i)                                            *
c                                                                      *
c     where                                                            *
c                                                                      *
c       i is an integer input variable set to 1, 2, or 3 which         *
c         selects the desired machine parameter. if the machine has    *
c         t base b digits and its smallest and largest exponents are   *
c         emin and emax, respectively, then these parameters are       *
c                                                                      *
c         dpmpar(1) = b**(1 - t), the machine precision,               *
c                                                                      *
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,           *
c                                                                      *
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.    *
c                                                                      *
c     argonne national laboratory. minpack project. june 1983.         *
c     burton s. garbow, kenneth e. hillstrom, jorge j. more            *
c                                                                      *
c***********************************************************************
      function dpmpar(i)
c
      integer*4 mcheps(2),minmag(2),maxmag(2),i
      real*8    dmach(3),dpmpar
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
c
c.....machine constants for the ibm 360/370 series,
c     the amdahl 470/v6, the icl 2900, the itel as/6,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
c     data mcheps(1),mcheps(2) / z34100000, z00000000 /
c     data minmag(1),minmag(2) / z00100000, z00000000 /
c     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
c
c.....machine constants for the prime 400.
c
c     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
c     data minmag(1),minmag(2) / :10000000000, :00000100000 /
c     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
c
c     machine constants for the data general mv/serie
c
c     data minmag /    4000000k,17777777777k/
c     data maxmag /17777777777k,37777777777k/
c     data mcheps / 6404000000k,    4000000k/
c
c.....mips unix
c
      data dmach(1) /1.0661664009683891d-16/
      data dmach(2) /2.2250738585072012d-300/
      data dmach(3) /1.7976931348623158d+300/
c
      dpmpar = dmach(i)
      return
      end
