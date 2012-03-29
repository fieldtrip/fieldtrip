c@**********************************************************************
c                                                                      *
c     function spmpar                                                  *
c                                                                      *
c     this function provides single precision machine parameters       *
c     when the appropriate set of data statements is activated (by     *
c     removing the c from column 1) and all other data statements are  *
c     rendered inactive. most of the parameter values were obtained    *
c     from the corresponding bell laboratories port library function.  *
c                                                                      *
c     the function statement is                                        *
c                                                                      *
c       function spmpar (i)                                            *
c                                                                      *
c     where                                                            *
c                                                                      *
c       i is an integer input variable set to 1, 2, or 3 which         *
c         selects the desired machine parameter. if the machine has    *
c         t base b digits and its smallest and largest exponents are   *
c         emin and emax, respectively, then these parameters are       *
c                                                                      *
c         spmpar(1) = b**(1 - t), the machine precision,               *
c                                                                      *
c         spmpar(2) = b**(emin - 1), the smallest magnitude,           *
c                                                                      *
c         spmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.    *
c                                                                      *
c     argonne national laboratory. minpack project. june 1983.         *
c     burton s. garbow, kenneth e. hillstrom, jorge j. more            *
c                                                                      *
c***********************************************************************
      function spmpar(i)
c
      integer*4 mcheps(1),minmag(1),maxmag(1),i
      real*4    smach(3),spmpar
      equivalence (smach(1),mcheps(1))
      equivalence (smach(2),minmag(1))
      equivalence (smach(3),maxmag(1))
c
c.....machine constants for the data general mv/serie
c
c     data minmag /    4000000k/
c     data maxmag /17777777777k/
c     data mcheps / 7404000000k/
c
c.....mips unix
c
      data smach(1) /5.75664494e-08/
      data smach(2) /3.40282356e+38/
      data smach(3) /1.17549429e-38/
c
      spmpar = smach(i)
      return
      end
