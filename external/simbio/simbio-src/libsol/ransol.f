c=======================================================================
c
      function ransol(idum)
c
c-----------------------------------------------------------------------
c---random number generator nach donald e. knuth
c---aus: numerical recipes,b.p.flannery,w.h.press/cambridge university 
c---press
c-----------------------------------------------------------------------
      save ma,inext,inextp
     
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
      dimension ma(55)
      data iff /0/
      real ransol = 0.0
      real fac = 0.0
c      save ma,inext,inextp
c
      if (idum.lt.0.or.iff.eq.0) then
         iff=1
         mj=mseed-iabs(idum)
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 10 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if (mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 10      continue
         do 30 k=1,4
            do 20 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 20         continue
 30      continue
         inext=0
         inextp=31
         idum=1
      end if
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      ransol=mj*fac
      return
      end
c=======================================================================
