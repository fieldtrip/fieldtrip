C     glmnetMex.F
C
C     Lasso and elastic-net regularized generalized linear models
      
C     [a0,ca,ia,nin,rsq,alm,nlp,jerr] = ...
C        glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,w,ka)
C     [a0,ca,ia,nin,dev,alm,nlp,jerr] = ...
C        glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,nc,maxit,kopt)
C    
C     Extremely efficient procedures for fitting the entire lasso or
C     elastic-net regularization path for linear regression, logistic and
C     multinomial regression models. The algorithm uses cyclical coordinate
C     descent in a pathwise as described in the paper on the maintainer's
C     website.
C    
C     NOTES: This is a MEX-file wrapper of GLMnet.f for MATLAB. Should be called
C     only by glmnet.m. For details about input and output arguments, see
C     GLMnet.f.
C    
C     LICENSE: GPL-2
C    
C     DATE: 13 Jul 2009
C    
C     AUTHORS:
C        Algorithm designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
C        Fortran code written by Jerome Friedman 
C        MATLAB wrapper written and maintained by Hui Jiang, jiangh@stanford.edu 
C        Department of Statistics, Stanford University, Stanford, California, USA.
C    
C     REFERENCES:
C        Friedman, J., Hastie, T. and Tibshirani, R. (2009)
C        Regularization Paths for Generalized Linear Models via Coordinate Descent.
C        To appear, Journal of Statistical Software
C    
C     EXAMPLE:
C        parm = 1.0;
C        x = [1 1; 2 2; 3 3];
C        y = [1 3 2]';
C        jd = 0;
C        vp = [1 1]';
C        ne = 3;
C        nx = 2;
C        nlam = 100;
C        flmin = 0.0001;
C        ulam = 0;
C        thr = 1.0e-4;
C        isd = 1;
C        w = [1 1 1]';
C        ka = 2;
C        [a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,w,ka)
C     
C     DEVELOPMENT: 13 Jul 2009: Original version of glmnetMex.f written.
C
C-----------------------------------------------------------------------

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha
C     64-bit platform

      integer plhs(*), prhs(*)
      integer nlhs, nrhs
C-----------------------------------------------------------------------

C     Input
      real parm,flmin,thr
      integer ka,no,ni,nc,ne,nx,nlam,isd,maxit,kopt
      real, dimension (:), allocatable :: x,y,w,vp,ulam
      integer, dimension (:), allocatable :: ix,jx,jd

C     Output
      integer lmu,nlp,jerr
      real, dimension (:), allocatable :: a0,ca,alm,dev,rsq
      integer, dimension (:), allocatable :: ia,nin

C     Temporary      
      integer temp_pr
      integer temp_m, temp_n
      integer task
      
C     Check for proper number of arguments. 
      if(nrhs .eq. 14 .and. nlhs .eq. 8) then
         task = 1;
      elseif(nrhs .eq. 15 .and. nlhs .eq. 8) then
         task = 2;
      else
         call mexErrMsgTxt('Incorrect number of arguments.')
      endif

C     Get input

      temp_pr = mxGetPr(prhs(1))
      call getreal(temp_pr,parm,1)

      temp_pr = mxGetPr(prhs(2))
      no = mxGetM(prhs(2))
      ni = mxGetN(prhs(2))
      allocate(x(1:no*ni))
      call getreal(temp_pr,x,no*ni)     
   
      temp_pr = mxGetPr(prhs(4))
      temp_m = mxGetM(prhs(4))
      temp_n = mxGetN(prhs(4))
      allocate(jd(temp_m*temp_n))
      call getinteger(temp_pr,jd,temp_m*temp_n)
      
      temp_pr = mxGetPr(prhs(5))
      allocate(vp(1:ni))
      call getreal(temp_pr,vp,ni)

      temp_pr = mxGetPr(prhs(6))
      call getinteger(temp_pr,ne,1)

      temp_pr = mxGetPr(prhs(7))
      call getinteger(temp_pr,nx,1)
      
      temp_pr = mxGetPr(prhs(8))
      call getinteger(temp_pr,nlam,1)

      temp_pr = mxGetPr(prhs(9))
      call getreal(temp_pr,flmin,1)     
 
      temp_pr = mxGetPr(prhs(10))
      temp_m = mxGetM(prhs(10))
      temp_n = mxGetN(prhs(10))
      allocate(ulam(1:temp_m * temp_n))
      call getreal(temp_pr,ulam,temp_m * temp_n)
      
      temp_pr = mxGetPr(prhs(11))
      call getreal(temp_pr,thr,1)
      
      temp_pr = mxGetPr(prhs(12))
      call getinteger(temp_pr,isd,1)

      if (task .eq. 1) then
         temp_pr = mxGetPr(prhs(3))
         allocate(y(1:no))
         call getreal(temp_pr,y,no)

         temp_pr = mxGetPr(prhs(13))
         allocate(w(1:no))
         call getreal(temp_pr,w,no)      
         
         temp_pr = mxGetPr(prhs(14))
         call getinteger(temp_pr,ka,1)
      elseif (task .eq. 2) then
         temp_pr = mxGetPr(prhs(13))
         call getinteger(temp_pr,nc,1)
          
         temp_pr = mxGetPr(prhs(14))
         call getinteger(temp_pr,maxit,1)
          
         temp_pr = mxGetPr(prhs(15))
         call getinteger(temp_pr,kopt,1)        

         temp_pr = mxGetPr(prhs(3))
         allocate(y(1:no*(max(2,nc))))
         call getreal(temp_pr,y,no*(max(2,nc)))
      endif
      
C     Allocate memory for output
      allocate(ia(1:nx))
      call zerointeger(ia,nx)
      allocate(nin(1:nlam))
      call zerointeger(nin,nlam)
      allocate(alm(1:nlam))
      call zeroreal(alm,nlam)
      if (task .eq. 1) then
         allocate(a0(1:nlam))
         call zeroreal(a0,nlam)
         
         allocate(ca(1:nx*nlam))
         call zeroreal(ca,nx*nlam)
         
         allocate(rsq(1:nlam))
         call zeroreal(rsq,nlam)
      elseif (task .eq. 2) then
         allocate(a0(1:nc*nlam))
         call zeroreal(a0,nc*nlam)
         
         allocate(ca(1:nx*nc*nlam))
         call zeroreal(ca,nx*nc*nlam)

         allocate(dev(1:nlam))
         call zeroreal(dev,nlam)
      endif

C     Call glmnet
      lmu = 0
      nlp = 0
      jerr = 0
      if (task .eq. 1) then
          call elnet(ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam
     *,flmin,ulam,thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      elseif (task .eq. 2) then
          call lognet(parm,no,ni,nc,x,y,jd,vp,ne,nx,nlam,flmin
     *,ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      endif

C     Prepare output
      plhs(3) = mxCreateDoubleMatrix(nx,1,0)
      temp_pr = mxGetPr(plhs(3))
      call putinteger(ia,temp_pr,nx)
      
      plhs(4) = mxCreateDoubleMatrix(lmu,1,0)
      temp_pr = mxGetPr(plhs(4))
      call putinteger(nin,temp_pr,lmu)

      plhs(6) = mxCreateDoubleMatrix(lmu,1,0)
      temp_pr = mxGetPr(plhs(6))
      call putreal(alm,temp_pr,lmu)
      
      plhs(7) = mxCreateDoubleMatrix(1,1,0)
      temp_pr = mxGetPr(plhs(7))
      call putinteger(nlp,temp_pr,1)

      plhs(8) = mxCreateDoubleMatrix(1,1,0)
      temp_pr = mxGetPr(plhs(8))
      call putinteger(jerr,temp_pr,1)

      if (task .eq. 1) then
         plhs(1) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putreal(a0,temp_pr,lmu)   

         plhs(2) = mxCreateDoubleMatrix(nx,lmu,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(ca,temp_pr,nx*lmu)  
         
         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putreal(rsq,temp_pr,lmu)               
      elseif (task .eq. 2) then         
         plhs(1) = mxCreateDoubleMatrix(nc,lmu,0)
         temp_pr = mxGetPr(plhs(1))
         call putreal(a0,temp_pr,nc*lmu)

         plhs(2) = mxCreateDoubleMatrix(nx*nc,lmu,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(ca,temp_pr,nx*nc*lmu)

         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putreal(dev,temp_pr,lmu)
      endif

C     Deallocate memory      
      deallocate(x)
      deallocate(y)
      deallocate(jd)
      deallocate(vp)
      deallocate(ulam)
      deallocate(a0)
      deallocate(ca)
      deallocate(ia)
      deallocate(nin)
      deallocate(alm)
C     For logistic elastic net         
      if (task .eq. 1) then
          deallocate(w)
          deallocate(rsq)
      elseif (task .eq. 2) then   
          deallocate(dev)
      endif
      
      return
      end

C     End of subroutine mexFunction
      
      subroutine real8toreal(x, y, size)
      integer size
      real*8 x(size)
      real y(size)
      do 10 i=1,size
         y(i)= x(i)
 10   continue
      return
      end

      subroutine realtoreal8(x, y, size)
      integer size
      real x(size)
      real*8 y(size)
      do 20 i=1,size
         y(i)= x(i)
 20   continue
      return
      end
 
      subroutine real8tointeger(x, y, size)
      integer size
      real*8 x(size)
      integer y(size)
      do 30 i=1,size
         y(i)= x(i)
 30   continue
      return
      end
 
      subroutine integertoreal8(x, y, size)
      integer size
      integer x(size)
      real*8 y(size)
      do 40 i=1,size
         y(i)= x(i)
 40   continue
      return
      end
      
      subroutine getreal(pr,x,size)
      integer pr,size
      real x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call mxCopyPtrToReal8(pr,temp,size)
      call real8toreal(temp,x,size)
      deallocate(temp)      
      return
      end
      
      subroutine getinteger(pr,x,size)
      integer pr,size
      integer x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call mxCopyPtrToReal8(pr,temp,size)
      call real8tointeger(temp,x,size)
      deallocate(temp)      
      return
      end      
      
      subroutine putreal(x,pr,size)
      integer pr,size
      real x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call realtoreal8(x,temp,size)
      call mxCopyReal8ToPtr(temp,pr,size)
      deallocate(temp)      
      return
      end
      
      subroutine putinteger(x,pr,size)
      integer pr,size
      integer x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call integertoreal8(x,temp,size)
      call mxCopyReal8ToPtr(temp,pr,size)
      deallocate(temp)      
      return
      end            
      
      subroutine zeroreal(x,size)
      integer size
      real x(size)
      do 90 i=1,size
         x(i) = 0
 90   continue     
      return 
      end
      
      subroutine zerointeger(x,size)
      integer size
      integer x(size)
      do 100 i=1,size
         x(i) = 0
 100  continue
      return
      end