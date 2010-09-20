c
c                          GLMnet (5/17/08)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c            lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c             isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse column format
c
c
c other inputs:
c
c   ka = algorithm flag
c      ka=1 => covariance updating algorithm
c      ka=2 => naive algorithm
c   parm = family member index (0 <= parm <= 1)
c        = 0.0 => ridge
c        = 1.0 => lasso
c   no = number of observations
c   ni = number of predictor variables
c   y(no) = response vector
c   w(no)= observation weights
c   jd(jd(1)+1) = predictor variable deletion flag
c      jd(1) = 0  => use all variables
c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
c   vp(ni) = relative penalties for each predictor variable
c      vp(j) = 0 => jth variable unpenalized
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   nlam = (maximum) number of lamda values
c   flmin = user control of lamda values (>=0)
c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
c      flmin >= 1.0 => use supplied lamda values (see below)
c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
c   thr = convergence threshold for each lamda solution.
c      iterations stop when the maximum standardized coefficient
c      change from the previous iteration is less than thr
c      (suggested value, thr=1.0e-4)
c   isd = standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
c
c output:
c
c   lmu = actual number of lamda values (solutions)
c   a0(lmu) = intercept values for each solution
c   ca(nx,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(lmu) = number of compressed coefficients for each solution
c   rsq(lmu) = R**2 values for each solution
c   alm(lmu) = lamda values corresponding to each solution
c   nlp = total passes over the data summed over all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr != 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
c
c
c Note: x, y and w are overwritten by programs
c
c
c least-squares utility routines:
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(n) = model predictions
c
c
c
c
c          Symmetric binomial/multinomial logistic elastic net
c
c
c dense predictor matrix:
c
c call lognet (parm,no,ni,nc,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c              maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,jd,vp,ne,nx,nlam,flmin,
c             ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse column format
c
c
c other inputs:
c
c   parm, no, ni, jd, vp, ne, nx, nlam, flmin, ulam, thr, isd, same as above.
c
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point(overwritten)
c   maxit = maximum number of iterations allowed for any lamda value
c           (suggested value, maxit = 100)
c   kopt = optimization flag
c      kopt = 0 => Newton-Raphson
c      kpot = 1 => modified Newton-Raphson (recommended)
c
c
c output:
c
c   lmu, ia, nin, alm, nlp, same as above
c
c   a0(nc,lmu) = intercept values for each class at each solution
c   ca(nx,nc,lmu) = compressed coefficient values for each class at
c                each solution
c   dev(lmu) = fraction of explained devience for each solution
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8000 + k => null probability < 1.0e-5 for class k
c         jerr = 9000 + k => null probability for class k
c                            > 1.0 - 1.0e-5
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output returned
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations. Solutions for
c            larger lamdas returned
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value. Solutions for
c            larger lamdas returned
c
c
c
c logistic utilitity routines:
c
c uncompress coefficient vector for particular solution:
c
c call luncomp(ni,nx,nc,ca,ia,nin,a)
c
c input:
c
c    ni, nx, nc = same as above
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c    a(ni,nc) =  uncompressed coefficient vectors
c                 referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor vectors:
c
c call lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans);
c
c input:
c
c    nt = number of observations
c    x(nt,ni) = full (uncompressed) predictor vectors
c    nc, nx = same as above
c    a0(nc) = intercepts
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c ans(nc,nt) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nc, nx = same as above
c    a0(nc) = intercept
c    ca(nx,nc) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nc,n) = model predictions
c
c
c
c              
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    352 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          353
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          354
      integer jd(*),ia(nx),nin(nlam)                                        355
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     358
      jerr=10000                                                            358
      return                                                                358
10021 continue                                                              359
      allocate(vq(1:ni),stat=jerr)                                          359
      if(jerr.ne.0) return                                                  360
      vq=max(0.0,vp)                                                        360
      vq=vq*ni/sum(vq)                                                      361
      if(ka .ne. 1)goto 10041                                               362
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    365 
     *,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            366
10041 continue                                                              367
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    370 
     *  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              371
10031 continue                                                              371
      deallocate(vq)                                                        372
      return                                                                373
      end                                                                   374
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    377 
     *hr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           378
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         379
      integer jd(*),ia(nx),nin(nlam)                                        380
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           385
      allocate(xm(1:ni),stat=ierr)                                          385
      jerr=jerr+ierr                                                        386
      allocate(xs(1:ni),stat=ierr)                                          386
      jerr=jerr+ierr                                                        387
      allocate(ju(1:ni),stat=ierr)                                          387
      jerr=jerr+ierr                                                        388
      allocate(xv(1:ni),stat=ierr)                                          388
      jerr=jerr+ierr                                                        389
      allocate(vlam(1:nlam),stat=ierr)                                      389
      jerr=jerr+ierr                                                        390
      if(jerr.ne.0) return                                                  391
      call chkvars(no,ni,x,ju)                                              392
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  393
      if(maxval(ju) .gt. 0)goto 10071                                       393
      jerr=7777                                                             393
      return                                                                393
10071 continue                                                              394
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               395
      if(jerr.ne.0) return                                                  396
      if(flmin.ge.1.0) vlam=ulam/ys                                         397
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm    399 
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  400
10080 do 10081 k=1,lmu                                                      400
      alm(k)=ys*alm(k)                                                      400
      nk=nin(k)                                                             401
10090 do 10091 l=1,nk                                                       401
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          401
10091 continue                                                              402
10092 continue                                                              402
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         403
10081 continue                                                              404
10082 continue                                                              404
      deallocate(xm,xs,g,ju,xv,vlam)                                        405
      return                                                                406
      end                                                                   407
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        408
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  408
      integer ju(ni)                                                        409
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           412
      if(jerr.ne.0) return                                                  413
      w=w/sum(w)                                                            413
      v=sqrt(w)                                                             414
10100 do 10101 j=1,ni                                                       414
      if(ju(j).eq.0)goto 10101                                              415
      xm(j)=dot_product(w,x(:,j))                                           415
      x(:,j)=v*(x(:,j)-xm(j))                                               416
      xv(j)=dot_product(x(:,j),x(:,j))                                      417
10101 continue                                                              418
10102 continue                                                              418
      if(isd .ne. 0)goto 10121                                              418
      xs=1.0                                                                418
      goto 10131                                                            419
10121 continue                                                              420
10140 do 10141 j=1,ni                                                       420
      if(ju(j).eq.0)goto 10141                                              420
      xs(j)=sqrt(xv(j))                                                     420
      x(:,j)=x(:,j)/xs(j)                                                   420
10141 continue                                                              421
10142 continue                                                              421
      xv=1.0                                                                422
10131 continue                                                              423
10111 continue                                                              423
      ym=dot_product(w,y)                                                   423
      y=v*(y-ym)                                                            423
      ys=sqrt(dot_product(y,y))                                             423
      y=y/ys                                                                423
      g=0.0                                                                 424
10150 do 10151 j=1,ni                                                       424
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             424
10151 continue                                                              425
10152 continue                                                              425
      deallocate(v)                                                         426
      return                                                                427
      end                                                                   428
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    430 
     *xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    431 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    432 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       433
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           439
      jerr=jerr+ierr                                                        440
      allocate(mm(1:ni),stat=ierr)                                          440
      jerr=jerr+ierr                                                        441
      allocate(da(1:ni),stat=ierr)                                          441
      jerr=jerr+ierr                                                        442
      if(jerr.ne.0) return                                                  443
      bta=max(beta,1.0e-3)                                                  443
      omb=1.0-bta                                                           444
      if(flmin .ge. 1.0)goto 10171                                          444
      eqs=max(eps,flmin)                                                    444
      alf=eqs**(1.0/(nlam-1))                                               444
10171 continue                                                              445
      rsq=0.0                                                               445
      a=0.0                                                                 445
      mm=0                                                                  445
      nlp=0                                                                 445
      nin=nlp                                                               445
      iz=0                                                                  445
      mnl=min(mnlam,nlam)                                                   446
10180 do 10181 m=1,nlam                                                     447
      if(flmin .lt. 1.0)goto 10201                                          447
      alm=ulam(m)                                                           447
      goto 10191                                                            448
10201 if(m .le. 2)goto 10211                                                448
      alm=alm*alf                                                           448
      goto 10191                                                            449
10211 if(m .ne. 1)goto 10221                                                449
      alm=big                                                               449
      goto 10231                                                            450
10221 continue                                                              450
      alm=0.0                                                               451
10240 do 10241 j=1,ni                                                       451
      if(ju(j).eq.0)goto 10241                                              451
      if(vp(j).le.0.0)goto 10241                                            452
      alm=max(alm,abs(g(j))/vp(j))                                          453
10241 continue                                                              454
10242 continue                                                              454
      alm=alf*alm/bta                                                       455
10231 continue                                                              456
10191 continue                                                              456
      dem=alm*omb                                                           456
      ab=alm*bta                                                            456
      rsq0=rsq                                                              456
      jz=1                                                                  457
10250 continue                                                              457
10251 continue                                                              457
      if(iz*jz.ne.0) go to 10260                                            457
      nlp=nlp+1                                                             457
      dlx=0.0                                                               458
10270 do 10271 k=1,ni                                                       458
      if(ju(k).eq.0)goto 10271                                              459
      ak=a(k)                                                               459
      u=g(k)+ak*xv(k)                                                       459
      v=abs(u)-vp(k)*ab                                                     459
      a(k)=0.0                                                              460
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         461
      if(a(k).eq.ak)goto 10271                                              462
      if(mm(k) .ne. 0)goto 10291                                            462
      nin=nin+1                                                             462
      if(nin.gt.nx)goto 10272                                               463
10300 do 10301 j=1,ni                                                       463
      if(ju(j).eq.0)goto 10301                                              464
      if(mm(j) .eq. 0)goto 10321                                            464
      c(j,nin)=c(k,mm(j))                                                   464
      goto 10301                                                            464
10321 continue                                                              465
      if(j .ne. k)goto 10341                                                465
      c(j,nin)=xv(j)                                                        465
      goto 10301                                                            465
10341 continue                                                              466
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   467
10301 continue                                                              468
10302 continue                                                              468
      mm(k)=nin                                                             468
      ia(nin)=k                                                             469
10291 continue                                                              470
      del=a(k)-ak                                                           470
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      471
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     472
10350 do 10351 j=1,ni                                                       472
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               472
10351 continue                                                              473
10352 continue                                                              473
10271 continue                                                              474
10272 continue                                                              474
      if(dlx.lt.thr)goto 10252                                              474
      if(nin.gt.nx)goto 10252                                               475
10260 continue                                                              475
      iz=1                                                                  475
      da(1:nin)=a(ia(1:nin))                                                476
10360 continue                                                              476
10361 continue                                                              476
      nlp=nlp+1                                                             476
      dlx=0.0                                                               477
10370 do 10371 l=1,nin                                                      477
      k=ia(l)                                                               477
      ak=a(k)                                                               477
      u=g(k)+ak*xv(k)                                                       477
      v=abs(u)-vp(k)*ab                                                     478
      a(k)=0.0                                                              479
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         480
      if(a(k).eq.ak)goto 10371                                              481
      del=a(k)-ak                                                           481
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      482
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     483
10380 do 10381 j=1,nin                                                      483
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  483
10381 continue                                                              484
10382 continue                                                              484
10371 continue                                                              485
10372 continue                                                              485
      if(dlx.lt.thr)goto 10362                                              485
      goto 10361                                                            486
10362 continue                                                              486
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      487
10390 do 10391 j=1,ni                                                       487
      if(mm(j).ne.0)goto 10391                                              488
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            489
10391 continue                                                              490
10392 continue                                                              490
      jz=0                                                                  491
      goto 10251                                                            492
10252 continue                                                              492
      if(nin.gt.nx)goto 10182                                               493
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 493
      kin(m)=nin                                                            494
      rsqo(m)=rsq                                                           494
      almo(m)=alm                                                           494
      lmu=m                                                                 495
      if(m.lt.mnl)goto 10181                                                495
      if(flmin.ge.1.0)goto 10181                                            496
      me=0                                                                  496
10400 do 10401 j=1,nin                                                      496
      if(ao(j,m).ne.0.0) me=me+1                                            496
10401 continue                                                              496
10402 continue                                                              496
      if(me.gt.ne)goto 10182                                                497
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     497
      if(rsq.gt.rsqmax)goto 10182                                           498
10181 continue                                                              499
10182 continue                                                              499
      deallocate(a,mm,c,da)                                                 500
      return                                                                501
      end                                                                   502
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    504 
     *r,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           505
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         506
      integer jd(*),ia(nx),nin(nlam)                                        507
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          512
      allocate(xs(1:ni),stat=ierr)                                          512
      jerr=jerr+ierr                                                        513
      allocate(ju(1:ni),stat=ierr)                                          513
      jerr=jerr+ierr                                                        514
      allocate(xv(1:ni),stat=ierr)                                          514
      jerr=jerr+ierr                                                        515
      allocate(vlam(1:nlam),stat=ierr)                                      515
      jerr=jerr+ierr                                                        516
      if(jerr.ne.0) return                                                  517
      call chkvars(no,ni,x,ju)                                              518
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  519
      if(maxval(ju) .gt. 0)goto 10421                                       519
      jerr=7777                                                             519
      return                                                                519
10421 continue                                                              520
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                521
      if(jerr.ne.0) return                                                  522
      if(flmin.ge.1.0) vlam=ulam/ys                                         523
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm    525 
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  526
10430 do 10431 k=1,lmu                                                      526
      alm(k)=ys*alm(k)                                                      526
      nk=nin(k)                                                             527
10440 do 10441 l=1,nk                                                       527
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          527
10441 continue                                                              528
10442 continue                                                              528
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         529
10431 continue                                                              530
10432 continue                                                              530
      deallocate(xm,xs,ju,xv,vlam)                                          531
      return                                                                532
      end                                                                   533
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         534
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        534
      integer ju(ni)                                                        535
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           538
      if(jerr.ne.0) return                                                  539
      w=w/sum(w)                                                            539
      v=sqrt(w)                                                             540
10450 do 10451 j=1,ni                                                       540
      if(ju(j).eq.0)goto 10451                                              541
      xm(j)=dot_product(w,x(:,j))                                           541
      x(:,j)=v*(x(:,j)-xm(j))                                               542
      xv(j)=dot_product(x(:,j),x(:,j))                                      543
10451 continue                                                              544
10452 continue                                                              544
      if(isd .ne. 0)goto 10471                                              544
      xs=1.0                                                                544
      goto 10481                                                            545
10471 continue                                                              545
10490 do 10491 j=1,ni                                                       545
      if(ju(j).eq.0)goto 10491                                              545
      xs(j)=sqrt(xv(j))                                                     545
      x(:,j)=x(:,j)/xs(j)                                                   545
10491 continue                                                              546
10492 continue                                                              546
      xv=1.0                                                                547
10481 continue                                                              548
10461 continue                                                              548
      ym=dot_product(w,y)                                                   548
      y=v*(y-ym)                                                            548
      ys=sqrt(dot_product(y,y))                                             548
      y=y/ys                                                                549
      deallocate(v)                                                         550
      return                                                                551
      end                                                                   552
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,x    554 
     *v,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    555 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    556 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       557
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           562
      allocate(mm(1:ni),stat=ierr)                                          562
      jerr=jerr+ierr                                                        563
      if(jerr.ne.0) return                                                  564
      bta=max(beta,1.0e-3)                                                  564
      omb=1.0-bta                                                           565
      if(flmin .ge. 1.0)goto 10511                                          565
      eqs=max(eps,flmin)                                                    565
      alf=eqs**(1.0/(nlam-1))                                               565
10511 continue                                                              566
      rsq=0.0                                                               566
      a=0.0                                                                 566
      mm=0                                                                  566
      nlp=0                                                                 566
      nin=nlp                                                               566
      iz=0                                                                  566
      mnl=min(mnlam,nlam)                                                   567
10520 do 10521 m=1,nlam                                                     568
      if(flmin .lt. 1.0)goto 10541                                          568
      alm=ulam(m)                                                           568
      goto 10531                                                            569
10541 if(m .le. 2)goto 10551                                                569
      alm=alm*alf                                                           569
      goto 10531                                                            570
10551 if(m .ne. 1)goto 10561                                                570
      alm=big                                                               570
      goto 10571                                                            571
10561 continue                                                              571
      alm=0.0                                                               572
10580 do 10581 j=1,ni                                                       572
      if(ju(j).eq.0)goto 10581                                              572
      if(vp(j).le.0.0)goto 10581                                            573
      alm=max(alm,abs(dot_product(y,x(:,j)))/vp(j))                         574
10581 continue                                                              575
10582 continue                                                              575
      alm=alf*alm/bta                                                       576
10571 continue                                                              577
10531 continue                                                              577
      dem=alm*omb                                                           577
      ab=alm*bta                                                            577
      rsq0=rsq                                                              577
      jz=1                                                                  578
10590 continue                                                              578
10591 continue                                                              578
      if(iz*jz.ne.0) go to 10260                                            578
      nlp=nlp+1                                                             578
      dlx=0.0                                                               579
10600 do 10601 k=1,ni                                                       579
      if(ju(k).eq.0)goto 10601                                              579
      gk=dot_product(y,x(:,k))                                              580
      ak=a(k)                                                               580
      u=gk+ak*xv(k)                                                         580
      v=abs(u)-vp(k)*ab                                                     580
      a(k)=0.0                                                              581
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         582
      if(a(k).eq.ak)goto 10601                                              583
      if(mm(k) .ne. 0)goto 10621                                            583
      nin=nin+1                                                             583
      if(nin.gt.nx)goto 10602                                               584
      mm(k)=nin                                                             584
      ia(nin)=k                                                             585
10621 continue                                                              586
      del=a(k)-ak                                                           586
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        587
      y=y-del*x(:,k)                                                        587
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     588
10601 continue                                                              589
10602 continue                                                              589
      if(dlx.lt.thr)goto 10592                                              589
      if(nin.gt.nx)goto 10592                                               590
10260 continue                                                              590
      iz=1                                                                  591
10630 continue                                                              591
10631 continue                                                              591
      nlp=nlp+1                                                             591
      dlx=0.0                                                               592
10640 do 10641 l=1,nin                                                      592
      k=ia(l)                                                               592
      gk=dot_product(y,x(:,k))                                              593
      ak=a(k)                                                               593
      u=gk+ak*xv(k)                                                         593
      v=abs(u)-vp(k)*ab                                                     593
      a(k)=0.0                                                              594
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         595
      if(a(k).eq.ak)goto 10641                                              596
      del=a(k)-ak                                                           596
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        597
      y=y-del*x(:,k)                                                        597
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     598
10641 continue                                                              599
10642 continue                                                              599
      if(dlx.lt.thr)goto 10632                                              599
      goto 10631                                                            600
10632 continue                                                              600
      jz=0                                                                  601
      goto 10591                                                            602
10592 continue                                                              602
      if(nin.gt.nx)goto 10522                                               603
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 603
      kin(m)=nin                                                            604
      rsqo(m)=rsq                                                           604
      almo(m)=alm                                                           604
      lmu=m                                                                 605
      if(m.lt.mnl)goto 10521                                                605
      if(flmin.ge.1.0)goto 10521                                            606
      me=0                                                                  606
10650 do 10651 j=1,nin                                                      606
      if(ao(j,m).ne.0.0) me=me+1                                            606
10651 continue                                                              606
10652 continue                                                              606
      if(me.gt.ne)goto 10522                                                607
      if(rsq-rsq0.lt.sml*rsq)goto 10522                                     607
      if(rsq.gt.rsqmax)goto 10522                                           608
10521 continue                                                              609
10522 continue                                                              609
      deallocate(a,mm)                                                      610
      return                                                                611
      end                                                                   612
      subroutine chkvars(no,ni,x,ju)                                        613
      real x(no,ni)                                                         613
      integer ju(ni)                                                        614
10660 do 10661 j=1,ni                                                       614
      ju(j)=0                                                               614
      t=x(1,j)                                                              615
10670 do 10671 i=2,no                                                       615
      if(x(i,j).eq.t)goto 10671                                             615
      ju(j)=1                                                               615
      goto 10672                                                            615
10671 continue                                                              616
10672 continue                                                              616
10661 continue                                                              617
10662 continue                                                              617
      return                                                                618
      end                                                                   619
      subroutine uncomp(ni,ca,ia,nin,a)                                     620
      real ca(*),a(ni)                                                      620
      integer ia(*)                                                         621
      a=0.0                                                                 621
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   622
      return                                                                623
      end                                                                   624
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 625
      real ca(nin),x(n,*),f(n)                                              625
      integer ia(nin)                                                       626
      f=a0                                                                  626
      if(nin.le.0) return                                                   627
10680 do 10681 i=1,n                                                        627
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       627
10681 continue                                                              628
10682 continue                                                              628
      return                                                                629
      end                                                                   630
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl    633 
     *min,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               634
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         635
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            636
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10701                                     639
      jerr=10000                                                            639
      return                                                                639
10701 continue                                                              640
      allocate(vq(1:ni),stat=jerr)                                          640
      if(jerr.ne.0) return                                                  641
      vq=max(0.0,vp)                                                        641
      vq=vq*ni/sum(vq)                                                      642
      if(ka .ne. 1)goto 10721                                               643
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam    646 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10731                                                            647
10721 continue                                                              648
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,    651 
     *thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10731 continue                                                              652
10711 continue                                                              652
      deallocate(vq)                                                        653
      return                                                                654
      end                                                                   655
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi    658 
     *n,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               659
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         660
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            661
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           666
      allocate(xm(1:ni),stat=ierr)                                          666
      jerr=jerr+ierr                                                        667
      allocate(xs(1:ni),stat=ierr)                                          667
      jerr=jerr+ierr                                                        668
      allocate(ju(1:ni),stat=ierr)                                          668
      jerr=jerr+ierr                                                        669
      allocate(xv(1:ni),stat=ierr)                                          669
      jerr=jerr+ierr                                                        670
      allocate(vlam(1:nlam),stat=ierr)                                      670
      jerr=jerr+ierr                                                        671
      if(jerr.ne.0) return                                                  672
      call spchkvars(no,ni,x,ix,ju)                                         673
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  674
      if(maxval(ju) .gt. 0)goto 10751                                       674
      jerr=7777                                                             674
      return                                                                674
10751 continue                                                              675
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)       676
      if(jerr.ne.0) return                                                  677
      if(flmin.ge.1.0) vlam=ulam/ys                                         678
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t    680 
     *hr,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  681
10760 do 10761 k=1,lmu                                                      681
      alm(k)=ys*alm(k)                                                      681
      nk=nin(k)                                                             682
10770 do 10771 l=1,nk                                                       682
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          682
10771 continue                                                              683
10772 continue                                                              683
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         684
10761 continue                                                              685
10762 continue                                                              685
      deallocate(xm,xs,g,ju,xv,vlam)                                        686
      return                                                                687
      end                                                                   688
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j    689 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                      689
      integer ix(*),jx(*),ju(ni)                                            690
      w=w/sum(w)                                                            691
10780 do 10781 j=1,ni                                                       691
      if(ju(j).eq.0)goto 10781                                              692
      jb=ix(j)                                                              692
      je=ix(j+1)-1                                                          692
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                              693
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                  694
10781 continue                                                              695
10782 continue                                                              695
      if(isd .ne. 0)goto 10801                                              695
      xs=1.0                                                                695
      goto 10811                                                            696
10801 continue                                                              696
10820 do 10821 j=1,ni                                                       696
      if(ju(j).ne.0) xs(j)=sqrt(xv(j))                                      696
10821 continue                                                              696
10822 continue                                                              696
      xv=1.0                                                                696
10811 continue                                                              697
10791 continue                                                              697
      ym=dot_product(w,y)                                                   697
      y=y-ym                                                                697
      ys=sqrt(dot_product(w,y**2))                                          697
      y=y/ys                                                                697
      g=0.0                                                                 698
10830 do 10831 j=1,ni                                                       698
      if(ju(j).eq.0)goto 10831                                              698
      jb=ix(j)                                                              698
      je=ix(j+1)-1                                                          699
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)            700
10831 continue                                                              701
10832 continue                                                              701
      return                                                                702
      end                                                                   703
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,    705 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    706 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                               707
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)           708
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                           709
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           715
      jerr=jerr+ierr                                                        716
      allocate(mm(1:ni),stat=ierr)                                          716
      jerr=jerr+ierr                                                        717
      allocate(da(1:ni),stat=ierr)                                          717
      jerr=jerr+ierr                                                        718
      if(jerr.ne.0) return                                                  719
      bta=max(beta,1.0e-3)                                                  719
      omb=1.0-bta                                                           720
      if(flmin .ge. 1.0)goto 10851                                          720
      eqs=max(eps,flmin)                                                    720
      alf=eqs**(1.0/(nlam-1))                                               720
10851 continue                                                              721
      rsq=0.0                                                               721
      a=0.0                                                                 721
      mm=0                                                                  721
      nlp=0                                                                 721
      nin=nlp                                                               721
      iz=0                                                                  721
      mnl=min(mnlam,nlam)                                                   722
10860 do 10861 m=1,nlam                                                     723
      if(flmin .lt. 1.0)goto 10881                                          723
      alm=ulam(m)                                                           723
      goto 10871                                                            724
10881 if(m .le. 2)goto 10891                                                724
      alm=alm*alf                                                           724
      goto 10871                                                            725
10891 if(m .ne. 1)goto 10901                                                725
      alm=big                                                               725
      goto 10911                                                            726
10901 continue                                                              726
      alm=0.0                                                               727
10920 do 10921 j=1,ni                                                       727
      if(ju(j).eq.0)goto 10921                                              727
      if(vp(j).le.0.0)goto 10921                                            728
      alm=max(alm,abs(g(j))/vp(j))                                          729
10921 continue                                                              730
10922 continue                                                              730
      alm=alf*alm/bta                                                       731
10911 continue                                                              732
10871 continue                                                              732
      dem=alm*omb                                                           732
      ab=alm*bta                                                            732
      rsq0=rsq                                                              732
      jz=1                                                                  733
10930 continue                                                              733
10931 continue                                                              733
      if(iz*jz.ne.0) go to 10260                                            733
      nlp=nlp+1                                                             733
      dlx=0.0                                                               734
10940 do 10941 k=1,ni                                                       734
      if(ju(k).eq.0)goto 10941                                              735
      ak=a(k)                                                               735
      u=g(k)+ak*xv(k)                                                       735
      v=abs(u)-vp(k)*ab                                                     735
      a(k)=0.0                                                              736
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         737
      if(a(k).eq.ak)goto 10941                                              738
      if(mm(k) .ne. 0)goto 10961                                            738
      nin=nin+1                                                             738
      if(nin.gt.nx)goto 10942                                               739
10970 do 10971 j=1,ni                                                       739
      if(ju(j).eq.0)goto 10971                                              740
      if(mm(j) .eq. 0)goto 10991                                            740
      c(j,nin)=c(k,mm(j))                                                   740
      goto 10971                                                            740
10991 continue                                                              741
      if(j .ne. k)goto 11011                                                741
      c(j,nin)=xv(j)                                                        741
      goto 10971                                                            741
11011 continue                                                              742
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))        744
10971 continue                                                              745
10972 continue                                                              745
      mm(k)=nin                                                             745
      ia(nin)=k                                                             746
10961 continue                                                              747
      del=a(k)-ak                                                           747
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      748
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     749
11020 do 11021 j=1,ni                                                       749
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               749
11021 continue                                                              750
11022 continue                                                              750
10941 continue                                                              751
10942 continue                                                              751
      if(dlx.lt.thr)goto 10932                                              751
      if(nin.gt.nx)goto 10932                                               752
10260 continue                                                              752
      iz=1                                                                  752
      da(1:nin)=a(ia(1:nin))                                                753
11030 continue                                                              753
11031 continue                                                              753
      nlp=nlp+1                                                             753
      dlx=0.0                                                               754
11040 do 11041 l=1,nin                                                      754
      k=ia(l)                                                               755
      ak=a(k)                                                               755
      u=g(k)+ak*xv(k)                                                       755
      v=abs(u)-vp(k)*ab                                                     755
      a(k)=0.0                                                              756
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         757
      if(a(k).eq.ak)goto 11041                                              758
      del=a(k)-ak                                                           758
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      759
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     760
11050 do 11051 j=1,nin                                                      760
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  760
11051 continue                                                              761
11052 continue                                                              761
11041 continue                                                              762
11042 continue                                                              762
      if(dlx.lt.thr)goto 11032                                              762
      goto 11031                                                            763
11032 continue                                                              763
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      764
11060 do 11061 j=1,ni                                                       764
      if(mm(j).ne.0)goto 11061                                              765
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            766
11061 continue                                                              767
11062 continue                                                              767
      jz=0                                                                  768
      goto 10931                                                            769
10932 continue                                                              769
      if(nin.gt.nx)goto 10862                                               770
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 770
      kin(m)=nin                                                            771
      rsqo(m)=rsq                                                           771
      almo(m)=alm                                                           771
      lmu=m                                                                 772
      if(m.lt.mnl)goto 10861                                                772
      if(flmin.ge.1.0)goto 10861                                            773
      me=0                                                                  773
11070 do 11071 j=1,nin                                                      773
      if(ao(j,m).ne.0.0) me=me+1                                            773
11071 continue                                                              773
11072 continue                                                              773
      if(me.gt.ne)goto 10862                                                774
      if(rsq-rsq0.lt.sml*rsq)goto 10862                                     774
      if(rsq.gt.rsqmax)goto 10862                                           775
10861 continue                                                              776
10862 continue                                                              776
      deallocate(a,mm,c,da)                                                 777
      return                                                                778
      end                                                                   779
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,    781 
     *ulam,  thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                               782
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         783
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            784
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          789
      allocate(xs(1:ni),stat=ierr)                                          789
      jerr=jerr+ierr                                                        790
      allocate(ju(1:ni),stat=ierr)                                          790
      jerr=jerr+ierr                                                        791
      allocate(xv(1:ni),stat=ierr)                                          791
      jerr=jerr+ierr                                                        792
      allocate(vlam(1:nlam),stat=ierr)                                      792
      jerr=jerr+ierr                                                        793
      if(jerr.ne.0) return                                                  794
      call spchkvars(no,ni,x,ix,ju)                                         795
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  796
      if(maxval(ju) .gt. 0)goto 11091                                       796
      jerr=7777                                                             796
      return                                                                796
11091 continue                                                              797
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)        798
      if(jerr.ne.0) return                                                  799
      if(flmin.ge.1.0) vlam=ulam/ys                                         800
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t    802 
     *hr,xm,xs,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  803
11100 do 11101 k=1,lmu                                                      803
      alm(k)=ys*alm(k)                                                      803
      nk=nin(k)                                                             804
11110 do 11111 l=1,nk                                                       804
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          804
11111 continue                                                              805
11112 continue                                                              805
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         806
11101 continue                                                              807
11102 continue                                                              807
      deallocate(xm,xs,ju,xv,vlam)                                          808
      return                                                                809
      end                                                                   810
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je    811 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                            811
      integer ix(*),jx(*),ju(ni)                                            812
      w=w/sum(w)                                                            813
11120 do 11121 j=1,ni                                                       813
      if(ju(j).eq.0)goto 11121                                              814
      jb=ix(j)                                                              814
      je=ix(j+1)-1                                                          814
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                              815
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                  816
11121 continue                                                              817
11122 continue                                                              817
      if(isd .ne. 0)goto 11141                                              817
      xs=1.0                                                                817
      goto 11151                                                            818
11141 continue                                                              818
11160 do 11161 j=1,ni                                                       818
      if(ju(j).ne.0) xs(j)=sqrt(xv(j))                                      818
11161 continue                                                              818
11162 continue                                                              818
      xv=1.0                                                                818
11151 continue                                                              819
11131 continue                                                              819
      ym=dot_product(w,y)                                                   819
      y=y-ym                                                                819
      ys=sqrt(dot_product(w,y**2))                                          819
      y=y/ys                                                                820
      return                                                                821
      end                                                                   822
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,    824 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    825 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                               826
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)           827
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                           828
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           833
      allocate(mm(1:ni),stat=ierr)                                          833
      jerr=jerr+ierr                                                        834
      if(jerr.ne.0) return                                                  835
      bta=max(beta,1.0e-3)                                                  835
      omb=1.0-bta                                                           836
      if(flmin .ge. 1.0)goto 11181                                          836
      eqs=max(eps,flmin)                                                    836
      alf=eqs**(1.0/(nlam-1))                                               836
11181 continue                                                              837
      rsq=0.0                                                               837
      a=0.0                                                                 837
      mm=0                                                                  837
      o=0.0                                                                 837
      nlp=0                                                                 837
      nin=nlp                                                               837
      iz=0                                                                  837
      mnl=min(mnlam,nlam)                                                   838
11190 do 11191 m=1,nlam                                                     839
      if(flmin .lt. 1.0)goto 11211                                          839
      alm=ulam(m)                                                           839
      goto 11201                                                            840
11211 if(m .le. 2)goto 11221                                                840
      alm=alm*alf                                                           840
      goto 11201                                                            841
11221 if(m .ne. 1)goto 11231                                                841
      alm=big                                                               841
      goto 11241                                                            842
11231 continue                                                              842
      alm=0.0                                                               843
11250 do 11251 j=1,ni                                                       843
      if(ju(j).eq.0)goto 11251                                              843
      if(vp(j).le.0.0)goto 11251                                            844
      jb=ix(j)                                                              844
      je=ix(j+1)-1                                                          845
      alm=max(alm,abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))     847 
     * /(vp(j)*xs(j))))
11251 continue                                                              848
11252 continue                                                              848
      alm=alf*alm/bta                                                       849
11241 continue                                                              850
11201 continue                                                              850
      dem=alm*omb                                                           850
      ab=alm*bta                                                            850
      rsq0=rsq                                                              850
      jz=1                                                                  851
11260 continue                                                              851
11261 continue                                                              851
      if(iz*jz.ne.0) go to 10260                                            851
      nlp=nlp+1                                                             851
      dlx=0.0                                                               852
11270 do 11271 k=1,ni                                                       852
      if(ju(k).eq.0)goto 11271                                              852
      jb=ix(k)                                                              852
      je=ix(k+1)-1                                                          853
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)            854
      ak=a(k)                                                               854
      u=gk+ak*xv(k)                                                         854
      v=abs(u)-vp(k)*ab                                                     854
      a(k)=0.0                                                              855
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         856
      if(a(k).eq.ak)goto 11271                                              857
      if(mm(k) .ne. 0)goto 11291                                            857
      nin=nin+1                                                             857
      if(nin.gt.nx)goto 11272                                               858
      mm(k)=nin                                                             858
      ia(nin)=k                                                             859
11291 continue                                                              860
      del=a(k)-ak                                                           860
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        861
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                          862
      o=o+del*xm(k)/xs(k)                                                   862
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     863
11271 continue                                                              864
11272 continue                                                              864
      if(dlx.lt.thr)goto 11262                                              864
      if(nin.gt.nx)goto 11262                                               865
10260 continue                                                              865
      iz=1                                                                  866
11300 continue                                                              866
11301 continue                                                              866
      nlp=nlp+1                                                             866
      dlx=0.0                                                               867
11310 do 11311 l=1,nin                                                      867
      k=ia(l)                                                               867
      jb=ix(k)                                                              867
      je=ix(k+1)-1                                                          868
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)            869
      ak=a(k)                                                               869
      u=gk+ak*xv(k)                                                         869
      v=abs(u)-vp(k)*ab                                                     869
      a(k)=0.0                                                              870
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         871
      if(a(k).eq.ak)goto 11311                                              872
      del=a(k)-ak                                                           872
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        873
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                          874
      o=o+del*xm(k)/xs(k)                                                   874
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     875
11311 continue                                                              876
11312 continue                                                              876
      if(dlx.lt.thr)goto 11302                                              876
      goto 11301                                                            877
11302 continue                                                              877
      jz=0                                                                  878
      goto 11261                                                            879
11262 continue                                                              879
      if(nin.gt.nx)goto 11192                                               880
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 880
      kin(m)=nin                                                            881
      rsqo(m)=rsq                                                           881
      almo(m)=alm                                                           881
      lmu=m                                                                 882
      if(m.lt.mnl)goto 11191                                                882
      if(flmin.ge.1.0)goto 11191                                            883
      me=0                                                                  883
11320 do 11321 j=1,nin                                                      883
      if(ao(j,m).ne.0.0) me=me+1                                            883
11321 continue                                                              883
11322 continue                                                              883
      if(me.gt.ne)goto 11192                                                884
      if(rsq-rsq0.lt.sml*rsq)goto 11192                                     884
      if(rsq.gt.rsqmax)goto 11192                                           885
11191 continue                                                              886
11192 continue                                                              886
      deallocate(a,mm)                                                      887
      return                                                                888
      end                                                                   889
      subroutine spchkvars(no,ni,x,ix,ju)                                   890
      real x(*)                                                             890
      integer ix(*),ju(ni)                                                  891
11330 do 11331 j=1,ni                                                       891
      ju(j)=0                                                               891
      jb=ix(j)                                                              891
      nj=ix(j+1)-jb                                                         891
      if(nj.eq.0)goto 11331                                                 892
      je=ix(j+1)-1                                                          893
      if(nj .ge. no)goto 11351                                              893
11360 do 11361 i=jb,je                                                      893
      if(x(i).eq.0.0)goto 11361                                             893
      ju(j)=1                                                               893
      goto 11362                                                            893
11361 continue                                                              893
11362 continue                                                              893
      goto 11371                                                            894
11351 continue                                                              894
      t=x(jb)                                                               894
11380 do 11381 i=jb+1,je                                                    894
      if(x(i).eq.t)goto 11381                                               894
      ju(j)=1                                                               894
      goto 11382                                                            894
11381 continue                                                              894
11382 continue                                                              894
11371 continue                                                              895
11341 continue                                                              895
11331 continue                                                              896
11332 continue                                                              896
      return                                                                897
      end                                                                   898
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                          899
      real ca(*),x(*),f(n)                                                  899
      integer ia(*),ix(*),jx(*)                                             900
      f=a0                                                                  901
11390 do 11391 j=1,nin                                                      901
      k=ia(j)                                                               901
      kb=ix(k)                                                              901
      ke=ix(k+1)-1                                                          902
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                              903
11391 continue                                                              904
11392 continue                                                              904
      return                                                                905
      end                                                                   906
      function row_prod(i,j,ia,ja,ra,w)                                     907
      integer ia(*),ja(*)                                                   907
      real ra(*),w(*)                                                       908
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(    910 
     *i),ia(j+1)-ia(j),w)
      return                                                                911
      end                                                                   912
      function dot(x,y,mx,my,nx,ny,w)                                       913
      real x(*),y(*),w(*)                                                   913
      integer mx(*),my(*)                                                   914
      i=1                                                                   914
      j=i                                                                   914
      s=0.0                                                                 915
11400 continue                                                              915
11401 continue                                                              915
11410 continue                                                              916
11411 if(mx(i).ge.my(j))goto 11412                                          916
      i=i+1                                                                 916
      if(i.gt.nx) go to 11420                                               916
      goto 11411                                                            917
11412 continue                                                              917
      if(mx(i).eq.my(j)) go to 11430                                        918
11440 continue                                                              918
11441 if(my(j).ge.mx(i))goto 11442                                          918
      j=j+1                                                                 918
      if(j.gt.ny) go to 11420                                               918
      goto 11441                                                            919
11442 continue                                                              919
      if(mx(i).eq.my(j)) go to 11430                                        919
      goto 11401                                                            920
11430 continue                                                              920
      s=s+w(mx(i))*x(i)*y(j)                                                921
      i=i+1                                                                 921
      if(i.gt.nx)goto 11402                                                 921
      j=j+1                                                                 921
      if(j.gt.ny)goto 11402                                                 922
      goto 11401                                                            923
11402 continue                                                              923
11420 continue                                                              923
      dot=s                                                                 924
      return                                                                925
      end                                                                   926
      subroutine lognet (parm,no,ni,nc,x,y,jd,vp,ne,nx,nlam,flmin,ulam,t    928 
     *hr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),vp(ni),ulam(nlam)                       929
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   930
      integer jd(*),ia(nx),nin(nlam)                                        931
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11461                                     935
      jerr=10000                                                            935
      return                                                                935
11461 continue                                                              936
      allocate(ww(1:no),stat=jerr)                                          937
      allocate(ju(1:ni),stat=ierr)                                          937
      jerr=jerr+ierr                                                        938
      allocate(vq(1:ni),stat=ierr)                                          938
      jerr=jerr+ierr                                                        939
      allocate(xm(1:ni),stat=ierr)                                          939
      jerr=jerr+ierr                                                        940
      if(isd .le. 0)goto 11481                                              940
      allocate(xs(1:ni),stat=ierr)                                          940
      jerr=jerr+ierr                                                        940
11481 continue                                                              941
      if(jerr.ne.0) return                                                  942
      call chkvars(no,ni,x,ju)                                              943
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  944
      if(maxval(ju) .gt. 0)goto 11501                                       944
      jerr=7777                                                             944
      return                                                                944
11501 continue                                                              945
      vq=max(0.0,vp)                                                        945
      vq=vq*ni/sum(vq)                                                      946
11510 do 11511 i=1,no                                                       946
      ww(i)=sum(y(i,:))                                                     946
      y(i,:)=y(i,:)/ww(i)                                                   946
11511 continue                                                              946
11512 continue                                                              946
      ww=ww/sum(ww)                                                         947
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                              948
      if(nc .ne. 1)goto 11531                                               949
      call lognet2n(parm,no,ni,x,y(:,1),ww,ju,vq,ne,nx,nlam,flmin,ulam,t    951 
     *hr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      goto 11541                                                            952
11531 continue                                                              953
      call lognetn(parm,no,ni,nc,x,y,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,    955 
     *  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
11541 continue                                                              956
11521 continue                                                              956
      if(jerr.gt.0) return                                                  957
11550 do 11551 k=1,lmu                                                      957
      nk=nin(k)                                                             958
11560 do 11561 ic=1,nc                                                      958
      if(isd .le. 0)goto 11581                                              958
11590 do 11591 l=1,nk                                                       958
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                       958
11591 continue                                                              958
11592 continue                                                              958
11581 continue                                                              959
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))             960
11561 continue                                                              961
11562 continue                                                              961
11551 continue                                                              962
11552 continue                                                              962
      deallocate(ww,ju,vq,xm)                                               962
      if(isd.gt.0) deallocate(xs)                                           963
      return                                                                964
      end                                                                   965
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                        966
      real x(no,ni),w(no),xm(ni),xs(ni)                                     966
      integer ju(ni)                                                        967
11600 do 11601 j=1,ni                                                       967
      if(ju(j).eq.0)goto 11601                                              968
      xm(j)=dot_product(w,x(:,j))                                           968
      x(1:no,j)=x(1:no,j)-xm(j)                                             969
      if(isd .le. 0)goto 11621                                              970
      xs(j)=sqrt(dot_product(w*x(:,j),x(:,j)))                              971
      x(1:no,j)=x(1:no,j)/xs(j)                                             972
11621 continue                                                              973
11601 continue                                                              974
11602 continue                                                              974
      return                                                                975
      end                                                                   976
      subroutine lognet2n(parm,no,ni,x,y,w,ju,vp,ne,nx,nlam,flmin,ulam,s    978 
     *hr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=    980 
     *5, devmax=0.999)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           981
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                          982
      integer ju(ni),m(nx),kin(nlam)                                        983
      real, dimension (:), allocatable :: b,bs,v,r,xv,q                         
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                           988
      allocate(xv(1:ni),stat=ierr)                                          988
      jerr=jerr+ierr                                                        989
      allocate(bs(0:ni),stat=ierr)                                          989
      jerr=jerr+ierr                                                        990
      allocate(mm(1:ni),stat=ierr)                                          990
      jerr=jerr+ierr                                                        991
      allocate(r(1:no),stat=ierr)                                           991
      jerr=jerr+ierr                                                        992
      allocate(v(1:no),stat=ierr)                                           992
      jerr=jerr+ierr                                                        993
      allocate(q(1:no),stat=ierr)                                           993
      jerr=jerr+ierr                                                        994
      if(jerr.ne.0) return                                                  995
      fmax=log(1.0/pmin-1.0)                                                995
      fmin=-fmax                                                            995
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                       996
      bta=max(parm,1.0e-3)                                                  996
      omb=1.0-bta                                                           997
      q0=dot_product(w,y)                                                   997
      if(q0 .gt. pmin)goto 11641                                            997
      jerr=8001                                                             997
      return                                                                997
11641 continue                                                              998
      if(q0 .lt. 1.0-pmin)goto 11661                                        998
      jerr=9001                                                             998
      return                                                                998
11661 continue                                                              999
      vi=q0*(1.0-q0)                                                        999
      b(1:ni)=0.0                                                           999
      b(0)=log(q0/(1.0-q0))                                                 999
      v=vi*w                                                               1000
      r=w*(y-q0)                                                           1000
      dev1=-(b(0)*q0+log(1.0-q0))                                          1000
      q=q0                                                                 1001
      if(isd .le. 0)goto 11681                                             1001
      xv=0.25                                                              1001
      goto 11691                                                           1002
11681 continue                                                             1002
11700 do 11701 j=1,ni                                                      1002
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1002
11701 continue                                                             1002
11702 continue                                                             1002
11691 continue                                                             1003
11671 continue                                                             1003
      xmz=vi                                                               1003
      bs=0.0                                                               1004
      if(flmin .ge. 1.0)goto 11721                                         1004
      eqs=max(eps,flmin)                                                   1004
      alf=eqs**(1.0/(nlam-1))                                              1004
11721 continue                                                             1005
      m=0                                                                  1005
      mm=0                                                                 1005
      nlp=0                                                                1005
      nin=nlp                                                              1005
      mnl=min(mnlam,nlam)                                                  1006
11730 do 11731 ilm=1,nlam                                                  1007
      if(flmin .lt. 1.0)goto 11751                                         1007
      al=ulam(ilm)                                                         1007
      goto 11741                                                           1008
11751 if(ilm .le. 2)goto 11761                                             1008
      al=al*alf                                                            1008
      goto 11741                                                           1009
11761 if(ilm .ne. 1)goto 11771                                             1009
      al=big                                                               1009
      goto 11781                                                           1010
11771 continue                                                             1010
      al=0.0                                                               1011
11790 do 11791 j=1,ni                                                      1011
      if(ju(j).eq.0)goto 11791                                             1011
      if(vp(j).le.0.0)goto 11791                                           1012
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1013
11791 continue                                                             1014
11792 continue                                                             1014
      al=alf*al/bta                                                        1015
11781 continue                                                             1016
11741 continue                                                             1016
      al2=al*omb                                                           1016
      al1=al*bta                                                           1016
      nit=0                                                                1017
11800 continue                                                             1017
11801 continue                                                             1017
      bs(0)=b(0)                                                           1017
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1018
11810 continue                                                             1018
11811 continue                                                             1018
      nlp=nlp+1                                                            1018
      dlx=0.0                                                              1019
11820 do 11821 k=1,ni                                                      1019
      if(ju(k).eq.0)goto 11821                                             1020
      bk=b(k)                                                              1020
      gk=dot_product(r,x(:,k))                                             1021
      u=gk+xv(k)*b(k)                                                      1021
      au=abs(u)-vp(k)*al1                                                  1022
      if(au .gt. 0.0)goto 11841                                            1022
      b(k)=0.0                                                             1022
      goto 11851                                                           1023
11841 continue                                                             1023
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1023
11851 continue                                                             1024
11831 continue                                                             1024
      d=b(k)-bk                                                            1024
      if(abs(d).le.0.0)goto 11821                                          1024
      dlx=max(dlx,abs(d))                                                  1025
      r=r-d*v*x(:,k)                                                       1026
      if(mm(k) .ne. 0)goto 11871                                           1026
      nin=nin+1                                                            1026
      if(nin.gt.nx)goto 11822                                              1027
      mm(k)=nin                                                            1027
      m(nin)=k                                                             1028
11871 continue                                                             1029
11821 continue                                                             1030
11822 continue                                                             1030
      if(nin.gt.nx)goto 11812                                              1031
      d=sum(r)/xmz                                                         1032
      if(d .eq. 0.0)goto 11891                                             1032
      b(0)=b(0)+d                                                          1032
      dlx=max(dlx,abs(d))                                                  1032
      r=r-d*v                                                              1032
11891 continue                                                             1033
      if(dlx.lt.shr)goto 11812                                             1034
11900 continue                                                             1034
11901 continue                                                             1034
      nlp=nlp+1                                                            1034
      dlx=0.0                                                              1035
11910 do 11911 l=1,nin                                                     1035
      k=m(l)                                                               1035
      bk=b(k)                                                              1035
      gk=dot_product(r,x(:,k))                                             1036
      u=gk+xv(k)*b(k)                                                      1036
      au=abs(u)-vp(k)*al1                                                  1037
      if(au .gt. 0.0)goto 11931                                            1037
      b(k)=0.0                                                             1037
      goto 11941                                                           1038
11931 continue                                                             1038
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1038
11941 continue                                                             1039
11921 continue                                                             1039
      d=b(k)-bk                                                            1039
      if(abs(d).le.0.0)goto 11911                                          1039
      dlx=max(dlx,abs(d))                                                  1040
      r=r-d*v*x(:,k)                                                       1041
11911 continue                                                             1042
11912 continue                                                             1042
      d=sum(r)/xmz                                                         1043
      if(d .eq. 0.0)goto 11961                                             1043
      b(0)=b(0)+d                                                          1043
      dlx=max(dlx,abs(d))                                                  1043
      r=r-d*v                                                              1043
11961 continue                                                             1045
      if(dlx.lt.shr)goto 11902                                             1045
      goto 11901                                                           1046
11902 continue                                                             1046
      goto 11811                                                           1047
11812 continue                                                             1047
      if(nin.gt.nx)goto 11802                                              1048
      if(abs(b(0)-bs(0)) .ge. shr)goto 11981                               1048
      ix=0                                                                 1049
11990 do 11991 j=1,nin                                                     1049
      if(abs(b(m(j))-bs(m(j))).lt.shr)goto 11991                           1049
      ix=1                                                                 1049
      goto 11992                                                           1049
11991 continue                                                             1050
11992 continue                                                             1050
      if(ix.eq.0)goto 11802                                                1051
11981 continue                                                             1052
12000 do 12001 i=1,no                                                      1052
      fi=b(0)                                                              1053
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1054
      if(fi .ge. fmin)goto 12021                                           1054
      q(i)=0.0                                                             1054
      goto 12011                                                           1054
12021 if(fi .le. fmax)goto 12031                                           1054
      q(i)=1.0                                                             1054
      goto 12041                                                           1055
12031 continue                                                             1055
      q(i)=1.0/(1.0+exp(-fi))                                              1055
12041 continue                                                             1056
12011 continue                                                             1056
12001 continue                                                             1057
12002 continue                                                             1057
      v=w*q*(1.0-q)                                                        1057
      xmz=sum(v)                                                           1057
      if(xmz.le.vmin)goto 11802                                            1057
      r=w*(y-q)                                                            1058
      if(kopt .ne. 0)goto 12061                                            1059
12070 do 12071 j=1,nin                                                     1059
      xv(m(j))=dot_product(v,x(:,m(j))**2)                                 1059
12071 continue                                                             1060
12072 continue                                                             1060
12061 continue                                                             1061
      nit=nit+1                                                            1061
      if(nit .le. maxit)goto 12091                                         1061
      jerr=-ilm                                                            1061
      return                                                               1061
12091 continue                                                             1062
      goto 11801                                                           1063
11802 continue                                                             1063
      if(nin .le. nx)goto 12111                                            1063
      jerr=-10000-ilm                                                      1063
      goto 11732                                                           1063
12111 continue                                                             1064
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1064
      kin(ilm)=nin                                                         1065
      a0(ilm)=b(0)                                                         1065
      alm(ilm)=al                                                          1065
      lmu=ilm                                                              1066
      devi=dev2(no,w,y,q,pmin)                                             1067
      dev(ilm)=(dev1-devi)/dev1                                            1067
      if(xmz.le.vmin)goto 11732                                            1068
      if(ilm.lt.mnl)goto 11731                                             1068
      if(flmin.ge.1.0)goto 11731                                           1069
      me=0                                                                 1069
12120 do 12121 j=1,nin                                                     1069
      if(a(j,ilm).ne.0.0) me=me+1                                          1069
12121 continue                                                             1069
12122 continue                                                             1069
      if(me.gt.ne)goto 11732                                               1070
      if(dev(ilm).gt.devmax)goto 11732                                     1070
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 11732                             1071
11731 continue                                                             1072
11732 continue                                                             1072
      deallocate(b,bs,v,r,xv,q,mm)                                         1073
      return                                                               1074
      end                                                                  1075
      function dev2(n,w,y,p,pmin)                                          1076
      real w(n),y(n),p(n)                                                  1077
      pmax=1.0-pmin                                                        1077
      s=0.0                                                                1078
12130 do 12131 i=1,n                                                       1078
      pi=min(max(pmin,p(i)),pmax)                                          1079
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1080
12131 continue                                                             1081
12132 continue                                                             1081
      dev2=s                                                               1082
      return                                                               1083
      end                                                                  1084
      subroutine lognetn(parm,no,ni,nc,x,y,w,ju,vp,ne,nx,nlam,flmin,ulam   1086 
     *,shr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1088 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),w(no),vp(ni),ulam(nlam)                       1089
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1090
      integer ju(ni),m(nx),kin(nlam)                                       1091
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: di,v,r                                
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1102
      jerr=jerr+ierr                                                       1103
      allocate(v(1:no),stat=ierr)                                          1103
      jerr=jerr+ierr                                                       1104
      allocate(mm(1:ni),stat=ierr)                                         1104
      jerr=jerr+ierr                                                       1105
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1105
      jerr=jerr+ierr                                                       1106
      allocate(sxp(1:no),stat=ierr)                                        1106
      jerr=jerr+ierr                                                       1107
      allocate(di(1:no),stat=ierr)                                         1107
      jerr=jerr+ierr                                                       1108
      if(jerr.ne.0) return                                                 1109
      pmax=1.0-pmin                                                        1109
      emin=pmin/pmax                                                       1109
      emax=1.0/emin                                                        1110
      pfm=(1.0+pmin)*pmin                                                  1110
      pfx=(1.0-pmin)*pmax                                                  1110
      vmin=pfm*pmax                                                        1111
      bta=max(parm,1.0e-3)                                                 1111
      omb=1.0-bta                                                          1111
      dev1=0.0                                                             1112
12140 do 12141 ic=1,nc                                                     1112
      q0=dot_product(w,y(:,ic))                                            1113
      if(q0 .gt. pmin)goto 12161                                           1113
      jerr =8000+ic                                                        1113
      return                                                               1113
12161 continue                                                             1114
      if(q0 .lt. 1.0-pmin)goto 12181                                       1114
      jerr =9000+ic                                                        1114
      return                                                               1114
12181 continue                                                             1115
      vi=q0*(1.0-q0)                                                       1115
      v=vi*w                                                               1115
      b(1:ni,ic)=0.0                                                       1116
      b(0,ic)=log(q0)                                                      1116
      dev1=dev1-q0*b(0,ic)                                                 1117
12141 continue                                                             1118
12142 continue                                                             1118
      if(isd .le. 0)goto 12201                                             1118
      xv=0.25                                                              1118
      goto 12211                                                           1119
12201 continue                                                             1119
12220 do 12221 j=1,ni                                                      1119
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1119
12221 continue                                                             1119
12222 continue                                                             1119
12211 continue                                                             1120
12191 continue                                                             1120
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1120
      sxp=0.0                                                              1121
12230 do 12231 ic=1,nc                                                     1121
      q(:,ic)=exp(b(0,ic))                                                 1121
      sxp=sxp+q(:,ic)                                                      1121
12231 continue                                                             1122
12232 continue                                                             1122
      if(flmin .ge. 1.0)goto 12251                                         1122
      eqs=max(eps,flmin)                                                   1122
      alf=eqs**(1.0/(nlam-1))                                              1122
12251 continue                                                             1123
      m=0                                                                  1123
      mm=0                                                                 1123
      nin=0                                                                1123
      nlp=0                                                                1123
      mnl=min(mnlam,nlam)                                                  1123
      bs=0.0                                                               1124
12260 do 12261 ilm=1,nlam                                                  1125
      if(flmin .lt. 1.0)goto 12281                                         1125
      al=ulam(ilm)                                                         1125
      goto 12271                                                           1126
12281 if(ilm .le. 2)goto 12291                                             1126
      al=al*alf                                                            1126
      goto 12271                                                           1127
12291 if(ilm .ne. 1)goto 12301                                             1127
      al=big                                                               1127
      goto 12311                                                           1128
12301 continue                                                             1128
      al=0.0                                                               1129
12320 do 12321 ic=1,nc                                                     1129
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1130
12330 do 12331 j=1,ni                                                      1130
      if(ju(j).eq.0)goto 12331                                             1130
      if(vp(j).le.0.0)goto 12331                                           1131
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1132
12331 continue                                                             1133
12332 continue                                                             1133
12321 continue                                                             1134
12322 continue                                                             1134
      al=alf*al/bta                                                        1135
12311 continue                                                             1136
12271 continue                                                             1136
      al2=al*omb                                                           1136
      al1=al*bta                                                           1136
      nit=0                                                                1137
12340 continue                                                             1137
12341 continue                                                             1137
      ix=0                                                                 1137
      jx=ix                                                                1137
      ig=0                                                                 1138
12350 do 12351 ic=1,nc                                                     1138
      bs(0,ic)=b(0,ic)                                                     1139
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1140
      xmz=0.0                                                              1141
12360 do 12361 i=1,no                                                      1141
      pic=q(i,ic)/sxp(i)                                                   1142
      if(pic .ge. pfm)goto 12381                                           1142
      pic=0.0                                                              1142
      v(i)=0.0                                                             1142
      goto 12371                                                           1143
12381 if(pic .le. pfx)goto 12391                                           1143
      pic=1.0                                                              1143
      v(i)=0.0                                                             1143
      goto 12401                                                           1144
12391 continue                                                             1144
      v(i)=w(i)*pic*(1.0-pic)                                              1144
      xmz=xmz+v(i)                                                         1144
12401 continue                                                             1145
12371 continue                                                             1145
      r(i)=w(i)*(y(i,ic)-pic)                                              1146
12361 continue                                                             1147
12362 continue                                                             1147
      if(xmz.le.vmin)goto 12351                                            1147
      ig=1                                                                 1148
      if(kopt .ne. 0)goto 12421                                            1149
12430 do 12431 j=1,nin                                                     1149
      xv(m(j),ic)=dot_product(v,x(:,m(j))**2)                              1149
12431 continue                                                             1150
12432 continue                                                             1150
12421 continue                                                             1151
12440 continue                                                             1151
12441 continue                                                             1151
      nlp=nlp+1                                                            1151
      dlx=0.0                                                              1152
12450 do 12451 k=1,ni                                                      1152
      if(ju(k).eq.0)goto 12451                                             1153
      bk=b(k,ic)                                                           1153
      gk=dot_product(r,x(:,k))                                             1154
      u=gk+xv(k,ic)*b(k,ic)                                                1154
      au=abs(u)-vp(k)*al1                                                  1155
      if(au .gt. 0.0)goto 12471                                            1155
      b(k,ic)=0.0                                                          1155
      goto 12481                                                           1156
12471 continue                                                             1156
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1156
12481 continue                                                             1157
12461 continue                                                             1157
      d=b(k,ic)-bk                                                         1157
      if(abs(d).le.0.0)goto 12451                                          1157
      dlx=max(dlx,abs(d))                                                  1158
      r=r-d*v*x(:,k)                                                       1159
      if(mm(k) .ne. 0)goto 12501                                           1159
      nin=nin+1                                                            1160
      if(nin .le. nx)goto 12521                                            1160
      jx=1                                                                 1160
      goto 12452                                                           1160
12521 continue                                                             1161
      mm(k)=nin                                                            1161
      m(nin)=k                                                             1162
12501 continue                                                             1163
12451 continue                                                             1164
12452 continue                                                             1164
      if(jx.gt.0)goto 12442                                                1165
      d=sum(r)/xmz                                                         1166
      if(d .eq. 0.0)goto 12541                                             1166
      b(0,ic)=b(0,ic)+d                                                    1166
      dlx=max(dlx,abs(d))                                                  1166
      r=r-d*v                                                              1166
12541 continue                                                             1167
      if(dlx.lt.shr)goto 12442                                             1168
12550 continue                                                             1168
12551 continue                                                             1168
      nlp=nlp+1                                                            1168
      dlx=0.0                                                              1169
12560 do 12561 l=1,nin                                                     1169
      k=m(l)                                                               1169
      bk=b(k,ic)                                                           1170
      gk=dot_product(r,x(:,k))                                             1171
      u=gk+xv(k,ic)*b(k,ic)                                                1171
      au=abs(u)-vp(k)*al1                                                  1172
      if(au .gt. 0.0)goto 12581                                            1172
      b(k,ic)=0.0                                                          1172
      goto 12591                                                           1173
12581 continue                                                             1173
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1173
12591 continue                                                             1174
12571 continue                                                             1174
      d=b(k,ic)-bk                                                         1174
      if(abs(d).le.0.0)goto 12561                                          1175
      dlx=max(dlx,abs(d))                                                  1175
      r=r-d*v*x(:,k)                                                       1176
12561 continue                                                             1177
12562 continue                                                             1177
      d=sum(r)/xmz                                                         1178
      if(d .eq. 0.0)goto 12611                                             1178
      b(0,ic)=b(0,ic)+d                                                    1179
      dlx=max(dlx,abs(d))                                                  1179
      r=r-d*v                                                              1180
12611 continue                                                             1181
      if(dlx.lt.shr)goto 12552                                             1181
      goto 12551                                                           1182
12552 continue                                                             1182
      goto 12441                                                           1183
12442 continue                                                             1183
      if(jx.gt.0)goto 12352                                                1184
      if(abs(b(0,ic)-bs(0,ic)).gt.shr) ix=1                                1185
      if(ix .ne. 0)goto 12631                                              1186
12640 do 12641 j=1,nin                                                     1187
      if(abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 12661                   1187
      ix=1                                                                 1187
      goto 12642                                                           1187
12661 continue                                                             1188
12641 continue                                                             1189
12642 continue                                                             1189
12631 continue                                                             1190
12670 do 12671 i=1,no                                                      1190
      fi=b(0,ic)                                                           1192
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1193
      fi=min(max(exmn,fi),exmx)                                            1193
      sxp(i)=sxp(i)-q(i,ic)                                                1194
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))              1195
      sxp(i)=sxp(i)+q(i,ic)                                                1196
12671 continue                                                             1197
12672 continue                                                             1197
12351 continue                                                             1198
12352 continue                                                             1198
      if(jx.gt.0)goto 12342                                                1198
      if(ix.eq.0)goto 12342                                                1198
      if(ig.eq.0)goto 12342                                                1199
      s=-sum(b(0,:))/nc                                                    1199
      b(0,:)=b(0,:)+s                                                      1199
      di=s                                                                 1200
12680 do 12681 j=1,nin                                                     1200
      l=m(j)                                                               1201
      if(vp(l) .gt. 0.0)goto 12701                                         1201
      s=sum(b(l,:))/nc                                                     1201
      goto 12711                                                           1202
12701 continue                                                             1202
      s=elc(parm,nc,b(l,:),is)                                             1202
12711 continue                                                             1203
12691 continue                                                             1203
      b(l,:)=b(l,:)-s                                                      1203
      di=di-s*x(:,l)                                                       1204
12681 continue                                                             1205
12682 continue                                                             1205
      di=exp(di)                                                           1205
      sxp=sxp*di                                                           1205
12720 do 12721 ic=1,nc                                                     1205
      q(:,ic)=q(:,ic)*di                                                   1205
12721 continue                                                             1206
12722 continue                                                             1206
      nit=nit+1                                                            1206
      if(nit .le. maxit)goto 12741                                         1206
      jerr=-ilm                                                            1206
      return                                                               1206
12741 continue                                                             1207
      goto 12341                                                           1208
12342 continue                                                             1208
      if(jx .le. 0)goto 12761                                              1208
      jerr=-10000-ilm                                                      1208
      goto 12262                                                           1208
12761 continue                                                             1208
      devi=0.0                                                             1209
12770 do 12771 ic=1,nc                                                     1210
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1210
      a0(ic,ilm)=b(0,ic)                                                   1211
12780 do 12781 i=1,no                                                      1211
      if(y(i,ic).le.0.0)goto 12781                                         1212
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1213
12781 continue                                                             1214
12782 continue                                                             1214
12771 continue                                                             1215
12772 continue                                                             1215
      kin(ilm)=nin                                                         1215
      alm(ilm)=al                                                          1215
      lmu=ilm                                                              1216
      dev(ilm)=(dev1-devi)/dev1                                            1216
      if(ig.eq.0)goto 12262                                                1217
      if(ilm.lt.mnl)goto 12261                                             1217
      if(flmin.ge.1.0)goto 12261                                           1218
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12262             1219
      if(dev(ilm).gt.devmax)goto 12262                                     1219
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12262                             1220
12261 continue                                                             1221
12262 continue                                                             1221
      deallocate(sxp,b,bs,v,r,xv,q,mm,is)                                  1222
      return                                                               1223
      end                                                                  1224
      function elc(parm,n,a,m)                                             1225
      real a(n)                                                            1225
      integer m(n)                                                         1226
      fn=n                                                                 1226
      am=sum(a)/fn                                                         1227
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 12801                       1227
      elc=am                                                               1227
      return                                                               1227
12801 continue                                                             1228
12810 do 12811 i=1,n                                                       1228
      m(i)=i                                                               1228
12811 continue                                                             1228
12812 continue                                                             1228
      call psort7(a,m,1,n)                                                 1229
      if(a(m(1)) .ne. a(m(n)))goto 12831                                   1229
      elc=a(1)                                                             1229
      return                                                               1229
12831 continue                                                             1230
      if(mod(n,2) .ne. 1)goto 12851                                        1230
      ad=a(m(n/2+1))                                                       1230
      goto 12861                                                           1231
12851 continue                                                             1231
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1231
12861 continue                                                             1232
12841 continue                                                             1232
      if(parm .ne. 1.0)goto 12881                                          1232
      elc=ad                                                               1232
      return                                                               1232
12881 continue                                                             1233
      b1=min(am,ad)                                                        1233
      b2=max(am,ad)                                                        1233
      k2=1                                                                 1234
12890 continue                                                             1234
12891 if(a(m(k2)).gt.b1)goto 12892                                         1234
      k2=k2+1                                                              1234
      goto 12891                                                           1234
12892 continue                                                             1234
      k1=k2-1                                                              1235
12900 continue                                                             1235
12901 if(a(m(k2)).ge.b2)goto 12902                                         1235
      k2=k2+1                                                              1235
      goto 12901                                                           1236
12902 continue                                                             1236
      r=parm/((1.0-parm)*fn)                                               1236
      is=0                                                                 1236
      sm=n-2*(k1-1)                                                        1237
12910 do 12911 k=k1,k2-1                                                   1237
      sm=sm-2.0                                                            1237
      s=r*sm+am                                                            1238
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 12931                   1238
      is=k                                                                 1238
      goto 12912                                                           1238
12931 continue                                                             1239
12911 continue                                                             1240
12912 continue                                                             1240
      if(is .eq. 0)goto 12951                                              1240
      elc=s                                                                1240
      return                                                               1240
12951 continue                                                             1240
      r2=2.0*r                                                             1240
      s1=a(m(k1))                                                          1240
      am2=2.0*am                                                           1241
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1241
      elc=s1                                                               1242
12960 do 12961 k=k1+1,k2                                                   1242
      s=a(m(k))                                                            1242
      if(s.eq.s1)goto 12961                                                1243
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1244
      if(c .ge. cri)goto 12981                                             1244
      cri=c                                                                1244
      elc=s                                                                1244
12981 continue                                                             1244
      s1=s                                                                 1245
12961 continue                                                             1246
12962 continue                                                             1246
      return                                                               1247
      end                                                                  1248
      function nintot(ni,nx,nc,a,m,nin,is)                                 1249
      real a(nx,nc)                                                        1249
      integer m(nx),is(ni)                                                 1250
      is=0                                                                 1250
      nintot=0                                                             1251
12990 do 12991 ic=1,nc                                                     1251
13000 do 13001 j=1,nin                                                     1251
      k=m(j)                                                               1251
      if(is(k).ne.0)goto 13001                                             1252
      if(a(j,ic).eq.0.0)goto 13001                                         1252
      is(k)=k                                                              1252
      nintot=nintot+1                                                      1253
13001 continue                                                             1253
13002 continue                                                             1253
12991 continue                                                             1254
12992 continue                                                             1254
      return                                                               1255
      end                                                                  1256
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1257
      real ca(nx,nc),a(ni,nc)                                              1257
      integer ia(nx)                                                       1258
      a=0.0                                                                1259
13010 do 13011 ic=1,nc                                                     1259
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1259
13011 continue                                                             1260
13012 continue                                                             1260
      return                                                               1261
      end                                                                  1262
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1263
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1263
      integer ia(nx)                                                       1264
13020 do 13021 i=1,nt                                                      1264
13030 do 13031 ic=1,nc                                                     1264
      ans(ic,i)=a0(ic)                                                     1266
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1267 
     *:nin)))
13031 continue                                                             1267
13032 continue                                                             1267
13021 continue                                                             1268
13022 continue                                                             1268
      return                                                               1269
      end                                                                  1270
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,jd,vp,ne,nx,nlam,flmi   1272 
     *n,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
      real x(*),y(no,max(2,nc)),vp(ni),ulam(nlam)                          1273
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1274
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1275
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13051                                    1279
      jerr=10000                                                           1279
      return                                                               1279
13051 continue                                                             1280
      allocate(ww(1:no),stat=jerr)                                         1281
      allocate(ju(1:ni),stat=ierr)                                         1281
      jerr=jerr+ierr                                                       1282
      allocate(vq(1:ni),stat=ierr)                                         1282
      jerr=jerr+ierr                                                       1283
      allocate(xm(1:ni),stat=ierr)                                         1283
      jerr=jerr+ierr                                                       1284
      allocate(xs(1:ni),stat=ierr)                                         1284
      jerr=jerr+ierr                                                       1285
      if(jerr.ne.0) return                                                 1286
      call spchkvars(no,ni,x,ix,ju)                                        1287
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1288
      if(maxval(ju) .gt. 0)goto 13071                                      1288
      jerr=7777                                                            1288
      return                                                               1288
13071 continue                                                             1289
      vq=max(0.0,vp)                                                       1289
      vq=vq*ni/sum(vq)                                                     1290
13080 do 13081 i=1,no                                                      1290
      ww(i)=sum(y(i,:))                                                    1290
      y(i,:)=y(i,:)/ww(i)                                                  1290
13081 continue                                                             1290
13082 continue                                                             1290
      ww=ww/sum(ww)                                                        1291
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1292
      if(nc .ne. 1)goto 13101                                              1293
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),ww,ju,vq,ne,nx,nlam,flm   1295 
     *in,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev,alm,nlp,je
     *rr)
      goto 13111                                                           1296
13101 continue                                                             1297
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,ww,ju,vq,ne,nx,nlam,flmin,   1299 
     *ulam,thr,  isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev,alm,nlp,jerr)
13111 continue                                                             1300
13091 continue                                                             1300
      if(jerr.gt.0) return                                                 1301
13120 do 13121 k=1,lmu                                                     1301
      nk=nin(k)                                                            1302
13130 do 13131 ic=1,nc                                                     1302
      if(isd .le. 0)goto 13151                                             1302
13160 do 13161 l=1,nk                                                      1302
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1302
13161 continue                                                             1302
13162 continue                                                             1302
13151 continue                                                             1303
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1304
13131 continue                                                             1305
13132 continue                                                             1305
13121 continue                                                             1306
13122 continue                                                             1306
      deallocate(ww,ju,vq,xm,xs)                                           1307
      return                                                               1308
      end                                                                  1309
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1310
      real x(*),w(no),xm(ni),xs(ni)                                        1310
      integer ix(*),jx(*),ju(ni)                                           1311
13170 do 13171 j=1,ni                                                      1311
      if(ju(j).eq.0)goto 13171                                             1311
      jb=ix(j)                                                             1311
      je=ix(j+1)-1                                                         1312
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1313
      if(isd.gt.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1314 
     *)**2)
13171 continue                                                             1315
13172 continue                                                             1315
      if(isd.eq.0) xs=1.0                                                  1316
      return                                                               1317
      end                                                                  1318
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,w,ju,vp,ne,nx,nlam,     1320 
     *flmin,ulam,shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev,alm,nlp,jer
     *r)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1322 
     *5, devmax=0.999)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                              1323
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1324
      real xb(ni),xs(ni)                                                   1324
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1325
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1330
      allocate(xm(0:ni),stat=ierr)                                         1330
      jerr=jerr+ierr                                                       1331
      allocate(xv(1:ni),stat=ierr)                                         1331
      jerr=jerr+ierr                                                       1332
      allocate(bs(0:ni),stat=ierr)                                         1332
      jerr=jerr+ierr                                                       1333
      allocate(mm(1:ni),stat=ierr)                                         1333
      jerr=jerr+ierr                                                       1334
      allocate(q(1:no),stat=ierr)                                          1334
      jerr=jerr+ierr                                                       1335
      allocate(r(1:no),stat=ierr)                                          1335
      jerr=jerr+ierr                                                       1336
      allocate(v(1:no),stat=ierr)                                          1336
      jerr=jerr+ierr                                                       1337
      allocate(sc(1:no),stat=ierr)                                         1337
      jerr=jerr+ierr                                                       1338
      if(jerr.ne.0) return                                                 1339
      fmax=log(1.0/pmin-1.0)                                               1339
      fmin=-fmax                                                           1339
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1340
      bta=max(parm,1.0e-3)                                                 1340
      omb=1.0-bta                                                          1341
      q0=dot_product(w,y)                                                  1341
      if(q0 .gt. pmin)goto 13191                                           1341
      jerr=8001                                                            1341
      return                                                               1341
13191 continue                                                             1342
      if(q0 .lt. 1.0-pmin)goto 13211                                       1342
      jerr=9001                                                            1342
      return                                                               1342
13211 continue                                                             1343
      vi=q0*(1.0-q0)                                                       1343
      b(1:ni)=0.0                                                          1343
      b(0)=log(q0/(1.0-q0))                                                1343
      v=vi*w                                                               1344
      r=w*(y-q0)                                                           1344
      dev1=-(b(0)*q0+log(1.0-q0))                                          1344
      q=q0                                                                 1345
      if(isd .le. 0)goto 13231                                             1345
      xv=0.25                                                              1345
      goto 13241                                                           1346
13231 continue                                                             1347
13250 do 13251 j=1,ni                                                      1347
      if(ju(j).eq.0)goto 13251                                             1347
      jb=ix(j)                                                             1347
      je=ix(j+1)-1                                                         1348
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1349
13251 continue                                                             1350
13252 continue                                                             1350
13241 continue                                                             1351
13221 continue                                                             1351
      xm(0)=vi                                                             1351
      nlp=0                                                                1351
      nin=nlp                                                              1352
      if(flmin .ge. 1.0)goto 13271                                         1352
      eqs=max(eps,flmin)                                                   1352
      alf=eqs**(1.0/(nlam-1))                                              1352
13271 continue                                                             1353
      m=0                                                                  1353
      mm=0                                                                 1353
      nin=0                                                                1353
      o=0.0                                                                1353
      svr=o                                                                1353
      mnl=min(mnlam,nlam)                                                  1353
      bs=0.0                                                               1354
13280 do 13281 ilm=1,nlam                                                  1355
      if(flmin .lt. 1.0)goto 13301                                         1355
      al=ulam(ilm)                                                         1355
      goto 13291                                                           1356
13301 if(ilm .le. 2)goto 13311                                             1356
      al=al*alf                                                            1356
      goto 13291                                                           1357
13311 if(ilm .ne. 1)goto 13321                                             1357
      al=big                                                               1357
      goto 13331                                                           1358
13321 continue                                                             1358
      al=0.0                                                               1359
13340 do 13341 j=1,ni                                                      1359
      if(ju(j).eq.0)goto 13341                                             1359
      if(vp(j).le.0.0)goto 13341                                           1360
      jb=ix(j)                                                             1360
      je=ix(j+1)-1                                                         1360
      jn=ix(j+1)-ix(j)                                                     1361
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1362
      gj=dot_product(sc(1:jn),x(jb:je))                                    1363
      gj=(gj-svr*xb(j))/xs(j)                                              1364
      al=max(al,abs(gj)/vp(j))                                             1365
13341 continue                                                             1366
13342 continue                                                             1366
      al=alf*al/bta                                                        1367
13331 continue                                                             1368
13291 continue                                                             1368
      al2=al*omb                                                           1368
      al1=al*bta                                                           1368
      nit=0                                                                1369
13350 continue                                                             1369
13351 continue                                                             1369
      bs(0)=b(0)                                                           1369
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1370
13360 continue                                                             1370
13361 continue                                                             1370
      nlp=nlp+1                                                            1370
      dlx=0.0                                                              1371
13370 do 13371 k=1,ni                                                      1371
      if(ju(k).eq.0)goto 13371                                             1372
      jb=ix(k)                                                             1372
      je=ix(k+1)-1                                                         1372
      jn=ix(k+1)-ix(k)                                                     1372
      bk=b(k)                                                              1373
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1374
      gk=dot_product(sc(1:jn),x(jb:je))                                    1375
      gk=(gk-svr*xb(k))/xs(k)                                              1376
      u=gk+xv(k)*b(k)                                                      1376
      au=abs(u)-vp(k)*al1                                                  1377
      if(au .gt. 0.0)goto 13391                                            1377
      b(k)=0.0                                                             1377
      goto 13401                                                           1378
13391 continue                                                             1378
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1378
13401 continue                                                             1379
13381 continue                                                             1379
      d=b(k)-bk                                                            1379
      if(abs(d).le.0.0)goto 13371                                          1379
      dlx=max(dlx,abs(d))                                                  1380
      if(mm(k) .ne. 0)goto 13421                                           1380
      nin=nin+1                                                            1380
      if(nin.gt.nx)goto 13372                                              1381
      mm(k)=nin                                                            1381
      m(nin)=k                                                             1381
      sc(1:jn)=v(jx(jb:je))                                                1382
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1383
13421 continue                                                             1384
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1385
      o=o+d*(xb(k)/xs(k))                                                  1386
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1387
13371 continue                                                             1388
13372 continue                                                             1388
      if(nin.gt.nx)goto 13362                                              1389
      d=svr/xm(0)                                                          1390
      if(d .eq. 0.0)goto 13441                                             1390
      b(0)=b(0)+d                                                          1390
      dlx=max(dlx,abs(d))                                                  1390
      r=r-d*v                                                              1390
13441 continue                                                             1391
      svr=svr-d*xm(0)                                                      1391
      if(dlx.lt.shr)goto 13362                                             1392
13450 continue                                                             1392
13451 continue                                                             1392
      nlp=nlp+1                                                            1392
      dlx=0.0                                                              1393
13460 do 13461 l=1,nin                                                     1393
      k=m(l)                                                               1393
      jb=ix(k)                                                             1393
      je=ix(k+1)-1                                                         1394
      jn=ix(k+1)-ix(k)                                                     1394
      bk=b(k)                                                              1395
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1396
      gk=dot_product(sc(1:jn),x(jb:je))                                    1397
      gk=(gk-svr*xb(k))/xs(k)                                              1398
      u=gk+xv(k)*b(k)                                                      1398
      au=abs(u)-vp(k)*al1                                                  1399
      if(au .gt. 0.0)goto 13481                                            1399
      b(k)=0.0                                                             1399
      goto 13491                                                           1400
13481 continue                                                             1400
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1400
13491 continue                                                             1401
13471 continue                                                             1401
      d=b(k)-bk                                                            1401
      if(abs(d).le.0.0)goto 13461                                          1401
      dlx=max(dlx,abs(d))                                                  1402
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1403
      o=o+d*(xb(k)/xs(k))                                                  1404
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1405
13461 continue                                                             1406
13462 continue                                                             1406
      d=svr/xm(0)                                                          1407
      if(d .eq. 0.0)goto 13511                                             1407
      b(0)=b(0)+d                                                          1407
      dlx=max(dlx,abs(d))                                                  1407
      r=r-d*v                                                              1407
13511 continue                                                             1408
      svr=svr-d*xm(0)                                                      1409
      if(dlx.lt.shr)goto 13452                                             1409
      goto 13451                                                           1410
13452 continue                                                             1410
      goto 13361                                                           1411
13362 continue                                                             1411
      if(nin.gt.nx)goto 13352                                              1412
      if(abs(b(0)-bs(0)) .ge. shr)goto 13531                               1412
      kx=0                                                                 1413
13540 do 13541 j=1,nin                                                     1413
      if(abs(b(m(j))-bs(m(j))).lt.shr)goto 13541                           1413
      kx=1                                                                 1413
      goto 13542                                                           1413
13541 continue                                                             1414
13542 continue                                                             1414
      if(kx.eq.0)goto 13352                                                1415
13531 continue                                                             1416
      sc=b(0)                                                              1416
      b0=0.0                                                               1417
13550 do 13551 j=1,nin                                                     1417
      l=m(j)                                                               1417
      jb=ix(l)                                                             1417
      je=ix(l+1)-1                                                         1418
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1419
      b0=b0-b(l)*xb(l)/xs(l)                                               1420
13551 continue                                                             1421
13552 continue                                                             1421
      sc=sc+b0                                                             1422
13560 do 13561 i=1,no                                                      1422
      fi=sc(i)                                                             1423
      if(fi .ge. fmin)goto 13581                                           1423
      q(i)=0.0                                                             1423
      goto 13571                                                           1423
13581 if(fi .le. fmax)goto 13591                                           1423
      q(i)=1.0                                                             1423
      goto 13601                                                           1424
13591 continue                                                             1424
      q(i)=1.0/(1.0+exp(-fi))                                              1424
13601 continue                                                             1425
13571 continue                                                             1425
13561 continue                                                             1426
13562 continue                                                             1426
      v=w*q*(1.0-q)                                                        1426
      r=w*(y-q)                                                            1427
      xm(0)=sum(v)                                                         1427
      if(xm(0).lt.vmin)goto 13352                                          1427
      svr=sum(r)                                                           1427
      o=0.0                                                                1428
13610 do 13611 l=1,nin                                                     1428
      j=m(l)                                                               1429
      jb=ix(j)                                                             1429
      je=ix(j+1)-1                                                         1429
      jn=ix(j+1)-ix(j)                                                     1430
      sc(1:jn)=v(jx(jb:je))                                                1431
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1432
      if(kopt .ne. 0)goto 13631                                            1433
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1434
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1436
13631 continue                                                             1437
13611 continue                                                             1438
13612 continue                                                             1438
      nit=nit+1                                                            1438
      if(nit .le. maxit)goto 13651                                         1438
      jerr=-ilm                                                            1438
      return                                                               1438
13651 continue                                                             1439
      goto 13351                                                           1440
13352 continue                                                             1440
      if(nin .le. nx)goto 13671                                            1440
      jerr=-10000-ilm                                                      1440
      goto 13282                                                           1440
13671 continue                                                             1441
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1441
      kin(ilm)=nin                                                         1442
      a0(ilm)=b(0)                                                         1442
      alm(ilm)=al                                                          1442
      lmu=ilm                                                              1443
      devi=dev2(no,w,y,q,pmin)                                             1444
      dev(ilm)=(dev1-devi)/dev1                                            1445
      if(ilm.lt.mnl)goto 13281                                             1445
      if(flmin.ge.1.0)goto 13281                                           1446
      me=0                                                                 1446
13680 do 13681 j=1,nin                                                     1446
      if(a(j,ilm).ne.0.0) me=me+1                                          1446
13681 continue                                                             1446
13682 continue                                                             1446
      if(me.gt.ne)goto 13282                                               1447
      if(dev(ilm).gt.devmax)goto 13282                                     1447
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13282                             1448
      if(xm(0).lt.vmin)goto 13282                                          1449
13281 continue                                                             1450
13282 continue                                                             1450
      deallocate(xm,b,bs,v,r,sc,xv,q,mm)                                   1451
      return                                                               1452
      end                                                                  1453
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,w,ju,vp,ne,nx,nlam,f   1455 
     *lmin,ulam,  shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev,alm,nlp,je
     *rr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1457 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)             1458
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1459
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1460
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: sc,xm,v,r                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1471
      jerr=jerr+ierr                                                       1472
      allocate(r(1:no),stat=ierr)                                          1472
      jerr=jerr+ierr                                                       1473
      allocate(v(1:no),stat=ierr)                                          1473
      jerr=jerr+ierr                                                       1474
      allocate(mm(1:ni),stat=ierr)                                         1474
      jerr=jerr+ierr                                                       1475
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1475
      jerr=jerr+ierr                                                       1476
      allocate(sxp(1:no),stat=ierr)                                        1476
      jerr=jerr+ierr                                                       1477
      allocate(sc(1:no),stat=ierr)                                         1477
      jerr=jerr+ierr                                                       1478
      if(jerr.ne.0) return                                                 1479
      pmax=1.0-pmin                                                        1479
      emin=pmin/pmax                                                       1479
      emax=1.0/emin                                                        1480
      pfm=(1.0+pmin)*pmin                                                  1480
      pfx=(1.0-pmin)*pmax                                                  1480
      vmin=pfm*pmax                                                        1481
      bta=max(parm,1.0e-3)                                                 1481
      omb=1.0-bta                                                          1481
      dev1=0.0                                                             1482
13690 do 13691 ic=1,nc                                                     1482
      q0=dot_product(w,y(:,ic))                                            1483
      if(q0 .gt. pmin)goto 13711                                           1483
      jerr =8000+ic                                                        1483
      return                                                               1483
13711 continue                                                             1484
      if(q0 .lt. 1.0-pmin)goto 13731                                       1484
      jerr =9000+ic                                                        1484
      return                                                               1484
13731 continue                                                             1485
      vi=q0*(1.0-q0)                                                       1485
      v=vi*w                                                               1485
      b(1:ni,ic)=0.0                                                       1486
      b(0,ic)=log(q0)                                                      1486
      xm(0)=vi                                                             1486
      dev1=dev1-q0*b(0,ic)                                                 1487
13691 continue                                                             1488
13692 continue                                                             1488
      if(isd .le. 0)goto 13751                                             1488
      xv=0.25                                                              1488
      goto 13761                                                           1489
13751 continue                                                             1490
13770 do 13771 j=1,ni                                                      1490
      if(ju(j).eq.0)goto 13771                                             1490
      jb=ix(j)                                                             1490
      je=ix(j+1)-1                                                         1491
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1492
13771 continue                                                             1493
13772 continue                                                             1493
13761 continue                                                             1494
13741 continue                                                             1494
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1494
      sxp=0.0                                                              1495
13780 do 13781 ic=1,nc                                                     1495
      q(:,ic)=exp(b(0,ic))                                                 1495
      sxp=sxp+q(:,ic)                                                      1495
13781 continue                                                             1496
13782 continue                                                             1496
      if(flmin .ge. 1.0)goto 13801                                         1496
      eqs=max(eps,flmin)                                                   1496
      alf=eqs**(1.0/(nlam-1))                                              1496
13801 continue                                                             1497
      m=0                                                                  1497
      mm=0                                                                 1497
      nin=0                                                                1497
      nlp=0                                                                1497
      mnl=min(mnlam,nlam)                                                  1497
      bs=0.0                                                               1497
      svr=0.0                                                              1497
      o=0.0                                                                1498
13810 do 13811 ilm=1,nlam                                                  1499
      if(flmin .lt. 1.0)goto 13831                                         1499
      al=ulam(ilm)                                                         1499
      goto 13821                                                           1500
13831 if(ilm .le. 2)goto 13841                                             1500
      al=al*alf                                                            1500
      goto 13821                                                           1501
13841 if(ilm .ne. 1)goto 13851                                             1501
      al=big                                                               1501
      goto 13861                                                           1502
13851 continue                                                             1502
      al=0.0                                                               1503
13870 do 13871 ic=1,nc                                                     1503
      v=q(:,ic)/sxp                                                        1503
      r=w*(y(:,ic)-v)                                                      1503
      v=w*v*(1.0-v)                                                        1504
13880 do 13881 j=1,ni                                                      1504
      if(ju(j).eq.0)goto 13881                                             1504
      if(vp(j).le.0.0)goto 13881                                           1505
      jb=ix(j)                                                             1505
      je=ix(j+1)-1                                                         1505
      jn=ix(j+1)-ix(j)                                                     1506
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1507
      gj=dot_product(sc(1:jn),x(jb:je))                                    1508
      gj=(gj-svr*xb(j))/xs(j)                                              1509
      al=max(al,abs(gj)/vp(j))                                             1510
13881 continue                                                             1511
13882 continue                                                             1511
13871 continue                                                             1512
13872 continue                                                             1512
      al=alf*al/bta                                                        1513
13861 continue                                                             1514
13821 continue                                                             1514
      al2=al*omb                                                           1514
      al1=al*bta                                                           1514
      nit=0                                                                1515
13890 continue                                                             1515
13891 continue                                                             1515
      ixx=0                                                                1515
      jxx=ixx                                                              1515
      ig=0                                                                 1516
13900 do 13901 ic=1,nc                                                     1516
      bs(0,ic)=b(0,ic)                                                     1517
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1518
      xm(0)=0.0                                                            1518
      svr=0.0                                                              1518
      o=0.0                                                                1519
13910 do 13911 i=1,no                                                      1519
      pic=q(i,ic)/sxp(i)                                                   1520
      if(pic .ge. pfm)goto 13931                                           1520
      pic=0.0                                                              1520
      v(i)=0.0                                                             1520
      goto 13921                                                           1521
13931 if(pic .le. pfx)goto 13941                                           1521
      pic=1.0                                                              1521
      v(i)=0.0                                                             1521
      goto 13951                                                           1522
13941 continue                                                             1522
      v(i)=w(i)*pic*(1.0-pic)                                              1522
      xm(0)=xm(0)+v(i)                                                     1522
13951 continue                                                             1523
13921 continue                                                             1523
      r(i)=w(i)*(y(i,ic)-pic)                                              1523
      svr=svr+r(i)                                                         1524
13911 continue                                                             1525
13912 continue                                                             1525
      if(xm(0).le.vmin)goto 13901                                          1525
      ig=1                                                                 1526
13960 do 13961 l=1,nin                                                     1526
      j=m(l)                                                               1527
      jb=ix(j)                                                             1527
      je=ix(j+1)-1                                                         1528
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1529
      if(kopt .ne. 0)goto 13981                                            1530
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1531
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1532
13981 continue                                                             1533
13961 continue                                                             1534
13962 continue                                                             1534
13990 continue                                                             1534
13991 continue                                                             1534
      nlp=nlp+1                                                            1534
      dlx=0.0                                                              1535
14000 do 14001 k=1,ni                                                      1535
      if(ju(k).eq.0)goto 14001                                             1536
      jb=ix(k)                                                             1536
      je=ix(k+1)-1                                                         1536
      jn=ix(k+1)-ix(k)                                                     1536
      bk=b(k,ic)                                                           1537
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1538
      gk=dot_product(sc(1:jn),x(jb:je))                                    1539
      gk=(gk-svr*xb(k))/xs(k)                                              1540
      u=gk+xv(k,ic)*b(k,ic)                                                1540
      au=abs(u)-vp(k)*al1                                                  1541
      if(au .gt. 0.0)goto 14021                                            1541
      b(k,ic)=0.0                                                          1541
      goto 14031                                                           1542
14021 continue                                                             1542
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1542
14031 continue                                                             1543
14011 continue                                                             1543
      d=b(k,ic)-bk                                                         1543
      if(abs(d).le.0.0)goto 14001                                          1543
      dlx=max(dlx,abs(d))                                                  1544
      if(mm(k) .ne. 0)goto 14051                                           1544
      nin=nin+1                                                            1545
      if(nin .le. nx)goto 14071                                            1545
      jxx=1                                                                1545
      goto 14002                                                           1545
14071 continue                                                             1546
      mm(k)=nin                                                            1546
      m(nin)=k                                                             1547
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             1548
14051 continue                                                             1549
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1550
      o=o+d*(xb(k)/xs(k))                                                  1551
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1552
14001 continue                                                             1553
14002 continue                                                             1553
      if(jxx.gt.0)goto 13992                                               1554
      d=svr/xm(0)                                                          1555
      if(d .eq. 0.0)goto 14091                                             1555
      b(0,ic)=b(0,ic)+d                                                    1555
      dlx=max(dlx,abs(d))                                                  1556
      r=r-d*v                                                              1556
      svr=svr-d*xm(0)                                                      1557
14091 continue                                                             1558
      if(dlx.lt.shr)goto 13992                                             1559
14100 continue                                                             1559
14101 continue                                                             1559
      nlp=nlp+1                                                            1559
      dlx=0.0                                                              1560
14110 do 14111 l=1,nin                                                     1560
      k=m(l)                                                               1560
      jb=ix(k)                                                             1560
      je=ix(k+1)-1                                                         1561
      jn=ix(k+1)-ix(k)                                                     1561
      bk=b(k,ic)                                                           1562
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1563
      gk=dot_product(sc(1:jn),x(jb:je))                                    1564
      gk=(gk-svr*xb(k))/xs(k)                                              1565
      u=gk+xv(k,ic)*b(k,ic)                                                1565
      au=abs(u)-vp(k)*al1                                                  1566
      if(au .gt. 0.0)goto 14131                                            1566
      b(k,ic)=0.0                                                          1566
      goto 14141                                                           1567
14131 continue                                                             1567
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1567
14141 continue                                                             1568
14121 continue                                                             1568
      d=b(k,ic)-bk                                                         1568
      if(abs(d).le.0.0)goto 14111                                          1569
      dlx=max(dlx,abs(d))                                                  1570
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1571
      o=o+d*(xb(k)/xs(k))                                                  1572
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1573
14111 continue                                                             1574
14112 continue                                                             1574
      d=svr/xm(0)                                                          1575
      if(d .eq. 0.0)goto 14161                                             1575
      b(0,ic)=b(0,ic)+d                                                    1575
      dlx=max(dlx,abs(d))                                                  1576
      r=r-d*v                                                              1576
      svr=svr-d*xm(0)                                                      1577
14161 continue                                                             1578
      if(dlx.lt.shr)goto 14102                                             1578
      goto 14101                                                           1579
14102 continue                                                             1579
      goto 13991                                                           1580
13992 continue                                                             1580
      if(jxx.gt.0)goto 13902                                               1581
      if(abs(b(0,ic)-bs(0,ic)).gt.shr) ixx=1                               1582
      if(ixx .ne. 0)goto 14181                                             1583
14190 do 14191 j=1,nin                                                     1584
      if(abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 14211                   1584
      ixx=1                                                                1584
      goto 14192                                                           1584
14211 continue                                                             1585
14191 continue                                                             1586
14192 continue                                                             1586
14181 continue                                                             1587
      sc=b(0,ic)                                                           1587
      b0=0.0                                                               1588
14220 do 14221 j=1,nin                                                     1588
      l=m(j)                                                               1588
      jb=ix(l)                                                             1588
      je=ix(l+1)-1                                                         1589
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   1590
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            1591
14221 continue                                                             1592
14222 continue                                                             1592
      sc=min(max(exmn,sc+b0),exmx)                                         1593
      sxp=sxp-q(:,ic)                                                      1594
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)                    1595
      sxp=sxp+q(:,ic)                                                      1596
13901 continue                                                             1597
13902 continue                                                             1597
      if(jxx.gt.0)goto 13892                                               1597
      if(ixx.eq.0)goto 13892                                               1597
      if(ig.eq.0)goto 13892                                                1598
      s=-sum(b(0,:))/nc                                                    1598
      b(0,:)=b(0,:)+s                                                      1598
      sc=s                                                                 1598
      b0=0.0                                                               1599
14230 do 14231 j=1,nin                                                     1599
      l=m(j)                                                               1600
      if(vp(l) .gt. 0.0)goto 14251                                         1600
      s=sum(b(l,:))/nc                                                     1600
      goto 14261                                                           1601
14251 continue                                                             1601
      s=elc(parm,nc,b(l,:),is)                                             1601
14261 continue                                                             1602
14241 continue                                                             1602
      b(l,:)=b(l,:)-s                                                      1603
      jb=ix(l)                                                             1603
      je=ix(l+1)-1                                                         1604
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         1605
      b0=b0+s*xb(l)/xs(l)                                                  1606
14231 continue                                                             1607
14232 continue                                                             1607
      sc=sc+b0                                                             1607
      sc=exp(sc)                                                           1607
      sxp=sxp*sc                                                           1607
14270 do 14271 ic=1,nc                                                     1607
      q(:,ic)=q(:,ic)*sc                                                   1607
14271 continue                                                             1608
14272 continue                                                             1608
      nit=nit+1                                                            1608
      if(nit .le. maxit)goto 14291                                         1608
      jerr=-ilm                                                            1608
      return                                                               1608
14291 continue                                                             1609
      goto 13891                                                           1610
13892 continue                                                             1610
      if(jxx .le. 0)goto 14311                                             1610
      jerr=-10000-ilm                                                      1610
      goto 13812                                                           1610
14311 continue                                                             1610
      devi=0.0                                                             1611
14320 do 14321 ic=1,nc                                                     1612
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1612
      a0(ic,ilm)=b(0,ic)                                                   1613
14330 do 14331 i=1,no                                                      1613
      if(y(i,ic).le.0.0)goto 14331                                         1614
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1615
14331 continue                                                             1616
14332 continue                                                             1616
14321 continue                                                             1617
14322 continue                                                             1617
      kin(ilm)=nin                                                         1617
      alm(ilm)=al                                                          1617
      lmu=ilm                                                              1618
      dev(ilm)=(dev1-devi)/dev1                                            1618
      if(ig.eq.0)goto 13812                                                1619
      if(ilm.lt.mnl)goto 13811                                             1619
      if(flmin.ge.1.0)goto 13811                                           1620
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13812             1621
      if(dev(ilm).gt.devmax)goto 13812                                     1621
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13812                             1622
13811 continue                                                             1623
13812 continue                                                             1623
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc)                            1624
      return                                                               1625
      end                                                                  1626
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  1627
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   1627
      integer ia(*),ix(*),jx(*)                                            1628
14340 do 14341 ic=1,nc                                                     1628
      f(ic,:)=a0(ic)                                                       1628
14341 continue                                                             1629
14342 continue                                                             1629
14350 do 14351 j=1,nin                                                     1629
      k=ia(j)                                                              1629
      kb=ix(k)                                                             1629
      ke=ix(k+1)-1                                                         1630
14360 do 14361 ic=1,nc                                                     1630
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    1630
14361 continue                                                             1631
14362 continue                                                             1631
14351 continue                                                             1632
14352 continue                                                             1632
      return                                                               1633
      end                                                                  1634
      subroutine psort7 (v,a,ii,jj)                                             
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      real v                                                                    
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end                                                                       