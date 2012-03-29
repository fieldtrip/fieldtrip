!  $6 21/03/2003  Anwander A.  released for SIMBIO
!  $5 29/01/03    aa mxkn1d changed to be used with meg sensors
!  $4 26/07/00    aa mxelcs added as const of a large estimation for the max. 
!                    nb. of elements per linefor the symmetric compact storage 
!  $3 18/07/00    aa qfsize removed as global variable
!  $2 12/07/00    aa new variables for neurofem extention
!  $1 30/06/00    aa new const for size of the compact storrage

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem.inc :      	 include parameters in the neurofem files
!                              -------------------
!     begin                : Mai 2000
!
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!    NeuroFEM license:
!    =================
!    Copyright (c) 2007 by 
!    Dr.Carsten Wolters, Dr.Alfred Anwander,  
!    Prof.Dr.H.Buchner, Prof.Dr.G.Knoll, Dr.Adrian Rienaecker, 
!    Rainer Beckmann. 
!
!    Permission is hereby granted, free of charge, to any person 
!    obtaining a copy of the NeuroFEM software and associated 
!    documentation files (the "Software"), to deal in the Software 
!    without restrictions, including without limitations the rights to 
!    use, copy, modify, merge, publish, distribute, sublicense, 
!    and/or sell copies of the Software, and to permit persons to 
!    whom the Software is furnished to do so, subject to the following 
!    conditions:
!
!    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
!    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
!    OR OTHER DEALINGS IN THE SOFTWARE. THE ENTIRE RISK AS TO THE 
!    QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH YOU.
!
!    The above copyright notice, permission and warranty notices 
!    shall be included in all copies of the Software. A copyright 
!    notice of the contributing authors and the above permission 
!    and warranty notices shall be included in all source code 
!    files of the Software. 
!
!---->---1---------2---------3---------4---------5---------6---------7--<

      implicit double precision(a-h,o-z)
      real ran2
!---for linear    elements:  parameter (mxkn1d=2,mxkn2d= 4,mxkn3d= 8)
!---for quadratic elements:  parameter (mxkn1d=3,mxkn2d= 8,mxkn3d=20)
!---for cubic     elements:  parameter (mxkn1d=4,mxkn2d=12,mxkn3d=32)
!     isoparametric meg sensor elements uses 3 nodes
      parameter (mxkn1d=3,mxkn2d= 4,mxkn3d=8)
!---for compact storrage(cs), max. nb. of elements per line in symmetric cs
      parameter (mxelcs=30)
      parameter (maxdim=3,mxmxkn=32)
!%%% maxstr determines the start position of keywords
!%%% not yet active 
      parameter (maxfil=50,maxfm1=maxfil-1,maxstr=30)
      parameter (maxlei=10,mvagam=6,mxloop=4)
      parameter (pi=3.141592653589793d0,pi4mue=1.d-10)
      parameter (ftconv=1.d9,flconv=1.d6)
      parameter ( z0= 0.d0, z1= 1.d0, z2= 2.d0, z3= 3.d0, z4= 4.d0,     &
     &            z5= 5.d0, z6= 6.d0, z7= 7.d0, z8= 8.d0, z9= 9.d0,     &
     &           z10=10.d0,z11=11.d0,z12=12.d0,z13=13.d0,z14=14.d0,     &
     &         z50=50.d0,z100=100.d0,z180=180.d0 )
      parameter ( f2=0.5d0,f4=0.25d0,f5=0.2d0,f12=0.08333333333333333d0)
      parameter ( tol=1.d-3,wuzeps=3.d-8,bignum=1.d100 )
      parameter ( tolaus=1.d-15,tolmeg=1.d-7)
!%%% tmpfkm: minimum value of the decreasing factor (Simulated Annealing)
      parameter ( tmpfkm=0.8d0 )
!
      character*80 zeil  ,filnam,  filakt,chifrm
      character*20 cabb
      character*1  cartax,chavol,chasur,chasvu,chasvv
      character date*8,time*8,lognam*80
      character*6 chapre,chapos,chaset
      integer qfreeu,qfrecl,qclen
!%%%18/07/00aa      integer qfsize      
      logical qfilda
      logical lrelin,lreint,lmetro,ldummy
      logical linmat,lrefsl,linver,lsimsr,ldipol,                       &
     &        lgobac,lalpha,lsurko,lkonve,                              &
     &        ldipan,lnodlo,linfck,ldirsl,lnorco,                       &
     &        ltrans,ldevsc,lsiman,lostop,lmncor,                       &
     &        lperor,lnorch,lsurvl,leigen,lreduk,                       &
     &        lmofwd,lmocgs,liomem,lfocus,llanzp,                       &
     &        lconsi,lcholy,lqrfak,lqrfax,lsngvl,                       &
     &        linsvd,ldepsc,lcovar,lgeorg,logmeg,                       &
     &        lcglin,lentro,ll1nrm,ll2nrm,llucut,                       &
     &        lognoi,llamit,logneu,lrango,logsur,                       &
     &        laniso,lchisl,lfsolv,lwrani
!
      common /file1/ filnam(maxfil)
      common /imelb/ cabb,labb,nttast,ntterm
      common /foflb/ lrelin,lreint
!
      common /cha1/ cartax(3),chavol,chasur
      common /geo1/ npogeo,nelgeo,ndmgeo,lengeo,numper,                 &
     &              maxnei,lenasy,ndmflx,ndimmt,ndirhs,                 &
     &              ndimat,ndmsrc,norder,numana,lenflf,                 &
     &              naniso
      common /geo2/ xyzmax(maxdim),xyzmin(maxdim),                      &
     &              xyzdst(maxdim),xyzcgr(maxdim),                      &
     &              xyzmom(maxdim),dipmax(maxdim),                      &
     &              dstmax
      common /geo3/ surfac,volume,weight
      common /fla1/ npoflf,nelflf,ndmflf,numalf,nelsur
      common /meg1/ npomeg,nelmeg,ndmmeg,nummeg,numeeg,                 &
     &              mgagam,numtyp,numsys,melmeg
      common /tol1/ isolvr,intgrd,iglatt,                               &
     &              igladi
      common /tol2/ tolsol,tolinv,resiun
      common /inpu/ permeb,epsilo,densty
      common /log1/ linmat,lrefsl,linver,lsimsr,ldipol,                 &
     &              lgobac,lalpha,lsurko,lkonve,                        &
     &              ldipan,lnodlo,linfck,ldirsl,lnorco,                 &
     &              ltrans,ldevsc,lsiman,lostop,lmncor,                 &
     &              lperor,lnorch,lsurvl,leigen,lreduk,                 &
     &              lmofwd,lmocgs,liomem,lfocus,llanzp,                 &
     &              lconsi,lcholy,lqrfak,lqrfax,lsngvl,                 &
     &              linsvd,ldepsc,lcovar,lgeorg,logmeg,                 &
     &              lcglin,lentro,ll1nrm,ll2nrm,llucut,                 &
     &              lognoi,llamit,logneu,lrango,logsur,                 &
     &              laniso,lchisl,lwrani
      common /inf1/ ninfpo,nevlpo,ncorno
      common /inv1/ leninv,npoinv,nelinv,mxknin,ndminv,                 &
     &              invway,numtim,numinc,incsta,incend,                 &
     &              nblopt,numpvt,numsar,lawork,iunco
      common /inv2/ xlmbda,xlmbfi,xlmbin,diplam,                        &
     &              trasmo,tradep,tradat, xlmdn, xlmup 
      common /inv3/ potsqr,potave,potcor
      common /dip1/ disdip,radius,dismax,scadip,thrdip
      common /anli/ iseed ,numdip,nmanmt,nmanrs,maxtry,                 &
     &              lsqway,numsol,itim
      common /anld/ tmpsta,temper,tmpfak,avepot,                        &
     &              epsnrm,potnrm,defect,signoi
!
      common /laso/ nval,nfig,nperm,nmval,nmvec,nblock,                 &
     &              maxop,maxj,nvala,nuwrk,nmodes,iunio,                &
     &              nstore,nustei,nwriln
      common /lnzc/ lenval,lenlan,lenrow,lencol,                        &
     &              ipar(20),ldylan,inerti,nstpsh,                      &
     &              numvec,numval
      common /dlan/ rpar(5),dshift,relacc
      common /nsvd/ nwrksv
      common /loui/ rpostl,funcut
      common /fuxx/ chisqr,chimax,chimin,chivar,chival
      common /leii/ ianfle(maxlei),iendle(maxlei),numlei
      common /leid/ wertle(maxlei,6),volule(maxlei),powele(maxlei)
      common /fsol/ lfsolv,ifsolv

