!  $6 21/10/2005  Anwander A.  changed ndmgeo line 1310/1314 to 3
!  $5 21/03/2003  Anwander A.  released for SIMBIO
!  $4 12/07/00    aa  writing out equation system for fast solver methods
!  $3 30/06/00    aa  dwork and ndrest not used, removed from simtim
!  $2 30/06/00    aa  multiple return in do loops removed
!  $1 30/06/00    aa  the field xlk is allocated with size 3*32

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem1.f : subroutines for preprocessing
!                              -------------------
!     begin                : Mai 2000
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

      subroutine reafil
      include 'neurofem.inc'
	 
!     
!---reads filenames from filename-file=filnam(1)
      if(.not.qfilda(filnam( 1))) write(*,910) 'filnam( 1)= ',          &
     &     filnam( 1)(1:qclen(filnam( 1)))
!
      iunit=qfreeu()
      call qfopse(iunit,filnam(1),'un','fo',ierr)
      if (ierr.gt.0) then
         write (*,900) 'fehler filnam(1) in reafil'
         stop
      end if
 10   read(iunit,900,end=50) zeil
!
!---files 1-7 are always used
      call chsuch(zeil,'geometriefile'       ,13,30,filnam( 2),n,*10)
      call chsuch(zeil,'geometry file'       ,13,30,filnam( 2),n,*10)
      call chsuch(zeil,'cauchy-inputfile'    ,16,30,filnam( 3),n,*10)
      call chsuch(zeil,'cauchy input file'   ,17,30,filnam( 3),n,*10)
      call chsuch(zeil,'toleranzfile'        ,12,30,filnam( 4),n,*10)
      call chsuch(zeil,'cauchy tolerance file',21,30,filnam( 4),n,*10)
      call chsuch(zeil,'leitfaehigkeiten'    ,16,30,filnam( 5),n,*10)
      call chsuch(zeil,'conductivity'        ,12,30,filnam( 5),n,*10)
      call chsuch(zeil,'randbedingungsfile'  ,18,30,filnam( 6),n,*10)
      call chsuch(zeil,'boundary condition'  ,18,30,filnam( 6),n,*10)
      call chsuch(zeil,'cortexgeometrie'     ,15,30,filnam( 7),n,*10)
      call chsuch(zeil,'cortex geometry'     ,15,30,filnam( 7),n,*10)
!
!---file 18 in case of heat transfer coefficients
      call chsuch(zeil,'alpha-theta-u-file'  ,18,30,filnam(18),n,*10)
      call chsuch(zeil,'alfa teta u file'    ,16,30,filnam(18),n,*10)
!
!---for reference solution (forward)
      call chsuch(zeil,'quellinputfile'      ,14,30,filnam( 8),n,*10)
      call chsuch(zeil,'source input file'   ,17,30,filnam( 8),n,*10)
      call chsuch(zeil,'eeg-messwertfile'    ,16,30,filnam( 9),n,*10)
      call chsuch(zeil,'eeg measurement'     ,15,30,filnam( 9),n,*10)
      call chsuch(zeil,'referenzergebnis'    ,16,30,filnam(10),n,*10)
      call chsuch(zeil,'reference solution'  ,18,30,filnam(10),n,*10)
!
!
!---for influence matrix
      call chsuch(zeil,'einflussknotenfile'  ,18,30,filnam(11),n,*10)
      call chsuch(zeil,'influence nodes'     ,15,30,filnam(11),n,*10)
      call chsuch(zeil,'auswerteknotenfile'  ,18,30,filnam(12),n,*10)
      call chsuch(zeil,'evaluation nodes'    ,16,30,filnam(12),n,*10)
!
!---for influence matrix and inverse solution
      call chsuch(zeil,'einflussmatrixfile'  ,18,30,filnam(13),n,*10)
      call chsuch(zeil,'influence matrix'    ,16,30,filnam(13),n,*10)
!---for eeg inverse reconstruction
      call chsuch(zeil,'eeg-rauschen'        ,12,30,filnam(14),n,*10)
      call chsuch(zeil,'eeg noise file'      ,14,30,filnam(14),n,*10)
!
!---for inverse solution
      call chsuch(zeil,'inverses-ergebnis'   ,17,30,filnam(15),n,*10)
      call chsuch(zeil,'inverse result'      ,14,30,filnam(15),n,*10)
      call chsuch(zeil,'quelloutputfile'     ,15,30,filnam(16),n,*10)
      call chsuch(zeil,'source output file'  ,18,30,filnam(16),n,*10)
!
!---for source simulation
      call chsuch(zeil,'simulationsinputfile',20,30,filnam(17),n,*10)
      call chsuch(zeil,'source simulation'   ,17,30,filnam(17),n,*10)
!
!---convergence output cg-solution
      call chsuch(zeil,'konvergenzoutputfile',20,30,filnam(19),n,*10)
      call chsuch(zeil,'convergence output'  ,18,30,filnam(19),n,*10)
!
!---transient time increments
      call chsuch(zeil,'zeitinkremente-file' ,19,30,filnam(22),n,*10)
      call chsuch(zeil,'time increments'     ,15,30,filnam(22),n,*10)
!
!---nodal order file
      call chsuch(zeil,'knoten-reihenfolge'  ,18,30,filnam(23),n,*10)
      call chsuch(zeil,'order of nodes'      ,14,30,filnam(23),n,*10)
!
!---curry--dipole
      call chsuch(zeil,'curry-dipole-format' ,19,30,filnam(24),n,*10)
      call chsuch(zeil,'curry dipole file'   ,17,30,filnam(24),n,*10)
!
!---memory management
      call chsuch(zeil,'speicherbelegung'    ,16,30,filnam(25),n,*10)
      call chsuch(zeil,'memory allocation'   ,17,30,filnam(25),n,*10)
!
!---annealing information
      call chsuch(zeil,'annealing-info'      ,14,30,filnam(26),n,*10)
      call chsuch(zeil,'annealing info'      ,14,30,filnam(26),n,*10)
!
!---surface file for correlation between potentials
      call chsuch(zeil,'korrelationsknoten'  ,18,30,filnam(27),n,*10)
      call chsuch(zeil,'correlation nodes'   ,17,30,filnam(27),n,*10)
!
!---reference surface results
      call chsuch(zeil,'oberflaechen-ergebnis-ref'  ,25,30,             &
     &     filnam(28),n,*10)
      call chsuch(zeil,'ref surface results'        ,19,30,             &
     &     filnam(28),n,*10)
      call chsuch(zeil,'oberflaechen-ergebnis-inv'  ,25,30,             &
     &     filnam(29),n,*10)
      call chsuch(zeil,'inv surface results'        ,19,30,             &
     &     filnam(29),n,*10)
      call chsuch(zeil,'normierte-eigenvektoren'    ,23,30,             &
     &     filnam(30),n,*10)
      call chsuch(zeil,'normed eigenvectors'        ,19,30,             &
     &     filnam(30),n,*10)
      call chsuch(zeil,'information-eigenwerte'     ,22,30,             &
     &     filnam(31),n,*10)
      call chsuch(zeil,'eigenvalue information'     ,22,30,             &
     &     filnam(31),n,*10)
      call chsuch(zeil,'eigenvektoren-binaer'       ,20,30,             &
     &     filnam(32),n,*10)
      call chsuch(zeil,'binary eigenvectors'        ,19,30,             &
     &     filnam(32),n,*10)
      call chsuch(zeil,'singulaerwertzerlegung'     ,22,30,             &
     &     filnam(33),n,*10)
      call chsuch(zeil,'singular value decomposition',28,30,            &
     &     filnam(33),n,*10)
      call chsuch(zeil,'svd-rechte-vektoren'        ,19,30,             &
     &     filnam(34),n,*10)
      call chsuch(zeil,'svd right side vectors'     ,22,30,             &
     &     filnam(34),n,*10)
!
!---meg-files
      call chsuch(zeil,'meg-geometrie'              ,13,30,             &
     &     filnam(35),n,*10)
      call chsuch(zeil,'meg geometry'               ,12,30,             &
     &     filnam(35),n,*10)
      call chsuch(zeil,'meg-system'                 ,10,30,             &
     &     filnam(36),n,*10)
      call chsuch(zeil,'meg system'                 ,10,30,             &
     &     filnam(36),n,*10)
      call chsuch(zeil,'meg-kanaele'                ,11,30,             &
     &     filnam(37),n,*10)
      call chsuch(zeil,'meg channels'               ,12,30,             &
     &     filnam(37),n,*10)
      call chsuch(zeil,'meg-messwerte'              ,13,30,             &
     &     filnam(38),n,*10)
      call chsuch(zeil,'meg measurement'            ,15,30,             &
     &     filnam(38),n,*10)
      call chsuch(zeil,'meg-rauschen'               ,12,30,             &
     &     filnam(39),n,*10)
      call chsuch(zeil,'meg noise file'             ,14,30,             &
     &     filnam(39),n,*10)
      call chsuch(zeil,'meg-secondary-flux'         ,18,30,             &
     &     filnam(40),n,*10)
      call chsuch(zeil,'meg secondary flux'         ,18,30,             &
     &     filnam(40),n,*10)
      call chsuch(zeil,'meg-sys-transform'          ,17,30,             &
     &     filnam(41),n,*10)
      call chsuch(zeil,'meg sys transform'          ,17,30,             &
     &     filnam(41),n,*10)
!---surface-file
      call chsuch(zeil,'oberflaechenfile'           ,16,30,             &
     &     filnam(42),n,*10)
      call chsuch(zeil,'file of surfaces'           ,16,30,             &
     &     filnam(42),n,*10)
!
!---C.Wolters, 8.01.2001
!---Writing out FE-equation system for PILUTS-solver
      call chsuch(zeil,'piluts gleichungssystem'    ,23,30,             &
     &     filnam(43),n,*10)
      call chsuch(zeil,'piluts equation system'     ,22,30,             &
     &     filnam(43),n,*10)
!
!---C.Wolters, 8.01.2001
!---PEBBLES input-file for AMG-CG-control
      call chsuch(zeil,'pebbles inp-file'           ,16,30,             &
     &     filnam(44),n,*10)
!
!---C.Wolters, 7.7.2000
!---Writing out tensor-valued conductivity file for sphere simulation 
!---studies 
      call chsuch(zeil,'anisotrope Kugelsimulation' ,26,30,             &
     &     filnam(45),n,*10)
      call chsuch(zeil,'aniso-sphere simulation'    ,23,30,             &
     &     filnam(45),n,*10)
      goto 10
   50 call qfclos(iunit,0)
!
      return
 900  format(a)
 910  format('Attention! File: ',a,a, ' is missing!')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reainp
      include 'neurofem.inc'
!
!---reads input file
      iunit=qfreeu()
      call qfopse(iunit,filnam(3),'un','fo',ierr)
      if (ierr.ne.0) stop 'error reainp'
!
   10 read (iunit,900,end=20) zeil
      call qctolo(zeil)
!
      call lsuch(zeil,'dipolquellen'                ,12,30,ldipol,*10)
      call lsuch(zeil,'dipole sources'              ,14,30,ldipol,*10)
      call lsuch(zeil,'dipol-analytisch'            ,16,30,ldipan,*10)
      call lsuch(zeil,'analytical dipole'           ,17,30,ldipan,*10)
      call lsuch(zeil,'mittelwertkorrektur'         ,19,30,lmncor,*10)
      call lsuch(zeil,'average correction'          ,18,30,lmncor,*10)
      call lsuch(zeil,'einflussmatrizen'            ,16,30,linmat,*10)
      call lsuch(zeil,'lead field matrix'           ,17,30,linmat,*10)
      call lsuch(zeil,'singulaerwertzerlegung'      ,22,30,linsvd,*10)
      call lsuch(zeil,'singular value decomposition',28,30,linsvd,*10)
      call lsuch(zeil,'eigenwertbestimmung'         ,19,30,leigen,*10)
      call lsuch(zeil,'determine eigenvalue'        ,20,30,leigen,*10)
      call lsuch(zeil,'eigenwert-lanzpack'          ,18,30,llanzp,*10)
      call lsuch(zeil,'eigenvalue lanzpack'         ,19,30,llanzp,*10)
      call lsuch(zeil,'lanzpack-konsistent'         ,19,30,lconsi,*10)
      call lsuch(zeil,'consistent lanzpack'         ,19,30,lconsi,*10)
      call lsuch(zeil,'eigenwert-io-im-memory'      ,22,30,liomem,*10)
      call lsuch(zeil,'eigenvalue in memory'        ,20,30,liomem,*10)
      call lsuch(zeil,'eigenwertreduktion'          ,18,30,lreduk,*10)
      call lsuch(zeil,'eigenvalue reduction'        ,20,30,lreduk,*10)
      call lsuch(zeil,'modale-loesung'              ,14,30,lmofwd,*10)
      call lsuch(zeil,'modal solution'              ,14,30,lmofwd,*10)
      call lsuch(zeil,'modale-cg-startvektoren'     ,23,30,lmocgs,*10)
      call lsuch(zeil,'modal cg start vectors'      ,22,30,lmocgs,*10)
      call lsuch(zeil,'quellsimulation'             ,15,30,lsimsr,*10)
      call lsuch(zeil,'source simulation'           ,17,30,lsimsr,*10)
      call lsuch(zeil,'referenzloesung'             ,15,30,lrefsl,*10)
      call lsuch(zeil,'reference solution'          ,18,30,lrefsl,*10)
      call lsuch(zeil,'inverse loesung'             ,15,30,linver,*10)
      call lsuch(zeil,'inverse solution'            ,16,30,linver,*10)
      call isuch(zeil,'art des einflussraums'       ,21,30,invway,*10)
      call isuch(zeil,'type of influence space'     ,23,30,invway,*10)
      call isuch(zeil,'art-der-least-squares'       ,21,30,lsqway,*10)
      call isuch(zeil,'least squares method'        ,20,30,lsqway,*10)
      call lsuch(zeil,'focussierte-loesung'         ,19,30,lfocus,*10)
      call lsuch(zeil,'focused solution'            ,16,30,lfocus,*10)
      call lsuch(zeil,'focussiert-george'           ,17,30,lgeorg,*10)
      call lsuch(zeil,'focused george'              ,14,30,lgeorg,*10)
      call lsuch(zeil,'cut-george-louis'            ,16,30,llucut,*10)
      call lsuch(zeil,'cut george louis'            ,16,30,llucut,*10)
      call lsuch(zeil,'mit-tiefenwichtung'          ,18,30,ldepsc,*10)
      call lsuch(zeil,'use depth weigth'            ,16,30,ldepsc,*10)
      call lsuch(zeil,'mit-covarianzmatrix'         ,19,30,lcovar,*10)
      call lsuch(zeil,'use covariance matrix'       ,21,30,lcovar,*10)
      call lsuch(zeil,'einflussmatrix-check'        ,20,30,linfck,*10)
      call lsuch(zeil,'check influence matrix'      ,22,30,linfck,*10)
      call lsuch(zeil,'konvergenztest'              ,14,30,lkonve,*10)
      call lsuch(zeil,'convergence test'            ,16,30,lkonve,*10)
      call lsuch(zeil,'chi square loesung'          ,18,30,lchisl,*10)
      call lsuch(zeil,'chi square solution'         ,19,30,lchisl,*10)
      call lsuch(zeil,'anisotrope leitfaehigkeit'   ,25,30,laniso,*10)
      call lsuch(zeil,'anisotropic conductivities'  ,26,30,laniso,*10)
      call lsuch(zeil,'mit waermeuebergangszahl'    ,24,30,lalpha,*10)
      call lsuch(zeil,'heat transfer coefficient'   ,25,30,lalpha,*10)
      call lsuch(zeil,'knotenbezogene lasten'       ,21,30,lnodlo,*10)
      call lsuch(zeil,'nodal loads'                 ,11,30,lnodlo,*10)
      call lsuch(zeil,'direkter-inverser-loeser'    ,24,30,ldirsl,*10)
      call lsuch(zeil,'direct inverse solver'       ,21,30,ldirsl,*10)
      call lsuch(zeil,'invers-linear-cg-loeser'     ,23,30,lcglin,*10)
      call lsuch(zeil,'linear cg inverse solver'    ,24,30,lcglin,*10)
      call lsuch(zeil,'entropie-regularisierung'    ,24,30,lentro,*10)
      call lsuch(zeil,'entropy regularization'      ,22,30,lentro,*10)
      call lsuch(zeil,'l1-norm-regularisierung'     ,23,30,ll1nrm,*10)
      call lsuch(zeil,'l1 norm regularization'      ,22,30,ll1nrm,*10)
      call lsuch(zeil,'l2-norm-regularisierung'     ,23,30,ll2nrm,*10)
      call lsuch(zeil,'l2 norm regularization'      ,22,30,ll2nrm,*10)
      call lsuch(zeil,'normalen-constraint'         ,19,30,lnorco,*10)
      call lsuch(zeil,'normal constraint'           ,17,30,lnorco,*10)
      call lsuch(zeil,'transiente-loesung'          ,18,30,ltrans,*10)
      call lsuch(zeil,'transient solution'          ,18,30,ltrans,*10)
      call lsuch(zeil,'deviation scan'              ,14,30,ldevsc,*10)
      call lsuch(zeil,'deviation scan'              ,14,30,ldevsc,*10)
      call lsuch(zeil,'simulated annealing'         ,19,30,lsiman,*10)
      call lsuch(zeil,'simulated annealing'         ,19,30,lsiman,*10)
      call lsuch(zeil,'leitfaehigkeiten geordnet'   ,25,30,lperor,*10)
      call lsuch(zeil,'sorted conductivity'         ,19,30,lperor,*10)
      call lsuch(zeil,'flaechennormalen-check'      ,22,30,lnorch,*10)
      call lsuch(zeil,'check surface normals'       ,21,30,lnorch,*10)
      call lsuch(zeil,'magnetische-loesung'         ,19,30,logmeg,*10)
      call lsuch(zeil,'meg solution'                ,12,30,logmeg,*10)
      call lsuch(zeil,'messdaten-verrauschen'       ,21,30,lognoi,*10)
      call lsuch(zeil,'add noise to measurements'   ,25,30,lognoi,*10)
      call lsuch(zeil,'lambda-iteration'            ,16,30,llamit,*10)
      call lsuch(zeil,'lambda iteration'            ,16,30,llamit,*10)
!
!
!---C.Wolters, 7.07.2000
!---Generation of an anisotropic permeability file for a 
!---sphere simulation
      call lsuch(zeil,'aniso simulation'            ,16,30,lwrani,*10)
!
!%%% auskommentiert am 24.3.97 Robert Pohlmeier, wegen unvollstaendiger Implementierung
!      call lsuch(zeil,'neumann-korrektur'           ,17,30,logneu,*10)
!      call lsuch(zeil,'neumann correction'          ,18,30,logneu,*10)
!      call lsuch(zeil,'fundamentalloesung'          ,18,30,lrango,*10)
!      call lsuch(zeil,'fundamental solution'        ,20,30,lrango,*10)
!      call lsuch(zeil,'volumenint-surface'          ,18,30,logsur,*10)
!      call lsuch(zeil,'volumenint surface'          ,18,30,logsur,*10)
!
      goto 10
 20   call qfclos(iunit,0) 
!
!
      if ((invway.eq.1).and.lnorco) then
         lnorco=.false.
         write(*,*)' WARNING: normal constraint is not possible ',      &
     &             ' with a volume as influence space!'
      end if
!
      ndmflx=1
      if (ldipol.and..not.lnorco) ndmflx=ndmgeo
      numper=npogeo
!
!---element related properties
      numper=nelgeo
      if (.not.llanzp) lconsi=.false.
!     
      if (lmncor) then
         if (linver.or.linmat) then
            lmncor=.false.
            write(*,*)                                                  &
     &  ' WARNING: Correction to zero average potential will not be '
            write(*,*)                                                  &
     &  ' performed for inverse solutions or influence matrices!'
         end if
      end if
!
      if (lsqway.eq.0) lcholy=.true.
      if (lsqway.eq.1) lqrfak=.true.
      if (lsqway.eq.2) lqrfax=.true.
      if (lsqway.eq.3) lsngvl=.true.
      if (lsqway.lt.0.or.lsqway.gt.3) then
         write (*,*) 'choose lsqway appropriately'
         stop
      end if
!
      if (invway.ne.1.and.invway.ne.2) then
         write (*,*) 'invway not defined'
         stop
      end if
!
      if (linver) then
         idum=0
         if (ldirsl) idum=idum+1
         if (lcglin) idum=idum+1
         if (lentro) idum=idum+1
         if (ll1nrm) idum=idum+1
         if (ll2nrm) idum=idum+1
         if (ldevsc) idum=idum+1
         if (lsiman) idum=idum+1
         if (idum.ne.1) then
            write(*,*)'exactly one inverse method has to be chosen!'
            stop
         end if
      end if
!
      if (linver.and.lsngvl.and.(.not.lcovar)) then
         write(*,*) ' TSVD requires the use of the covariance matrix'
         stop
      end if
!
      if (linver.and.lgeorg.and.lfocus) then
         write(*,*) ' GEORGE and FOCUSS algorithm are not compatible!'
         stop
      end if
!
      if (linver.and.llamit.and.lkonve) then
         write(*,*) ' convergence test and lambda iteration are not ',  &
     &        'compatible!'
         stop
      end if
!
      if (laniso) then
         naniso=6
         if (ndmgeo.ne.3) then
            write(*,*)' Anisotropic conductivity values require'
            write(*,*)' a threedimensionale structure!'
         end if
      end if
!
      return
  900 format(a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reatol
!
!---reads tolerances
      include 'neurofem.inc'
!
      incsta=0
      incend=0
      iunit=qfreeu()
      call qfopse(iunit,filnam(4),'un','fo',ierr)
      if ( ierr.ne.0 ) stop 'error reatol'
   10 read(iunit,900,end=800) zeil(1:80)
      call qctolo(zeil)
!
      call dsuch(zeil,'toleranz in der loesung'      ,23,30,tolsol,*10)
      call dsuch(zeil,'tolerance forward solution'   ,26,30,tolsol,*10)
      call dsuch(zeil,'toleranz inverse loesung'     ,24,30,tolinv,*10)
      call dsuch(zeil,'tolerance inverse solution'   ,26,30,tolinv,*10)
      call isuch(zeil,'gleichungsloeser'             ,16,30,isolvr,*10)
      call isuch(zeil,'solver method'                ,13,30,isolvr,*10)
      call isuch(zeil,'anzahl der zeitinkremente'    ,25,30,numinc,*10)
      call isuch(zeil,'time increments'              ,15,30,numinc,*10)
      call isuch(zeil,'analytische loesung'          ,19,30,numana,*10)
      call isuch(zeil,'analytical solution'          ,19,30,numana,*10)
      call isuch(zeil,'startzeitpunkt'               ,14,30,incsta,*10)
      call isuch(zeil,'begin at time point'          ,19,30,incsta,*10)
      call isuch(zeil,'endzeitpunkt'                 ,12,30,incend,*10)
      call isuch(zeil,'end at time point'            ,17,30,incend,*10)
      call isuch(zeil,'integrationsgrad'             ,16,30,intgrd,*10)
      call isuch(zeil,'degree of integration'        ,21,30,intgrd,*10)
      call isuch(zeil,'glaettungsgrad'               ,14,30,iglatt,*10)
      call isuch(zeil,'method of smoothness'         ,20,30,iglatt,*10)
      call isuch(zeil,'dipol-glaette'                ,13,30,igladi,*10)
      call isuch(zeil,'dipole smoothness'            ,17,30,igladi,*10)
      call isuch(zeil,'dipol-ordnung'                ,13,30,norder,*10)
      call isuch(zeil,'dipole order'                 ,12,30,norder,*10)
      call dsuch(zeil,'dipol-skalierung'             ,16,30,scadip,*10)
      call dsuch(zeil,'dipole scale'                 ,12,30,scadip,*10)
      call dsuch(zeil,'dipol-lambda'                 ,12,30,diplam,*10)
      call dsuch(zeil,'dipole lambda'                ,13,30,diplam,*10)
      call dsuch(zeil,'dipol-distanz'                ,13,30,disdip,*10)
      call dsuch(zeil,'dipole distance'              ,15,30,disdip,*10)
      call dsuch(zeil,'lambda-invers'                ,13,30,xlmbda,*10)
      call dsuch(zeil,'lambda inverse'               ,14,30,xlmbda,*10)
      call dsuch(zeil,'lambda-ende'                  ,11,30,xlmbfi,*10)
      call dsuch(zeil,'lambda end'                   ,10,30,xlmbfi,*10)
      call dsuch(zeil,'lambda-inkrement'             ,16,30,xlmbin,*10)
      call dsuch(zeil,'lambda increment'             ,16,30,xlmbin,*10)
      call dsuch(zeil,'chi-square wert'              ,15,30,chival,*10)
      call dsuch(zeil,'desired chi square value'     ,24,30,chival,*10)
      call dsuch(zeil,'epsilon-dirac'                ,13,30,epsilo,*10)
      call dsuch(zeil,'epsilon dirac'                ,13,30,epsilo,*10)
      call isuch(zeil,'anzahl der dipole'            ,17,30,numdip,*10)
      call isuch(zeil,'number of dipoles'            ,17,30,numdip,*10)
      call dsuch(zeil,'temperatur-start'             ,16,30,tmpsta,*10)
      call dsuch(zeil,'start temperature'            ,17,30,tmpsta,*10)
      call dsuch(zeil,'temperatur-faktor'            ,17,30,tmpfak,*10)
      call dsuch(zeil,'decreasing factor of temp'    ,25,30,tmpfak,*10)
      call isuch(zeil,'maximale anzahl anneal'       ,22,30,maxtry,*10)
      call isuch(zeil,'maximum annealing steps'      ,23,30,maxtry,*10)
      call isuch(zeil,'eigenwertanzahl'              ,15,30,nval  ,*10)
      call isuch(zeil,'number of eigenvalues'        ,21,30,nval  ,*10)
      call isuch(zeil,'eigenwertstellen'             ,16,30,nfig  ,*10)
      call isuch(zeil,'digits of eigenvalue'         ,20,30,nfig  ,*10)
      call isuch(zeil,'eigenwert-maxop'              ,15,30,maxop ,*10)
      call isuch(zeil,'eigenvalue maxop'             ,16,30,maxop ,*10)
      call isuch(zeil,'eigenwert-maxj'               ,14,30,maxj  ,*10)
      call isuch(zeil,'eigenvalue maxj'              ,15,30,maxj  ,*10)
      call isuch(zeil,'eigenwert-block'              ,15,30,nblock,*10)
      call isuch(zeil,'eigenvalue block'             ,16,30,nblock,*10)
      call isuch(zeil,'modenanzahl'                  ,11,30,nmodes,*10)
      call isuch(zeil,'number of modes'              ,15,30,nmodes,*10)
      call dsuch(zeil,'dynamic-shift'                ,13,30,dshift,*10)
      call dsuch(zeil,'dynamic shift'                ,13,30,dshift,*10)
      call dsuch(zeil,'relative-accuracy'            ,17,30,relacc,*10)
      call dsuch(zeil,'relative accuracy'            ,17,30,relacc,*10)
      call isuch(zeil,'eigenwert-store'              ,15,30,nstore,*10)
      call isuch(zeil,'eigenvalue store'             ,16,30,nstore,*10)
      call isuch(zeil,'lanzpack-inertia'             ,16,30,inerti,*10)
      call isuch(zeil,'lanzpack inertia'             ,16,30,inerti,*10)
      call isuch(zeil,'lanzpack-steps-per-shift'     ,24,30,nstpsh,*10)
      call isuch(zeil,'lanzpack steps per shift'     ,24,30,nstpsh,*10)
      call isuch(zeil,'svd-workarray'                ,13,30,nwrksv,*10)
      call isuch(zeil,'svd workarray'                ,13,30,nwrksv,*10)
      call dsuch(zeil,'louis-r-parameter'            ,17,30,rpostl,*10)
      call dsuch(zeil,'r parameter'                  ,11,30,rpostl,*10)
      call dsuch(zeil,'dipole-threshold'             ,16,30,thrdip,*10)
      call dsuch(zeil,'dipole threshold'             ,16,30,thrdip,*10)
!
      goto 10
  800 call qfclos(iunit,0)
!
      mgagam=intgrd
!
      if (abs(scadip).lt.wuzeps) then
         scadip=z1
         write(*,910) scadip
      end if
!
      if (incend.gt.numinc+1) then
         incend=numinc+1
         write(*,920) incend
      end if
      numtim=incend-incsta+1
      if (.not.ltrans) numtim=1
!
      if (linver.and.ldevsc.and.numdip.ne.1) then
         write(*,940) 
         numdip=1
      end if
!
      if (rpostl.le.z1) then
         write(*,930) rpostl
      end if
!
      if (lkonve.and.((xlmbin.eq.z1).or.(xlmbin.le.z0))) then
         write (*,970)
         lkonve=.false.
      end if
!
      if (tmpfak.lt.tmpfkm) then
         write(*,950) tmpfkm
         tmpfak=tmpfkm
      else if(tmpfak.ge.z1) then
         write(*,960)
         stop
      end if
!
      if (nmodes.gt.abs(nval)) then
         write(*,900) 'nmodes > abs(nval)'
         stop
      end if
!
      if ((isolvr.ge.5).or.(isolvr.le.0)) then
         isolvr = 2
         write(*,980) isolvr
      end if

      return
  900 format (a)
  910 format (1x,' scadip changed to ',g12.5,' (!!!)')
  920 format (' incend changed to ',i6)
  930 format (' Warning: rpostl should be larger than 1 ',f12.5)
  940 format (' Warning: numdip changed to 1 for deviation scan')
  950 format (' Warning: decreasing factor is set to: ',f12.5)
  960 format (' DECREASING FAKTOR MUST BE LESS THAN ONE')
  970 format (' Warning: lambda increment must be positive ',           &
     &        'and not equal 1, no convergence test is performed')
  980 format (' isolvr changed to ',i6)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaflx(nodevl,nodinf,flxmat,vecref)
      include 'neurofem.inc'
      dimension nodevl(nevlpo),nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),vecref(npogeo)
!
      call dset(npogeo,z0,vecref,1)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam(13),'un','un',ierr)
      if (ierr.gt.0) then
         write (*,900) 'error filnam(13) in reaflx'
         stop
      end if
!
      read(iunit) npogte
      if (npogte.lt.0) then
         npogte=-npogte
         read(iunit) numete
         read(iunit) nevlte
         if (npogte.eq.npogeo) then
            write(*,926) 'OK'
         else 
            write(*,926) 'MISMATCH'
            write(*,905) 'structure: ', npogeo,                         &
     &                   '   matrix: ', npogte
            ierr=1
         end if
         if (numete.eq.numeeg) then
            write(*,927) 'OK'
         else 
            write(*,927) 'MISMATCH'
            write(*,905) 'structure: ', numeeg,                         &
     &                   '   matrix: ', numete
            ierr=1
         end if
      else
         nevlte=npogte
      end if
      read(iunit) ninfte
      read(iunit) ndmfte
      read(iunit) chksum
      ierr=0
      if (nevlte.eq.nevlpo) then
         write(*,920) 'OK'
      else 
         write(*,920) 'MISMATCH'
         write(*,905) 'structure: ', nevlpo,                            &
     &                '   matrix: ', nevlte
         ierr=1
      end if
      if (ninfte.eq.ninfpo) then
         write(*,910) 'OK'
      else
         write(*,910) 'MISMATCH'
         write(*,905) 'structure: ', ninfpo,                            &
     &                '   matrix: ', ninfte
         ierr=2
      end if
      if (ndmfte.eq.ndmflx) then
         write(*,925) 'OK'
      else
         write(*,925) 'MISMATCH'
         write(*,905) 'structure: ', ndmflx,                            &
     &                '   matrix: ', ndmfte
         ierr=3
      end if
      if (ierr.eq.0) then
         write(*,930) chksum
      else
         write(*,900) 'Influence-matrix and given structure mismatch'
         stop
      end if
!
      read(iunit) (nodinf(i),i=1,ninfpo)
      read(iunit) (nodevl(i),i=1,numeeg)
      read(iunit) (vecref(i),i=1,npogeo)
      read(iunit) (((flxmat(i,j,k),i=1,nevlpo),j=1,ninfpo),k=1,ndmflx)
      call qfclos(iunit,0)
!
      numtot=ninfpo*nevlpo*ndmflx
      tessum=dsum(numtot,flxmat,1)/dble(numtot)
      write(*,940) chksum-tessum
!
      return
  900 format(a)
  905 format(27x,a,i8,a,i8)
  910 format('Number of Influence  nodes ',a)
  920 format('Number of Evaluation nodes ',a)
  925 format('Dimension of the matrix    ',a)
  926 format('Number of Geometry nodes   ',a)
  927 format('Number of EEG nodes        ',a)
  930 format('Checksum of influence-matrix is: ',e23.16)
  940 format('Comparison of checksums: ',e23.16)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine messin(nodinf,invbnd,bndinv)
      include 'neurofem.inc'
      dimension nodinf(ninfpo),invbnd(npoinv)
      dimension bndinv(npoinv)
!
      call iset(npoinv, 1,invbnd,1)
      call dset(npoinv,z0,bndinv,1)
!
      do 20 i=1,ninfpo
            invbnd(nodinf(i))=0
 20   continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine mesmod(nodevl,dkondg,potmes,vecref,covinv,dprodg,      &
     &           potmeg,covmeg)
      include 'neurofem.inc'
      dimension nodevl(nevlpo)
      dimension dkondg(npogeo),potmes(nevlpo),vecref(npogeo),           &
     &          covinv(nevlpo),dprodg(npogeo),potmeg(nummeg),           &
     &          covmeg(nummeg)
!
!      open (53,file='pot.gnu')
!      write(53,930)
      iseed=-2
      ran=z0
      signoi=z0
      potave=z0
      valnoi=z0
      if (numeeg.gt.0) then
         do 10 i=1,numeeg
            knoglo=nodevl(i)
            potmes(i)=(dkondg(knoglo)-vecref(knoglo))
            dummy=dprodg(knoglo)*dprodg(knoglo)
            if (dummy.gt.z0) then
               covinv(i)=z1/dummy
            else
               write(*,*) 'if using the covariance matrix, the ',       &
     &         'noise of each measurement channel must not be zero!'
               stop 'check covariance matrix'
            end if
            if (lognoi) then
               ran=sign(abs(dprodg(knoglo)),z2*dble(ran2(iseed))-z1)
               potmes(i)=potmes(i)+ran
            end if
!     write(53,920) i,nodevl(i),potmes(i)-ran,potmes(i)
            signoi=signoi+sqrt(potmes(i)*potmes(i)*covinv(i))
            potave=potave+potmes(i)*potmes(i)
            valnoi=valnoi+dummy
   10    continue
         raueeg=potave/valnoi
      else
         raueeg=z0
      end if
      valnoi=z0
      avemeg=z0
      if (logmeg) then
         do 20 j=1,nummeg
            jpos=numeeg+j
            potmes(jpos)=potmeg(j)
            dummy=covmeg(j)*covmeg(j)
            if (dummy.gt.z0) then
               covinv(jpos)=z1/dummy
            else
               write(*,*) 'if using the covariance matrix, the ',       &
     &         'noise of each measurement channel must not be zero!'
               stop 'check covariance matrix'
            end if
            if (lognoi) then
               ran=sign(abs(covmeg(j)),z2*dble(ran2(iseed))-z1)
               potmes(jpos)=potmes(jpos)+ran
            end if
!            write(53,920) jpos,jpos,potmes(jpos)-ran,potmes(jpos)
            signoi=signoi+sqrt(potmeg(j)*potmeg(j)*covinv(jpos))
            avemeg=avemeg+potmes(jpos)*potmes(jpos)
            valnoi=valnoi+dummy
   20    continue
         raumeg=avemeg/valnoi
      else
         raumeg=z0
      end if
      signo2=(dble(numeeg)*raueeg+dble(nummeg)*raumeg)/dble(nevlpo)
      signo2=sqrt(signo2)
      signoi=signoi/dble(nevlpo)
!
      potsqr=potave
      potave=sqrt(potave)
      if (lcovar) write(*,910) signoi
      if (lcovar) write(*,900) signo2
!
!      close (53)
      return
  900 format(' Theoretical     Signal-to-noise ratio is: ',f19.5)
  910 format(' Single averaged Signal-to-noise ratio is: ',f19.5)
  920 format(1x,i6,1x,i6,1x,e12.5,1x,e12.5)
  930 format('# i node pot-ori pot-perturbed')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wrimat(ievlor,inflor,flxmat,vecref, nb_electrodes,     &
     &                  nb_nodes, nb_sources, nb_source_dir)
      include 'neurofem.inc'
      dimension ievlor(nb_electrodes),inflor(nb_sources)
      dimension flxmat(nb_electrodes,nb_sources,nb_source_dir)
	  dimension vecref(nb_nodes)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam(13),'un','un',ierr)
      if (ierr.gt.0) then
         write (*,900) 'error filnam(13) in wrimat'
         stop
      end if
!
!--determine checksum
      numtot=nb_sources*nb_electrodes*nb_source_dir
      chksum=dsum(numtot,flxmat,1)/dble(numtot)
      write(*,930) chksum
!
      write(iunit) -nb_nodes
      write(iunit) nb_electrodes
      write(iunit) nb_electrodes
      write(iunit) nb_sources
      write(iunit) nb_source_dir
      write(iunit) chksum
!
      write(iunit) (inflor(i),i=1,nb_sources)
      write(iunit) (ievlor(i),i=1,numeeg)
      write(iunit) (vecref(i),i=1,nb_nodes)
      write(iunit) (((flxmat(i,j,k),i=1,nb_electrodes),                 &
     &                         j=1,nb_sources),k=1,nb_source_dir)
      call qfclos(iunit,0)
!
      return
 900  format(a)
 930  format('Checksum of influence-matrix is: ',e23.16)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine calmat(knogeo,itygeo,necgeo,ipogeo,indgeo,             &
     &                  sysmat,volmat,xyzgeo,pergeo,                    &
     &                  knoflf,ityflf,necflf,ipoflf,indflf,surmat,      &
     &                  sysrng,knoalf,ityalf,necalf,                    &
     &                  alfmat,alfval,tmpamb,vcbamb,xyzflf)
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo),    &
     &                 ipogeo(npogeo),indgeo(lengeo),ipoflf(npoflf),    &
     &                 indflf(lenflf),                                  &
     &          knoflf(mxkn2d,nelflf),ityflf(nelflf),necflf(nelflf),    &
     &          knoalf(mxkn2d,numalf),ityalf(numalf),necalf(numalf)
      dimension sysmat(lengeo),volmat(lengeo),xyzgeo(npogeo,ndmgeo),    &
     &          pergeo(numper,naniso),surmat(lenflf),alfmat(lengeo),    &
     &          alfval(numalf),tmpamb(numalf),vcbamb(npogeo),           &
     &          xyzflf(npoflf,ndmgeo),sysrng(lengeo)
      dimension syslok(mxkn3d,mxkn3d),vollok(mxkn3d,mxkn3d),            &
     &          surlok(mxkn2d,mxkn2d),alflok(mxkn2d,mxkn2d),            &
     &          vcblok(mxkn2d),rnglok(mxkn3d,mxkn3d)
!
!---zero system matrices
      call dset(lengeo,z0,sysmat,1)
      call dset(lengeo,z0,volmat,1)
      call dset(mxkn3d*mxkn3d,z0,syslok,1)
      call dset(mxkn3d*mxkn3d,z0,vollok,1)
      call dset(mxkn2d*mxkn2d,z0,alflok,1)
      call dset(mxkn2d,z0,vcblok,1)
!
!---compile volume matrices
      cpuvol=qscput(9,0,ierr)
      call percen(0,nelgeo)
      do 10 iel=1,nelgeo
         call volele(knogeo,itygeo,necgeo,xyzgeo,pergeo,                &
     &        syslok,vollok,rnglok,numper,npogeo,   iel )
         do 20 i=1,necgeo(iel)
            iglo=knogeo(i,iel)
            do 30 j=1,necgeo(iel)
               jglo=knogeo(j,iel)
               if (jglo.le.iglo) then
                  if (iglo.gt.1) then
                     ia=ipogeo(iglo-1)+1
                     in=ipogeo(iglo)-ipogeo(iglo-1)
                     iposi=ia+ipoj(jglo,indgeo(ia),in)
                  else
                     iposi=1
                  end if
                  sysmat(iposi) = sysmat(iposi) + syslok(i,j)
                  volmat(iposi) = volmat(iposi) + vollok(i,j)
                  if (lrango) then
                     sysrng(iposi) = sysrng(iposi) + rnglok(i,j)
                  end if
               end if
   30       continue
   20    continue
         call percen(iel,nelgeo)
   10 continue
      write(*,900) qscput(9,1,ierr)
      dum=qscput(9,0,ierr)
!
      if (invway.eq.1) goto 800
!     
!---  compile surface matrix in global numbering
      call dset(lenflf,z0,surmat,1)
      call dset(mxkn2d*mxkn2d,z0,surlok,1)
      call percen(0,nelflf)
      do 40 iel=1,nelflf
         call surele(knoflf,ityflf,necflf,xyzflf,surlok, iel )
         do 50 i=1,necflf(iel)
            iglo=knoflf(i,iel)
            do 60 j=1,necflf(iel)
               jglo=knoflf(j,iel)
               if (jglo.le.iglo) then
                  if (iglo.gt.1) then
                     ia=ipoflf(iglo-1)+1
                     in=ipoflf(iglo)-ipoflf(iglo-1)
                     iposi=ia+ipoj(jglo,indflf(ia),in)
                  else
                     iposi=1
                  end if
                  surmat(iposi) = surmat(iposi) + surlok(i,j)
               end if
   60       continue
   50    continue
         call percen(iel,nelflf)
   40 continue
      write(*,910) qscput(9,1,ierr)
!
  800 if (.not.lalpha) goto 1000
!
!---compile alpha matrix in global numbering
      call dset(lengeo,z0,alfmat,1)
      call dset(npogeo,z0,vcbamb,1)
      call percen(0,numalf)
      do 70 iel=1,numalf
         call alfele(knoalf,ityalf,necalf,xyzgeo,alfval,                &
     &        tmpamb,alflok,vcblok, iel )
         do 80 i=1,necalf(iel)
            iglo=knoalf(i,iel)
            vcbamb(iglo) = vcbamb(iglo) + vcblok(i)
            do 90 j=1,necalf(iel)
               jglo=knoalf(j,iel)
               if (jglo.le.iglo) then
                  if (iglo.gt.1) then
                     ia=ipogeo(iglo-1)+1
                     in=ipogeo(iglo)-ipogeo(iglo-1)
                     iposi=ia+ipoj(jglo,indgeo(ia),in)
                  else
                     iposi=1
                  end if
                  alfmat(iposi) = alfmat(iposi) + alflok(i,j)
               end if
   90       continue
   80    continue
         call percen(iel,numalf)
   70 continue
!
 1000 return
  900 format (' CPU--time for  volume integration: ',g12.5)
  910 format (' CPU--time for surface integration: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine calsti(knogeo,itygeo,necgeo,ipogeo,indgeo,             &
     &                  sysmat,volmat,xyzgeo,pergeo,sysrng)
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo),    &
     &                 ipogeo(npogeo),indgeo(lengeo)
      dimension sysmat(lengeo),volmat(lengeo),xyzgeo(npogeo,ndmgeo),    &
     &          pergeo(numper,naniso),sysrng(lengeo)
      dimension syslok(mxkn3d,mxkn3d),vollok(mxkn3d,mxkn3d),            &
     &          rnglok(mxkn3d,mxkn3d)
!
!---zero system matrices
      call dset(lengeo,z0,sysmat,1)
      call dset(lengeo,z0,volmat,1)
      if (lrango) then
         call dset(lengeo,z0,sysrng,1)
      end if
      call dset(mxkn3d*mxkn3d,z0,syslok,1)
      call dset(mxkn3d*mxkn3d,z0,vollok,1)
      call dset(mxkn3d*mxkn3d,z0,rnglok,1)
!
!---compile matrices
!      cpuvol=qscput(9,0,ierr)
      call percen(0,nelgeo)
      do 10 iel=1,nelgeo
         call volele(knogeo,itygeo,necgeo,xyzgeo,pergeo,                &
     &        syslok,vollok,rnglok,numper,npogeo,   iel )
         do 20 i=1,necgeo(iel)
            iglo=knogeo(i,iel)
            do 30 j=1,necgeo(iel)
               jglo=knogeo(j,iel)
               if (jglo.le.iglo) then
                  if (iglo.gt.1) then
                     ia=ipogeo(iglo-1)+1
                     in=ipogeo(iglo)-ipogeo(iglo-1)
                     iposi=ia+ipoj(jglo,indgeo(ia),in)
                  else
                     iposi=1
                  end if
                  sysmat(iposi) = sysmat(iposi) + syslok(i,j)
                  volmat(iposi) = volmat(iposi) + vollok(i,j)
                  if (lrango) then
                     sysrng(iposi) = sysrng(iposi) + rnglok(i,j)
                  end if
               end if
   30       continue
   20    continue
!         call percen(iel,nelgeo)
   10 continue
!      if (lrango) then
!         write(*,900) qscput(9,1,ierr)
!      else
!         write(*,910) qscput(9,1,ierr)
!      end if
!      dum=qscput(9,0,ierr)
!
! 1000 return
!  900 format (' CPU--time for stiffness/volume/rng integration: ',g12.5)
!  910 format (' CPU--time for stiffness/volume integration: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine calvol(knogeo,itygeo,necgeo,ipogeo,indgeo,             &
     &                  volmat,xyzgeo,pergeo,sysrng)
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo),    &
     &                 ipogeo(npogeo),indgeo(lengeo)
      dimension volmat(lengeo),sysrng(lengeo),xyzgeo(npogeo,ndmgeo),    &
     &          pergeo(numper,naniso)
      dimension vollok(mxkn3d,mxkn3d),rnglok(mxkn3d,mxkn3d),            &
     &          syslok(mxkn3d,mxkn3d)
!
!---zero the matrices
      call dset(lengeo,z0,volmat,1)
      if (lrango) then
         call dset(lengeo,z0,sysrng,1)
      end if
      call dset(mxkn3d*mxkn3d,z0,vollok,1)
      call dset(mxkn3d*mxkn3d,z0,syslok,1)
      call dset(mxkn3d*mxkn3d,z0,rnglok,1)
!
!---compile volume matrix
      cpuvol=qscput(9,0,ierr)
      call percen(0,nelgeo)
      do 10 iel=1,nelgeo
         call volele(knogeo,itygeo,necgeo,xyzgeo,pergeo,                &
     &        syslok,vollok,rnglok,numper,npogeo,   iel )
         do 20 i=1,necgeo(iel)
            iglo=knogeo(i,iel)
            do 30 j=1,necgeo(iel)
               jglo=knogeo(j,iel)
               if (jglo.le.iglo) then
                  if (iglo.gt.1) then
                     ia=ipogeo(iglo-1)+1
                     in=ipogeo(iglo)-ipogeo(iglo-1)
                     iposi=ia+ipoj(jglo,indgeo(ia),in)
                  else
                     iposi=1
                  end if
                  volmat(iposi) = volmat(iposi) + vollok(i,j)
                  if (lrango) then
                     sysrng(iposi) = sysrng(iposi) + rnglok(i,j)
                  end if
               end if
   30       continue
   20    continue
         call percen(iel,nelgeo)
   10 continue
      if (lrango) then
         write(*,900) qscput(9,1,ierr)
      else
         write(*,910) qscput(9,1,ierr)
      end if
      dum=qscput(9,0,ierr)
!
 1000 return
  900 format (' CPU--time for volume/rng integration: ',g12.5)
  910 format (' CPU--time for volume integration: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine calsur(knoflf,ityflf,necflf,ipoflf,indflf,surmat)
      include 'neurofem.inc'
      dimension ipoflf(npoflf),indflf(lenflf),                          &
     &          knoflf(mxkn2d,nelflf),ityflf(nelflf),necflf(nelflf)
      dimension surmat(lenflf), xyzflf(npoflf,ndmgeo)
      dimension surlok(mxkn2d,mxkn2d)

!
!---  compile surface matrix in global numbering
      call dset(lenflf,z0,surmat,1)
      call dset(mxkn2d*mxkn2d,z0,surlok,1)
      call percen(0,nelflf)
      do 40 iel=1,nelflf
         call surele(knoflf,ityflf,necflf,xyzflf,surlok, iel )
         do 50 i=1,necflf(iel)
            iglo=knoflf(i,iel)
            do 60 j=1,necflf(iel)
               jglo=knoflf(j,iel)
               if (jglo.le.iglo) then
                  if (iglo.gt.1) then
                     ia=ipoflf(iglo-1)+1
                     in=ipoflf(iglo)-ipoflf(iglo-1)
                     iposi=ia+ipoj(jglo,indflf(ia),in)
                  else
                     iposi=1
                  end if
                  surmat(iposi) = surmat(iposi) + surlok(i,j)
               end if
   60       continue
   50    continue
         call percen(iel,nelflf)
   40 continue
      write(*,910) qscput(9,1,ierr)
!
 1000 return
  910 format (' CPU--time for surface integration: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine volele(knogeo,itygeo,necgeo,xyzgeo,pergeo,             &
     &                  syslok,vollok,rnglok,inumel,inumpo,             &
     &                  iel)
!
!---calculation of element matrices
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,*),itygeo(*),necgeo(*),                   &
     &          loknod(32)
      dimension pergeo(inumel,*),xyzgeo(inumpo,*)
      dimension syslok(mxkn3d,mxkn3d),vollok(mxkn3d,mxkn3d),            &
     &          perlok(mxkn3d,     6),rnglok(mxkn3d,mxkn3d)

!$$$30.06.00aa bug reported :
!the field xlk is allocated with size 3*20 (3*mxmx3d)
!and used in the funktion fofagn with the size 3*32 *** error !
!allocation with (3*mxmxkn) resolves the problem
!      dimension xlk(3,mxmx3d),dinpkt(4),pergau(6)
      dimension xlk(3,mxmxkn),dinpkt(4),pergau(6)
      dimension fag(3,mxmxkn),fwe(mxmxkn),fagcpy(3),fagrng(3)
!
!---turn to local node numbers
      call iset(32,0,loknod(1),1)
      do 10 i=1,necgeo(iel)
         kno=knogeo(i,iel)
         loknod(i)=kno
         do 12 idim=1,ndmgeo
            xlk(idim,i) = xyzgeo(kno,idim)
   12    continue
!
!
!---element related properties
         idum=iel
!
         do 15 j=1,naniso
            perlok(i,j) = pergeo(idum,j)
   15    continue
!
!---zero element stiffness matrix and right hand side
         fag(1,i) = z0
         fag(2,i) = z0
         fag(3,i) = z0
         do 10 j=1,necgeo(iel)
            syslok(j,i) = z0
            vollok(j,i) = z0
            rnglok(j,i) = z0
   10 continue
!
!---determine number of gaussian points
      call itpanz(itygeo(iel),necgeo(iel),intgrd,ninpkt,ierr)
      if (ierr.ne.0) then 
         write(*,*)'STOP volele 1; problem with element ', iel
         stop
      end if
!
!---loop over all gaussian points
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
         call itpdat(itygeo(iel),intgrd,iint,dinpkt,ierr)
         if (ierr.ne.0) then
            write(*,*)'STOP volele 2; problem with element ', iel
            stop
         end if
         call fofagn(itygeo(iel),dinpkt(1),dinpkt(2),dinpkt(3),         &
     &               loknod,xlk,fwe,fag,detj,ierr)
         if (ierr.ne.0) then
            write(*,*)'STOP volele 3; problem with element/detj ',      &
     &           iel, detj
            STOP
         end if
!
!---determine values at gaussian points
         call dset(naniso,z0,pergau,1)
         do 40 i = 1, necgeo(iel)
            do 45 j=1,naniso
               pergau(j) = pergau(j) + perlok(i,j) * fwe(i)
   45       continue
   40    continue
         xfak   =  detj   * dinpkt(4)
!     perfak =  pergau * xfak
!     perfa2 = (pergau-permeb) * xfak
!
!---loop over all nodal points and summation of contributions
      do 30 i=1,necgeo(iel)
!
!---multiply with electric conductivity
         if (laniso) then
            fagcpy(1) = pergau(1)*fag(1,i)                              &
     &                + pergau(4)*fag(2,i)+pergau(6)*fag(3,i)
            fagcpy(2) = pergau(2)*fag(2,i)                              &
     &                + pergau(5)*fag(3,i)+pergau(4)*fag(1,i) 
            fagcpy(3) = pergau(3)*fag(3,i)                              &
     &                + pergau(6)*fag(1,i)+pergau(5)*fag(2,i) 
            fagrng(1) = (pergau(1)-permeb)*fag(1,i)                     &
     &                + pergau(4)*fag(2,i)+pergau(6)*fag(3,i)
            fagrng(2) = (pergau(2)-permeb)*fag(2,i)                     &
     &                + pergau(5)*fag(3,i)+pergau(4)*fag(1,i) 
            fagrng(3) = (pergau(3)-permeb)*fag(3,i)                     &
     &                + pergau(6)*fag(1,i)+pergau(5)*fag(2,i) 
         else
            fagcpy(1) = pergau(1)*fag(1,i)
            fagcpy(2) = pergau(1)*fag(2,i)
            fagcpy(3) = pergau(1)*fag(3,i)
            fagrng(1) = (pergau(1)-permeb)*fag(1,i)
            fagrng(2) = (pergau(1)-permeb)*fag(2,i)
            fagrng(3) = (pergau(1)-permeb)*fag(3,i)
         end if
!$$$30/06/00aa multible return not possible in f90
!         do 30 j=1,necgeo(iel)
         do 35 j=1,necgeo(iel)
            sysdum=z0
            rngdum=z0
            do 50 k=1,ndmgeo
               sysdum = sysdum + fagcpy(k)*fag(k,j)
               rngdum = rngdum + fagrng(k)*fag(k,j)
   50       continue
!            syslok(j,i) = syslok(j,i) + sysdum * perfak
!            rnglok(j,i) = rnglok(j,i) + sysdum * perfa2
            syslok(j,i) = syslok(j,i) + sysdum * xfak
            rnglok(j,i) = rnglok(j,i) + rngdum * xfak
            vollok(j,i) = vollok(j,i) + fwe(j)*fwe(i)*xfak
   35    continue
   30    continue
   20 continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine surele(knoflf,ityflf,necflf,xyzflf,surlok,   iel)
!
!---calculation of surface element matrix
      include 'neurofem.inc'
      dimension knoflf(mxkn2d,nelflf),ityflf(nelflf),necflf(nelflf),    &
     &          lokflf(12)
      dimension xyzflf(npoflf,ndmflf)
      dimension surlok(mxkn2d,mxkn2d)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn)
!
      call iset(12,0,lokflf(1),1)
!
!---turn to local node numbers
      do 10 i=1,necflf(iel)
        kno=knoflf(i,iel)
        lokflf(i)=kno
        do 15 idim=1,ndmflf
           xlk(idim,i) = xyzflf(kno,idim)
 15     continue
!
!---zero element stiffness matrix and right hand side
         fag(2,i) = z0
         fag(3,i) = z0
         do 10 j=1,necflf(iel)
            surlok(j,i) = z0
   10 continue
!
!---determine number of gaussian points
      call itpanz(ityflf(iel),necflf(iel),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'surele 1'
!
!---loop over all gaussian points
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
      call itpdat(ityflf(iel),intgrd,iint,dinpkt,ierr)
      if (ierr.ne.0) stop 'surele 2'
      call fofagn(ityflf(iel),dinpkt(1),dinpkt(2),dinpkt(3),            &
     &            lokflf(1),xlk,fwe,fag,detj,ierr)
      if (ierr.ne.0) stop 'surele 3'
!
       xfak =   detj * dinpkt(4)
!
!---loop over all nodal points and summation of contributions
      do 30 i=1,necflf(iel)
!$$$30/06/00aa multible return not possible in f90
!         do 30 j=1,necflf(iel)
         do 35 j=1,necflf(iel)
            surlok(j,i) = surlok(j,i) + fwe(j)*fwe(i)*xfak

 35      continue
 30      continue
 20   continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine alfele(knoalf,ityalf,necalf,xyzgeo,alfval,             &
     &     tmpamb,alflok,vcblok,   iel)
!
!---calculation of alpha-surface element matrix
      include 'neurofem.inc'
      dimension knoalf(mxkn2d,numalf),ityalf(numalf),necalf(numalf),    &
     &          lokflf(12)
      dimension xyzgeo(npogeo,ndmgeo),alfval(numalf),tmpamb(numalf)
      dimension alflok(mxkn2d,mxkn2d),vcblok(mxkn2d)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn)
!
      call iset(12,0,lokflf(1),1)
!
!---turn to local node numbers
      do 10 i=1,necalf(iel)
        kno=knoalf(i,iel)
        lokflf(i)=kno
        do 15 idim=1,ndmgeo
           xlk(idim,i) = xyzgeo(kno,idim)
 15     continue
!
!---zero element stiffness matrix and right hand side
         fag(2,i) = z0
         fag(3,i) = z0
         vcblok(i)= z0
         do 10 j=1,necalf(iel)
            alflok(j,i) = z0
   10 continue
!
!---determine number of gaussian points
      call itpanz(ityalf(iel),necalf(iel),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'alfele 1'
!
!---loop over all gaussian points
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
      call itpdat(ityalf(iel),intgrd,iint,dinpkt,ierr)
      if (ierr.eq.1) stop 'surele 2'
      call fofagn(ityalf(iel),dinpkt(1),dinpkt(2),dinpkt(3),            &
     &            lokflf(1),xlk,fwe,fag,detj,ierr)
      if (ierr.ne.0) stop 'alfele 3'
!
       xfak =   detj * dinpkt(4) * alfval(iel)
       tmpf =   tmpamb(iel)
!
!---loop over all nodal point and summation of contributions
      do 30 i=1,necalf(iel)
            vcblok(i) = vcblok(i) + fwe(i)*xfak*tmpf
!$$$30/06/00aa multible return not possible in f90
!         do 30 j=1,necalf(iel)
         do 35 j=1,necalf(iel)
            alflok(j,i) = alflok(j,i) + fwe(j)*fwe(i)*xfak
 35      continue
 30      continue
 20   continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine extgeo(xyzgeo)
      include 'neurofem.inc'
      dimension xyzgeo(npogeo,3)
!
      dismax=z0
      radius=z0
      do 100 i=1,3
         xyzmax(i)=xyzgeo(idmax(npogeo,xyzgeo(1,i),1),i)
         xyzmin(i)=xyzgeo(idmin(npogeo,xyzgeo(1,i),1),i)
         xyzdst(i)=xyzmax(i)-xyzmin(i)
         radius=max(radius,xyzdst(i)/z2)
         dismax=dismax+xyzdst(i)*xyzdst(i)
 100  continue
      dismax=sqrt(dismax)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine simsrc(mrksim,mrkdip,inflno,                           &
     &                  xyzflf,dkondg,dipsim,dipole,                    &
     &                  vcnflf)
!
!---source simulation
      include 'neurofem.inc'
      dimension mrksim(npoflf,ndmflx),inflno(npoflf),                   &
     &          mrkdip(npoflf,ndmflf)
      dimension xyzflf(npoflf,ndmflf),dkondg(npoflf),                   &
     &          dipsim(npoflf,ndmflx),dipole(npoflf,ndmflf),            &
     &          vcnflf(npoflf,ndmflf)
      dimension diskno(maxdim)
!
      ndmdip=1
      if (ldipol) ndmdip=ndmflf
      call iset(npoflf*ndmdip, 0,mrkdip,1)
      call dset(npoflf*ndmdip,z0,dipole,1)
      call dset(       maxdim,z0,diskno,1)
!
!---how many master sources are available?
      numfhr=isum(npoflf*ndmflx,mrksim,1)
      numwrk=numfhr
      write(*,'(a)') ' Source simulation '
!
!---loop over all nodes for finding the master sources
      do 10 i=1,npoflf
         ifrcou=0
         do 300 k=1,ndmflx
            if (mrksim(i,k).gt.0) ifrcou=ifrcou+1
  300    continue
         numwrk=numwrk-ifrcou
!
!---if no master source, go to next point
         if (ifrcou.eq.0) goto 10
!
         do 20 j=1,npoflf
!
!---if no source allowed here, goto next point
            if (inflno(j).lt.1) goto 20
!
            dif=z0
            do 30 k=1,ndmflf
               diskno(k)=abs(xyzflf(j,k)-xyzflf(i,k))/dismax
               dis=diskno(k)
               dif=dif+dis*dis
   30       continue
!
!---normal operation
            if (lnorco.and.ldipol) then
               streng=dipsim(i,ndmflx)
               do 410 k=1,ndmflf
                  mrkdip(j,k)=1
                  dira=exp(-dif/(z2*epsilo*epsilo)/sqrt(z2*pi))
                  dipole(j,k)=dipole(j,k)+dira*streng*vcnflf(j,k)
  410          continue
            else
               do 40 k=1,ndmflx
                  mrkdip(j,k)=1
                  dira=exp(-dif/(z2*epsilo*epsilo)/sqrt(z2*pi))
                  dipole(j,k)=dipole(j,k)+dira*dipsim(i,k)
   40          continue
            end if
   20    continue
         write(*,900) numwrk+1
         if (numwrk.eq.0) goto 50
   10 continue
!
!---sources are distributed
   50 write(*,910) isum(npoflf*ndmflf,mrkdip,1)
!
!---determine sensible tolerances
      dipsup=z0
      do 222 k=1,ndmdip
         dipmax(k)=abs(dipole(idamax(npoflf,dipole(1,k),1),k))
         dipsup=max(dipmax(k),dipsup)
  222 continue
      toldip=max(thrdip*dipsup,wuzeps)
!
      do 60 i=1,npoflf
         if (lnorco) then
            amount=z0
            do 600 k=1,ndmdip
               amount=amount+dipole(i,k)*dipole(i,k)
  600       continue
            amount=sqrt(amount)
            if (amount.lt.toldip) then
               do 610 k=1,ndmflf
                  mrkdip(i,k)=0
                  dipole(i,k)=z0
  610          continue
            end if
         else
            do 70 k=1,ndmflx
               if (mrkdip(i,k).ne.0) then
                  if (abs(dipole(i,k)).lt.toldip) then
                     mrkdip(i,k)=0
                     dipole(i,k)=z0
                  end if
               end if
   70       continue
         end if
   60 continue
!
      write(*,920) isum(npoflf*ndmdip,mrkdip,1)
!
!      open (56,file='focus.ref')
!      do 800 i=1,npoflf
!         dum=z0
!         do 810 k=1,ndmflf
!            dum=dum+dipole(i,k)*dipole(i,k)
!  810    continue
!         dum=sqrt(dum)
!         write(56,'(3(1x,e12.5))') -z1,dble(i),dum
!  800 continue
!      close (56)
!
      return
  900 format ('DoF number ',i6,' completed in source simulation')
  910 format ('Number of sources after source simulation: ',i6)
  920 format ('Number of sources above         threshold: ',i6)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine mesmrk(nodevl,potmes,vecref,xyzgeo)
      include 'neurofem.inc'
      dimension nodevl(nevlpo)
      dimension potmes(nevlpo),vecref(npogeo),xyzgeo(npogeo,ndmgeo)
!
      iunit=qfreeu()
      call qfopse(iunit,'mes-inv.gnu','un','fo',ierr)
      if (ierr.gt.0) stop 'mes-inv.gnu'
!
      do 10 i=1,numeeg
         iglo  = nodevl(i)
         value = potmes(i)+vecref(i)
         write(iunit,900) (xyzgeo(iglo,j),j=1,ndmgeo),value
 10   continue
!
      call qfclos(iunit,0)
!
      return
 900  format(g12.5,1x,g12.5,1x,g12.5,1x,g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine srcint(ipogeo,indgeo,iposrc,                           &
     &     srcmat,dkondg,source,dprodg)
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),iposrc(npogeo)
      dimension srcmat(lengeo),dkondg(npogeo),source(npogeo),           &
     &     dprodg(npogeo)
!
      call dset(npogeo,z0,dkondg,1)
!
      icount=0
      do 10 i=1,npogeo
         if (iposrc(i).gt.0) then
            icount=icount+1
            dkondg(i)=source(i)
            call matpro(srcmat,dkondg,dprodg,                           &
     &           ipogeo,indgeo,lengeo,npogeo)
            write(*,900) icount,i,dsum(npogeo,dprodg,1)
            dkondg(i)=z0
         end if
 10   continue
!
      call matpro(srcmat,source,dprodg,                                 &
     &     ipogeo,indgeo,lengeo,npogeo)
      write(*,910) dsum(npogeo,dprodg,1)
!
      return
  900 format('Source # ',i6,' node ',i6,' integral strength ',g12.5,    &
     &       ' [A]')
  910 format('Integral of sources (zero?) ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine namgen (filnam,number)
      character*80 filnam
      integer qclen
!
      numend=qclen(filnam)
      filnam(numend-3:numend)='-000'
      if (number.lt.10) then
         write(filnam(numend:numend),'(i1)') number
      else if (number.lt.100) then
         write(filnam(numend-1:numend),'(i2)') number
      else if (number.lt.1000) then
         write(filnam(numend-2:numend),'(i3)') number
      else 
         stop 'file-name generator overflow'
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$aa 30/06/00 dwork and ndrest not used, removed from simtim
!      subroutine simtim(mrksim,mrkdip,knoord,                           &
!     &                  dipsim,dipole,vcnflf,timinc,dwork,ndrest)
      subroutine simtim(mrksim,mrkdip,knoord,                           &
     &                  dipsim,dipole,vcnflf,timinc)
!
!---source simulation for transient problems
      include 'neurofem.inc'
      dimension mrksim(npogeo,ndmflx),mrkdip(npogeo,ndmgeo),            &
     &          knoord(npogeo)
      dimension dipsim(npogeo,ndmflx),dipole(npogeo,ndmgeo),            &
     &          vcnflf(npogeo,ndmgeo),timinc(numinc)
      if (.not.ldipol) stop 'simtim can only be done for dipoles'
!
!---how many master nodes are available?
      master=0
      do 10 i=1,npogeo
         idum=0
         do 20 k=1,ndmflx
            if (mrksim(i,k).gt.0) idum=idum+1
   20    continue
         if(idum.gt.0) master=master+1
   10 continue
      write(*,900) numinc+1,master
!
!---how many order nodes are available?
      norder=0
      do 100 i=1,npogeo
         if (knoord(i).gt.0) norder=norder+1
  100 continue
      write(*,920) norder
!
!---what time range needs to be considered?
      timtot=z0
      do 30 i=1,numinc
         timtot=timtot+timinc(i)
   30 continue
!
!---what time lag will be appropriate
      timlag=timtot/dble(norder)
      epstim=epsilo*timtot
!
!---loop over all time steps
      timact=z0
      iorder=0
      do 40 itim=0,numinc
         write(*,910) itim,timact
         call iset(npogeo*ndmgeo, 0,mrkdip,1)
         call dset(npogeo*ndmgeo,z0,dipole,1)
         iorder=mod(iorder+1,norder)
!
!---loop over all nodes for finding the master sources
            do 50 jorder=1,norder
               knoakt=knoord(jorder)
               dist=dble(jorder-iorder)
               diftim=timtot*sqrt(dist*dist)
               dirac=exp(-diftim/(z2*epstim))
               if (dirac.gt.wuzeps) then
                  if (lnorco) then
                     streng=dipsim(knoakt,ndmflx)*dirac
                     do 60 k=1,ndmgeo
                        mrkdip(knoakt,k)=1
                        dipole(knoakt,k)=streng*vcnflf(knoakt,k)
   60                continue
                  else
                     do 70 k=1,ndmflx
                        mrkdip(knoakt,k)=1
                        dipole(knoakt,k)=dipsim(knoakt,k)               &
     &                       * dirac*vcnflf(knoakt,k)
   70                continue
                  end if
               end if
   50       continue
            call namgen(filnam(8)(1:qclen(filnam(8))),itim)
            call wriknw(mrkdip,dipole,npogeo,ndmgeo,                    &
     &           filnam( 8)(1:qclen(filnam( 8))))
            if(itim.ne.numinc) timact=timact+timinc(itim+1)
   40 continue
!
      return
  900 format(1x,' Time steps ',i3,' number of master sources ',i6)
  910 format(1x,' Working on time step ',i3,' time is ',g12.5)
  920 format(1x,' Sequence nodes for transient source simulation',i6)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaord(knoord,nummax,nfound,filnam)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character*80 filnam,zeil
      dimension knoord(nummax),kno(12)
!
      nfound=0
      call iset(nummax,0,knoord,1)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         write(*,'(a)') 'error while opening ',filnam,' in -reaord-'
         stop
      end if
!
 20   read(iunit,'(a)',end=9998) zeil
      call qctolo(zeil)
      ikart=index(zeil,'boi - knotennummernkarte')
      if (ikart.le.0) goto 20
!
!---wenn richtige karte gefunden, karte bis 'eoi -' abarbeiten
 30   read(iunit,'(a)') zeil
      call qctolo(zeil)
      if (index(zeil,'eoi -').gt.0) goto 300
      read (zeil,'(bn,7x,12i6)',end=9998) (kno(i),i=1,12)
      do 40 i=1,12
         node=kno(i)
         if (node.ne.0) then
            nfound=nfound+1
            knoord(nfound)=node
         end if
 40   continue
      goto 30
 300  call qfclos(iunit,0)
!
!
      return
 9998 stop 'fehler reaord, file corrupted'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine curdip(xyzflf,dipole,npoflf,ndmflf,                    &
     &                  filnam,thresh,lsiman)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character*80 filnam
      logical lsiman
      parameter (z0=0.d0,z1=1.d0,z2=2.d0,f2=0.5d0,f4=0.25d0,f5=0.2d0)
      dimension xyzflf(npoflf,ndmflf),dipole(npoflf,ndmflf)
!
      valmax=z0
      do 110 i=1,npoflf
         valakt=z0
         do 120 k=1,ndmflf
            valakt=valakt+dipole(i,k)*dipole(i,k)
  120    continue
         valmax=max(valmax,valakt)
  110 continue
      valmax=sqrt(valmax)
      if( lsiman ) valmax = 0.d0
!
      icount=0
      do 130 i=1,npoflf
         valakt=z0
         do 140 k=1,ndmflf
            valakt=valakt+dipole(i,k)*dipole(i,k)
  140    continue
         valakt=sqrt(valakt)
         if (valakt.gt.thresh*valmax) icount=icount+1
  130 continue
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) stop 'curdip: qfopse'
      write(iunit,900) 1,z1,z2
      write(iunit,910) icount,0,0,0,f5
!
      do 30 i=1,npoflf
         valakt=z0
         do 150 k=1,ndmflf
            valakt=valakt+dipole(i,k)*dipole(i,k)
  150    continue
         valakt=sqrt(valakt)
         if (valakt.gt.thresh*valmax) then
            if (ndmflf.eq.1) then
               write(iunit,912) (xyzflf(i,k),k=1,ndmflf),               &
     &              (dipole(i,k),k=1,ndmflf)
            else if (ndmflf.eq.2) then
               write(iunit,915) (xyzflf(i,k),k=1,ndmflf),               &
     &              (dipole(i,k),k=1,ndmflf)
            else
               write(iunit,920) (xyzflf(i,k),k=1,ndmflf),               &
     &              (dipole(i,k),k=1,ndmflf)
            end if
         end if
   30 continue
      write(iunit,930) f2,f2
      call qfclos(iunit,0)
!
      return
  900 format(i3,1x,f6.1,1x,f6.1)
  910 format(i4,1x,i3,1x,i3,1x,i3,1x,f6.1)
  912 format(1(f6.1,1x),1(e12.5,1x))
  915 format(2(f6.1,1x),2(e12.5,1x))
  920 format(3(f6.1,1x),3(e12.5,1x))
  930 format(2(f6.1,1x))
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaper(filakt,ifrgeo,pergeo,nknowe)
      include 'neurofem.inc'
      dimension ifrgeo(numper)
      dimension pergeo(numper,naniso)
!
      iunit=qfreeu()
      call qfopse(iunit,filakt,'un','fo',ierr)
      call iset(numper,0,ifrgeo,1)
      idum=0
      nknowe=0
   10 read(iunit,900,end=100) zeil
      if (index(zeil,'#').eq.1) goto 10
      if (zeil.eq.' ') goto 10
      idum=idum+1
      if (laniso) then
         read(zeil,910)ianfle(idum),iendle(idum),(wertle(idum,k),k=1,3)
         read(iunit,900,end=100) zeil
         read(zeil,915)                          (wertle(idum,k),k=4,6)
      else
         read(zeil,910)ianfle(idum),iendle(idum),wertle(idum,1)
      end if
      inum=iendle(idum)-ianfle(idum)+1
      nknowe=nknowe+inum
      call iset(inum,1,ifrgeo(ianfle(idum)),1)
      do 20 k=1,naniso
         call dset(inum,wertle(idum,k),pergeo(ianfle(idum),k),1)
   20 continue
      write(*,920) ianfle(idum),iendle(idum),(wertle(idum,k),k=1,naniso)
      goto 10
  100 continue
      call qfclos(iunit,0)
      numlei=idum
!
!      if (ianfle(1).ne.1) write(*,930) 1
!      do 30 i=1,idum-1
!         if (iendle(i)+1.ne.ianfle(i+1)) write(*,930) i
!   30 continue
!      if (iendle(idum).ne.nelgeo) write(*,930) idum
!
      return
  900 format(a)
  910 format(2i10,3e13.5)
  915 format( 20x,3e13.5)
  920 format(' Elements ',i6,' to ',i6,' conductivity',6g13.5)
  930 format(' Warning: bad definition of conductivity ',i3)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!      subroutine modnum(ipoint,iordnu,npoint,nordnu)
!      dimension ipoint(npoint),iordnu(nordnu)
!
!      icount=0
!      do 10 i=1,npoint
!         if (ipoint(i).ne.0) then
!            icount=icount+1
!            iordnu(icount)=i
!         end if
!   10 continue
!
!      return
!      end
!---->---1---------2---------3---------4---------5---------6---------7--<
! aa 21/01/2002 the order of the channels are not kept in the old convention
! we change the label in ipoint to the channel number
! now nodes in iordnu are organised in the same order as the channels
! node numbering in fortran convention : first node has number 1!!
      subroutine modnum(ipoint,iordnu,npoint,nordnu)
      dimension ipoint(npoint),iordnu(nordnu)
!
      do 10 i=1,npoint
         if (ipoint(i).ne.0) then
            iordnu(ipoint(i))=i
         end if
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function chifrm(i)
      character*80 chifrm
!
      chifrm(1:2)='(i'
      idum=int(log10(real(i)))+1
      write(chifrm(3:8),'(i6)') idum
      chifrm(9:9)=')'
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine percen(iakt,iges)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0,z1=1.d0,z2=2.d0,z100=100.d0)
      save idone
!
      if (iakt.eq.0) then
         write(*,890)
         write(*,900)
         call qdtext(' >',1)
         idone=0
      else if (iakt.lt.iges) then
         numpkt=50*iakt/iges-idone
         do 10 i=1,numpkt
            call qdtext('.',1)
   10    continue
         idone=idone+numpkt
      else if (iakt.eq.iges) then
         numpkt=50*iakt/iges-idone
         do 20 i=1,numpkt
            call qdtext('.',1)
   20    continue
         call qdtext('<',1)
         write(*,910) ' OK'
      end if
!
      return
  890 format(' Completed percentage of this task indicated by dots')
  900 format(' 0---10---20---30---40---50---60---70---80---90--100')
  910 format(a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine rodchk(knogeo,itygeo,necgeo,xyzgeo)
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo)
      dimension xyzgeo(npogeo,ndmgeo)
      dimension difvec(3)
!
!---check rod elements
      icount=0
      do 10 iel=1,nelgeo
         if (itygeo(iel).eq.201.or.itygeo(iel).eq.301) then
            kn1=knogeo(1,iel)
            kn2=knogeo(2,iel)
            difges=z0
            do 20 i=1,ndmgeo
               difvec(i)=xyzgeo(kn2,i)-xyzgeo(kn1,i)
               difges=difges+difvec(i)*difvec(i)
   20       continue
            difges=sqrt(difges)
            if (difges.lt.wuzeps) icount=icount+1
         end if
   10 continue
!
      if (icount.gt.0) then
         write(*,900) icount
         stop 'volchk'
      end if
!
      return
  900 format(' Number of improper beam elements ',i8)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wrieig(valeig,veceig,vallnz,veclnz,lanrel)
      include 'neurofem.inc'
      dimension valeig(nvala,nustei),veceig(npogeo,nvala)
      dimension lanrel(lenval)
      dimension vallnz(lenval,nustei),veclnz(npogeo,lenval)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam(31),'un','fo',ierr)
      if (ierr.ne.0) then
         write(*,'(a)') 'error filnam(31) '
         stop 'wrieig'
      end if
      write(iunit,900)
      write(    *,900)
      if (llanzp) then
         do 5 i=1,nwriln
            parsum=parsum+vallnz(i,1)
            write(iunit,910) i,(vallnz(i,j), j=1,nustei),parsum
            write(    *,910) i,(vallnz(i,j), j=1,nustei),parsum
    5    continue
      else
         parsum=z0
         do 10 i=1,nwriln
            parsum=parsum+valeig(i,1)
            write(iunit,910) i,(valeig(i,j), j=1,nustei),parsum
            write(    *,910) i,(valeig(i,j), j=1,nustei),parsum
   10    continue
      end if
      call qfclos(iunit,0)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam(32),'un','un',ierr)
      if (ierr.ne.0) then
         write(*,'(a)') 'error filnam(32) '
         stop 'wrieig'
      end if
      write(*,920) 'writing file: ',filnam(32)(1:qclen(filnam(32)))
      if (llanzp) then
         write(iunit)  1
         write(iunit)  nwriln
         write(iunit)  npogeo
         write(iunit)((vallnz(i,j),i=1, nwriln),j=1,2    )
         write(iunit)((veclnz(i,lanrel(j)),i=1,npogeo),j=1,nwriln)
      else
         write(iunit)  0
         write(iunit)  nwriln
         write(iunit)  npogeo
         write(iunit)((valeig(i,j),i=1,nvala),j=1,2)
         write(iunit)((veceig(i,j),i=1,npogeo),j=1,nvala)
      end if
      write(*,920) 'wrote   file: ',filnam(32)(1:qclen(filnam(32)))
      call qfclos(iunit,0)
!
      return
  900 format('# eigenwert res-nrm acc-est err-est')
  910 format(1x,i5,5(1x,g12.5,:))
  920 format(a,a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaeig(valeig,veceig,vallnz,veclnz,icall)
      include 'neurofem.inc'
      dimension valeig(nvala,nustei),veceig(npogeo,nvala)
      dimension vallnz(numval,nustei),veclnz(npogeo,numvec)
!
      if (.not.qfilda(filnam(32))) return
      iunit=qfreeu()
      call qfopse(iunit,filnam(32),'un','un',ierr)
      if (ierr.ne.0) then
         write(*,'(a)') 'error filnam(32) '
         stop 'reaeig'
      end if
      if (icall.eq.0) then
         read(iunit)  ndum
         if (llanzp)  then
            nsol=1
         else
            nsol=0
         end if
         if (ndum.ne.nsol.and..not.leigen) then
            write(*,950)
            if(nsol.eq.1) llanzp=.true.
            if(nsol.eq.0) llanzp=.false.
         end if
         read(iunit)  nwriln
      else
         write(*,920) 'reading file: ',filnam(32)(1:qclen(filnam(32)))
         read(iunit)  ldummy
         read(iunit)  nwriln
         if (nwriln.gt.nvala) then
            write(*,930)
            stop 'reaeig'
         end if
         read(iunit) ndummy
         if (ndummy.ne.npogeo) stop 'npogeo >< ndummy in reaeig'
         if (llanzp) then
            read(iunit)((vallnz(i,j),i=1,nwriln),j=1,2)
            read(iunit)((veclnz(i,j),i=1,npogeo),j=1,nwriln)
         else
            read(iunit)((valeig(i,j),i=1,nwriln),j=1,2)
            read(iunit)((veceig(i,j),i=1,npogeo),j=1,nwriln)
         end if
         write(*,920) 'read    file: ',filnam(32)(1:qclen(filnam(32)))
      end if
      call qfclos(iunit,0)
      if (nwriln.lt.nmodes) then
         write(*,940) nwriln
         nmodes=nwriln
      end if
!
      return
  900 format('# eigenwert res-nrm acc-est err-est')
  910 format(1x,i5,4(1x,g12.5))
  920 format(a,a)
  930 format(' ERROR (!!!) nwriln > nvala in reaeig')
  940 format(' WARNING: nmodes redefined to ',i6)
  950 format(' WARNING: redefined name of the eigenvector package')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine sellod(ipoasy,indasy,mrkvec,neigeo,isello,             &
     &           xyzgeo,xyzflf)
      include 'neurofem.inc'
      dimension ipoasy(npogeo+1),indasy(lenasy),mrkvec(maxnei),         &
     &          neigeo(npoflf),isello(npogeo)
      dimension xyzgeo(npogeo,ndmgeo),xyzflf(npoflf,ndmflf)
!
      call iset(npogeo,0,isello,1)
      do 10 iknflf=1,npoflf
         call neibor(ipoasy,indasy,mrkvec,xyzgeo,xyzflf,                &
     &        neigeo(iknflf),iknflf,numnei)
         do 20 j=1,numnei
            isello(mrkvec(j))=1
   20    continue
   10 continue
!
      numlod=isum(npogeo,isello,1)
      rel=dble(numlod)/dble(npogeo)*z100
      write(*,900) numlod,npogeo,rel
!
      return
  900 format(' Number of possible load nodes in volume structure ',i8,  &
     &     /,' compared to total number of nodes ',i8,' makes ',g12.5,  &
     &       ' [%]')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine setmeg(ifrmeg,nummeg)
      implicit double precision (a-h,o-z)
      dimension ifrmeg(nummeg)
!
      do 10 i=1,nummeg
         ifrmeg(i)=i
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine curtim(ichopt,nodinf,xyzflf,vcnflf,diptim,             &
     &                  filakt,thresh)
      include 'neurofem.inc'
      dimension ichopt(numdip),nodinf(ninfpo)
      dimension xyzflf(npoflf,ndmflf),vcnflf(npoflf,ndmflf),            &
     &          diptim(numsol,numtim)
      dimension dipdum(3)
!
! -- set itim = 1; for further developing eliminate dimension diptim(numsol,numtim)
!                  to dimension(numsol)  
      itim = 1
!
! -- determine maximum value of strength
      valmax=z0
      do 110 i=1,numdip
         valakt=z0
         do 120 k=1,ndmflx
            j=(i-1)*ndmflx+k
            valakt=valakt+diptim(j,itim)*diptim(j,itim)
  120    continue
         valmax=max(valmax,valakt)
  110 continue
      valmax=sqrt(valmax)
!
! -- count values above threshhold
      icount=0
      do 130 i=1,numdip
         valakt=z0
         do 140 k=1,ndmflx
            j=(i-1)*ndmflx+k
            valakt=valakt+diptim(j,itim)*diptim(j,itim)
  140    continue
         valakt=sqrt(valakt)
         if (valakt.gt.thresh*valmax) icount=icount+1
  130 continue
!
      iunit=qfreeu()
      call qfopse(iunit,filakt,'un','fo',ierr)
      if (ierr.gt.0) stop 'curtim: qfopse'
      write(iunit,900) 1,z1,z2
      write(iunit,910) icount,0,0,0,f5
!
      do 30 i=1,numdip
         kno=nodinf(ichopt(i))
         valakt=z0
         do 150 k=1,ndmflx
            j=(i-1)*ndmflx+k
            valakt=valakt+diptim(j,itim)*diptim(j,itim)
  150    continue
         valakt=sqrt(valakt)
         if (valakt.gt.thresh*valmax) then
            j=(i-1)*ndmflx
            call dset(3,z0,dipdum,1)
            if (lnorco) then
               do 400 k=1,ndmflf
                  dipdum(k)=diptim(i,itim)*vcnflf(kno,k)
  400          continue
            else
               do 410 k=1,ndmflf
                  dipdum(k)=diptim(j+k,itim)
  410          continue
            end if
            if( .not. ldipol ) then
               dipdum(2)=z0
               dipdum(3)=z0
            end if               
            if (ndmflf.eq.1) then
               write(iunit,912) (xyzflf(kno,k),k=1,ndmflf),             &
     &              (dipdum(k),k=1,ndmflf)
            else if (ndmflf.eq.2) then
               write(iunit,915) (xyzflf(kno,k),k=1,ndmflf),             &
     &              (dipdum(k),k=1,ndmflf)
            else
               write(iunit,920) (xyzflf(kno,k),k=1,ndmflf),             &
     &              (dipdum(k),k=1,ndmflf)
            end if
         end if
   30 continue
      write(iunit,930) f2,f2
      call qfclos(iunit,0)
!
      return
  900 format(i3,1x,f6.1,1x,f6.1)
  910 format(i4,1x,i3,1x,i3,1x,i3,1x,f6.1)
  912 format(1(f6.1,1x),1(e12.5,1x))
  915 format(2(f6.1,1x),2(e12.5,1x))
  920 format(3(f6.1,1x),3(e12.5,1x))
  930 format(2(f6.1,1x))
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function gasdev(iniran,cutval,jcut)
!
      implicit double precision(a-h,o-z)
      real ran2
      save iset
      data iset/0/
!
  200 if (iset.eq.0) then
  100    v1=2.d0*dble(ran2(iniran))-1.d0
         v2=2.d0*dble(ran2(iniran))-1.d0
         r=v1*v1+v2*v2
         if (r.ge.1.d0.or.r.eq.0.d0) goto 100
         fac=sqrt(-2.d0*log(r)/r)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      end if
!
!---  falls gausscut: abschneiden
      if ((jcut.eq.1).and.(gasdev.gt.cutval)) goto 200
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine chknum(intvek,numvek,minent,maxent,ierr)
      implicit double precision(a-h,o-z)
      dimension intvek(numvek)
!
      ierr=0
      do 10 i=1,numvek
         if (intvek(i).lt.minent.or.intvek(i).gt.maxent) ierr=1
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine cmpvek(iduknw,idukn2,npogeo,ierr)
      implicit double precision(a-h,o-z)
      dimension iduknw(npogeo),idukn2(npogeo)
!
      ierr=0
      do 10 i=1,npogeo
         if (iduknw(i).ne.idukn2(i)) ierr=1
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine chkflf(ityflf,nelflf,ierr)
      implicit double precision(a-h,o-z)
      dimension ityflf(nelflf)
!
      ierr=0
      do 10 i=1,nelflf
         if ((ityflf(i).lt.302).or.(mod(ityflf(i),10).ne.2)) ierr=1
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine chkfil(icall)
      include 'neurofem.inc'
!
      ierr=0
!     
!---checks before reading the command file
!
      if (icall.eq.1) then
!
!---files that are always neccesary
!
         if(.not.qfilda(filnam( 2))) then
            write(*,910) 'filnam( 2)= ',filnam( 2)(1:qclen(filnam( 2)))
            ierr=ierr+1
         end if
         if(.not.qfilda(filnam( 3))) then
            write(*,910) 'filnam( 3)= ',filnam( 3)(1:qclen(filnam( 3)))
            ierr=ierr+1
         end if
         if(.not.qfilda(filnam( 4))) then
            write(*,910) 'filnam( 4)= ',filnam( 4)(1:qclen(filnam( 4)))
            ierr=ierr+1
         end if
         if(.not.qfilda(filnam( 5))) then
            write(*,910) 'filnam( 5)= ',filnam( 5)(1:qclen(filnam( 5)))
            ierr=ierr+1
         end if
         if(.not.qfilda(filnam( 6))) then
            write(*,910) 'filnam( 6)= ',filnam( 6)(1:qclen(filnam( 6)))
            ierr=ierr+1
         end if
         if(.not.qfilda(filnam(12))) then
            write(*,910) 'filnam(12)= ',filnam(12)(1:qclen(filnam(12)))
            ierr=ierr+1
         end if
         if (ierr.gt.0) then
            write(*,920) ierr
            stop 'FILES ARE MISSING!'
         end if
!
!---checks after reading the command file
!
      else if (icall.eq.2) then
!
! files depending on input switches
!
         if (invway.eq.2) then
            if(.not.qfilda(filnam( 7))) then
               write(*,910) 'filnam( 7)= ',                             &
     &              filnam( 7)(1:qclen(filnam( 7)))
               ierr=ierr+1
            end if
         end if
         if (lalpha) then
            if(.not.qfilda(filnam(18))) then
               write(*,910) 'filnam(18)= ',                             &
     &              filnam(18)(1:qclen(filnam(18)))
               ierr=ierr+1
            end if
         end if
         if (logmeg) then
            if(.not.qfilda(filnam(35))) then
               write(*,910) 'filnam(35)= ',                             &
     &              filnam(35)(1:qclen(filnam(35)))
               ierr=ierr+1
            end if
            if(.not.qfilda(filnam(36))) then
               write(*,910) 'filnam(36)= ',                             &
     &              filnam(36)(1:qclen(filnam(36)))
               ierr=ierr+1
            end if
            if(.not.qfilda(filnam(37))) then
               write(*,910) 'filnam(37)= ',                             &
     &              filnam(37)(1:qclen(filnam(37)))
               ierr=ierr+1
            end if
            if(.not.qfilda(filnam(41))) then
               write(*,910) 'filnam(41)= ',                             &
     &              filnam(41)(1:qclen(filnam(41)))
               ierr=ierr+1
            end if
         end if
         if (ldipan) then
            if(.not.qfilda(filnam(27))) then
               write(*,910) 'filnam(27)= ',                             &
     &              filnam(27)(1:qclen(filnam(27)))
               ierr=ierr+1
            end if
         end if
         if (lrango) then
            if(.not.qfilda(filnam(42))) then
               write(*,910) 'filnam(42)= ',                             &
     &              filnam(42)(1:qclen(filnam(42)))
               ierr=ierr+1
            end if
         end if
!
! files depending on desired computing jobs
!
         if (lsimsr) then
            if(.not.qfilda(filnam(17))) then
               write(*,910) 'filnam(17)= ',                             &
     &              filnam(17)(1:qclen(filnam(17)))
               ierr=ierr+1
            end if
         end if
         if (lsimsr.and.ltrans) then
            if(.not.qfilda(filnam(22))) then
               write(*,910) 'filnam(22)= ',                             &
     &              filnam(22)(1:qclen(filnam(22)))
               ierr=ierr+1
            end if
            if(.not.qfilda(filnam(23))) then
               write(*,910) 'filnam(23)= ',                             &
     &              filnam(23)(1:qclen(filnam(23)))
               ierr=ierr+1
            end if
         end if
         if (lrefsl.and.(.not.lsimsr).and.(.not.lalpha)) then
            if(.not.qfilda(filnam( 8))) then
               write(*,910) 'filnam( 8)= ',                             &
     &              filnam( 8)(1:qclen(filnam( 8)))
               ierr=ierr+1
            end if
         end if
         if (linmat.or.linver.or.lsimsr.or.linsvd) then
            if(.not.qfilda(filnam(11))) then
               write(*,910) 'filnam(11)= ',                             &
     &              filnam(11)(1:qclen(filnam(11)))
               ierr=ierr+1
            end if
         end if
         if (linver.and.(.not.linmat)) then
            if(.not.qfilda(filnam(13))) then
               write(*,910) 'filnam(13)= ',                             &
     &              filnam(13)(1:qclen(filnam(13)))
               ierr=ierr+1
            end if
         end if
         if (linver.and.(.not.lrefsl)) then
            if(.not.qfilda(filnam( 9))) then
               write(*,910) 'filnam( 9)= ',                             &
     &              filnam( 9)(1:qclen(filnam( 9)))
               ierr=ierr+1
            end if
         end if
         if (linver.and.(.not.lrefsl).and.logmeg) then
            if(.not.qfilda(filnam(38))) then
               write(*,910) 'filnam(38)= ',                             &
     &              filnam(38)(1:qclen(filnam(38)))
               ierr=ierr+1
            end if
         end if
         if (linver.and.lcovar) then
            if(.not.qfilda(filnam(14))) then
               write(*,910) 'filnam(14)= ',                             &
     &              filnam(14)(1:qclen(filnam(14)))
               ierr=ierr+1
            end if
         end if
         if (linver.and.lcovar.and.logmeg) then
            if(.not.qfilda(filnam(39))) then
               write(*,910) 'filnam(39)= ',                             &
     &              filnam(39)(1:qclen(filnam(39)))
               ierr=ierr+1
            end if
         end if

! Files, depending on the chosen solver-method
! CW: 8.1.2001: 

! CW: 8.1.2001: AMG-CG needs input-file
         if (isolvr.eq.3) then
            if(.not.qfilda(filnam(44))) then
               write(*,910) 'filnam(44)= ',                             &
     &              filnam(44)(1:qclen(filnam(44)))
               ierr=ierr+1
            end if
         end if
         if (ierr.gt.0) then
            write(*,920) ierr
            stop 'FILES ARE MISSING!'
         end if
      end if
!
      return
  900 format(a)
  910 format('Attention! File: ',a,a, ' is missing!')
  915 format('Attention! Solver-method not well defined!')
  920 format(1x,i3,' essential input files are missing!')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaten(filnam,markkn,werkno,npoin,iknowe,numfre)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      integer qfreeu,qclen
      character*80 zeil, filnam
      dimension markkn(npoin)
      dimension werkno(npoin,numfre),wread(6)
!
!---check implemented tensor sizes
      if (numfre.ne.6) then
         write(*,*)' until now only symmetric 3*3 tensors implemented'
         goto 99
      end if
!
!--initialize arrays
      call iset(npoin       , 0,markkn,1)
      call dset(npoin*numfre,z0,werkno,1)
      markfl=0
      iknowe=0
!
!---open file
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         write(*,*)' error while opening file: ',filnam(1:qclen(filnam))
         goto 99
      end if
!
!---read contents of file
   10 read (iunit,'(a)',end=30) zeil
      call qctolo(zeil)
      if(index(zeil,'boi - tensorvaluefile').gt.0) then
         markfl=1
         read (iunit,900,end=30) zeil
      end if
      if(index(zeil,'eoi - tensorvaluefile').gt.0) goto 100
!
      if ( index(zeil,'boi - tensor').gt.0 ) then
   20    read (iunit,900,end=30) zeil
         call qctolo(zeil)
         if ( index(zeil,'eoi - tensor').gt.0 ) goto 10
         read ( zeil,910)        iread,(wread(i),i=1,3)
         read (iunit,900,end=30) zeil
         read ( zeil,920)              (wread(i),i=4,6)
         iknowe=iknowe+1
         markkn(iread  ) = 1
         call dcopy(6,wread,1,werkno(iread,1),npoin)
         goto 20
      end if
      goto 10
   30 write(*,*)' tensor value file is corrupt'
      goto 99
!
 100  call qfclos(iunit,0)
!
      if (markfl.eq.0) then
         write(*,*)' Label boi - tensorvaluefile is missing!' 
         goto 99
      end if
      if (markfl.gt.1) then
         write(*,*)' Label boi - tensorvaluefile more than one time!' 
         goto 99
      end if
      if (isum(npoin,markkn,1).ne.iknowe) then
         write(*,*)' Warning: at least one node is assigned twice to!'
      end if
!
      return
!
   99 write(*,*)' stop in reaten'
      stop
!
  900 format(a)
  910 format(bn,i10,3(f13.0))
  920 format(bn,10x,3(f13.0))
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!     Call SELSRC: select sources according to -thrdip-
      subroutine selsrc( mrkinv, dipinv, thrdip, npoflf, ndmflf )
      implicit double precision (a-h,o-z)
      character*80 filnam
      parameter (z0=0.d0)
      dimension mrkinv(npoflf,ndmflf), dipinv(npoflf,ndmflf)
!
! select maximum value
      valmax=z0
      do 110 i=1,npoflf
         valakt=z0
         do 120 k=1,ndmflf
            valakt=valakt+dipinv(i,k)*dipinv(i,k)
  120    continue
         valmax=max(valmax,valakt)
  110 continue
      valmax=sqrt(valmax)
!
! cancel pointer entries for values less then -thrdip*valmax-
      icount=0
      do 130 i=1,npoflf
         valakt=z0
         do 140 k=1,ndmflf
            valakt=valakt+dipinv(i,k)*dipinv(i,k)
  140    continue
         valakt=sqrt(valakt)
         if (valakt.lt.thrdip*valmax) then
            do k = 1,ndmflf
               mrkinv(i,k) = 0
            end do
         end if
  130 continue
!
      RETURN
      END
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine version
      write(*,900)
 900  format(1x,' Program CAUCHY',/,/,                                  &
     &       1x,'        VERSION: 1.8.4   -Oktober 1997-',//,           &
     &       1x,' Warning: This is an experimental code',/,             &
     &       1x,'          it may contain errors',//,                   &
     &       1x,'          PLEASE READ THE DOCUMENTATION ! ',/)
!
      RETURN
      END
!---->---1---------2---------3---------4---------5---------6---------7--<
! If i is contained in vector ivec, nzero will be true and the position 
! of i in ivec will be assigned to ipos.
! INPUT: 
!     ivec(n): vector with column indices of non-zero elements
!     i:       index to be tested
! Return value:
!      0: i is not contained in ivec
!      1: i is contained in ivec
! OUTPUT: 
!     ipos: position of i in the column vector ivec, 0 if not contained
! Version 1.0, C.Wolters (24.5.2000)
!=======================================================================
      function nzero(i,ivec,n,ipos)
      dimension ivec(n)
      integer ipos,i
!
!---is i contained in ivec?
      nzero=0
      ipos=0
      do 10 k=1,n
      if ( ivec(k).eq.i ) then
        ipos=k
        nzero=1
        goto 100
      end if
   10 continue
  100 return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
! Transform symmetric compact storage format to unsymmetric compact 
! storage format for fast solvers. Because "subroutine stous" was 
! already called before to determine "maxnei", idiasy and indasy are
! already known, but will be tested for correctness.
! INPUT: 
!     ipodia(npogeo): point. to the diag. el. of each row, symm.format
!     indexj(lengeo): column indices of non-zero el., symm.format
!     sysmat(lengeo): stiffness matrix in symm. compact storage
!     idiasy(npogeo): point. to the diag. el. of each row, unsymm.format
!     indasy(lenasy): column indices of non-zero el., unsymm.format
! OUTPUT: 
!     sysasy(lenasy): stiffness matrix in unsymm. compact storage
! Version 1.0, C.Wolters (24.5.2000)
!=======================================================================
      subroutine stouma(ipodia, indexj, sysmat, idiasy, indasy, sysasy)
      include 'neurofem.inc'
      dimension ipodia(npogeo),indexj(lengeo),sysmat(lengeo)
      dimension idiasy(npogeo),indasy(lenasy),sysasy(lenasy)
!
!-----row by row, the full matrix will be written out
      iposiu=0
      cpuvol=qscput(9,0,ierr)
      call percen(0,npogeo)
      do 100 iglo=1,npogeo
!-----storage of diagonal element first
         iposiu=iposiu+1
         if (idiasy(iglo).ne.iposiu) then
            write(*,*) 'Error: idiasy(iglo).ne.iposiu'
            stop 'stouma: error'
         end if
         if (indasy(iposiu).ne.iglo) then
            write(*,*) 'Error: indasy(iposiu).ne.iglo'
            stop 'stouma: error'
         end if
         sysasy(iposiu)=sysmat(ipodia(iglo))
!-----storage of elements at the left side of the diagonal element
         if (iglo.ne.1) then         
            ia=ipodia(iglo-1)+1
            in=ipodia(iglo)-ipodia(iglo-1)-1
            do 150 icol=ia,(ia+in-1)
               iposiu=iposiu+1
               if (indasy(iposiu).ne.indexj(icol)) then
                  write(*,*) 'Error: indasy(iposiu).ne.indexj(icol)'
                  stop 'stouma: error'
               end if
               sysasy(iposiu)=sysmat(icol)
  150       continue
         end if
!-----storage of elements at the right side of the diagonal element
         do 200 irow=(iglo+1),npogeo
            if (indasy(iposiu+1).eq.irow) then
               ia=ipodia(irow-1)+1
               in=ipodia(irow)-ipodia(irow-1)
               if ((nzero(iglo,indexj(ia),in,ipos)).eq.1) then
                  iposis=ia+ipos-1
                  iposiu=iposiu+1
                  if (indasy(iposiu).ne.irow) then
                     write(*,*) 'Error: indasy(iposiu).ne.irow'
                     stop 'stouma: error'
                  end if
                  sysasy(iposiu)=sysmat(iposis)
               end if
            end if
  200    continue
         call percen(iglo,npogeo)
  100 continue     
      write(*,900) qscput(9,1,ierr)
  900 format (' CPU--time for  stouma: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
! Transform symmetric compact storage format to unsymmetric compact 
! storage format for fast solvers. Counts the nonzero-elements for each 
! block of a node-wise partitioning (for the PILUTS-solver). 
! Because "subroutine stous" was 
! already called before to determine "maxnei", idiasy and indasy are
! already known, but will be tested for correctness.
! INPUT: 
!   ipodia(npogeo): point. to the diag. el. of each row, symm.format
!   indexj(lengeo): column indices of non-zero el., symm.format
!   sysmat(lengeo): stiffness matrix in symm. compact storage
!   idiasy(npogeo): point. to the diag. el. of each row, unsymm.format
!   indasy(lenasy): column indices of non-zero el., unsymm.format
!   iblock(1:2*nprocs): index of first and last row per processor
! OUTPUT: 
!     sysasy(lenasy): stiffness matrix in unsymm. compact storage
!     isumnz(1:num_procs): Number of nonzeros of the blocks
!     
! Version 1.0, C.Wolters (3.7.2001) For PILUTS-Coupling
!=======================================================================
      subroutine stoupi(ipodia, indexj, sysmat, idiasy, indasy, sysasy,
     &           nprocs,iblock,isumnz)
      include 'neurofem.inc'
      dimension ipodia(npogeo),indexj(lengeo),sysmat(lengeo)
!      dimension idiasy(npogeo),indasy(lenasy),sysasy(lenasy)
!      dimension iblock(2*nprocs),isumnz(nprocs)
      dimension idiasy(npogeo),indasy(*),sysasy(*)
      dimension iblock(*),isumnz(*)
!
!-----row by row, the full matrix will be written out
      iposiu=0
      cpuvol=qscput(9,0,ierr)
      call percen(0,npogeo)
      do 50 i=1,nprocs
         isumnz(i)=0
         do 100 iglo=iblock(2*i-1),iblock(2*i)
!-----storage of diagonal element first
            iposiu=iposiu+1
            if (idiasy(iglo).ne.iposiu) then
!               write(*,*) 'Error: idiasy(iglo).ne.iposiu'
               stop 'stoupi: error'
            end if
            if (indasy(iposiu).ne.iglo) then
               write(*,*) 'Error: indasy(iposiu).ne.iglo'
               stop 'stoupi: error'
            end if
            sysasy(iposiu)=sysmat(ipodia(iglo))
            isumnz(i)=isumnz(i)+1
!-----storage of elements at the left side of the diagonal element
            if (iglo.ne.1) then         
               ia=ipodia(iglo-1)+1
               in=ipodia(iglo)-ipodia(iglo-1)-1
               do 150 icol=ia,(ia+in-1)
                  iposiu=iposiu+1
                  if (indasy(iposiu).ne.indexj(icol)) then
                     write(*,*) 'Error: indasy(iposiu).ne.indexj(icol)'
                     stop 'stoupi: error'
                  end if
                  sysasy(iposiu)=sysmat(icol)
                  isumnz(i)=isumnz(i)+1
  150          continue
            end if
!-----storage of elements at the right side of the diagonal element
            do 200 irow=(iglo+1),npogeo
               if (indasy(iposiu+1).eq.irow) then
                  ia=ipodia(irow-1)+1
                  in=ipodia(irow)-ipodia(irow-1)
                  if ((nzero(iglo,indexj(ia),in,ipos)).eq.1) then
                     iposis=ia+ipos-1
                     iposiu=iposiu+1
                     if (indasy(iposiu).ne.irow) then
                        write(*,*) 'Error: indasy(iposiu).ne.irow'
                        stop 'stoupi: error'
                     end if
                     sysasy(iposiu)=sysmat(iposis)
                     isumnz(i)=isumnz(i)+1
                  end if
               end if
  200       continue
            call percen(iglo,npogeo)
  100    continue  
   50 continue   
      write(*,900) qscput(9,1,ierr)
  900 format (' CPU--time for  stoupi: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
! Write out equation system in AMG-format.
! Version 1.0, C.Wolters (24.5.2000)
!=======================================================================
      subroutine amgsav(fequa,fsolv,idiasy,indasy,sysasy,vecbgo,vecsol)
      include 'neurofem.inc'
      dimension vecsol(npogeo),vecbgo(npogeo)
      dimension idiasy(npogeo),indasy(lenasy),sysasy(lenasy)
      character*(*) fequa,fsolv
!      integer qfreeu,qfrecl
!---write out equation system in AMG-format
!
!---write out CAUCHY's solution vector in AMG-format
      lrerpl=qfrecl(1,1,0)
!.....Anzahl der Records im Ergebnisfile
      nrerpl=1+npogeo
!.....open binary save-file
      nterpl=qfreeu()
      call qfsize(nrerpl)
      call qfcrdi(nterpl,fsolv,'un',lrerpl,ierror)
      if (ierror.ne.0)  then
         write(*,8001) 'Error when writing',fsolv(1:qclen(fsolv))
         stop 'amgsav: error'
      endif
!.....write number of unknowns
      write(nterpl,rec=1) npogeo
!.....write out solution vector in AMG-format: row, value
      do 200 irow=1,npogeo
         write(nterpl,rec=irow+1,err=9990) irow,vecsol(irow)
 200  continue
      call qfclos(ntlese,0)
!
!----write out sysasy and (-vecbgo) in AMG-format
      lrerpl=qfrecl(2,1,0)
!.....Anzahl der Records im Ergebnisfile
      nrerpl=2+lenasy+npogeo
!.....open binary save-file
      nterpl=qfreeu()
      call qfsize(nrerpl)
      call qfcrdi(nterpl,fequa,'un',lrerpl,ierror)
      if (ierror.ne.0)  then
         write(*,8002) fequa
         stop 'amgsav: error'
      endif
!.....write number of unknowns
      irec=1
      write(nterpl,rec=irec) npogeo
!.....write number of non-zero matrix entries
      irec=irec+1
      write(nterpl,rec=irec) lenasy
!.....write out stiffness-matrix in AMG-format: row, column, value
      ilength=1
      do 300 irow=1,npogeo
         if (irow.eq.npogeo) then
            nonzero=lenasy-idiasy(irow)+1
         else
            nonzero=idiasy(irow+1)-idiasy(irow)
         end if
   	 do 400 icol=ilength,(ilength+nonzero-1)
            write(nterpl,rec=(icol+2),err=9992) irow,indasy(icol),      &
     &            sysasy(icol)
 400     continue
         ilength=ilength+nonzero
 300  continue
!----write out right hand side: row, value (assuming sysasy*x=vecbgo)
      ilength=ilength+1
      do 500 irow=1,npogeo
         write(nterpl,rec=(ilength+irow),err=9992) irow,                &
     &                                             vecbgo(irow)
 500     continue
      call qfclos(ntlese,0)
!---End of equation system output for AMG
!
      return
 8002 format(1x,'error while writing the result file ',a)
 9992 stop 'wriamg: error while writing of a result record '
 8001 format(a,a)
 9990 stop 'wrisol: error while writing of a result record for AMG-test'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
! Write out equation system in format for parallelized Krylov-solvers.
! Version 1.0, U.Hartmann, C.Wolters (25.5.2000)
!=======================================================================
!aa 25.01.2002 write out not needed ... difficult to integrate in linux and sgi...
! can be integrated in a final version 
      subroutine krysav(fequa,idiasy,indasy,sysasy,vecb,vecsol)
      include 'neurofem.inc'
      dimension vecsol(npogeo),vecb(npogeo)
      dimension idiasy(npogeo),indasy(lenasy),sysasy(lenasy)
      character*(*) fequa
!      integer qfreeu,qfrecl
!
! fd: file descriptor (C stream I/O)            
      integer fd
      integer sizeb
      integer pos
      integer num
      integer np1
      parameter (isizeb = 4)
      parameter (dsizeb = 8)
!
      double precision sta(npogeo)
!
! setup of test right hand side and start vector
      do i=1,npogeo
         sta(i)=0.0
      end do
!
! Opening file
!aa 25.01.2002 write out not needed ... difficult to integrate in linux and sgi...
! can be integrated in a final version 
!      call mopen(fd,fequa)  
!
! Writing the order of the matrix
      sizeb = isizeb
      num = 1
      pos = 20
!      call i_write(npogeo, sizeb, num, pos, fd)
!
! Writing maximal number of nonzeroes/row. This has been calculated 
! in subroutine neichk in cauch6.f and stored in 
! maxnei (common /geo1/ in neurofem.inc) in subroutine prepro.
      sizeb = isizeb
      num = 1
      pos = 24
!      call i_write(maxnei, sizeb, num, pos, fd)
!
! Writing diag positions                
      sizeb = isizeb
      pos = 28
      num = npogeo
!      call i_write(idiasy, sizeb, num, pos, fd)
!
! Writing the number of the non-zeros + 1
      sizeb = isizeb
      num = 1
      pos = 28 + npogeo * 4
      np1 = lenasy + 1
!      call i_write(np1, sizeb, num, pos, fd)
!
! Writing the values of the non-zeros
      sizeb = dsizeb
      pos = 28 + ((npogeo+1) * 4)
      num = lenasy
!      call d_write(sysasy, sizeb, num, pos, fd)
!
! Writing the column indices
      sizeb = isizeb
      pos = 28 + ((npogeo+1)*4) + (lenasy*8)
      num = lenasy
!      call i_write(indasy, sizeb, num, pos, fd)
!
! Writing the right hand side           
      sizeb = dsizeb
      num = npogeo
      pos = 32 + npogeo * 4 + lenasy * 12
!      call d_write(vecb, sizeb, num, pos, fd)
!
! Writing the start vector
      sizeb = dsizeb
      num = npogeo
      pos = 32 + npogeo * 12 + lenasy * 12
!      call d_write(sta, sizeb, num, pos, fd)            
!
      end             
!---->---1---------2---------3---------4---------5---------6---------7--<
! Write out equation system in Fast-solver-format. Right hand side of 
! the equation system is equal to -vecbgo, because within CAUCHY, 
! sysasy*vecsol+vecbgo=0 has been solved. 
! Version 1.0, C.Wolters (25.5.2000)
!=======================================================================
      subroutine wriequ(fequa,fsolv,ipodia,indexj,sysmat,vecsol,        &
     &                  vecbgo,vecb,idiasy,indasy,sysasy)
      include 'neurofem.inc'
      dimension ipodia(npogeo),indexj(lengeo),sysmat(lengeo)
      dimension vecsol(npogeo),vecbgo(npogeo),vecb(npogeo)
      dimension idiasy(npogeo),indasy(lenasy),sysasy(lenasy)
!
      character*(*) fequa,fsolv
!
!-----Transform the geometry matrix in symmetric compact storage format to 
!-----unsymmetric compact storage format.
      write(*,*) 'Transform stiff.-matrix to unsym. comp. stor. format:'
      call stouma(ipodia, indexj, sysmat, idiasy, indasy, sysasy)
!
!-----Right hand side in vecbgo is equal to -vecb, because within CAUCHY, 
!-----sysasy*vecsol+vecbgo=0 has been solved
      do 200 irow=1,npogeo
         vecb(irow)=-vecbgo(irow)
 200  continue
!
      if (ifsolv.eq.1) then
         write(*,*) 'Write out in AMG-format'
!        A and b in fequa, x in fsolv
         call amgsav(fequa,fsolv,idiasy,indasy,sysasy,vecb,vecsol)
      else if (ifsolv.eq.2) then
         write(*,*) 'Write out in Krylov-format'
!        A, b and x in fequa
         call krysav(fequa,idiasy,indasy,sysasy,vecb,vecsol)
      end if
!
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!---C.Wolters: Generation of an anisotropic sphere model from an isotropic
!---           one
!---5.7.2000
      subroutine genani(knogeo,xyzgeo,pergeo,necgeo,filani)
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo)
      dimension xyzgeo(npogeo,ndmgeo)
      dimension pergeo(nelgeo,1)
      dimension necgeo(nelgeo)
      character*80 filani
      dimension vcenter(3),vradi(3),vtang1(3),vtang2(3),vperm(3)
      dimension vec(3)
!
!---open output tensor-file and write header
      iunit=qfreeu()
      call qfopse(iunit,filani,'un','fo',ierr)
!
      if (ierr.gt.0) then
         write(*,930) ierr
         stop
      endif
!---write header
      write (iunit,940) 
      write (iunit,955)
      write (iunit,955) 
      write (iunit,950)
!
!     read coordinates of the center of the sphere
      write(*,*) 'Please enter the x-coordinate of the sphere-center-poi&
     &nt:'
      read*, vcenter(1)
      write(*,*) 'Please enter the y-coordinate of the sphere-center-poi&
     &nt:'
      read*, vcenter(2)
      write(*,*) 'Please enter the z-coordinate of the sphere-center-poi&
     &nt:'
      read*, vcenter(3)
!
!     read radial and tangential conductivities of the skull
!     sphere
      write(*,*) 'Please enter the radial conduct. of the skull-sphere'
      read*, vperm(1)
      write(*,*) 'Please enter the tang1-conduct. of the skull-sphere'
      read*, vperm(2)
      write(*,*) 'Please enter the tang2-conduct. of the skull-sphere'
      read*, vperm(3)
!
!     threshold for skull conduc., each element with a smaller 
!     cond. is defined to be located in the skull and will be assigned
!     a skull-sphere conduct.-tensor
      write(*,*) 'Please enter the threshold cond., under which all'
      write(*,*) 'elements will be taken as skull-elements:'
      read*, vskull
!
!     initialize local variables
      iknowe=0
      iskull=0
      vtest=0.0
!     determination of maximal and minimal "skull" radii
      vmaxsk=0.0
      vminsk=bignum
!
      do 100 ielem=1,nelgeo
         if (pergeo(ielem,1).lt.vskull) then
!           the element is a skull-element:
            iskull=iskull+1
!
            if (necgeo(ielem).eq.8) then
!              its a cube element:
!              Determine inner and outer skull surface extremal points
               do 110 j=1,necgeo(ielem)
                  vnorm=0.0
                  do 120 i=1,3
                     vec(i)=xyzgeo(knogeo(j,ielem),i)-vcenter(i)
                     vnorm=vnorm+vec(i)*vec(i)
  120             continue
                  vnorm=sqrt(vnorm)
                  if (vnorm.gt.vmaxsk) then
                     vmaxsk=vnorm
                  else if (vnorm.lt.vminsk) then
                     vminsk=vnorm
                  end if
  110          continue
!
!              determine radial direction: 
!              Midpoint between node 1 and 7 is taken as an 
!              approximation to the centerpoint
               vnorm=0.0
               do 130 i=1,3
                  vradi(i)=(xyzgeo(knogeo(1,ielem),i)-vcenter(i))+      &
     &                      0.5*(xyzgeo(knogeo(7,ielem),i)-             &
     &                           xyzgeo(knogeo(1,ielem),i))
                  vnorm=vnorm+vradi(i)*vradi(i)
  130          continue
               vnorm=sqrt(vnorm)
               imax=0
               radmax=0.0
               do 140 i=1,3
                  vradi(i)=vradi(i)/vnorm
                  if (abs(vradi(i)).gt.abs(radmax)) then
                     radmax=vradi(i)
                     imax=i
                  end if
  140          continue
            else
!              not yet implemented
               stop
            end if
!
!           determine first tangential direction so that
!           <vradi,vtang1>=0
            if (imax.eq.3) then
               vtang1(1)=0.0
               vtang1(2)=1.0
               vtang1(3)=-vradi(2)/vradi(3)
            else if (imax.eq.2) then
               vtang1(1)=0.0
               vtang1(3)=1.0
               vtang1(2)=-vradi(3)/vradi(2)
            else 
               vtang1(2)=0.0
               vtang1(3)=1.0
               vtang1(1)=-vradi(3)/vradi(1)
            end if
            vnorm=sqrt(1.0+vtang1(imax)*vtang1(imax))       
            do 150 i=1,3
               vtang1(i)=vtang1(i)/vnorm
  150       continue
!
!           determine second tangential direction through vector
!           product vtang2 = vradi x vtang1
            vtang2(1)=vradi(2)*vtang1(3)-vradi(3)*vtang1(2)
            vtang2(2)=vradi(3)*vtang1(1)-vradi(1)*vtang1(3)
            vtang2(3)=vradi(1)*vtang1(2)-vradi(2)*vtang1(1)
            vnorm=sqrt(vtang2(1)*vtang2(1)+vtang2(2)*vtang2(2)+         &
     &                 vtang2(3)*vtang2(3))       
            do 160 i=1,3
               vtang2(i)=vtang2(i)/vnorm
  160       continue
!
!    test, if the vectors built an orthonormal basis of R^3
            vtest=vtang1(1)*vtang2(1)+vtang1(2)*vtang2(2)+              &
     &            vtang1(3)*vtang2(3)+                                  &
     &            vtang1(1)*vradi(1)+vtang1(2)*vradi(2)+                &
     &            vtang1(3)*vradi(3)+                                   &
     &            vtang1(1)*vtang2(1)+vtang1(2)*vtang2(2)+              &
     &            vtang1(3)*vtang2(3)+                                  &
     &            vradi(1)*vtang2(1)+vradi(2)*vtang2(2)+                &
     &            vradi(3)*vtang2(3)
            vtest=vtest+                                                &
     &            vtang1(1)*vtang1(1)+vtang1(2)*vtang1(2)+              &
     &            vtang1(3)*vtang1(3)+                                  &
     &            vtang2(1)*vtang2(1)+vtang2(2)*vtang2(2)+              &
     &            vtang2(3)*vtang2(3)+                                  &
     &            vradi(1)*vradi(1)+vradi(2)*vradi(2)+                  &
     &            vradi(3)*vradi(3) 
            vtest=vtest-3.0      
            if (abs(vtest).gt.wuzeps) then       
               write(*,995) 'Tensor eigenvectors are not orthogonal: Ele&
     &t ', ielem
            end if
!
!
!           The eigenvector matrix is: R=(vradi,vtang1,vtang2)
!           The conductivity tensor is then: sigma=R*D*R^t
            sig_xx = vradi(1)*vperm(1)*vradi(1) +                       &
     &               vtang1(1)*vperm(2)*vtang1(1)+                      &
     &               vtang2(1)*vperm(3)*vtang2(1)
            sig_yy = vradi(2)*vperm(1)*vradi(2) +                       &
     &               vtang1(2)*vperm(2)*vtang1(2)+                      &
     &               vtang2(2)*vperm(3)*vtang2(2)
            sig_zz = vradi(3)*vperm(1)*vradi(3) +                       &
     &               vtang1(3)*vperm(2)*vtang1(3)+                      &
     &               vtang2(3)*vperm(3)*vtang2(3)
            sig_xy = vradi(2)*vperm(1)*vradi(1) +                       &
     &               vtang1(2)*vperm(2)*vtang1(1)+                      &
     &               vtang2(2)*vperm(3)*vtang2(1)
            sig_xz = vradi(3)*vperm(1)*vradi(1) +                       &
     &               vtang1(3)*vperm(2)*vtang1(1)+                      &
     &               vtang2(3)*vperm(3)*vtang2(1)
            sig_yz = vradi(3)*vperm(1)*vradi(2) +                       &
     &               vtang1(3)*vperm(2)*vtang1(2)+                      &
     &               vtang2(3)*vperm(3)*vtang2(2)
!
!    test, if sigma*eigenvector-vperm*eigenvector=0
            vtest=sig_xx*vradi(1)+sig_xy*vradi(2)+                      &
     &            sig_xz*vradi(3)-vperm(1)*vradi(1)+                    &
     &            sig_xy*vradi(1)+sig_yy*vradi(2)+                      &
     &            sig_yz*vradi(3)-vperm(1)*vradi(2)+                    &
     &            sig_xz*vradi(1)+sig_yz*vradi(2)+                      &
     &            sig_zz*vradi(3)-vperm(1)*vradi(3)+                    &
     &            sig_xx*vtang1(1)+sig_xy*vtang1(2)+                    &
     &            sig_xz*vtang1(3)-vperm(2)*vtang1(1)+                  &
     &            sig_xy*vtang1(1)+sig_yy*vtang1(2)+                    &
     &            sig_yz*vtang1(3)-vperm(2)*vtang1(2)+                  &
     &            sig_xz*vtang1(1)+sig_yz*vtang1(2)+                    &
     &            sig_zz*vtang1(3)-vperm(2)*vtang1(3)+                  &
     &            sig_xx*vtang2(1)+sig_xy*vtang2(2)+                    &
     &            sig_xz*vtang2(3)-vperm(3)*vtang2(1)+                  &
     &            sig_xy*vtang2(1)+sig_yy*vtang2(2)+                    &
     &            sig_yz*vtang2(3)-vperm(3)*vtang2(2)+                  &
     &            sig_xz*vtang2(1)+sig_yz*vtang2(2)+                    &
     &            sig_zz*vtang2(3)-vperm(3)*vtang2(3)       
            if (abs(vtest).gt.wuzeps) then       
               write(*,995) 'Tensor wrong: Element ', ielem
            end if
!
!---        following page 130 of CAUCHY's user guide, 
!---        the tensor is written out
! ORDER HAS TO BE VARIFIED!!
            write(iunit,970) ielem,sig_xx,sig_yy,sig_zz
            write(iunit,980) sig_xy,sig_yz,sig_xz
            iknowe=iknowe+1
         else 
            write(iunit,970) ielem,pergeo(ielem,1),pergeo(ielem,1),     &
     &                        pergeo(ielem,1)
            write(iunit,980) 0.0,0.0,0.0
            iknowe=iknowe+1
         end if
  100 continue
!
!---write ending in iunit and close file
      write(*,995) 'Written tensorvalued conduct.:', iknowe
      write(*,995) 'Written skull sphere elements:', iskull
      write(*,996) 'Maximum skull node norm:', vmaxsk
      write(*,996) 'Minimum skull node norm:', vminsk
      write (iunit,991)
      write (iunit,955) 
      write (iunit,955)
      write (iunit,992)
      call qfclos(iunit,0)
!
  930 format (1x,' Error when opening file: IOSTAT= ',i3)
  940 format('BOI - TENSORVALUEFILE')
  950 format('BOI - TENSOR')
  955 format('========================================================')
  960 format(bn,3(1x,i7,2x,f12.0))
  970 format(bn,3x,i7,3(1x,g12.5))
  980 format(bn,10x,3(1x,g12.5))
  990 format('keyword for node-value file is missing!')
  991 format('EOI - TENSOR')
  992 format('EOI - TENSORVALUEFILE')
  995 format(a,i7)
  996 format(a,f12.0)
  998 format (1x,' Error when closing file: IOSTAT= ',i3)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
