#include "fintrf.h"
C	calc_stiff_matrix_val(nd,ele,cond,mele)
C	size(nd) = [#nodes,3]
C	size(ele) = [#nodes per ele,#ele]
C	size(cond) = [#ele,#aniso] #aniso = 1 / 6
C	size(mele) = [1,1]
C	class(nd) = double
C	class(ele) = int32
C	class(cond) = double
C	class(mele) = int32
      SUBROUTINE mexfunction(nlhs, plhs, nrhs, prhs)
      INCLUDE 'neurofem.inc'
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER :: plhs(*), prhs(*)

C-----------------------------------------------------------------------
C
      INTEGER :: nlhs, nrhs
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER :: mxcreatedoublematrix, mxgetpr, mxisint32,
     &mxcreatenumericarray, mxisdouble
      INTEGER :: mxclassidfromclassname

C-----------------------------------------------------------------------
C
      MWSIZE :: mxgetm, mxgetn
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------

      MWSIZE :: ierr,ninres,ndbres,maxele3d,ngrinod,ngriele,
     &allstat,nconduc,mconduc,mgriele,dims(2),nele,mgrinod
      INTEGER,ALLOCATABLE :: indexj(:)
      DOUBLE PRECISION,ALLOCATABLE :: ivor(:)
      MWPOINTER :: p_nodes,p_elem,p_conduc,
     &p_diinsy,p_symind,
     &p_sysmat,p_volmat,p_grinod,p_nele
      INTEGER,ALLOCATABLE :: elemty(:),elemnn(:)
      DOUBLE PRECISION,ALLOCATABLE :: volmat(:)
      INTEGER :: ntaszu,nterzu,nterro,complexflag,class,
     &ndim
      CHARACTER*200 :: line

      ierr = 0
      allstat = 0
      complexflag = 0
      ndim = 2

C     CHECK INPUT

      IF (nrhs>4) THEN
        CALL mexerrmsgtxt('more than FOUR arguments\n
     +Nodes, cells, conductivities, number of nodes')
      END IF
      IF (nrhs<4) THEN
        CALL mexerrmsgtxt('less than FOUR arguments\n
     +Nodes, cells, conductivities, number of nodes')
      END IF 

      if(mxisint32(prhs(2)) .eq. 0) then
         call mexerrmsgtxt ('Cells is wrong type\n
     +Input must be int32 array.')
      endif
C      if(mxisint32(prhs(4)) .eq. 0) then
C         call mexerrmsgtxt ('Number of nodes is wrong type\n
C     +Input must be int32.')
C      endif
      if(mxisint32(prhs(4)) .eq. 0) then
         call mexerrmsgtxt ('Size of elements is wrong type\n
     +Input must be int32.')
      endif
      if(mxisdouble(prhs(1)) .eq. 0) then
         call mexerrmsgtxt ('Nodes is wrong type\n
     +Input must be double array.')
      endif
      if(mxisdouble(prhs(3)) .eq. 0) then
         call mexerrmsgtxt ('Conductivities is wrong type\n
     +Input must be double array.')
      endif

C     ASSIGN INPUT POINTERS

      p_nodes = mxgetpr(prhs(1))
      p_elem = mxgetpr(prhs(2))
      p_conduc = mxgetpr(prhs(3))
C      p_grinod = mxgetpr(prhs(4))
      p_nele = mxgetpr(prhs(4))

C      CALL mxcopyptrtointeger4(p_grinod,ngrinod,1)
      CALL mxcopyptrtointeger4(p_nele,nele,1)
      ngrinod = mxgetm(prhs(1))
      mgrinod = mxgetn(prhs(1))
      ngriele = mxgetn(prhs(2))
      mgriele = mxgetm(prhs(2))
      nconduc = mxgetm(prhs(3))
      mconduc = mxgetn(prhs(3))
      maxele3d = 8
      ninres = ngrinod*mxelcs
      ndbres = ngrinod*mxelcs

C      write(line,*) mgrinod mconduc
C      CALL mexprintf(line//achar(10))

      if(mgrinod .ne. 3) then
         call mexerrmsgtxt ('Nodes has to be an array with
     +three columns.')
      endif
      if(ngriele .ne. nconduc) then
         call mexerrmsgtxt ('Conductivities and elements do not fit.')
      endif
      if((mconduc .ne. 6).and.(mconduc .ne. 1)) then
         call mexerrmsgtxt ('Conductivity has to be an array with
     +one or six columns.')
      endif

C     FILL COMMON BLOCKS

C     GridGeneratorVista::readgrid(), l.287

      numper = ngriele
      npogeo = ngrinod
      nelgeo = ngriele
      ndmgeo = 3
      npoflf = ngrinod
      nelflf = ngriele
      ndmflf = 3

C     GridGeneratorVista::readgrid(), l.471

      permeb = 0.00033

      IF (mconduc.eq.6) THEN
            naniso = 6
            laniso = .TRUE.
      ELSE
            naniso = 1
            laniso = .FALSE.
      END IF

C     ALLOCATE VARIABLES

      class = mxclassidfromclassname('int32')

      ALLOCATE(indexj(ninres), STAT = allstat)
      IF (allstat.ne.0) CALL mexerrmsgtxt('Not enough memory.')
      ALLOCATE(ivor(ndbres*2), STAT = allstat)
      IF (allstat.ne.0) CALL mexerrmsgtxt('Not enough memory.')
      ALLOCATE(elemty(ngriele), STAT = allstat)
      IF (allstat.ne.0) CALL mexerrmsgtxt('Not enough memory.')
      ALLOCATE(elemnn(ngriele), STAT = allstat)
      IF (allstat.ne.0) CALL mexerrmsgtxt('Not enough memory.')

      dims(1) = ngrinod
      dims(2) = 1
      plhs(1) = mxcreatenumericarray(ndim,dims,class,complexflag)

      p_diinsy = mxgetpr(plhs(1))

C     SET elemty
      IF (nele.eq.8) THEN
            elemty = 323
            elemnn = 8
      ELSE IF (nele.eq.4) THEN
            elemty = 303
            elemnn = 4
      ELSE
            CALL mexerrmsgtxt('Invalid number of nodes per element')
      END IF

C     INIT IME LIB AbstractSimulatorEEGMEGNeuroFEM, l.2998 / 105

      CALL qlbini(ntaszu,nterzu,nterro,cabb,labb)
      CALL fofini(lrelin,lreint)

C     SETUP SYMMETRIC STORAGE

      CALL cstoss(%val(p_elem),%val(p_diinsy),indexj,ivor,ngrinod,
     &ngriele,maxele3d,lengeo,ninres,ndbres,ierr)

      DEALLOCATE(ivor)

      dims(1) = lengeo
      dims(2) = 1
      plhs(2) = mxcreatenumericarray(2,dims,class,0)

      p_symind = mxgetpr(plhs(2))

      CALL icopy(lengeo,indexj,1,%val(p_symind),1)

      DEALLOCATE(indexj)

      lenflf = lengeo

C      maxnei = 27
C      ndimmt = 378
      nevlpo = 1

C      nevlpo = Number of Channels
      invway = 2
      numana = 1
      intgrd = 2
C      mgagam = intgrd
C      tolsol = .10000000D-9
C      thrdip = .10000000D+01
C      epsilo = .10000000D-08
C      igladi = 2
C      norder = 2
C      scadip = 20.000000
C      diplam = .10000000D-05
C      disdip = 1.0000000
C      ldipol = .TRUE.
C      lmncor = .FALSE.
      lrango = .FALSE.

C     END CONSTRUCTOR OF AbSimEEGMEGNFEM

C     AbSimEEGMEGNFEM, allocatecompilematrices / computevolumematrix

      ALLOCATE(volmat(lengeo))

      plhs(3) = mxcreatedoublematrix(lengeo,1,0)

      p_sysmat = mxgetpr(plhs(3))

      CALL calsti(%val(p_elem),elemty,elemnn,%val(p_diinsy),
     &%val(p_symind),%val(p_sysmat),volmat,%val(p_nodes),
     &%val(p_conduc),sysrng)

      DEALLOCATE(elemty)
      DEALLOCATE(elemnn)

      DEALLOCATE(volmat)

      RETURN
      END
