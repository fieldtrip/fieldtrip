      SUBROUTINE DTGEVC( JOB, SIDE, SELECT, N, A, LDA, B, LDB, VL, LDVL,
     $                   VR, LDVR, MM, M, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( N, * )
*     ..
*
*
*  Purpose
*  =======
*
*  DTGEVC computes selected left and/or right generalized eigenvectors
*  of a pair of real upper triangular matrices (A,B).  The
*  j-th generalized left and right eigenvectors are  y  and  x, resp.,
*  such that:
*       H                           H
*      y  (A - wB) = 0  or  (A - wB) y = 0   and    (A - wB)x = 0
*
*                                                                      H
*  Note: the left eigenvector is sometimes defined as the row vector  y
*        but DTGEVC computes the column vector y.
*  Reminder: the eigenvectors may be real or complex.  If complex, the
*        eigenvector for the eigenvalue w s.t. Im(w) > 0 is computed.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'A': compute All (left/right/left+right) generalized
*                 eigenvectors of (A,B);
*          = 'S': compute Selected (left/right/left+right) generalized
*                 eigenvectors of (A,B) -- see the description of the
*                 argument SELECT;
*          = 'B' or 'T': compute all (left/right/left+right) generalized
*                 eigenvectors of (A,B), and Back Transform them
*                 using the initial contents of VL/VR -- see the
*                 descriptions of the arguments VL and VR.
*
*  SIDE    (input) CHARACTER*1
*          Specifies for which side eigenvectors are to be computed:
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If JOB='S', then SELECT specifies the (generalized)
*          eigenvectors to be computed.  To get the eigenvector
*          corresponding to the j-th eigenvalue, set SELECT(j) to
*          .TRUE.  If the j-th and (j+1)-st eigenvalues are
*          conjugates, i.e., A(j+1,j) is nonzero, then only the
*          eigenvector for the first may be selected (the second being
*          just the conjugate of the first); this may be done by
*          setting either SELECT(j) or SELECT(j+1) to .TRUE.
*
*          If JOB='A', 'B', or 'T', SELECT is not referenced, and all
*          eigenvectors are selected.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          One of the pair of matrices whose generalized eigenvectors
*          are to be computed.  It must be block upper triangular, with
*          1-by-1 or 2-by-2 blocks on the diagonal, the 1-by-1 blocks
*          corresponding to real generalized eigenvalues and the 2-by-2
*          blocks corresponding to complex generalized eigenvalues.
*          The eigenvalues are computed from the diagonal blocks of A
*          and corresponding entries of B.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A.  LDA >= max(1, N).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The other of the pair of matrices whose generalized
*          eigenvectors are to be computed.  It must be upper
*          triangular, and if A has a 2-by-2 diagonal block in
*          rows/columns j,j+1, then the corresponding 2-by-2 block of B
*          must be diagonal with positive entries.
*
*  LDB     (input) INTEGER
*          The leading dimension of array B.  LDB >= max(1, N).
*
*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
*          On exit, the left eigenvectors (column vectors -- see the
*          note in "Purpose".)  Real eigenvectors take one column,
*          complex take two columns, the first for the real part and the
*          second for the imaginary part.
*          If JOB='A', then all left eigenvectors of (A,B) will be
*             computed and stored in VL.
*          If JOB='S', then only the eigenvectors selected by SELECT
*             will be computed, and they will be stored one right after
*             another in VL; the first selected eigenvector will go
*             in column 1 (and 2, if complex), the second in the next
*             column(s), etc.
*          If JOB='B' or 'T', then all left eigenvectors of (A,B)
*             will be computed and multiplied (on the left) by the
*             matrix found in VL on entry to DTGEVC.  Usually, this
*             will be the Q matrix computed by DGGHRD and DHGEQZ,
*             so that on exit, VL will contain the left eigenvectors
*             of the original matrix pair.
*          In any case, each eigenvector will be scaled so the largest
*          component of each vector has
*          abs(real part) + abs(imag. part)=1, *unless*  the diagonal
*          blocks in A and B corresponding to the eigenvector are both
*          zero (hence, 1-by-1), in which case the eigenvector will be
*          zero.
*          If SIDE = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of array VL.  LDVL >= 1; if SIDE = 'B'
*          or 'L', LDVL >= N.
*
*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
*          On exit, the right eigenvectors.  Real eigenvectors take one
*          column, complex take two columns, the first for the real
*          part and the second for the imaginary part.
*          If JOB='A', then all right eigenvectors of (A,B) will be
*             computed and stored in VR.
*          If JOB='S', then only the eigenvectors selected by SELECT
*             will be computed, and they will be stored one right after
*             another in VR; the first selected eigenvector will go
*             in column 1 (and 2, if complex), the second in the next
*             column(s), etc.
*          If JOB='B' or 'T', then all right eigenvectors of (A,B)
*             will be computed and multiplied (on the left) by the
*             matrix found in VR on entry to DTGEVC.  Usually, this
*             will be the Z matrix computed by DGGHRD and DHGEQZ,
*             so that on exit, VR will contain the right eigenvectors
*             of the original matrix pair.
*          In any case, each eigenvector will be scaled so the largest
*          component of each vector has
*          abs(real part) + abs(imag. part)=1, *unless*  the diagonal
*          blocks in A and B corresponding to the eigenvector are both
*          zero (hence, 1-by-1), in which case the eigenvector will be
*          zero.
*          If SIDE = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of array VR.  LDVR >= 1; if SIDE = 'B'
*          or 'R', LDVR >= N.
*
*  MM      (input) INTEGER
*          The number of columns in VL and/or VR.
*          If JOB='A', 'B', or 'T', then MM >= N.
*          If JOB='S', then MM must be at least the number of columns
*             required, as computed from SELECT.  Each .TRUE. value in
*             SELECT corresponding to a real eigenvalue (i.e., A(j+1,j)
*             and A(j,j-1) are zero) counts for one column, and each
*             .TRUE.  value corresponding to the first of a complex
*             conjugate pair (i.e., A(j+1,j) is not zero) counts for
*             two columns.  (.TRUE. values corresponding to the second
*             of a pair -- A(j,j-1) is not zero -- are ignored.)
*
*  M       (output) INTEGER
*          The number of columns in VL and/or VR actually
*          used to store the eigenvectors.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension ( N, 6 )
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex
*                eigenvalue.
*
*  Further Details
*  ===============
*
*  Allocation of workspace:
*  ---------- -- ---------
*
*     WORK( j, 1 ) = 1-norm of j-th column of A, above the diagonal
*     WORK( j, 2 ) = 1-norm of j-th column of B, above the diagonal
*     WORK( *, 3 ) = real part of eigenvector
*     WORK( *, 4 ) = imaginary part of eigenvector
*     WORK( *, 5 ) = real part of back-transformed eigenvector
*     WORK( *, 6 ) = imaginary part of back-transformed eigenvector
*
*
*  Rowwise vs. columnwise solution methods:
*  ------- --  ---------- -------- -------
*
*  Finding a generalized eigenvector consists basically of solving the
*  singular triangular system
*                                                    H
*     (A - w B) x = 0     (for right) or:   (A - w B) y = 0  (for left)
*
*  Consider finding the i-th right eigenvector (assume all eigenvalues
*  are real). The equation to be solved is:
*       n                   i
*  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1
*      k=j                 k=j
*
*  where  C = (A - w B)  (The components v(i+1:n) are 0.)
*
*  The "rowwise" method is:
*
*  (1)  v(i) := 1
*  for j = i-1,. . .,1:
*                          i
*      (2) compute  s = - sum C(j,k) v(k)   and
*                        k=j+1
*
*      (3) v(j) := s / C(j,j)
*
*  Step 2 is sometimes called the "dot product" step, since it is an
*  inner product between the j-th row and the portion of the eigenvector
*  that has been computed so far.
*
*  The "columnwise" method consists basically in doing the sums
*  for all the rows in parallel.  As each v(j) is computed, the
*  contribution of v(j) times the j-th column of C is added to the
*  partial sums.  Since FORTRAN arrays are stored columnwise, this has
*  the advantage that at each step, the entries of C that are accessed
*  are adjacent to one another, whereas with the rowwise method, the
*  entries accessed at a step are spaced LDA (and LDB) words apart.
*
*  When finding left eigenvectors, the matrix in question is the
*  transpose of the one in storage, so the rowwise method then
*  actually accesses columns of A and B at each step, and so is the
*  preferred method.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COMPL, COMPR, IL2BY2, ILABAD, ILALL, ILBACK,
     $                   ILBBAD, ILCOMP, ILCPLX, LSA, LSB
      INTEGER            I, IBEG, IEIG, IEND, IINFO, IJOB, IM, ISIDE, J,
     $                   JA, JC, JE, JR, JW, NA, NW
      DOUBLE PRECISION   ACOEF, ACOEFA, ANORM, ASCALE, BCOEFA, BCOEFI,
     $                   BCOEFR, BIG, BIGNUM, BNORM, BSCALE, CIM2A,
     $                   CIM2B, CIMAGA, CIMAGB, CRE2A, CRE2B, CREALA,
     $                   CREALB, DMIN, SAFMIN, SALFAR, SBETA, SCALE,
     $                   SMALL, TEMP, TEMP2, TEMP2I, TEMP2R, ULP, XMAX,
     $                   XSCALE
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   BDIAG( 2 ), SUM( 2, 2 ), SUMA( 2, 2 ),
     $                   SUMB( 2, 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DLABAD, DLACPY, DLAG2, DLALN2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      IF( LSAME( JOB, 'A' ) ) THEN
         IJOB = 1
         ILALL = .TRUE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         IJOB = 2
         ILALL = .FALSE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( JOB, 'B' ) .OR. LSAME( JOB, 'T' ) ) THEN
         IJOB = 3
         ILALL = .TRUE.
         ILBACK = .TRUE.
      ELSE
         IJOB = -1
         ILALL = .TRUE.
      END IF
*
      IF( LSAME( SIDE, 'R' ) ) THEN
         ISIDE = 1
         COMPL = .FALSE.
         COMPR = .TRUE.
      ELSE IF( LSAME( SIDE, 'L' ) ) THEN
         ISIDE = 2
         COMPL = .TRUE.
         COMPR = .FALSE.
      ELSE IF( LSAME( SIDE, 'B' ) ) THEN
         ISIDE = 3
         COMPL = .TRUE.
         COMPR = .TRUE.
      ELSE
         ISIDE = -1
      END IF
*
*     Count the number of eigenvectors to be computed
*
      IF( .NOT.ILALL ) THEN
         IM = 0
         ILCPLX = .FALSE.
         DO 10 J = 1, N
            IF( ILCPLX ) THEN
               ILCPLX = .FALSE.
               GO TO 10
            END IF
            IF( J.LT.N ) THEN
               IF( A( J+1, J ).NE.ZERO )
     $            ILCPLX = .TRUE.
            END IF
            IF( ILCPLX ) THEN
               IF( SELECT( J ) .OR. SELECT( J+1 ) )
     $            IM = IM + 2
            ELSE
               IF( SELECT( J ) )
     $            IM = IM + 1
            END IF
   10    CONTINUE
      ELSE
         IM = N
      END IF
*
*     Check 2-by-2 diagonal blocks of A, B
*
      ILABAD = .FALSE.
      ILBBAD = .FALSE.
      DO 20 J = 1, N - 1
         IF( A( J+1, J ).NE.ZERO ) THEN
            IF( B( J, J ).EQ.ZERO .OR. B( J+1, J+1 ).EQ.ZERO .OR.
     $          B( J, J+1 ).NE.ZERO )ILBBAD = .TRUE.
            IF( J.LT.N-1 ) THEN
               IF( A( J+2, J+1 ).NE.ZERO )
     $            ILABAD = .TRUE.
            END IF
         END IF
   20 CONTINUE
*
      INFO = 0
      IF( IJOB.LT.0 ) THEN
         INFO = -1
      ELSE IF( ISIDE.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILABAD ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( ILBBAD ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) THEN
         INFO = -10
      ELSE IF( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) THEN
         INFO = -12
      ELSE IF( MM.LT.IM ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTGEVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = IM
      IF( N.EQ.0 )
     $   RETURN
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      BIG = ONE / SAFMIN
      CALL DLABAD( SAFMIN, BIG )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SMALL = SAFMIN*N / ULP
      BIG = ONE / SMALL
      BIGNUM = ONE / ( SAFMIN*N )
*
*     Compute the 1-norm of each column of the strictly upper triangular
*     part (i.e., excluding all entries belonging to the diagonal
*     blocks) of A and B to check for possible overflow in the
*     triangular solver.
*
      ANORM = ABS( A( 1, 1 ) )
      IF( N.GT.1 )
     $   ANORM = ANORM + ABS( A( 2, 1 ) )
      BNORM = ABS( B( 1, 1 ) )
      WORK( 1, 1 ) = ZERO
      WORK( 1, 2 ) = ZERO
*
      DO 50 J = 2, N
         TEMP = ZERO
         TEMP2 = ZERO
         IF( A( J, J-1 ).EQ.ZERO ) THEN
            IEND = J - 1
         ELSE
            IEND = J - 2
         END IF
         DO 30 I = 1, IEND
            TEMP = TEMP + ABS( A( I, J ) )
            TEMP2 = TEMP2 + ABS( B( I, J ) )
   30    CONTINUE
         WORK( J, 1 ) = TEMP
         WORK( J, 2 ) = TEMP2
         DO 40 I = IEND + 1, MIN( J+1, N )
            TEMP = TEMP + ABS( A( I, J ) )
            TEMP2 = TEMP2 + ABS( B( I, J ) )
   40    CONTINUE
         ANORM = MAX( ANORM, TEMP )
         BNORM = MAX( BNORM, TEMP2 )
   50 CONTINUE
*
      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )
*
*     Left eigenvectors
*
      IF( COMPL ) THEN
         IEIG = 0
*
*        Main loop over eigenvalues
*
         ILCPLX = .FALSE.
         DO 220 JE = 1, N
*
*           Skip this iteration if (a) JOB='S' and SELECT=.FALSE., or
*           (b) this would be the second of a complex pair.
*           Check for complex eigenvalue, so as to be sure of which
*           entry(-ies) of SELECT to look at.
*
            IF( ILCPLX ) THEN
               ILCPLX = .FALSE.
               GO TO 220
            END IF
            NW = 1
            IF( JE.LT.N ) THEN
               IF( A( JE+1, JE ).NE.ZERO ) THEN
                  ILCPLX = .TRUE.
                  NW = 2
               END IF
            END IF
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE IF( ILCPLX ) THEN
               ILCOMP = SELECT( JE ) .OR. SELECT( JE+1 )
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( .NOT.ILCOMP )
     $         GO TO 220
*
*           Decide if (a) singular pencil, (b) real eigenvalue, or
*           (c) complex eigenvalue.
*
            IF( .NOT.ILCPLX ) THEN
               IF( ABS( A( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( B( JE, JE ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- zero eigenvector
*
                  IEIG = IEIG + 1
                  DO 60 JR = 1, N
                     VL( JR, IEIG ) = ZERO
   60             CONTINUE
                  GO TO 220
               END IF
            END IF
*
*           Clear vector
*
            DO 70 JR = 1, NW*N
               WORK( JR, 3 ) = ZERO
   70       CONTINUE
*                                                 T
*           Compute coefficients in  ( a A - b B )  y = 0
*              a  is  ACOEF
*              b  is  BCOEFR + i*BCOEFI
*
            IF( .NOT.ILCPLX ) THEN
*
*              Real eigenvalue
*
               TEMP = ONE / MAX( ABS( A( JE, JE ) )*ASCALE,
     $                ABS( B( JE, JE ) )*BSCALE, SAFMIN )
               SALFAR = ( TEMP*A( JE, JE ) )*ASCALE
               SBETA = ( TEMP*B( JE, JE ) )*BSCALE
               ACOEF = SBETA*ASCALE
               BCOEFR = SALFAR*BSCALE
               BCOEFI = ZERO
*
*              Scale to avoid underflow
*
               SCALE = ONE
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEF ).LT.SMALL
               LSB = ABS( SALFAR ).GE.SAFMIN .AND. ABS( BCOEFR ).LT.
     $               SMALL
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEF ),
     $                    ABS( BCOEFR ) ) ) )
                  IF( LSA ) THEN
                     ACOEF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEF = SCALE*ACOEF
                  END IF
                  IF( LSB ) THEN
                     BCOEFR = BSCALE*( SCALE*SALFAR )
                  ELSE
                     BCOEFR = SCALE*BCOEFR
                  END IF
               END IF
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR )
*
*              First component is 1
*
               WORK( JE, 3 ) = ONE
               XMAX = ONE
            ELSE
*
*              Complex eigenvalue
*
               CALL DLAG2( A( JE, JE ), LDA, B( JE, JE ), LDB, SAFMIN,
     $                     ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI )
               BCOEFI = -BCOEFI
               IF( BCOEFI.EQ.ZERO ) THEN
                  INFO = JE
                  RETURN
               END IF
*
*              Scale to avoid overflow
*
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               SCALE = ONE
               IF( SAFMIN*ACOEFA.GT.ASCALE )
     $            SCALE = ASCALE / ( SAFMIN*ACOEFA )
               IF( SAFMIN*BCOEFA.GT.BSCALE )
     $            SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
               IF( SCALE.LT.ONE ) THEN
                  ACOEF = SCALE*ACOEF
                  ACOEFA = ABS( ACOEF )
                  BCOEFR = SCALE*BCOEFR
                  BCOEFI = SCALE*BCOEFI
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               END IF
*
*              Compute first two components of eigenvector
*
               TEMP = ACOEF*A( JE+1, JE )
               TEMP2R = ACOEF*A( JE, JE ) - BCOEFR*B( JE, JE )
               TEMP2I = -BCOEFI*B( JE, JE )
               IF( ABS( TEMP ).GT.ABS( TEMP2R )+ABS( TEMP2I ) ) THEN
                  WORK( JE, 3 ) = ONE
                  WORK( JE, 4 ) = ZERO
                  WORK( JE+1, 3 ) = -TEMP2R / TEMP
                  WORK( JE+1, 4 ) = -TEMP2I / TEMP
               ELSE
                  WORK( JE+1, 3 ) = ONE
                  WORK( JE+1, 4 ) = ZERO
                  TEMP = ACOEF*A( JE, JE+1 )
                  WORK( JE, 3 ) = ( BCOEFR*B( JE+1, JE+1 )-ACOEF*
     $                            A( JE+1, JE+1 ) ) / TEMP
                  WORK( JE, 4 ) = BCOEFI*B( JE+1, JE+1 ) / TEMP
               END IF
               XMAX = MAX( ABS( WORK( JE, 3 ) )+ABS( WORK( JE, 4 ) ),
     $                ABS( WORK( JE+1, 3 ) )+ABS( WORK( JE+1, 4 ) ) )
            END IF
*
            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*                                           T
*           Triangular solve of  (a A - b B)  y = 0
*
*                                   T
*           (rowwise in  (a A - b B) , or columnwise in (a A - b B) )
*
            IL2BY2 = .FALSE.
*
            DO 160 J = JE + NW, N
               IF( IL2BY2 ) THEN
                  IL2BY2 = .FALSE.
                  GO TO 160
               END IF
*
               NA = 1
               BDIAG( 1 ) = B( J, J )
               IF( J.LT.N ) THEN
                  IF( A( J+1, J ).NE.ZERO ) THEN
                     IL2BY2 = .TRUE.
                     BDIAG( 2 ) = B( J+1, J+1 )
                     NA = 2
                  END IF
               END IF
*
*              Check whether scaling is necessary for dot products
*
               XSCALE = ONE / MAX( ONE, XMAX )
               TEMP = MAX( WORK( J, 1 ), WORK( J, 2 ),
     $                ACOEFA*WORK( J, 1 )+BCOEFA*WORK( J, 2 ) )
               IF( IL2BY2 )
     $            TEMP = MAX( TEMP, WORK( J+1, 1 ), WORK( J+1, 2 ),
     $                   ACOEFA*WORK( J+1, 1 )+BCOEFA*WORK( J+1, 2 ) )
               IF( TEMP.GT.BIGNUM*XSCALE ) THEN
                  DO 90 JW = 0, NW - 1
                     DO 80 JR = JE, J - 1
                        WORK( JR, JW+3 ) = XSCALE*WORK( JR, JW+3 )
   80                CONTINUE
   90             CONTINUE
                  XMAX = XMAX*XSCALE
               END IF
*
*              Compute dot products
*
*                    j-1
*              SUM = sum  conjg( a*A(k,j) - b*B(k,j) )*x(k)
*                    k=je
*
*              To reduce the op count, this is done as
*
*              _        j-1                  _        j-1
*              a*conjg( sum  A(k,j)*x(k) ) - b*conjg( sum  B(k,j)*x(k) )
*                       k=je                          k=je
*
*              which may cause underflow problems if A or B are close
*              to underflow.  (E.g., less than SMALL.)
*
*
*              A series of compiler directives to defeat vectorization
*              for the next loop
*
*$PL$ CMCHAR=' '
CDIR$          NEXTSCALAR
C$DIR          SCALAR
CDIR$          NEXT SCALAR
CVD$L          NOVECTOR
CDEC$          NOVECTOR
CVD$           NOVECTOR
*VDIR          NOVECTOR
*VOCL          LOOP,SCALAR
CIBM           PREFER SCALAR
*$PL$ CMCHAR='*'
*
               DO 120 JW = 1, NW
*
*$PL$ CMCHAR=' '
CDIR$             NEXTSCALAR
C$DIR             SCALAR
CDIR$             NEXT SCALAR
CVD$L             NOVECTOR
CDEC$             NOVECTOR
CVD$              NOVECTOR
*VDIR             NOVECTOR
*VOCL             LOOP,SCALAR
CIBM              PREFER SCALAR
*$PL$ CMCHAR='*'
*
                  DO 110 JA = 1, NA
                     SUMA( JA, JW ) = ZERO
                     SUMB( JA, JW ) = ZERO
*
                     DO 100 JR = JE, J - 1
                        SUMA( JA, JW ) = SUMA( JA, JW ) +
     $                                   A( JR, J+JA-1 )*
     $                                   WORK( JR, 2+JW )
                        SUMB( JA, JW ) = SUMB( JA, JW ) +
     $                                   B( JR, J+JA-1 )*
     $                                   WORK( JR, 2+JW )
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
*
*$PL$ CMCHAR=' '
CDIR$          NEXTSCALAR
C$DIR          SCALAR
CDIR$          NEXT SCALAR
CVD$L          NOVECTOR
CDEC$          NOVECTOR
CVD$           NOVECTOR
*VDIR          NOVECTOR
*VOCL          LOOP,SCALAR
CIBM           PREFER SCALAR
*$PL$ CMCHAR='*'
*
               DO 130 JA = 1, NA
                  IF( ILCPLX ) THEN
                     SUM( JA, 1 ) = -ACOEF*SUMA( JA, 1 ) +
     $                              BCOEFR*SUMB( JA, 1 ) -
     $                              BCOEFI*SUMB( JA, 2 )
                     SUM( JA, 2 ) = -ACOEF*SUMA( JA, 2 ) +
     $                              BCOEFR*SUMB( JA, 2 ) +
     $                              BCOEFI*SUMB( JA, 1 )
                  ELSE
                     SUM( JA, 1 ) = -ACOEF*SUMA( JA, 1 ) +
     $                              BCOEFR*SUMB( JA, 1 )
                  END IF
  130          CONTINUE
*
*                                  T
*              Solve  ( a A - b B )  y = SUM(,)
*              with scaling and perturbation of the denominator
*
               CALL DLALN2( .TRUE., NA, NW, DMIN, ACOEF, A( J, J ), LDA,
     $                      BDIAG( 1 ), BDIAG( 2 ), SUM, 2, BCOEFR,
     $                      BCOEFI, WORK( J, 3 ), N, SCALE, TEMP,
     $                      IINFO )
               IF( SCALE.LT.ONE ) THEN
                  DO 150 JW = 0, NW - 1
                     DO 140 JR = JE, J - 1
                        WORK( JR, 3+JW ) = SCALE*WORK( JR, 3+JW )
  140                CONTINUE
  150             CONTINUE
                  XMAX = SCALE*XMAX
               END IF
               XMAX = MAX( XMAX, TEMP )
  160       CONTINUE
*
*           Copy eigenvector to VL, back transforming if
*           JOB='B' or 'T'.
*
            IEIG = IEIG + 1
            IF( ILBACK ) THEN
               DO 170 JW = 0, NW - 1
                  CALL DGEMV( 'N', N, N+1-JE, ONE, VL( 1, JE ), LDVL,
     $                        WORK( JE, 3+JW ), 1, ZERO,
     $                        WORK( 1, 5+JW ), 1 )
  170          CONTINUE
               CALL DLACPY( ' ', N, NW, WORK( 1, 5 ), N, VL( 1, JE ),
     $                      LDVL )
               IBEG = 1
            ELSE
               CALL DLACPY( ' ', N, NW, WORK( 1, 3 ), N, VL( 1, IEIG ),
     $                      LDVL )
               IBEG = JE
            END IF
*
*           Scale eigenvector
*
            XMAX = ZERO
            IF( ILCPLX ) THEN
               DO 180 J = IBEG, N
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) )+
     $                   ABS( VL( J, IEIG+1 ) ) )
  180          CONTINUE
            ELSE
               DO 190 J = IBEG, N
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) ) )
  190          CONTINUE
            END IF
*
            IF( XMAX.GT.SAFMIN ) THEN
               XSCALE = ONE / XMAX
*
               DO 210 JW = 0, NW - 1
                  DO 200 JR = IBEG, N
                     VL( JR, IEIG+JW ) = XSCALE*VL( JR, IEIG+JW )
  200             CONTINUE
  210          CONTINUE
            END IF
            IEIG = IEIG + NW - 1
*
  220    CONTINUE
      END IF
*
*     Right eigenvectors
*
      IF( COMPR ) THEN
         IEIG = IM + 1
*
*        Main loop over eigenvalues
*
         ILCPLX = .FALSE.
         DO 500 JE = N, 1, -1
*
*           Skip this iteration if (a) JOB='S' and SELECT=.FALSE., or
*           (b) this would be the second of a complex pair.
*           Check for complex eigenvalue, so as to be sure of which
*           entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
*           or SELECT(JE-1).
*           If this is a complex pair, the 2-by-2 diagonal block
*           corresponding to the eigenvalue is in rows/columns JE-1:JE
*
            IF( ILCPLX ) THEN
               ILCPLX = .FALSE.
               GO TO 500
            END IF
            NW = 1
            IF( JE.GT.1 ) THEN
               IF( A( JE, JE-1 ).NE.ZERO ) THEN
                  ILCPLX = .TRUE.
                  NW = 2
               END IF
            END IF
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE IF( ILCPLX ) THEN
               ILCOMP = SELECT( JE ) .OR. SELECT( JE-1 )
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( .NOT.ILCOMP )
     $         GO TO 500
*
*           Decide if (a) singular pencil, (b) real eigenvalue, or
*           (c) complex eigenvalue.
*
            IF( .NOT.ILCPLX ) THEN
               IF( ABS( A( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( B( JE, JE ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- zero eigenvector
*
                  IEIG = IEIG - 1
                  DO 230 JR = 1, N
                     VR( JR, IEIG ) = ZERO
  230             CONTINUE
                  GO TO 500
               END IF
            END IF
*
*           Clear vector
*
            DO 250 JW = 0, NW - 1
               DO 240 JR = 1, N
                  WORK( JR, 3+JW ) = ZERO
  240          CONTINUE
  250       CONTINUE
*
*           Compute coefficients in  ( a A - b B ) x = 0
*              a  is  ACOEF
*              b  is  BCOEFR + i*BCOEFI
*
            IF( .NOT.ILCPLX ) THEN
*
*              Real eigenvalue
*
               TEMP = ONE / MAX( ABS( A( JE, JE ) )*ASCALE,
     $                ABS( B( JE, JE ) )*BSCALE, SAFMIN )
               SALFAR = ( TEMP*A( JE, JE ) )*ASCALE
               SBETA = ( TEMP*B( JE, JE ) )*BSCALE
               ACOEF = SBETA*ASCALE
               BCOEFR = SALFAR*BSCALE
               BCOEFI = ZERO
*
*              Scale to avoid underflow
*
               SCALE = ONE
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEF ).LT.SMALL
               LSB = ABS( SALFAR ).GE.SAFMIN .AND. ABS( BCOEFR ).LT.
     $               SMALL
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEF ),
     $                    ABS( BCOEFR ) ) ) )
                  IF( LSA ) THEN
                     ACOEF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEF = SCALE*ACOEF
                  END IF
                  IF( LSB ) THEN
                     BCOEFR = BSCALE*( SCALE*SALFAR )
                  ELSE
                     BCOEFR = SCALE*BCOEFR
                  END IF
               END IF
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR )
*
*              First component is 1
*
               WORK( JE, 3 ) = ONE
               XMAX = ONE
*
*              Compute contribution from column JE of A and B to sum
*              (See "Further Details", above.)
*
               DO 260 JR = 1, JE - 1
                  WORK( JR, 3 ) = BCOEFR*B( JR, JE ) - ACOEF*A( JR, JE )
  260          CONTINUE
            ELSE
*
*              Complex eigenvalue
*
               CALL DLAG2( A( JE-1, JE-1 ), LDA, B( JE-1, JE-1 ), LDB,
     $                     SAFMIN, ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI )
               IF( BCOEFI.EQ.ZERO ) THEN
                  INFO = JE - 1
                  RETURN
               END IF
*
*              Scale to avoid overflow
*
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               SCALE = ONE
               IF( SAFMIN*ACOEFA.GT.ASCALE )
     $            SCALE = ASCALE / ( SAFMIN*ACOEFA )
               IF( SAFMIN*BCOEFA.GT.BSCALE )
     $            SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
               IF( SCALE.LT.ONE ) THEN
                  ACOEF = SCALE*ACOEF
                  ACOEFA = ABS( ACOEF )
                  BCOEFR = SCALE*BCOEFR
                  BCOEFI = SCALE*BCOEFI
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               END IF
*
*              Compute first two components of eigenvector
*              and contribution to sums
*
               TEMP = ACOEF*A( JE, JE-1 )
               TEMP2R = ACOEF*A( JE, JE ) - BCOEFR*B( JE, JE )
               TEMP2I = -BCOEFI*B( JE, JE )
               IF( ABS( TEMP ).GE.ABS( TEMP2R )+ABS( TEMP2I ) ) THEN
                  WORK( JE, 3 ) = ONE
                  WORK( JE, 4 ) = ZERO
                  WORK( JE-1, 3 ) = -TEMP2R / TEMP
                  WORK( JE-1, 4 ) = -TEMP2I / TEMP
               ELSE
                  WORK( JE-1, 3 ) = ONE
                  WORK( JE-1, 4 ) = ZERO
                  TEMP = ACOEF*A( JE-1, JE )
                  WORK( JE, 3 ) = ( BCOEFR*B( JE-1, JE-1 )-ACOEF*
     $                            A( JE-1, JE-1 ) ) / TEMP
                  WORK( JE, 4 ) = BCOEFI*B( JE-1, JE-1 ) / TEMP
               END IF
*
               XMAX = MAX( ABS( WORK( JE, 3 ) )+ABS( WORK( JE, 4 ) ),
     $                ABS( WORK( JE-1, 3 ) )+ABS( WORK( JE-1, 4 ) ) )
*
*              Compute contribution from columns JE and JE-1
*              of A and B to the sums.
*
               CREALA = ACOEF*WORK( JE-1, 3 )
               CIMAGA = ACOEF*WORK( JE-1, 4 )
               CREALB = BCOEFR*WORK( JE-1, 3 ) - BCOEFI*WORK( JE-1, 4 )
               CIMAGB = BCOEFI*WORK( JE-1, 3 ) + BCOEFR*WORK( JE-1, 4 )
               CRE2A = ACOEF*WORK( JE, 3 )
               CIM2A = ACOEF*WORK( JE, 4 )
               CRE2B = BCOEFR*WORK( JE, 3 ) - BCOEFI*WORK( JE, 4 )
               CIM2B = BCOEFI*WORK( JE, 3 ) + BCOEFR*WORK( JE, 4 )
               DO 270 JR = 1, JE - 2
                  WORK( JR, 3 ) = -CREALA*A( JR, JE-1 ) +
     $                            CREALB*B( JR, JE-1 ) -
     $                            CRE2A*A( JR, JE ) + CRE2B*B( JR, JE )
                  WORK( JR, 4 ) = -CIMAGA*A( JR, JE-1 ) +
     $                            CIMAGB*B( JR, JE-1 ) -
     $                            CIM2A*A( JR, JE ) + CIM2B*B( JR, JE )
  270          CONTINUE
            END IF
*
            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*           Columnwise triangular solve of  (a A - b B)  x = 0
*
            IL2BY2 = .FALSE.
            DO 370 J = JE - NW, 1, -1
*
*              If a 2-by-2 block, is in position j-1:j, wait until
*              next iteration to process it (when it will be j:j+1)
*
               IF( .NOT.IL2BY2 .AND. J.GT.1 ) THEN
                  IF( A( J, J-1 ).NE.ZERO ) THEN
                     IL2BY2 = .TRUE.
                     GO TO 370
                  END IF
               END IF
               BDIAG( 1 ) = B( J, J )
               IF( IL2BY2 ) THEN
                  NA = 2
                  BDIAG( 2 ) = B( J+1, J+1 )
               ELSE
                  NA = 1
               END IF
*
*              Compute x(j) (and x(j+1), if 2-by-2 block)
*
               CALL DLALN2( .FALSE., NA, NW, DMIN, ACOEF, A( J, J ),
     $                      LDA, BDIAG( 1 ), BDIAG( 2 ), WORK( J, 3 ),
     $                      N, BCOEFR, BCOEFI, SUM, 2, SCALE, TEMP,
     $                      IINFO )
               IF( SCALE.LT.ONE ) THEN
*
                  DO 290 JW = 0, NW - 1
                     DO 280 JR = 1, JE
                        WORK( JR, 3+JW ) = SCALE*WORK( JR, 3+JW )
  280                CONTINUE
  290             CONTINUE
                  XMAX = MAX( SCALE*XMAX, TEMP )
               END IF
*
               DO 310 JW = 1, NW
                  DO 300 JA = 1, NA
                     WORK( J+JA-1, 2+JW ) = SUM( JA, JW )
  300             CONTINUE
  310          CONTINUE
*
*              w = w + x(j)*(a A(*,j) - b B(*,j) ) with scaling
*
               IF( J.GT.1 ) THEN
*
*                 Check whether scaling is necessary for sum.
*
                  XSCALE = ONE / MAX( ONE, XMAX )
                  TEMP = ACOEFA*WORK( J, 1 ) + BCOEFA*WORK( J, 2 )
                  IF( IL2BY2 )
     $               TEMP = MAX( TEMP, ACOEFA*WORK( J+1, 1 )+BCOEFA*
     $                      WORK( J+1, 2 ) )
                  IF( TEMP.GT.BIGNUM*XSCALE ) THEN
*
                     DO 330 JW = 0, NW - 1
                        DO 320 JR = 1, JE
                           WORK( JR, 3+JW ) = SCALE*WORK( JR, 3+JW )
  320                   CONTINUE
  330                CONTINUE
                     XMAX = XMAX*XSCALE
                  END IF
*
*                 Compute the contributions of the off-diagonals of
*                 column j (and j+1, if 2-by-2 block) of A and B to the
*                 sums.
*
*
                  DO 360 JA = 1, NA
                     IF( ILCPLX ) THEN
                        CREALA = ACOEF*WORK( J+JA-1, 3 )
                        CIMAGA = ACOEF*WORK( J+JA-1, 4 )
                        CREALB = BCOEFR*WORK( J+JA-1, 3 ) -
     $                           BCOEFI*WORK( J+JA-1, 4 )
                        CIMAGB = BCOEFI*WORK( J+JA-1, 3 ) +
     $                           BCOEFR*WORK( J+JA-1, 4 )
                        DO 340 JR = 1, J - 1
                           WORK( JR, 3 ) = WORK( JR, 3 ) -
     $                                     CREALA*A( JR, J+JA-1 ) +
     $                                     CREALB*B( JR, J+JA-1 )
                           WORK( JR, 4 ) = WORK( JR, 4 ) -
     $                                     CIMAGA*A( JR, J+JA-1 ) +
     $                                     CIMAGB*B( JR, J+JA-1 )
  340                   CONTINUE
                     ELSE
                        CREALA = ACOEF*WORK( J+JA-1, 3 )
                        CREALB = BCOEFR*WORK( J+JA-1, 3 )
                        DO 350 JR = 1, J - 1
                           WORK( JR, 3 ) = WORK( JR, 3 ) -
     $                                     CREALA*A( JR, J+JA-1 ) +
     $                                     CREALB*B( JR, J+JA-1 )
  350                   CONTINUE
                     END IF
  360             CONTINUE
               END IF
*
               IL2BY2 = .FALSE.
  370       CONTINUE
*
*           Copy eigenvector to VR, back transforming if
*           JOB='B' or 'T'.
*
            IEIG = IEIG - NW
            IF( ILBACK ) THEN
*
               DO 410 JW = 0, NW - 1
                  DO 380 JR = 1, N
                     WORK( JR, 5+JW ) = WORK( 1, 3+JW )*VR( JR, 1 )
  380             CONTINUE
*
*                 A series of compiler directives to defeat
*                 vectorization for the next loop
*
*
                  DO 400 JC = 2, JE
                     DO 390 JR = 1, N
                        WORK( JR, 5+JW ) = WORK( JR, 5+JW ) +
     $                                     WORK( JC, 3+JW )*VR( JR, JC )
  390                CONTINUE
  400             CONTINUE
  410          CONTINUE
*
               DO 430 JW = 0, NW - 1
                  DO 420 JR = 1, N
                     VR( JR, IEIG+JW ) = WORK( JR, 5+JW )
  420             CONTINUE
  430          CONTINUE
*
               IEND = N
            ELSE
               DO 450 JW = 0, NW - 1
                  DO 440 JR = 1, N
                     VR( JR, IEIG+JW ) = WORK( JR, 3+JW )
  440             CONTINUE
  450          CONTINUE
*
               IEND = JE
            END IF
*
*           Scale eigenvector
*
            XMAX = ZERO
            IF( ILCPLX ) THEN
               DO 460 J = 1, IEND
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) )+
     $                   ABS( VR( J, IEIG+1 ) ) )
  460          CONTINUE
            ELSE
               DO 470 J = 1, IEND
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) ) )
  470          CONTINUE
            END IF
*
            IF( XMAX.GT.SAFMIN ) THEN
               XSCALE = ONE / XMAX
               DO 490 JW = 0, NW - 1
                  DO 480 JR = 1, IEND
                     VR( JR, IEIG+JW ) = XSCALE*VR( JR, IEIG+JW )
  480             CONTINUE
  490          CONTINUE
            END IF
  500    CONTINUE
      END IF
*
      RETURN
*
*     End of DTGEVC
*
      END
