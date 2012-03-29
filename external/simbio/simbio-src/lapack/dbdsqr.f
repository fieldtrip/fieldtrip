      SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
     $                   LDU, C, LDC, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DBDSQR computes the singular value decomposition (SVD) of a real
*  N-by-N (upper or lower) bidiagonal matrix B:  B = Q * S * P' (P'
*  denotes the transpose of P), where S is a diagonal matrix with
*  non-negative diagonal elements (the singular values of B), and Q
*  and P are orthogonal matrices.
*
*  The routine computes S, and optionally computes U * Q, P' * VT,
*  or Q' * C, for given real input matrices U, VT, and C.
*
*  See "Computing  Small Singular Values of Bidiagonal Matrices With
*  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
*  LAPACK Working Note #3, for a detailed description of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  B is upper bidiagonal;
*          = 'L':  B is lower bidiagonal.
*
*  N       (input) INTEGER
*          The order of the matrix B.  N >= 0.
*
*  NCVT    (input) INTEGER
*          The number of columns of the matrix VT. NCVT >= 0.
*
*  NRU     (input) INTEGER
*          The number of rows of the matrix U. NRU >= 0.
*
*  NCC     (input) INTEGER
*          The number of columns of the matrix C. NCC >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the bidiagonal matrix B.
*          On exit, if INFO=0, the singular values of B in decreasing
*          order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) off-diagonal elements of the bidiagonal
*          matrix B.
*          On normal exit, E is destroyed.
*
*  VT      (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
*          On entry, an N-by-NCVT matrix VT.
*          On exit, VT is overwritten by P' * VT.
*          VT is not referenced if NCVT = 0.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.
*          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
*
*  U       (input/output) DOUBLE PRECISION array, dimension (LDU, N)
*          On entry, an NRU-by-N matrix U.
*          On exit, U is overwritten by U * Q.
*          U is not referenced if NRU = 0.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= max(1,NRU).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
*          On entry, an N-by-NCC matrix C.
*          On exit, C is overwritten by Q' * C.
*          C is not referenced if NCC = 0.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.
*          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                      (MAX( 1, 4*N-4 ))
*          WORK is not referenced if NCVT = NRU = NCC = 0.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  If INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm did not converge; D and E contain the
*                elements of a bidiagonal matrix which is orthogonally
*                similar to the input matrix B;  if INFO = i, i
*                elements of E have not converged to zero.
*
*  Internal Parameters
*  ===================
*
*  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
*          TOLMUL controls the convergence criterion of the QR loop.
*          If it is positive, TOLMUL*EPS is the desired relative
*             precision in the computed singular values.
*          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
*             desired absolute accuracy in the computed singular
*             values (corresponds to relative accuracy
*             abs(TOLMUL*EPS) in the largest singular value.
*          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
*             between 10 (for fast convergence) and .1/EPS
*             (for there to be some accuracy in the results).
*          Default is to lose at either one eighth or 2 of the
*             available decimal digits in each computed singular value
*             (whichever is smaller).
*
*  MAXITR  INTEGER, default = 6
*          MAXITR controls the maximum number of passes of the
*          algorithm through its inner loop. The algorithms stops
*          (and so fails to converge) if the number of passes
*          through the inner loop exceeds MAXITR*N**2.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
      DOUBLE PRECISION   MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ROTATE
      INTEGER            I, IDIR, IROT, ISUB, ITER, IUPLO, J, JOB, LL,
     $                   LLL, M, MAXIT, NM1, NM12, NM13, OLDLL, OLDM
      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, GAP,
     $                   GMAX, H, MU, OLDCS, OLDSN, R, SHIFT, SIGMN,
     $                   SIGMX, SINL, SINR, SLL, SMAX, SMIN, SMINL,
     $                   SMINLO, SMINOA, SN, THRESH, TOL, TOLMUL, UNFL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLAS2, DLASR, DLASV2, DROT, DSCAL,
     $                   DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IUPLO = 0
      IF( LSAME( UPLO, 'U' ) )
     $   IUPLO = 1
      IF( LSAME( UPLO, 'L' ) )
     $   IUPLO = 2
      IF( IUPLO.EQ.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NCVT.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NCC.LT.0 ) THEN
         INFO = -5
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR.
     $         ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
         INFO = -11
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR.
     $         ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSQR', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 )
     $   GO TO 190
*
*     ROTATE is true if any singular vectors desired, false otherwise
*
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
*
*     Get machine constants
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
*
*     If matrix lower bidiagonal, rotate to be upper bidiagonal
*     by applying Givens rotations on the left
*
      IF( IUPLO.EQ.2 ) THEN
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            IF( ROTATE ) THEN
               WORK( I ) = CS
               WORK( NM1+I ) = SN
            END IF
   10    CONTINUE
*
*        Update singular vectors if desired
*
         IF( NRU.GT.0 )
     $      CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U,
     $                  LDU )
         IF( NCC.GT.0 )
     $      CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C,
     $                  LDC )
      END IF
*
*     Compute approximate maximum, minimum singular values
*
      SMAX = ABS( D( N ) )
      DO 20 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( D( I ) ), ABS( E( I ) ) )
   20 CONTINUE
      SMINL = ZERO
      IF( TOL.GE.ZERO ) THEN
         SMINOA = ABS( D( 1 ) )
         IF( SMINOA.EQ.ZERO )
     $      GO TO 40
         MU = SMINOA
         DO 30 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            IF( SMINOA.EQ.ZERO )
     $         GO TO 40
   30    CONTINUE
   40    CONTINUE
         SMINOA = SMINOA / SQRT( DBLE( N ) )
      END IF
*
*     Prepare for main iteration loop for the singular values
*
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
      IF( NCC.EQ.0 .AND. NRU.EQ.0 .AND. NCVT.EQ.0 ) THEN
*
*        No singular vectors desired
*
         JOB = 0
      ELSE
*
*        Singular vectors desired
*
         JOB = 1
      END IF
      IF( TOL.GE.ZERO ) THEN
*
*        Relative accuracy desired
*
         THRESH = MAX( TOL*SMINOA, MAXIT*UNFL )
      ELSE
*
*        Absolute accuracy desired
*
         THRESH = MAX( ABS( TOL )*SMAX, MAXIT*UNFL )
      END IF
*
*     M points to last entry of unconverged part of matrix
*
      M = N
*
*     Begin main iteration loop
*
   50 CONTINUE
*
*     Check for convergence or exceeding iteration count
*
      IF( M.LE.1 )
     $   GO TO 190
      IF( ITER.GT.MAXIT )
     $   GO TO 230
*
*     Find diagonal block of matrix to work on
*
      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH )
     $   D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO 60 LLL = 1, M
         LL = M - LLL
         IF( LL.EQ.0 )
     $      GO TO 80
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH )
     $      D( LL ) = ZERO
         IF( ABSE.LE.THRESH )
     $      GO TO 70
         SMIN = MIN( SMIN, ABSS )
         SMAX = MAX( SMAX, ABSS, ABSE )
   60 CONTINUE
   70 CONTINUE
      E( LL ) = ZERO
*
*     Matrix splits since E(LL) = 0
*
      IF( LL.EQ.M-1 ) THEN
*
*        Convergence of bottom singular value, return to top of loop
*
         M = M - 1
         GO TO 50
      END IF
   80 CONTINUE
      LL = LL + 1
*
*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
*
      IF( LL.EQ.M-1 ) THEN
*
*        2 by 2 block, handle separately
*
         CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR,
     $                COSR, SINL, COSL )
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN
*
*        Compute singular vectors, if desired
*
         IF( NCVT.GT.0 )
     $      CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR,
     $                 SINR )
         IF( NRU.GT.0 )
     $      CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
         IF( NCC.GT.0 )
     $      CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL,
     $                 SINL )
         M = M - 2
         GO TO 50
      END IF
*
*     If working on new submatrix, choose shift direction
*     (from larger end diagonal entry towards smaller)
*
      IF( LL.GT.OLDM .OR. M.LT.OLDLL ) THEN
         IF( ABS( D( LL ) ).GE.ABS( D( M ) ) ) THEN
*
*           Chase bulge from top (big end) to bottom (small end)
*
            IDIR = 1
         ELSE
*
*           Chase bulge from bottom (big end) to top (small end)
*
            IDIR = 2
         END IF
      END IF
*
*     Apply convergence tests
*
      IF( IDIR.EQ.1 ) THEN
*
*        Run convergence test in forward direction
*        First apply standard test to bottom of matrix
*
         IF( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) THEN
            E( M-1 ) = ZERO
            GO TO 50
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           If relative accuracy desired,
*           apply convergence criterion forward
*
            MU = ABS( D( LL ) )
            SMINL = MU
            DO 90 LLL = LL, M - 1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 50
               END IF
               SMINLO = SMINL
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
   90       CONTINUE
*
*           If singular values only wanted, apply gap test to bottom
*           end of matrix
*
            IF( JOB.EQ.0 ) THEN
               GAP = SMINLO / SQRT( DBLE( M-LL ) ) - ABS( D( M ) )
               IF( GAP.GT.ZERO ) THEN
                  ABSS = ABS( D( M ) )
                  ABSE = ABS( E( M-1 ) )
                  GMAX = MAX( GAP, ABSS, ABSE )
                  IF( ( ABSE / GMAX )**2.LE.TOL*( GAP / GMAX )*
     $                ( ABSS / GMAX ) ) THEN
                     E( M-1 ) = ZERO
                     GO TO 50
                  END IF
               END IF
            END IF
         END IF
      ELSE
*
*        Run convergence test in backward direction
*        First apply standard test to top of matrix
*
         IF( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) THEN
            E( LL ) = ZERO
            GO TO 50
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           If relative accuracy desired,
*           apply convergence criterion backward
*
            MU = ABS( D( M ) )
            SMINL = MU
            DO 100 LLL = M - 1, LL, -1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 50
               END IF
               SMINLO = SMINL
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  100       CONTINUE
*
*           If singular values only wanted, apply gap test to top
*           end of matrix
*
            IF( JOB.EQ.0 ) THEN
               GAP = SMINLO / SQRT( DBLE( M-LL ) ) - ABS( D( LL ) )
               IF( GAP.GT.ZERO ) THEN
                  ABSS = ABS( D( LL ) )
                  ABSE = ABS( E( LL ) )
                  GMAX = MAX( GAP, ABSS, ABSE )
                  IF( ( ABSE / GMAX )**2.LE.TOL*( GAP / GMAX )*
     $                ( ABSS / GMAX ) ) THEN
                     E( LL ) = ZERO
                     GO TO 50
                  END IF
               END IF
            END IF
         END IF
      END IF
      OLDLL = LL
      OLDM = M
*
*     Compute shift.  First, test if shifting would ruin relative
*     accuracy, and if so set the shift to zero.
*
      IF( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE.
     $    MAX( EPS, HNDRTH*TOL ) ) THEN
*
*        Use a zero shift to avoid loss of relative accuracy
*
         SHIFT = ZERO
      ELSE
*
*        Compute the shift from 2-by-2 block at end of matrix
*
         IF( IDIR.EQ.1 ) THEN
            SLL = ABS( D( LL ) )
            CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
         ELSE
            SLL = ABS( D( M ) )
            CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
         END IF
*
*        Test if shift negligible, and if so set to zero
*
         IF( SLL.GT.ZERO ) THEN
            IF( ( SHIFT / SLL )**2.LT.EPS )
     $         SHIFT = ZERO
         END IF
      END IF
*
*     Increment iteration count
*
      ITER = ITER + M - LL
*
*     If SHIFT = 0, do simplified QR iteration
*
      IF( SHIFT.EQ.ZERO ) THEN
         IF( IDIR.EQ.1 ) THEN
*
*           Chase bulge from top to bottom
*
            CS = ONE
            OLDCS = ONE
*
*           Save cosines and sines if singular vectors desired
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( D( LL )*CS, E( LL ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( LL+1 )*SN, OLDCS, OLDSN,
     $                      D( LL ) )
               WORK( 1 ) = CS
               WORK( 1+NM1 ) = SN
               WORK( 1+NM12 ) = OLDCS
               WORK( 1+NM13 ) = OLDSN
               IROT = 1
               DO 110 I = LL + 1, M - 1
                  CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
                  E( I-1 ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
                  IROT = IROT + 1
                  WORK( IROT ) = CS
                  WORK( IROT+NM1 ) = SN
                  WORK( IROT+NM12 ) = OLDCS
                  WORK( IROT+NM13 ) = OLDSN
  110          CONTINUE
               H = D( M )*CS
               D( M ) = H*OLDCS
               E( M-1 ) = H*OLDSN
*
*              Update singular vectors
*
               IF( NCVT.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),
     $                        WORK( N ), VT( LL, 1 ), LDVT )
               IF( NRU.GT.0 )
     $            CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        U( 1, LL ), LDU )
               IF( NCC.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        C( LL, 1 ), LDC )
*
            ELSE
*
               CALL DLARTG( D( LL )*CS, E( LL ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( LL+1 )*SN, OLDCS, OLDSN,
     $                      D( LL ) )
               DO 120 I = LL + 1, M - 1
                  CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
                  E( I-1 ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
  120          CONTINUE
               H = D( M )*CS
               D( M ) = H*OLDCS
               E( M-1 ) = H*OLDSN
*
            END IF
*
*           Test convergence
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           Chase bulge from bottom to top
*
            CS = ONE
            OLDCS = ONE
*
*           Save cosines and sines if singular vectors desired
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( D( M )*CS, E( M-1 ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( M-1 )*SN, OLDCS, OLDSN, D( M ) )
               WORK( M-LL ) = CS
               WORK( M-LL+NM1 ) = -SN
               WORK( M-LL+NM12 ) = OLDCS
               WORK( M-LL+NM13 ) = -OLDSN
               IROT = M - LL
               DO 130 I = M - 1, LL + 1, -1
                  CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
                  E( I ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
                  IROT = IROT - 1
                  WORK( IROT ) = CS
                  WORK( IROT+NM1 ) = -SN
                  WORK( IROT+NM12 ) = OLDCS
                  WORK( IROT+NM13 ) = -OLDSN
  130          CONTINUE
               H = D( LL )*CS
               D( LL ) = H*OLDCS
               E( LL ) = H*OLDSN
*
*              Update singular vectors
*
               IF( NCVT.GT.0 )
     $            CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        VT( LL, 1 ), LDVT )
               IF( NRU.GT.0 )
     $            CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),
     $                        WORK( N ), U( 1, LL ), LDU )
               IF( NCC.GT.0 )
     $            CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),
     $                        WORK( N ), C( LL, 1 ), LDC )
*
            ELSE
*
               CALL DLARTG( D( M )*CS, E( M-1 ), CS, SN, R )
               CALL DLARTG( OLDCS*R, D( M-1 )*SN, OLDCS, OLDSN, D( M ) )
               DO 140 I = M - 1, LL + 1, -1
                  CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
                  E( I ) = OLDSN*R
                  CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN,
     $                         D( I ) )
  140          CONTINUE
               H = D( LL )*CS
               D( LL ) = H*OLDCS
               E( LL ) = H*OLDSN
*
            END IF
*
*           Test convergence
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
         END IF
      ELSE
*
*        Use nonzero shift
*
         IF( IDIR.EQ.1 ) THEN
*
*           Chase bulge from top to bottom
*
            F = ( ABS( D( LL ) )-SHIFT )*
     $          ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
*
*           Save cosines and sines if singular vectors desired
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( LL ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL )
               G = SINR*D( LL+1 )
               D( LL+1 ) = COSR*D( LL+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL ) = R
               F = COSL*E( LL ) + SINL*D( LL+1 )
               D( LL+1 ) = COSL*D( LL+1 ) - SINL*E( LL )
               G = SINL*E( LL+1 )
               E( LL+1 ) = COSL*E( LL+1 )
               WORK( 1 ) = COSR
               WORK( 1+NM1 ) = SINR
               WORK( 1+NM12 ) = COSL
               WORK( 1+NM13 ) = SINL
               IROT = 1
               DO 150 I = LL + 1, M - 2
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I-1 ) = R
                  F = COSR*D( I ) + SINR*E( I )
                  E( I ) = COSR*E( I ) - SINR*D( I )
                  G = SINR*D( I+1 )
                  D( I+1 ) = COSR*D( I+1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I ) + SINL*D( I+1 )
                  D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
                  IROT = IROT + 1
                  WORK( IROT ) = COSR
                  WORK( IROT+NM1 ) = SINR
                  WORK( IROT+NM12 ) = COSL
                  WORK( IROT+NM13 ) = SINL
  150          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( M-2 ) = R
               F = COSR*D( M-1 ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M-1 )
               G = SINR*D( M )
               D( M ) = COSR*D( M )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M-1 ) = R
               F = COSL*E( M-1 ) + SINL*D( M )
               D( M ) = COSL*D( M ) - SINL*E( M-1 )
               IROT = IROT + 1
               WORK( IROT ) = COSR
               WORK( IROT+NM1 ) = SINR
               WORK( IROT+NM12 ) = COSL
               WORK( IROT+NM13 ) = SINL
               E( M-1 ) = F
*
*              Update singular vectors
*
               IF( NCVT.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),
     $                        WORK( N ), VT( LL, 1 ), LDVT )
               IF( NRU.GT.0 )
     $            CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        U( 1, LL ), LDU )
               IF( NCC.GT.0 )
     $            CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC,
     $                        WORK( NM12+1 ), WORK( NM13+1 ),
     $                        C( LL, 1 ), LDC )
*
            ELSE
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( LL ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL )
               G = SINR*D( LL+1 )
               D( LL+1 ) = COSR*D( LL+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL ) = R
               F = COSL*E( LL ) + SINL*D( LL+1 )
               D( LL+1 ) = COSL*D( LL+1 ) - SINL*E( LL )
               G = SINL*E( LL+1 )
               E( LL+1 ) = COSL*E( LL+1 )
               DO 160 I = LL + 1, M - 2
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I-1 ) = R
                  F = COSR*D( I ) + SINR*E( I )
                  E( I ) = COSR*E( I ) - SINR*D( I )
                  G = SINR*D( I+1 )
                  D( I+1 ) = COSR*D( I+1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I ) + SINL*D( I+1 )
                  D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
  160          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( M-2 ) = R
               F = COSR*D( M-1 ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M-1 )
               G = SINR*D( M )
               D( M ) = COSR*D( M )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M-1 ) = R
               F = COSL*E( M-1 ) + SINL*D( M )
               D( M ) = COSL*D( M ) - SINL*E( M-1 )
               E( M-1 ) = F
*
            END IF
*
*           Test convergence
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           Chase bulge from bottom to top
*
            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT /
     $          D( M ) )
            G = E( M-1 )
*
*           Save cosines and sines if singular vectors desired
*
            IF( ROTATE ) THEN
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( M ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M )
               G = SINR*D( M-1 )
               D( M-1 ) = COSR*D( M-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M ) = R
               F = COSL*E( M-1 ) + SINL*D( M-1 )
               D( M-1 ) = COSL*D( M-1 ) - SINL*E( M-1 )
               G = SINL*E( M-2 )
               E( M-2 ) = COSL*E( M-2 )
               WORK( M-LL ) = COSR
               WORK( M-LL+NM1 ) = -SINR
               WORK( M-LL+NM12 ) = COSL
               WORK( M-LL+NM13 ) = -SINL
               IROT = M - LL
               DO 170 I = M - 1, LL + 2, -1
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I ) = R
                  F = COSR*D( I ) + SINR*E( I-1 )
                  E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
                  G = SINR*D( I-1 )
                  D( I-1 ) = COSR*D( I-1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I-1 ) + SINL*D( I-1 )
                  D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
                  IROT = IROT - 1
                  WORK( IROT ) = COSR
                  WORK( IROT+NM1 ) = -SINR
                  WORK( IROT+NM12 ) = COSL
                  WORK( IROT+NM13 ) = -SINL
  170          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( LL+1 ) = R
               F = COSR*D( LL+1 ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL+1 )
               G = SINR*D( LL )
               D( LL ) = COSR*D( LL )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL+1 ) = R
               F = COSL*E( LL ) + SINL*D( LL )
               D( LL ) = COSL*D( LL ) - SINL*E( LL )
               IROT = IROT - 1
               WORK( IROT ) = COSR
               WORK( IROT+NM1 ) = -SINR
               WORK( IROT+NM12 ) = COSL
               WORK( IROT+NM13 ) = -SINL
               E( LL ) = F
*
            ELSE
*
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( M ) + SINR*E( M-1 )
               E( M-1 ) = COSR*E( M-1 ) - SINR*D( M )
               G = SINR*D( M-1 )
               D( M-1 ) = COSR*D( M-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( M ) = R
               F = COSL*E( M-1 ) + SINL*D( M-1 )
               D( M-1 ) = COSL*D( M-1 ) - SINL*E( M-1 )
               G = SINL*E( M-2 )
               E( M-2 ) = COSL*E( M-2 )
               DO 180 I = M - 1, LL + 2, -1
                  CALL DLARTG( F, G, COSR, SINR, R )
                  E( I ) = R
                  F = COSR*D( I ) + SINR*E( I-1 )
                  E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
                  G = SINR*D( I-1 )
                  D( I-1 ) = COSR*D( I-1 )
                  CALL DLARTG( F, G, COSL, SINL, R )
                  D( I ) = R
                  F = COSL*E( I-1 ) + SINL*D( I-1 )
                  D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
  180          CONTINUE
               CALL DLARTG( F, G, COSR, SINR, R )
               E( LL+1 ) = R
               F = COSR*D( LL+1 ) + SINR*E( LL )
               E( LL ) = COSR*E( LL ) - SINR*D( LL+1 )
               G = SINR*D( LL )
               D( LL ) = COSR*D( LL )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( LL+1 ) = R
               F = COSL*E( LL ) + SINL*D( LL )
               D( LL ) = COSL*D( LL ) - SINL*E( LL )
               E( LL ) = F
*
            END IF
*
*           Test convergence
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
*
*           Update singular vectors if desired
*
            IF( NCVT.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),
     $                     WORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),
     $                     WORK( N ), C( LL, 1 ), LDC )
         END IF
      END IF
*
*     QR iteration finished, go back and check convergence
*
      GO TO 50
*
*     All singular values converged, so make them positive
*
  190 CONTINUE
      DO 200 I = 1, N
         IF( D( I ).LT.ZERO ) THEN
            D( I ) = -D( I )
*
*           Change sign of singular vectors, if desired
*
            IF( NCVT.GT.0 )
     $         CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         END IF
  200 CONTINUE
*
*     Sort the singular values into decreasing order (insertion sort on
*     singular values, but only one transposition per singular vector)
*
      DO 220 I = 1, N - 1
*
*        Scan for smallest D(I)
*
         ISUB = 1
         SMIN = D( 1 )
         DO 210 J = 2, N + 1 - I
            IF( D( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = D( J )
            END IF
  210    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
*
*           Swap singular values and vectors
*
            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            IF( NCVT.GT.0 )
     $         CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ),
     $                     LDVT )
            IF( NRU.GT.0 )
     $         CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
            IF( NCC.GT.0 )
     $         CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         END IF
  220 CONTINUE
      GO TO 250
*
*     Maximum number of iterations exceeded, failure to converge
*
  230 CONTINUE
      INFO = 0
      DO 240 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  240 CONTINUE
  250 CONTINUE
      RETURN
*
*     End of DBDSQR
*
      END
