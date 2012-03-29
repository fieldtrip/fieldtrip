      SUBROUTINE DGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), WORK( * ),
     $                   X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGGLM solves a generalized linear regression model (GLM) problem:
*
*          minimize y'*y     subject to    d = A*x + B*y
*            x,y
*
*  using a generalized QR factorization of A and B, where A is an N-by-M
*  matrix, B is a given N-by-P matrix, and d is a given N vector.  It
*  is also assumed that M <= N <= M+P and
*
*             rank(A) = M    and    rank([ A B ]) = N.
*
*  Under these assumptions, the constrained equation is always
*  consistent, and there is a unique solution x and a minimal 2-norm
*  solution y.
*
*  In particular, if matrix B is square nonsingular, then the problem
*  GLM is equivalent to the following weighted linear least squares
*  problem
*               minimize || inv(B)*(b-A*x) ||
*                  x
*  where ||.|| is vector 2-norm, and inv(B) denotes the inverse of
*  matrix B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of the matrices A and B.  N >= 0.
*
*  M       (input) INTEGER
*          The number of columns of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of columns of the matrix B.  P >= 0.
*          Assume that M <= N <= M+P.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
*          On entry, the N-by-M matrix A.
*          On exit, A is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max( 1,N ).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,P)
*          On entry, the N-by-P matrix B.
*          On exit, B is destroyed.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max( 1,N ).
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          On entry, D is the left hand side of the GLM equation.
*          On exit, D is destroyed.
*
*  X       (output) DOUBLE PRECISION array, dimension (M)
*  Y       (output) DOUBLE PRECISION array, dimension (P)
*          On exit, X and Y are the solutions of the GLM problem.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension ( LWORK )
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= M+P+max(N,M,P).
*          For optimum performance, LWORK >=
*          M+P+max(N,M,P)*max(NB1,NB2), where NB1 is the optimal
*          blocksize for the QR factorization of an N-by-M matrix A.
*          NB2 is the optimal blocksize for the RQ factorization of an
*          N-by-P matrix B.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  ===================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, NP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DGGQRF, DORMQR, DORMRQ, DTRSV,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.M .OR. N.GT.M+P ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGGLM', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. MAX( M, P ).EQ.0 )
     $   RETURN
*
*     GQR factorization of matrices A and B:
*
*            Q'*A = ( R11 ) M,    Q'*B*Z' = ( T11   T12 ) M
*                   (  0  ) N-M             (  0    T22 ) N-M
*                      M                     M+P-N  N-M
*
      NP = MIN( N, P )
      CALL DGGQRF( N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ),
     $             WORK( M+NP+1 ), LWORK-M-NP, INFO )
*
*     Solve the GLM problem
*
*     Update left-hand-side vector d = Q'*d = ( d1 ) M
*                                             ( d2 ) N-M
*
      CALL DORMQR( 'Left', 'Transpose', N, 1, M, A, LDA, WORK, D,
     $             MAX( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
*
*     Solve T22*y2 = d2 for y2
*
      CALL DTRSV( 'Upper', 'No transpose', 'Non unit', N-M,
     $            B( M+1, M+P-N+1 ), LDB, D( M+1 ), 1 )
      CALL DCOPY( N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 )
*
*     Set up y1 = 0
*
      DO 10 I = 1, M + P - N
         Y( I ) = ZERO
   10 CONTINUE
*
*     Update d1 = d1 - T12*y2
*
      CALL DGEMV( 'No transpose', M, N-M, -ONE, B( 1, M+P-N+1 ), LDB,
     $            Y( M+P-N+1 ), 1, ONE, D, 1 )
*
*     Solve triangular system: R11*x = d1
*
      CALL DTRSV( 'Upper', 'No Transpose', 'Non unit', M, A, LDA, D, 1 )
*
*     Copy D to X
*
      CALL DCOPY( M, D, 1, X, 1 )
*
*     y = Z'*y
*
      CALL DORMRQ( 'Left', 'Transpose', P, 1, MIN( P, N ),
     $             B( MAX( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y,
     $             MAX( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
*
      RETURN
*
*     End of DGGGLM
*
      END
