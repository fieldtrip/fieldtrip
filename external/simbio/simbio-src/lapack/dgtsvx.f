      SUBROUTINE DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,
     $                   DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          FACT, TRANS
      INTEGER            INFO, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), DF( * ),
     $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ),
     $                   FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DGTSVX uses the LU factorization to compute the solution to a real
*  system of linear equations A * X = B or A**T * X = B,
*  where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
*  matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*  Description
*  ===========
*
*  The following steps are performed:
*
*  1. If FACT = 'N', the LU decomposition is used to factor the matrix A
*     as A = L * U, where L is a product of permutation and unit lower
*     bidiagonal matrices and U is upper triangular with nonzeros in
*     only the main diagonal and first two superdiagonals.
*
*  2. The factored form of A is used to estimate the condition number
*     of the matrix A.  If the reciprocal of the condition number is
*     less than machine precision, steps 3 and 4 are skipped.
*
*  3. The system of equations is solved for X using the factored form
*     of A.
*
*  4. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  Arguments
*  =========
*
*  FACT    (input) CHARACTER*1
*          Specifies whether or not the factored form of A has been
*          supplied on entry.
*          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored
*                  form of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV
*                  will not be modified.
*          = 'N':  The matrix will be copied to DLF, DF, and DUF
*                  and factored.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of A.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of A.
*
*  DU      (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) superdiagonal elements of A.
*
*  DLF     (input or output) DOUBLE PRECISION array, dimension (N-1)
*          If FACT = 'F', then DLF is an input argument and on entry
*          contains the (n-1) multipliers that define the matrix L from
*          the LU factorization of A as computed by DGTTRF.
*
*          If FACT = 'N', then DLF is an output argument and on exit
*          contains the (n-1) multipliers that define the matrix L from
*          the LU factorization of A.
*
*  DF      (input or output) DOUBLE PRECISION array, dimension (N)
*          If FACT = 'F', then DF is an input argument and on entry
*          contains the n diagonal elements of the upper triangular
*          matrix U from the LU factorization of A.
*
*          If FACT = 'N', then DF is an output argument and on exit
*          contains the n diagonal elements of the upper triangular
*          matrix U from the LU factorization of A.
*
*  DUF     (input or output) DOUBLE PRECISION array, dimension (N-1)
*          If FACT = 'F', then DUF is an input argument and on entry
*          contains the (n-1) elements of the first superdiagonal of U.
*
*          If FACT = 'N', then DUF is an output argument and on exit
*          contains the (n-1) elements of the first superdiagonal of U.
*
*  DU2     (input or output) DOUBLE PRECISION array, dimension (N-2)
*          If FACT = 'F', then DU2 is an input argument and on entry
*          contains the (n-2) elements of the second superdiagonal of
*          U.
*
*          If FACT = 'N', then DU2 is an output argument and on exit
*          contains the (n-2) elements of the second superdiagonal of
*          U.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          If FACT = 'F', then IPIV is an input argument and on entry
*          contains the pivot indices from the LU factorization of A as
*          computed by DGTTRF.
*
*          If FACT = 'N', then IPIV is an output argument and on exit
*          contains the pivot indices from the LU factorization of A;
*          row i of the matrix was interchanged with row IPIV(i).
*          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates
*          a row interchange was not required.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          The N-by-NRHS right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          If INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  RCOND   (output) DOUBLE PRECISION
*          The estimate of the reciprocal condition number of the matrix
*          A.  If RCOND is less than the machine precision (in
*          particular, if RCOND = 0), the matrix is singular to working
*          precision.  This condition is indicated by a return code of
*          INFO > 0, and the solution and error bounds are not computed.
*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution, FERR(j) bounds the magnitude
*          of the largest entry in (X(j) - XTRUE) divided by the
*          magnitude of the largest entry in X(j).  The quality of the
*          error bound depends on the quality of the estimate of
*          norm(inv(A)) computed in the code; if the estimate of
*          norm(inv(A)) is accurate, the error bound is guaranteed.
*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in any
*          entry of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, and i is
*                <= N:  U(i,i) is exactly zero.  The factorization
*                       has not been completed unless i = N, but the
*                       factor U is exactly singular, so the solution
*                       and error bounds could not be computed.
*               = N+1:  RCOND is less than machine precision.  The
*                       factorization has been completed, but the
*                       matrix is singular to working precision, and
*                       the solution and error bounds have not been
*                       computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOFACT, NOTRAN
      CHARACTER          NORM
      DOUBLE PRECISION   ANORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGT
      EXTERNAL           LSAME, DLAMCH, DLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGTCON, DGTRFS, DGTTRF, DGTTRS, DLACPY,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGTSVX', -INFO )
         RETURN
      END IF
*
      IF( NOFACT ) THEN
*
*        Compute the LU factorization of A.
*
         CALL DCOPY( N, D, 1, DF, 1 )
         IF( N.GT.1 ) THEN
            CALL DCOPY( N-1, DL, 1, DLF, 1 )
            CALL DCOPY( N-1, DU, 1, DUF, 1 )
         END IF
         CALL DGTTRF( N, DLF, DF, DUF, DU2, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.NE.0 ) THEN
            IF( INFO.GT.0 )
     $         RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A.
*
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = DLANGT( NORM, N, DL, D, DU )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL DGTCON( NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK,
     $             IWORK, INFO )
*
*     Return if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
         INFO = N + 1
         RETURN
      END IF
*
*     Compute the solution vectors X.
*
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DGTTRS( TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX,
     $             INFO )
*
*     Use iterative refinement to improve the computed solutions and
*     compute error bounds and backward error estimates for them.
*
      CALL DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV,
     $             B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
      RETURN
*
*     End of DGTSVX
*
      END
