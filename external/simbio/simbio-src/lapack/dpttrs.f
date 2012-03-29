      SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DPTTRS solves a system of linear equations A * X = B with a
*  symmetric positive definite tridiagonal matrix A using the
*  factorization A = L*D*L**T or A = U**T*D*U computed by DPTTRF.
*  (The two forms are equivalent if A is real.)
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the diagonal matrix D from the
*          factorization computed by DPTTRF.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) off-diagonal elements of the unit bidiagonal factor
*          U or L from the factorization computed by DPTTRF.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Solve A * X = B using the factorization A = L*D*L',
*     overwriting each right hand side vector with its solution.
*
      DO 30 J = 1, NRHS
*
*        Solve L * x = b.
*
         DO 10 I = 2, N
            B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   10    CONTINUE
*
*        Solve D * L' * x = b.
*
         B( N, J ) = B( N, J ) / D( N )
         DO 20 I = N - 1, 1, -1
            B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   20    CONTINUE
   30 CONTINUE
*
      RETURN
*
*     End of DPTTRS
*
      END
