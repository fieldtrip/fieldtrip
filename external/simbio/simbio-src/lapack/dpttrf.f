      SUBROUTINE DPTTRF( N, D, E, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DPTTRF computes the factorization of a real symmetric positive
*  definite tridiagonal matrix A.
*
*  If the subdiagonal elements of A are supplied in the array E, the
*  factorization has the form A = L*D*L**T, where D is diagonal and L
*  is unit lower bidiagonal; if the superdiagonal elements of A are
*  supplied, it has the form A = U**T*D*U, where U is unit upper
*  bidiagonal.  (The two forms are equivalent if A is real.)
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          A.  On exit, the n diagonal elements of the diagonal matrix
*          D from the L*D*L**T factorization of A.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) off-diagonal elements of the tridiagonal
*          matrix A.
*          On exit, the (n-1) off-diagonal elements of the unit
*          bidiagonal factor L or U from the factorization of A.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite; if i < N, the factorization could
*                not be completed, while if i = N, the factorization was
*                completed, but D(N) = 0.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   DI, EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DPTTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute the L*D*L' (or U'*D*U) factorization of A.
*
      DO 10 I = 1, N - 1
*
*        Drop out of the loop if d(i) <= 0: the matrix is not positive
*        definite.
*
         DI = D( I )
         IF( DI.LE.ZERO )
     $      GO TO 20
*
*        Solve for e(i) and d(i+1).
*
         EI = E( I )
         E( I ) = EI / DI
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE
*
*     Check d(n) for positive definiteness.
*
      I = N
      IF( D( I ).GT.ZERO )
     $   GO TO 30
*
   20 CONTINUE
      INFO = I
*
   30 CONTINUE
      RETURN
*
*     End of DPTTRF
*
      END
