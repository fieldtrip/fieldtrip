      SUBROUTINE DLAZRO( M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAZRO initializes a 2-D array A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) DOUBLE PRECISION
*          The constant to which the diagonal elements are to be set.
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the leading m by n submatrix of A is set such that
*             A(i,j) = ALPHA,  1 <= i <= m, 1 <= j <= n, i <> j
*             A(i,i) = BETA,   1 <= i <= min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = ALPHA
   10    CONTINUE
   20 CONTINUE
*
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = BETA
   30 CONTINUE
*
      RETURN
*
*     End of DLAZRO
*
      END
