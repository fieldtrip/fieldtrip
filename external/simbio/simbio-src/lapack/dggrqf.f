      SUBROUTINE DGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK,
     $                   LWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
*
*  Purpose
*  =======
*
*  DGGRQF computes a generalized RQ factorization of an M-by-N matrix A
*  and a P-by-N matrix B:
*
*              A = R*Q,        B = Z*T*Q,
*
*  where Q is an M-by-M orthogonal matrix, and Z is an N-by-N orthogonal
*  matrix.  R and T assume one of the forms:
*
*  if M <= N,                   or if M > N
*
*       R = ( 0  R12 ) M,            R = ( R11 ) M-N
*            N-M  M                      ( R21 ) N
*                                           N
*  where R12 or R21 is upper triangular, and
*
*  if P >= N,                   or if P < N
*
*       T = ( T11 ) N                T = ( T11  T12 ) P
*           (  0  ) P-N                     P   N-P
*              N
*
*  where T11 is an upper triangular matrix.
*
*  In particular, if B is square and nonsingular, the GRQ factorization
*  of A and B implicitly gives the RQ factorization of matrix A*inv(B):
*
*                      A*inv(B) = (R*inv(T))*Z'
*
*  where inv(B) denotes the inverse of the matrix B, and Z' denotes the
*  transpose of the matrix Z.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, if M <= N, the upper triangle of the subarray
*          A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R;
*          if M > N, the elements on and above the (M-N)-th subdiagonal
*          contain the M-by-N upper trapezoidal matrix R; the remaining
*          elements, with the array TAUA, represent the orthogonal
*          matrix Q as a product of elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  TAUA    (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, the elements on and above the diagonal of the array
*          contain the MIN(P,N)-by-N upper trapezoidal matrix T (T is
*          upper triangular if P >= N); the elements below the diagonal,
*          with the array TAUB, represent the orthogonal matrix Z as a
*          product of elementary reflectors (see Further Details).
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= MAX(1,P).
*
*  TAUB    (output) DOUBLE PRECISION array, dimension (MIN(P,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= MAX(1,N,M,P).
*          For optimum performance LWORK >=
*          MAX(1,N,M,P)*MAX(NB1,NB2,NB3), where NB1 is the optimal
*          blocksize for the RQ factorization of the M-by-N matrix A.
*          NB2 is the optimal blocksize for the QR factorization of
*          the P-by-N matrix B.  NB3 is the optimal blocksize for
*          DORMRQ.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INF0= -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(M,N).
*
*  Each H(i) has the form
*
*     H(i) = I - taub * v * v'
*
*  where taub is a real scalar, and v is a real vector with
*  v(N-k+i+1:N) = 0 and v(N-k+i) = 1; v(1:N-k+i-1) is stored on exit in
*  A(M-k+i,1:n-k+i-1), and taub in TAUA(i).
*  To form Q explicitly, use LAPACK subroutine DORGRQ.
*  To use Q to update another matrix, use LAPACK subroutine DORMRQ.
*
*  The matrix Z is represented as a product of elementary reflectors
*
*     Z = H(1) H(2) . . . H(k), where k = min(P,N).
*
*  Each H(i) has the form
*
*     H(i) = I - taua * v * v'
*
*  where taua is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:M) is stored on exit in B(i+1:P,i),
*  and taua in TAUB(i).
*  To form Z explicitly, use LAPACK subroutine DORGQR.
*  To use Z to update another matrix, use LAPACK subroutine DORMQR.
*
*  =====================================================================
*
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ),
     $                   WORK( * )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGERQF, DORMRQ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( P.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGRQF', -INFO )
         RETURN
      END IF
*
*     RQ factorization of M-by-N matrix A: A = R*Q
*
      CALL DGERQF( M, N, A, LDA, TAUA, WORK, LWORK, INFO )
*
*     Update B := B*Q'
*
      CALL DORMRQ( 'Right', 'Transpose', P, N, MIN( M, N ),
     $             A( MAX( 1, M-N+1 ), 1 ), LDA, TAUA, B, LDB, WORK,
     $             LWORK, INFO )
*
*     QR factorization of P-by-N matrix B: B = Z*T
*
      CALL DGEQRF( P, N, B, LDB, TAUB, WORK, LWORK, INFO )
*
      RETURN
*
*     End of DGGRQF
      END
