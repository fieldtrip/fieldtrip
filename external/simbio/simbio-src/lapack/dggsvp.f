      SUBROUTINE DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,
     $                   TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,
     $                   IWORK, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P
      DOUBLE PRECISION   TOLA, TOLB
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGSVP computes orthogonal matrices U, V and Q such that
*
*   U'*A*Q = ( 0    A12  A13 ) K    , V'*B*Q = ( 0     0   B13 ) L
*            ( 0     0   A23 ) L               ( 0     0    0  ) P-L
*            ( 0     0    0  ) M-K-L            N-K-L  K    L
*             N-K-L  K    L
*
*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
*  upper triangular. If  M-K-L > 0, A23 is upper triangular, otherwise
*  A23 is upper trapezoidal.  K+L = the effective rank of (M+P)-by-N
*  matrix (A',B')'.  Z' denotes the transpose of Z.
*
*  This decomposition is the preprocessing step for computing the
*  Generalized Singular Value Decomposition (GSVD), see subroutine
*  DGGSVD.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  Orthogonal matrix U is computed;
*          = 'N':  U is not computed.
*
*  JOBV    (input) CHARACTER*1
*          = 'V':  Orthogonal matrix V is computed;
*          = 'N':  V is not computed.
*
*  JOBQ    (input) CHARACTER*1
*          = 'Q':  Orthogonal matrix Q is computed;
*          = 'N':  Q is not computed.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A contains the triangular (or trapezoidal) matrix
*          described in the Purpose section.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B contains the triangular matrix described in
*          the Purpose section.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  TOLA    (input) DOUBLE PRECISION
*  TOLB    (input) DOUBLE PRECISION
*          TOLA and TOLB are the thresholds to determine the effective
*          rank of matrix B and a subblock of A. Generally, they are set
*          to TOLA = MAX(M,N)*norm(A)*MAZHEPS,
*             TOLB = MAX(P,N)*norm(B)*MAZHEPS.
*          The size of TOLA and TOLB may affect the size of backward
*          errors of the decomposition.
*
*  K       (output) INTEGER
*  L       (output) INTEGER
*          On exit, K and L specify the dimension of the subblocks
*          described in Purpose.
*          K + L = effective numerical rank of (A',B')'.
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,M)
*          If JOBU = 'U', U contains the orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M).
*
*  V       (output) DOUBLE PRECISION array, dimension (LDV,M)
*          If JOBV = 'V', V contains the orthogonal matrix V.
*          If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P).
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
*          If JOBQ = 'Q', Q contains the orthogonal matrix Q.
*          If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N).
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  TAU     (workspace) DOUBLE PRECISION array, dimension (N)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(3*N,M,P))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*
*  Further Details
*  ===============
*
*  The subroutine uses LAPACK subroutine DGEQPF for the QR factorization
*  with column pivoting to detect the effective numerical rank of the
*  a matrix. It may be replaced by a better rank determination strategy.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FORWRD, WANTQ, WANTU, WANTV
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQPF, DGEQR2, DGERQ2, DLACPY, DLAPMT, DLAZRO,
     $                   DORG2R, DORM2R, DORMR2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      FORWRD = .TRUE.
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGSVP', -INFO )
         RETURN
      END IF
*
*     QR with column pivoting of B: B*P = V*( S11 S12 )
*                                           (  0   0  )
*
      DO 10 I = 1, N
         IWORK( I ) = 0
   10 CONTINUE
      CALL DGEQPF( P, N, B, LDB, IWORK, TAU, WORK, INFO )
*
*     Update A := A*P
*
      CALL DLAPMT( FORWRD, M, N, A, LDA, IWORK )
*
*     Determine the effective rank of matrix B.
*
      L = 0
      DO 20 I = 1, MIN( P, N )
         IF( ABS( B( I, I ) ).GT.TOLB )
     $      L = L + 1
   20 CONTINUE
*
      IF( WANTV ) THEN
*
*        Copy the details of V, and form V.
*
         CALL DLAZRO( P, P, ZERO, ZERO, V, LDV )
         IF( P.GT.1 )
     $      CALL DLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ),
     $                   LDV )
         CALL DORG2R( P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO )
      END IF
*
*     Clean up B
*
      DO 40 J = 1, L - 1
         DO 30 I = J + 1, L
            B( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      IF( P.GT.L )
     $   CALL DLAZRO( P-L, N, ZERO, ZERO, B( L+1, 1 ), LDB )
*
      IF( WANTQ ) THEN
*
*        Set Q = I and Update Q := Q*P
*
         CALL DLAZRO( N, N, ZERO, ONE, Q, LDQ )
         CALL DLAPMT( FORWRD, N, N, Q, LDQ, IWORK )
      END IF
*
      IF( P.GE.L .AND. N.NE.L ) THEN
*
*        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z
*
         CALL DGERQ2( L, N, B, LDB, TAU, WORK, INFO )
*
*        Update A := A*Z'
*
         CALL DORMR2( 'Right', 'Transpose', M, N, L, B, LDB, TAU, A,
     $                LDA, WORK, INFO )
*
         IF( WANTQ ) THEN
*
*           Update Q := Q*Z'
*
            CALL DORMR2( 'Right', 'Transpose', N, N, L, B, LDB, TAU, Q,
     $                   LDQ, WORK, INFO )
         END IF
*
*        Clean up B
*
         CALL DLAZRO( L, N-L, ZERO, ZERO, B, LDB )
         DO 60 J = N - L + 1, N
            DO 50 I = J - N + L + 1, L
               B( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
*
      END IF
*
*     Let              N-L     L
*                A = ( A11    A12 ) M,
*
*     then the following does the complete QR decomposition of A11:
*
*              A11 = U*(  0  T12 )*P1'
*                      (  0   0  )
*
      DO 70 I = 1, N - L
         IWORK( I ) = 0
   70 CONTINUE
      CALL DGEQPF( M, N-L, A, LDA, IWORK, TAU, WORK, INFO )
*
*     Determine the effective rank of A11
*
      K = 0
      DO 80 I = 1, MIN( M, N-L )
         IF( ABS( A( I, I ) ).GT.TOLA )
     $      K = K + 1
   80 CONTINUE
*
*     Update A12 := U'*A12, where A12 = A( 1:M, N-L+1:N )
*
      CALL DORM2R( 'Left', 'Transpose', M, L, MIN( M, N-L ), A, LDA,
     $             TAU, A( 1, N-L+1 ), LDA, WORK, INFO )
*
      IF( WANTU ) THEN
*
*        Copy the details of U, and form U
*
         CALL DLAZRO( M, M, ZERO, ZERO, U, LDU )
         IF( M.GT.1 )
     $      CALL DLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ),
     $                   LDU )
         CALL DORG2R( M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO )
      END IF
*
      IF( WANTQ ) THEN
*
*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
*
         CALL DLAPMT( FORWRD, N, N-L, Q, LDQ, IWORK )
      END IF
*
*     Clean up A: set the strictly lower triangular part of
*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.
*
      DO 100 J = 1, K - 1
         DO 90 I = J + 1, K
            A( I, J ) = ZERO
   90    CONTINUE
  100 CONTINUE
      IF( M.GT.K )
     $   CALL DLAZRO( M-K, N-L, ZERO, ZERO, A( K+1, 1 ), LDA )
*
      IF( N-L.GT.K ) THEN
*
*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
*
         CALL DGERQ2( K, N-L, A, LDA, TAU, WORK, INFO )
*
         IF( WANTQ ) THEN
*
*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1'
*
            CALL DORMR2( 'Right', 'Transpose', N, N-L, K, A, LDA, TAU,
     $                   Q, LDQ, WORK, INFO )
         END IF
*
*        Clean up A
*
         CALL DLAZRO( K, N-L-K, ZERO, ZERO, A, LDA )
         DO 120 J = N - L - K + 1, N - L
            DO 110 I = J - N + L + K + 1, K
               A( I, J ) = ZERO
  110       CONTINUE
  120    CONTINUE
*
      END IF
*
      IF( M.GT.K ) THEN
*
*        QR factorization of A( K+1:M,N-L+1:N )
*
         CALL DGEQR2( M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO )
*
         IF( WANTU ) THEN
*
*           Update U(:,K+1:M) := U(:,K+1:M)*U1
*
            CALL DORM2R( 'Right', 'No transpose', M, M-K, MIN( M-K, L ),
     $                   A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU,
     $                   WORK, INFO )
         END IF
*
*        Clean up
*
         DO 140 J = N - L + 1, N
            DO 130 I = J - N + K + L + 1, M
               A( I, J ) = ZERO
  130       CONTINUE
  140    CONTINUE
*
      END IF
*
      RETURN
*
*     End of DGGSVP
*
      END
