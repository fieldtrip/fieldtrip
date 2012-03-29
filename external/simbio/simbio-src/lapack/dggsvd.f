      SUBROUTINE DGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
     $                   LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,
     $                   IWORK, INFO )
*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), U( LDU, * ),
     $                   V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGSVD computes the generalized singular value decomposition (GSVD)
*  of the M-by-N matrix A and P-by-N matrix B:
*
*      U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R )               (1)
*
*  where U, V and Q are orthogonal matrices, and Z' is the transpose
*  of Z.  Let K+L = the numerical effective rank of the matrix (A',B')',
*  then R is a K+L-by-K+L nonsingular upper tridiagonal matrix, D1 and
*  D2 are "diagonal" matrices, and of the following structures,
*  respectively:
*
*  If M-K-L >= 0,
*
*     U'*A*Q = D1*( 0 R )
*
*            = K     ( I  0 ) * (  0   R11  R12 ) K
*              L     ( 0  C )   (  0    0   R22 ) L
*              M-K-L ( 0  0 )    N-K-L  K    L
*                      K  L
*
*     V'*B*Q = D2*( 0 R )
*
*            = L     ( 0  S ) * (  0   R11  R12 ) K
*              P-L   ( 0  0 )   (  0    0   R22 ) L
*                      K  L      N-K-L  K    L
*  where
*
*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*    S = diag( BETA(K+1),  ... , BETA(K+L) ), C**2 + S**2 = I.
*    The nonsingular triangular matrix R = ( R11 R12 ) is stored
*                                          (  0  R22 )
*    in A(1:K+L,N-K-L+1:N) on exit.
*
*  If M-K-L < 0,
*
*     U'*A*Q = D1*( 0 R )
*
*            = K   ( I  0    0   ) * ( 0    R11  R12  R13  ) K
*              M-K ( 0  C    0   )   ( 0     0   R22  R23  ) M-K
*                    K M-K K+L-M     ( 0     0    0   R33  ) K+L-M
*                                     N-K-L  K   M-K  K+L-M
*
*     V'*B*Q = D2*( 0 R )
*
*            = M-K   ( 0  S    0  ) * ( 0    R11  R12  R13  ) K
*              K+L-M ( 0  0    I  )   ( 0     0   R22  R23  ) M-K
*              P-L   ( 0  0    0  )   ( 0     0    0   R33  ) K+L-M
*                      K M-K K+L-M     N-K-L  K   M-K  K+L-M
*  where
*
*    C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*    S = diag( BETA(K+1),  ... , BETA(M) ), C**2 + S**2 = I.
*    R = ( R11 R12 R13 ) is a nonsingular upper triangular matrix,
*        (  0  R22 R23 )
*        (  0   0  R33 )
*    (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
*    ( 0  R22 R23 )
*    in B(M-K+1:L,N+M-K-L+1:N) on exit.
*
*  The routine computes C, S, R, and optionally the orthogonal
*  transformation matrices U, V and Q.
*
*  In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
*  A and B implicitly gives the SVD of the matrix A*inv(B):
*                       A*inv(B) = U*(D1*inv(D2))*V'.
*  If ( A',B')' has orthnormal columns, then the GSVD of A and B is also
*  equal to the CS decomposition of A and B. Furthermore, the GSVD can
*  be used to derive the solution of the eigenvalue problem:
*                       A'*A x = lambda* B'*B x.
*  In some literature, the GSVD of A and B is presented in the form
*                   U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 )          (2)
*  where U and V are orthogonal and X is nonsingular, D1 and D2 are
*  ``diagonal''.  It is easy to see that the GSVD form (1) can be
*  converted to the form (2) by taking the nonsingular matrix X as
*
*                       X = Q*( I   0    )
*                             ( 0 inv(R) ).
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
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  K       (output) INTEGER
*  L       (output) INTEGER
*          On exit, K and L specify the dimension of the subblocks
*          described in the Purpose section.
*          K + L = effective numerical rank of (A',B')'.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A contains the triangular matrix R, or part of R.
*          See Purpose for details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= MAX(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B contains the triangular matrix R if necessary.
*          See Purpose for details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDA >= MAX(1,P).
*
*  ALPHA   (output) DOUBLE PRECISION arrays, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          On exit, ALPHA and BETA contain the generalized singular
*          value pairs of A and B;
*          if M-K-L >= 0,
*            ALPHA(1:K) = ONE,  ALPHA(K+1:K+L) = C,
*            BETA(1:K)  = ZERO, BETA(K+1:K+L)  = S,
*          or if M-K-L < 0,
*            ALPHA(1:K)=ONE,  ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=ZERO
*            BETA(1:K) =ZERO, BETA(K+1:M) =S, BETA(M+1:K+L) =ONE
*          and
*            ALPHA(K+L+1:N) = ZERO
*            BETA(K+L+1:N)  = ZERO
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,M)
*          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= MAX(1,M).
*
*  V       (output) DOUBLE PRECISION array, dimension (LDV,P)
*          If JOBV = 'V', V contains the P-by-P orthogonal matrix V.
*          If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDA >= MAX(1,P).
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
*          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.
*          If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= MAX(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array,
*                      dimension (MAX(3*N,M,P)+N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output)INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, the Jacobi-type procedure failed to
*                converge.  For further details, see subroutine DTGSJA.
*
*  Internal Parameters
*  ===================
*
*  TOLA    DOUBLE PRECISION
*  TOLB    DOUBLE PRECISION
*          TOLA and TOLB are the thresholds to determine the effective
*          rank of (A',B')'. Generally, they are set to
*                   TOLA = MAX(M,N)*norm(A)*MAZHEPS,
*                   TOLB = MAX(P,N)*norm(B)*MAZHEPS.
*          The size of TOLA and TOLB may affect the size of backward
*          errors of the decomposition.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTQ, WANTU, WANTV
      INTEGER            NCYCLE
      DOUBLE PRECISION   ANORM, BNORM, TOLA, TOLB, ULP, UNFL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGGSVP, DTGSJA, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
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
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGSVD', -INFO )
         RETURN
      END IF
*
*     Compute the Frobenius norm of matrices A and B
*
      ANORM = DLANGE( '1', M, N, A, LDA, WORK )
      BNORM = DLANGE( '1', P, N, B, LDB, WORK )
*
*     Get machine precision and set up threshold for determining
*     the effective numerical rank of the matrices A and B.
*
      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP
*
*     Preprocessing
*
      CALL DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
     $             TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, WORK,
     $             WORK( N+1 ), INFO )
*
*     Compute the GSVD of two upper "triangular" matrices
*
      CALL DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB,
     $             TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,
     $             WORK, NCYCLE, INFO )
*
      RETURN
*
*     End of DGGSVD
*
      END
