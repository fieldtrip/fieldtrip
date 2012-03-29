      SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, E,
     $                   LDE, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDE, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( LDE, * ), LSCALE( * ), RSCALE( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGBAK forms the right or left eigenvectors of the generalized
*  eigenvalue problem by backward transformation on the computed
*  eigenvectors of the balanced matrix output by DGGBAL.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the type of backward transformation required:
*          = 'N':  do nothing, return immediately;
*          = 'P':  do backward transformation for permutation only;
*          = 'S':  do backward transformation for scaling only;
*          = 'B':  do backward transformations for both permutation and
*                  scaling.
*          JOB must be the same as the argument JOB supplied to DGGBAL.
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  E contains right eigenvectors;
*          = 'L':  E contains left eigenvectors.
*
*  N       (input) INTEGER
*          The number of rows of the matrix E.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          The integers ILO and IHI determined by DGGBAL.
*
*  LSCALE  (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and/or scaling factors applied
*          to the left side of A and B, as returned by DGGBAL.
*
*  RSCALE  (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and/or scaling factors applied
*          to the right side of A and B, as returned by DGGBAL.
*
*  M       (input) INTEGER
*          The number of columns of the matrix E.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (LDE,M)
*          On entry, the matrix of right or left eigenvectors to be
*          transformed, as returned by DTGEVC.
*          On exit, E is overwritten by the transformed eigenvectors.
*
*  LDE     (input) INTEGER
*          The leading dimension of the matrix E. LDE >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  See R.C. Ward, Balancing the generalized eigenvalue problem,
*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.LT.ILO .OR. IHI.GT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGBAK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( M.EQ.0 )
     $   RETURN
      IF( LSAME( JOB, 'N' ) )
     $   RETURN
*
      IF( ILO.EQ.IHI )
     $   GO TO 30
*
*     Backward balance
*
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward transformation on right eigenvectors
*
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               CALL DSCAL( M, RSCALE( I ), E( I, 1 ), LDE )
   10       CONTINUE
         END IF
*
*        Backward transformation on left eigenvectors
*
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               CALL DSCAL( M, LSCALE( I ), E( I, 1 ), LDE )
   20       CONTINUE
         END IF
      END IF
*
*     Backward permutation
*
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward permutation on right eigenvectors
*
         IF( RIGHTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 50
*
            DO 40 I = ILO - 1, 1, -1
               K = RSCALE( I )
               IF( K.EQ.I )
     $            GO TO 40
               CALL DSWAP( M, E( I, 1 ), LDE, E( K, 1 ), LDE )
   40       CONTINUE
*
   50       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 70
            DO 60 I = IHI + 1, N
               K = RSCALE( I )
               IF( K.EQ.I )
     $            GO TO 60
               CALL DSWAP( M, E( I, 1 ), LDE, E( K, 1 ), LDE )
   60       CONTINUE
         END IF
*
*        Backward permutation on left eigenvectors
*
   70    CONTINUE
         IF( LEFTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 90
            DO 80 I = ILO - 1, 1, -1
               K = LSCALE( I )
               IF( K.EQ.I )
     $            GO TO 80
               CALL DSWAP( M, E( I, 1 ), LDE, E( K, 1 ), LDE )
   80       CONTINUE
*
   90       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 110
            DO 100 I = IHI + 1, N
               K = LSCALE( I )
               IF( K.EQ.I )
     $            GO TO 100
               CALL DSWAP( M, E( I, 1 ), LDE, E( K, 1 ), LDE )
  100       CONTINUE
         END IF
      END IF
*
  110 CONTINUE
*
      RETURN
*
*     End of DGGBAK
*
      END
