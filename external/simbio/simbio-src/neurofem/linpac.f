c***********************************************************************
c                                                                      *
c     Bibliothek LINPACK                                               *
c                                                                      *
c     subroutine dgbc0 (abd,lda,n,ml,mu,ipvt,rcond,z)                  *
c     subroutine dgbdi (abd,lda,n,ml,mu,ipvt,det)                      *
c     subroutine dgbfa (abd,lda,n,ml,mu,ipvt,info)                     *
c     subroutine dgbsl (abd,lda,n,ml,mu,ipvt,b,job)                    *
c     subroutine dgec0 (a,lda,n,ipvt,rcond,z)                          *
c     subroutine dgedi (a,lda,n,ipvt,det,work,job)                     *
c     subroutine dgefa (a,lda,n,ipvt,info)                             *
c     subroutine dgesl (a,lda,n,ipvt,b,job)                            *
c     subroutine dpbfa (abd,lda,n,m,info)                              *
c     subroutine dpbsl (abd,lda,n,m,b)                                 *
c     subroutine dptsl (n,d,e,b)                                       *
c                                                                      *
c***********************************************************************
c                                                                      *
c     typed by and fortran 77 changes by                               *
c                                                                      *
c        adrian rienaecker, reinhard lechtape-grueter                  *
c        institut fuer maschinenelemente und maschinengestaltung       *
c        (ime), aachen                                                 *
c                                                                      *
c@**********************************************************************
c                                                                      *
c     subroutine dgbc0 (abd,lda,n,ml,mu,ipvt,rcond,z)                  *
c                                                                      *
c                                                                      *
c     dgbc0 factors a real band matrix by gaussian                     *
c     elimination and estimates the condition of the matrix.           *
c                                                                      *
c     if rcond is not needed, dgbfa is slightly faster.                *
c     to solve  a*x = b , follow dgbc0 by dgbsl.                       *
c     to compute determinant(a) , follow dgbc0 by dgbdi.               *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        abd      real(lda , n)                                        *
c                 contains the matrix in band storage. the columns     *
c                 of the matrix are stored in the columns of abd and   *
c                 the diagonals of the matrix are stored in rows       *
c                 ml+1 through 2*ml+mu+1 of abd.                       *
c                 see the comments below for more details.             *
c                                                                      *
c        lda      integer                                              *
c                 the leading dimension of the array abd.              *
c                 lda must be .ge. 2*ml+mu+1 .                         *
c                                                                      *
c        n        integer                                              *
c                 the order of the original matrix                     *
c                                                                      *
c        ml       integer                                              *
c                 number of diagonals below the main diagonal.         *
c                 0 .le. ml .lt. n .                                   *
c                                                                      *
c        mu       integer                                              *
c                 number of diagonals above the main diagonal.         *
c                 0 .le. mu .lt. n .                                   *
c                 more efficient if ml .le. mu                         *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        abd      an upper triangular matrix in band storage and       *
c                 the multipliers which were used to obtain it.        *
c                 the factorization can be written  a = l*u  where     *
c                 l is a product of permutation and unit lower         *
c                 triangular matrices and  u  is upper triangular.     *
c                                                                      *
c        ipvt     integer(n)                                           *
c                 an integer vector of pivot indices.                  *
c                                                                      *
c        rcond    real                                                 *
c                 an estimate of the reciprocal condition of a .       *
c                 for the system  a*x = b , relative perturbations     *
c                 in a and b of size epsilon may cause                 *
c                 relative perturbations in x of size epsilon/rcond .  *
c                 if rcond is so small that the logical expression     *
c                          1.0 + rcond .eq. 1.0                        *
c                 is true, then a may be singular to working           *
c                 precision. in particular, rcond is zero if           *
c                 exact singularity is detected or the estimate        *
c                 underflows.                                          *
c                                                                      *
c        z        real(n)                                              *
c                 a work vector whose contents are usually unimportant.*
c                 if a is close to a singular matrix, then z is        *
c                 an approximate null vector in the sense that         *
c                 norm(a*z) = rcond*norm(a)*norm(z)                    *
c                                                                      *
c      and storage                                                     *
c                                                                      *
c           if a is a band matrix, the following program segment       *
c           will set up the input.                                     *
c                                                                      *
c                 ml = (band width below the diagonal)                 *
c                 mu = (band width above the diagonal)                 *
c                 m = ml + mu + 1                                      *
c                 do 20 j = 1,n                                        *
c                    i1 = max(1,j-mu)                                  *
c                    i2 = min(n,j+ml)                                  *
c                    do 10 i = i1,i2                                   *
c                       k = i - j + m                                  *
c                       abd(k,j) = a(i,j)                              *
c              10    continue                                          *
c              20 continue                                             *
c                                                                      *
c           this uses rows ml+1 through 2*ml+mu+1 of abd.              *
c           in addition, the first ml rows in abd are used for         *
c           elements generated during the triangularization.           *
c           the total number of rows needed in abd is 2*ml*mu+1 .      *
c           the ml+mu by ml+mu upper left triangle and the             *
c           ml by ml lower right triangle are not referenced.          *
c                                                                      *
c     example.. if the original matrix is                              *
c                                                                      *
c           11 12 13  0  0  0                                          *
c           21 22 23 24  0  0                                          *
c            0 32 33 34 35  0                                          *
c            0  0 43 44 45 46                                          *
c            0  0  0 54 55 56                                          *
c            0  0  0  0 65 66                                          *
c                                                                      *
c        then n=6, ml=1, mu=2, lda .ge. 5 and abd should contain       *
c                                                                      *
c            *  *  *  +  +  +         * = not used                     *
c            *  * 13 24 35 46         + = used for pivoting            *
c            * 12 23 34 45 56                                          *
c           11 22 33 44 55 66                                          *
c           21 32 43 54 65  *                                          *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c        blas    ... daxpy, ddot, dscal, dasum                         *
c        fortran ... abs,max,min,sign                                  *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c        real dddot,ek,t,wk,wkm                                        *
c        real anorm,s,dasum,sm,ynorm                                   *
c        integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm                 *
c                                                                      *
c***********************************************************************
      subroutine dgbc0 (abd,lda,n,ml,mu,ipvt,rcond,z)
c
      implicit double precision (a-h,o-z)
      dimension abd(lda,n),z(n),ipvt(n)
c
c.....compute 1-norm of a
c
      anorm = 0.d0
      l = ml + 1
      is = l + mu
      do j = 1, n
         anorm = max(anorm,dasum(l,abd(is,j),1))
         if (is .gt. ml + 1)  is =is-1
         if (j  .le. mu    )  l  = l+1
         if (j  .ge. n-ml  )  l  = l-1
      end do
c
c.....factor
c
      call dgbfa(abd,lda,n,ml,mu,ipvt,info)
c
c.....rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c.....estimate = norm(z)/norm(y) where a*z = y and trans(a)*y = e .
c.....trans(a) is the transpose of a . the components of e are
c.....chosen to cause maximum local growth in the elements of w where
c.....trans(u)*w = e. the vectors are frequently rescaled to avoid
c.....overflow.
c
c.....solve trans(u)*w = e
c
      ek = 1.d0
      call dset (n,0.d0,z,1)
      m = ml + mu + 1
      ju = 0
      do k = 1,n
         if (z(k) .ne. 0.d0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .gt. abs(abd(m,k))) then
            s = abs(abd(m,k))/abs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
         end if
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (abd(m,k) .eq. 0.d0) then
            wk = 1.d0
            wkm = 1.d0
         else
            wk = wk/abd(m,k)
            wkm = wkm/abd(m,k)
         end if
         kp1 = k + 1
         ju = min(max(ju,mu+ipvt(k)),n)
         mm = m
         if (kp1 .le. ju) then
            do j = kp1,ju
               mm = mm - 1
               sm = sm + abs(z(j)+wkm*abd(mm,j))
               z(j) = z(j) + wk*abd(mm,j)
               s = s + abs(z(j))
            end do
            if (s .lt. sm) then
               t = wkm - wk
               wk = wkm
               mm = m
               do j = kp1,ju
                  mm = mm -1
                  z(j) = z(j) + t*abd(mm,j)
               end do
            end if
         end if
         z(k) = wk
      end do
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c.....solve trans(l)*y = w
c
      do kb = 1,n
         k = n + 1 - kb
         lm = min(ml,n-k)
         if (k .lt. n) z(k) = z(k) + ddot(lm,abd(m+1,k),1,z(k+1),1)
         if (abs(z(k)) .gt. 1.d0) then
            s = 1.d0/abs(z(k))
            call dscal(n,s,z,1)
         end if
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
      end do
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = 1.d0
c
c.....solve l*v = y
c
      do k = 1,n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         lm = min(ml,n-k)
         if (k .lt. n) call daxpy(lm,t,abd(m+1,k),1,z(k+1),1)
         if (abs(z(k)) .gt. 1.d0) then
            s = 1.d0/abs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
         end if
      end do
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c.....solve u*z = w
c
      do kb = 1,n
         k = n + 1 - kb
         if (abs(z(k)) .gt. abs(abd(m,k))) then
            s = abs(abd(m,k))/abs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
         end if
         if (abd(m,k) .ne. 0.d0) z(k) = z(k)/abd(m,k)
         if (abd(m,k) .eq. 0.d0) z(k) = 1.d0
         lm = min(k,m) - 1
         la = m - lm
         lz = k - lm
         t = -z(k)
         call daxpy(lm,t,abd(la,k),1,z(lz),1)
      end do
c
c.....make znorm = 1.d0
c
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
      if (anorm .ne. 0.d0) rcond = ynorm/anorm
      if (anorm .eq. 0.d0) rcond = 0.d0
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgbdi (abd,lda,n,ml,mu,ipvt,det)                      *
c                                                                      *
c     dgbdi computes the determinant of a band matrix                  *
c     using the factors computed by dgbco or dgbfa.                    *
c     if the inverse is neededm, use dgbsl n times.                    *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        abd     real(lda,n)                                           *
c                the output from dgbco or dgbfa.                       *
c                                                                      *
c        lda     integer                                               *
c                the leading dimension of the array abd.               *
c                                                                      *
c        n       integer                                               *
c                the order of the original matrix.                     *
c                                                                      *
c        ml      integer                                               *
c                number of diagonals below the main diagonal.          *
c                                                                      *
c        mu      integer                                               *
c                number of diagonals above the main diagonal.          *
c                                                                      *
c        ipvt    integer(n)                                            *
c                the pivot vector from dgbco or dgbfa.                 *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        det     real(2)                                               *
c        determinant of original matrix.                               *
c        determinant = det(1) * 10.0**det(2)                           *
c        with 1.0 .le. abs(det(1)) .lt. 10.0                           *
c        or det(1) = 0.0                                               *
c                                                                      *
c     linpack. this version dated 08/14/78.                            *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     fortran abs                                                      *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c***********************************************************************
      subroutine dgbdi (abd,lda,n,ml,mu,ipvt,det)
c
      integer lda,n,ml,mu,ipvt(n),i
      double precision abd(lda,n),det(2),ten
c
      m = ml + mu + 1
      det(1) = 1.d0
      det(2) = 0.d0
      ten = 10.d0
      do i=1,n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = abd(m,i)*det(1)
c
c........exit
c
         if (det(1) .eq. 0.d0) goto 9000
   10    if (abs(det(1)) .ge. 1.d0) goto 20
            det(1) = ten*det(1)
            det(2) = det(2) - 1.d0
         goto 10
   20    if (abs(det(1)) .lt. ten) goto 30
            det(1) = det(1)/ten
            det(2) = det(2) + 1.d0
         goto 20
   30    continue
      end do
 9000 return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgbfa (abd,lda,n,ml,mu,ipvt,info)                     *
c                                                                      *
c                                                                      *
c     dgbfa factors a real band matrix by elimination.                 *
c                                                                      *
c     dgbfa is usually called by dgbco, but can be called              *
c     directly with a saving in time if rcond is not needed.           *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        abd   real(lda,n)                                             *
c              contains the matrix in band storage. the columns        *
c              of the matrix are stored in the columbs of abd and      *
c              the diagonals of the matrix are stored in rows          *
c              ml+1 through 2*ml+mu+1 of abd.                          *
c              see the comments below for details.                     *
c                                                                      *
c        lda   integer                                                 *
c              the leading dimension of the array abd.                 *
c              lda must be .ge. 2*ml+mu+1.                             *
c                                                                      *
c        n     integer                                                 *
c              the order of the original matrix.                       *
c                                                                      *
c        ml    integer                                                 *
c              number of diagonals below the main diagonal.            *
c              0 .le. ml .lt. n .                                      *
c                                                                      *
c        mu    integer                                                 *
c              number of diagonals above the main diagonal.            *
c              0 .le. ml .lt. n .                                      *
c              more efficient if ml .le. mu .                          *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        abd   an upper triangular matrix in band storage and          *
c              the multipliers which were used to obtain it.           *
c              the factorisation can be written  a = l*u  where        *
c              l is a product of permutation and unit lower            *
c              triangular matrices and u is upper triangular.          *
c                                                                      *
c        ipvt  integer(n)                                              *
c              = 0  normal value                                       *
c              = k  if u(k,k) .eq. 0.0 .  this is not an error         *
c                   condition for this subroutine, but it does         *
c                   indicate that dgbsl will divide by zero if         *
c                   called. use rcond in dgbco for a reliable          *
c                   indication of singularity.                         *
c                                                                      *
c     band storage                                                     *
c                                                                      *
c           if a is a band matrix, the following program segment       *
c           will set up the input.                                     *
c                                                                      *
c                 ml = (band width below the diagonal)                 *
c                 mu = (band width above the diagonal)                 *
c                 m = ml + mu + 1                                      *
c                 do 20 j = 1,n                                        *
c                    i1 = max(1,j-mu)                                  *
c                    i2 = min(n,j+ml)                                  *
c                    do 10 i = i1,i2                                   *
c                       k = i - j + m                                  *
c                       abd(k,j) = a(i,j)                              *
c              10    continue                                          *
c              20 continue                                             *
c                                                                      *
c           this uses rows ml+1 through 2*ml+mu+1 of abd.              *
c           in addition, the first ml rows in abd are used for         *
c           elements generated during the triangularization.           *
c           the total number of rows needed in abd is 2*ml*mu+1 .      *
c           the ml+mu by ml+mu upper left triangle and the             *
c           ml by ml lower right triangle are not referenced.          *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     blas  daxpy, dscal, idamax                                       *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c     real t                                                           *
c     integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1            *
c                                                                      *
c***********************************************************************
      subroutine dgbfa (abd,lda,n,ml,mu,ipvt,info)
c
      implicit double precision (a-h,o-z)
      integer lda,n,ml,mu,ipvt,info
      dimension ipvt(n),abd(lda,n)
c
      m = ml + mu + 1
      info = 0
c
c.....zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min(n,m) - 1
      if (j1 .ge. j0) then
         do jz = j0,j1
            i0 = m + 1 - jz
            do i = i0,ml
               abd(i,jz) = 0.d0
            end do
         end do
      end if
      jz = j1
      ju = 0
c
c.....gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .ge. 1) then
         do k = 1,nm1
            kp1 = k + 1
c
c...........zero next fill-in column
c
            jz = jz + 1
            if (jz .le. n) then
               if (ml .ge. 1) then
                  do i = 1,ml
                     abd(i,jz) = 0.d0
                  end do
               end if
            end if
c
c...........find l = pivot index
c
            lm = min(ml,n-k)
            l = idamax(lm+1,abd(m,k),1) + m - 1
            ipvt(k) = l + k - m
c
c...........zero pivot implies this column already triangularized
c
            if (abd(l,k) .ne. 0.d0) then
c
c..............interchange if necessary
c
               if (l .ne. m) then
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               end if
c
c..............compute multipliers
c
               t = -1.d0/abd(m,k)
               call dscal(lm,t,abd(m+1,k),1)
c
c..............row elimination with column indexing
c
               ju = min(max(ju,mu+ipvt(k)),n)
               mm = m
               if (ju .ge. kp1) then
                  do j = kp1,ju
                     l = l -1
                     mm = mm -1
                     t = abd(l,j)
                     if (l .ne. mm) then
                        abd(l,j) = abd(mm,j)
                        abd(mm,j) = t
                     end if
                     call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
                  end do
               end if
            else
               info = k
            end if
         end do
      end if
      ipvt(n) = n
      if (abd(m,n) .eq. 0.d0) info = n
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgbsl (abd,lda,n,ml,mu,ipvt,b,job)                    *
c                                                                      *
c     dgbsl solves the real band system                                *
c     a * x = b or trans(a) * x = b                                    *
c     using the factors computed by dgbco or dgbfa                     *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        abd     real(lda,n)                                           *
c                the output from dgbco or dgbfa.                       *
c                                                                      *
c        lda     integer                                               *
c                the leading dimension of the array abd.               *
c                                                                      *
c        n       integer                                               *
c                the order of the original matrix.                     *
c                                                                      *
c        ml      integer                                               *
c                number of diagonals below the main diagonal.          *
c                                                                      *
c        mu      integer                                               *
c                number of diagonals above the main diagonal.          *
c                                                                      *
c        ipvt    integer(n)                                            *
c                the pivot vector from dgbco or dgbfa.                 *
c                                                                      *
c        b       real(n)                                               *
c                the right hand side vector.                           *
c                                                                      *
c        job     integer                                               *
c                =0        to solve a*x = b ,                          *
c                =nonzero  to solve trans(a)*x = b , where             *
c                          trans(a) is the transpose.                  *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        b       the solution vector x .                               *
c                                                                      *
c     error condition                                                  *
c                                                                      *
c        a division by zero will occur if the input factor contains a  *
c        zero on the diagonal. technically this indicates singularity  *
c        but it is often caused by improper arguments or improper      *
c        setting of lda. it will not occur if the subroutines are      *
c        called correctly and if dgbco has set rcond .gt. 0.d0         *
c                                                                      *
c     to compute inverse(a) * c where c is a matrix                    *
c     with p columns                                                   *
c          call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)                    *
c          if (rcond is too small) goto ...                            *
c          do 10 j = 1,p                                               *
c             call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)                *
c       10 continue                                                    *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     blas daxpy,ddot                                                  *
c     fortran min                                                      *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c     real ddot,t                                                      *
c     integer k,kb,l,la,lb,lm,m,nm1                                    *
c                                                                      *
c***********************************************************************
      subroutine dgbsl (abd,lda,n,ml,mu,ipvt,b,job)
c
      implicit double precision (a-h,o-z)
      integer lda,n,ml,mu,ipvt,job
      dimension ipvt(n),abd(lda,n),b(n)
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .eq. 0) then
c
c........job = 0 , solve a * x = b
c........first solve l*y = b
c
         if (ml .ne. 0) then
            if (nm1 .ge. 1) then
               do k = 1,nm1
                  lm = min(ml,n-k)
                  l = ipvt(k)
                  t = b(l)
                  if (l .ne. k) then
                     b(l) = b(k)
                     b(k) = t
                  end if
                  call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
               end do
            end if
         end if
c
c........now solve u*x = y
c
         do kb = 1,n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min(k,m) - 1
            la = m - lm
            lb = k -lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
         end do
      else
c
c........job = nonzero, solve trans(a) * x = b
c........first solve trans(u)*y = b
c
         do k = 1,n
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
         end do
c
c........now solve trans(l)*x = y
c
         if (ml .ne. 0) then
            if (nm1 .ge. 1) then
               do kb = 1,nm1
                  k = n -kb
                  lm = min(ml,n-k)
                  b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
                  l = ipvt(k)
                  if (l .ne. k) then
                     t = b(l)
                     b(l) = b(k)
                     b(k) = t
                  end if
               end do
            end if
         end if
      end if
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgec0 (a,lda,n,ipvt,rcond,z)                          *
c                                                                      *
c     dgec0 factors a real matrix by gaussian                          *
c     elimination and estimates the condition of the matrix.           *
c                                                                      *
c     if rcond is not needed, dgefa is slightly faster.                *
c     to solve  a*x = b , follow dgec0 by dgesl.                       *
c     to compute determinant(a), follow dgec0 by dgedi.                *
c     to compute inverse(a), follow dgec0 by dgedi.                    *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        a        real(lda , n)                                        *
c                 the matrix to be factored.                           *
c                                                                      *
c        lda      integer                                              *
c                 the leading dimension of the array a.                *
c                                                                      *
c        n        integer                                              *
c                 the order of the original matrix a.                  *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        a        an upper triangular matrix and the multipliers       *
c                 which were used to obtain it.                        *
c                 the factorization can be written  a = l*u  where     *
c                 l is a product of permutation and unit lower         *
c                 triangular matrices and  u  is upper triangular.     *
c                                                                      *
c        ipvt     integer(n)                                           *
c                 an integer vector of pivot indices.                  *
c                                                                      *
c        rcond    real                                                 *
c                 an estimate of the reciprocal condition of a .       *
c                 for the system  a*x = b , relative perturbations     *
c                 in a and b of size epsilon may cause                 *
c                 relative perturbations in x of size epsilon/rcond.   *
c                 if rcond is so small that the logical expression     *
c                          1.0 + rcond .eq. 1.0                        *
c                 is true, then a may be singular to working           *
c                 precision. in particular, rcond is zero if           *
c                 exact singularity is detected or the estimate        *
c                 underflows.                                          *
c                                                                      *
c        z        real(n)                                              *
c                 a work vector whose contents are usually unimportant.*
c                 if a is close to a singular matrix, then z is        *
c                 an approximate null vector in the sense that         *
c                 norm(a*z) = rcond*norm(a)*norm(z)                    *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c        linpack dgefa                                                 *
c        blas daxpy, ddot, dscal, dasum                                *
c        fortran abs,amax1,max0,min0,sign                              *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c     real ddot,ek,t,wk,wkm                                            *
c     real anorm,s,dasum,sm,ynorm                                      *
c     integer info,j,k,kb,kp1,l                                        *
c                                                                      *
c***********************************************************************
      subroutine dgec0 (a,lda,n,ipvt,rcond,z)
c
      implicit double precision (a-h,o-z)
      dimension a(lda,n),z(n),ipvt(n)
c
c.....compute 1-norm of a
c
      anorm = 0.d0
      do j = 1,n
         anorm = max(anorm,dasum(n,a(1,j),1))
      end do
c
c.....factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c.....rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c.....estimate = norm(z)/norm(y) where a*z = y and trans(a)*y = e .
c.....trans(a) is the transpose of a . the components of e are
c.....chosen to cause maximum local growth in the elements of w where
c.....trans(u)*w = e. the vectors are frequently rescaled to avoid
c.....overflow.
c
c.....solve trans(u)*w = e
c
      ek = 1.d0
      do j = 1,n
         z(j)=0.d0
      end do
      do k = 1, n
         if (z(k) .ne. 0.d0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .gt. abs(a(k,k))) then
            s = abs(a(k,k))/abs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
         end if
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (a(k,k) .ne. 0.d0) then
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         else
            wk = 1.d0
            wkm = 1.d0
         end if
         kp1 = k + 1
         if (kp1 .le. n) then
            do j = kp1,n
               sm = sm + abs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + abs(z(j))
            end do
            if (s .lt. sm) then
               t = wkm - wk
               wk = wkm
               do j = kp1,n
                  z(j) = z(j) + t*a(k,j)
               end do
            end if
         end if
         z(k) = wk
      end do
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c.....solve trans(l)*y = w
c
      do kb = 1,n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .gt. 1.d0) then
            s = 1.d0/abs(z(k))
            call dscal(n,s,z,1)
         end if
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
      end do
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.d0
c
c.....solve l*v = y
c
      do k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .gt. 1.d0) then
            s = 1.d0/abs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
         end if
      end do
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c.....solve u*z = w
c
      do kb = 1, n
         k = n + 1 - kb
         if (abs(z(k)) .gt. abs(a(k,k))) then
            s = abs(a(k,k))/abs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
         end if
         if (a(k,k) .ne. 0.d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.d0) z(k) = 1.d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
      end do
c
c.....make znorm = 1.d0
c
      s = 1.d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.d0) rcond = ynorm/anorm
      if (anorm .eq. 0.d0) rcond = 0.d0
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgedi (a,lda,n,ipvt,det,work,job)                     *
c                                                                      *
c     dgedi computes the determinant and inverse of a matrix           *
c     using the factors computed by dgeco or dgefa.                    *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        a       real(lda,n)                                           *
c                the output from dgeco or dgefa.                       *
c                                                                      *
c        lda     integer                                               *
c                the leading dimension of the array a.                 *
c                                                                      *
c        n       integer                                               *
c                the order of the original matrix.                     *
c                                                                      *
c        ipvt    integer(n)                                            *
c                the pivot vector from dgeco or dgefa.                 *
c                                                                      *
c        work    real(n)                                               *
c                work vector. contents destroyed.                      *
c                                                                      *
c        job     integer                                               *
c                = 11   both determinant and inverse.                  *
c                = 01   inverse only.                                  *
c                = 10   determinant only.                              *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        a       inverse of original matrix if requested.              *
c                otherwise unchanged.                                  *
c                                                                      *
c        det     real(2)                                               *
c                determinant of original matrix if requested.          *
c                otherwise not referenced.                             *
c                determinant = det(1) * 10.0**det(2)                   *
c                with 1.0 .le. abs(det(1)) .lt. 10.0                   *
c                or det(1) = 0.0                                       *
c                                                                      *
c     error condition                                                  *
c                                                                      *
c        a division by zero will occur if the input factor contains    *
c        a zero on the diagonal and the inverse is requested.          *
c        it will not occur if the subroutines are called correctly     *
c        and if dgec0 has set rcond .gt. 0.0 or dgefa set              *
c        info .eq. 0 .                                                 *
c                                                                      *
c     linpack. this version dated 08/14/78.                            *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     fortran abs                                                      *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c***********************************************************************
      subroutine dgedi (a,lda,n,ipvt,det,work,job)
c
      implicit double precision (a-h,o-z)
      integer lda,n,ipvt,job
      integer i,j,k,kb,kp1,l,nm1
      dimension a(lda,n),det(2),ipvt(n),work(n)
c
c.....compute determinant
c
      if (job/10 .ne. 0) then
         det(1) = 1.d0
         det(2) = 0.d0
         ten = 10.d0
         do i=1,n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c
c...........exit
c
            if (det(1) .eq. 0.d0) goto 60
   10          if (abs(det(1)) .lt. 1.d0) then
                  det(1) = ten*det(1)
                  det(2) = det(2) - 1.d0
                  goto 10
               end if
   20          if (abs(det(1)) .ge. ten) then
                  det(1) = det(1)/ten
                  det(2) = det(2) + 1.d0
                  goto 20
               end if
         end do
      end if
   60 continue
c
c.....compute inverse(u)
c
      if (mod(job,10) .ne. 0) then
         do k = 1, n
            a(k,k) = 1.d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .ge. kp1) then
               do j = kp1, n
                  t = a(k,j)
                  a(k,j) = 0.d0
                  call daxpy(k,t,a(1,k),1,a(1,j),1)
               end do
            end if
         end do
c
c........form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .ge. 1) then
            do kb = 1, nm1
               k = n - kb
               kp1 = k + 1
               do i = kp1 , n
                  work(i) = a(i,k)
                  a(i,k) = 0.d0
               end do
               do j = kp1 , n
                  t = work(j)
                  call daxpy(n,t,a(1,j),1,a(1,k),1)
               end do
               l = ipvt(k)
               if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
            end do
         end if
      end if
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgefa (a,lda,n,ipvt,info)                             *
c                                                                      *
c     dgefa factors a real matrix by elimination.                      *
c                                                                      *
c     dgefa is usually called by dgeco, but can be called              *
c     directly with a saving in time if rcond is not needed.           *
c     (time for dgec0) = (1 + 9/n)*(time for dgefa)                    *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        a     real(lda,n)                                             *
c              the matrix to be factored.                              *
c                                                                      *
c        lda   integer                                                 *
c              the leading dimension of the array a.                   *
c                                                                      *
c        n     integer                                                 *
c              the order of the original matrix.                       *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        a     an upper triangular matrix in and the multipliers       *
c              which were used to obtain it.                           *
c              the factorisation can be written  a = l*u  where        *
c              l is a product of permutation and unit lower            *
c              triangular matrices and u is upper triangular.          *
c                                                                      *
c        ipvt  integer(n)                                              *
c              = 0  normal value                                       *
c              = k  if u(k,k) .eq. 0.0 .  this is not an error         *
c                   condition for this subroutine, but it does         *
c                   indicate that dgesl will divide by zero if         *
c                   called. use rcond in dgeco for a reliable          *
c                   indication of singularity.                         *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     blas  saxpy, sscal, idamax                                       *
c                                                                      *
c     internal variables                                               *
c***********************************************************************
      subroutine dgefa (a,lda,n,ipvt,info)
c
      implicit double precision (a-h,o-z)
      integer lda,n,ipvt,info
      dimension ipvt(n),a(lda,n)
c     real t
      integer idamax,j,k,kp1,l,nm1
c
c.....gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .ge. 1) then
         do k = 1, nm1
            kp1 = k + 1
c
c...........find l = pivot index
c
            l = idamax(n-k+1,a(k,k),1) + k - 1
            ipvt(k) = l
c
c...........zero pivot implies this column already triangularized
c
            if (a(l,k) .ne. 0.d0) then
c
c..............interchange if necessary
c
               if (l .ne. k) then
                  t = a(l,k)
                  a(l,k) = a(k,k)
                  a(k,k) = t
               end if
c
c..............compute multipliers
c
               t = -1.d0/a(k,k)
               call dscal(n-k,t,a(k+1,k),1)
c
c..............row elimination with column indexing
c
               do j = kp1,n
                  t = a(l,j)
                  if (l .ne. k) then
                     a(l,j) = a(k,j)
                     a(k,j) = t
                  end if
                  call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
               end do
            else
               info = k
            end if
         end do
      end if
      ipvt(n) = n
      if (a(n,n) .eq. 0.d0) info = n
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dgesl (a,lda,n,ipvt,b,job)                            *
c                                                                      *
c     dgesl solves the real system                                     *
c     a * x = b or trans(a) * x = b                                    *
c     using the factors computed by dgeco or dgefa                     *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        a     real(lda,n)                                             *
c                the output from dgeco or dgefa.                       *
c                                                                      *
c        lda     integer                                               *
c                the leading dimension of the array a.                 *
c                                                                      *
c        n       integer                                               *
c                the order of the original matrix.                     *
c                                                                      *
c        ipvt    integer(n)                                            *
c                the pivot vector from dgeco or dgefa.                 *
c                                                                      *
c        b       real(n)                                               *
c                the right hand side vector.                           *
c                                                                      *
c        job     integer                                               *
c                =0        to solve a*x = b ,                          *
c                =nonzero  to solve trans(a)*x = b , where             *
c                          trans(a) is the transpose.                  *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        b       the solution vector x .                               *
c                                                                      *
c     error condition                                                  *
c                                                                      *
c        a division by zero will occur if the input factor contains a  *
c        zero on the diagonal. technically this indicates singularity  *
c        but it is often caused by improper arguments or improper      *
c        setting of lda. it will not occur if the subroutines are      *
c        called correctly and if dgeco has set rcond .gt. 0.           *
c                                                                      *
c     to compute inverse(a) * c where c is a matrix                    *
c     with p columns                                                   *
c          call dgeco(a,lda,n,ipvt,rcond,z)                            *
c          if (rcond is too small) goto ...                            *
c          do 10 j = 1, p                                              *
c             call dgesl(a,lda,n,ipvt,c(1,j),0)                        *
c       10 continue                                                    *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     blas daxpy,ddot                                                  *
c     fortran min                                                      *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c***********************************************************************
      subroutine dgesl (a,lda,n,ipvt,b,job)
c
      implicit double precision (a-h,o-z)
      integer lda,n,ipvt,job
      dimension ipvt(n),a(lda,n),b(n)
c
      nm1 = n - 1
      if (job .eq. 0) then
c
c........job = 0 , solve a * x = b
c........first solve l*y = b
c
         if (nm1 .ge. 1) then
            do k = 1,nm1
               l = ipvt(k)
               t = b(l)
               if (l .ne. k) then
                  b(l) = b(k)
                  b(k) = t
               end if
               call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
            end do
         end if
c
c........now solve u*x = y
c
         do kb = 1,n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
         end do
      else
c
c........job = nonzero, solve trans(a) * x = b
c........first solve trans(u)*y = b
c
         do k = 1,n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
         end do
c
c........now solve trans(l)*x = y
c
         if (nm1 .ge. 1) then
            do kb = 1,nm1
               k = n - kb
               b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .ne. k) then
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
               end if
            end do
         end if
      end if
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dpbfa (abd,lda,n,m,info)                              *
c                                                                      *
c                                                                      *
c     dpbfa factors a real symmetric positive definite matrix          *
c     stored in band form                                              *
c                                                                      *
c     dpbfa is usually called by dpbco, but can be called              *
c     directly with a saving in time if rcond is not needed.           *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        abd   real(lda,n)                                             *
c              the matrix to be factored. the columns of the upper     *
c              triangle are stored in the columns of abd and           *
c              the diagonals of the upper triangle are stored in the   *
c              rows of abd. see the comments below for details.        *
c                                                                      *
c        lda   integer                                                 *
c              the leading dimension of the array abd.                 *
c              lda must be .ge. m + 1                                  *
c                                                                      *
c        n     integer                                                 *
c              the order of the original matrix.                       *
c                                                                      *
c        m     integer                                                 *
c              the number of diagonals below the main diagonal.        *
c              0 .le. m .lt. n .                                       *
c                                                                      *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        abd   an upper triangular matrix r in band storage,           *
c              so that  a = trans(r)*r                                 *
c                                                                      *
c        info  integer                                                 *
c              = 0 for normal return                                   *
c              = k if the leading minor of order k is not              *
c                  positiv definit                                     *
c                                                                      *
c     band storage                                                     *
c                                                                      *
c           if a is a symmetric positive definite band matrix,         *
c           the following program segment will set up the input.       *
c                                                                      *
c                 m = (band width above the diagonal)                  *
c                 m = ml + mu + 1                                      *
c                 do 20 j = 1,n                                        *
c                    i1 = max(1,j-m)                                   *
c                    do 10 i = i1,j                                    *
c                       k = i - j + m +1                               *
c                       abd(k,j) = a(i,j)                              *
c              10    continue                                          *
c              20 continue                                             *
c                                                                      *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     blas  ddot                                                       *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c     real t, ddot                                                     *
c     integer ik,j,jk,k,mu                                             *
c                                                                      *
c***********************************************************************
      subroutine dpbfa (abd,lda,n,m,info)
c
      implicit double precision (a-h,o-z)
      integer lda,n,m,mu,info
      dimension abd(lda,n)
c
      do j = 1, n
         info = j
         s = 0.d0
         ik = m + 1
         jk = max(j-m,1)
         mu = max(m+2-j,1)
         if (m .lt. mu) goto 10
         do k = mu, m
            t = abd(k,j) - ddot(k-mu,abd(ik,jk),1,abd(mu,j),1)
            t = t/abd(m+1,jk)
            abd(k,j) = t
            s = s + t*t
            ik = ik - 1
            jk = jk + 1
         end do
   10    s = abd(m+1,j) - s
c
c........exit
c
         if (s .le. 0.d0) return
         abd(m+1,j) = sqrt(s)
      end do
      info = 0
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dpbsl (abd,lda,n,m,b)                                 *
c                                                                      *
c     dpbsl solves the real symmetric positiv definite band            *
c     system  a * x = b                                                *
c     using the factors computed by dpbco or dpbfa                     *
c                                                                      *
c     on entry                                                         *
c                                                                      *
c        abd     real(lda,n)                                           *
c                the output from dpbco or dpbfa.                       *
c                                                                      *
c        lda     integer                                               *
c                the leading dimension of the array abd.               *
c                                                                      *
c        n       integer                                               *
c                the order of the original matrix.                     *
c                                                                      *
c        m       integer                                               *
c                number of diagonals above the main diagonal.          *
c                                                                      *
c        b       real(n)                                               *
c                the right hand side vector.                           *
c                                                                      *
c     on return                                                        *
c                                                                      *
c        b       the solution vector x .                               *
c                                                                      *
c     error condition                                                  *
c                                                                      *
c        a division by zero will occur if the input factor contains a  *
c        zero on the diagonal. technically this indicates singularity  *
c        but it is often caused by improper arguments or improper      *
c        setting of lda. it will not occur if the subroutines are      *
c        called correctly and if dpbco has set rcond .gt. 0.d0         *
c                                                                      *
c     to compute inverse(a) * c where c is a matrix                    *
c     with p columns                                                   *
c          call dpbco(abd,lda,n,m,rcond,z,info)                        *
c          if (rcond is too small .or. info .ne. 0) goto ...           *
c          do 10 j = 1,p                                               *
c             call dpbsl(abd,lda,n,m,c(1,j))                           *
c       10 continue                                                    *
c                                                                      *
c     linpack. this version dated 08/14/78 .                           *
c     cleve moler, university of new mexico, argonne national lab.     *
c                                                                      *
c     subroutines and functions                                        *
c                                                                      *
c     blas saxpy,sdot                                                  *
c     fortran min                                                      *
c                                                                      *
c     internal variables                                               *
c                                                                      *
c     real ddot,t                                                      *
c     integer k,kb,la,lb,lm                                            *
c                                                                      *
c***********************************************************************
      subroutine dpbsl (abd,lda,n,m,b)
c
      implicit double precision (a-h,o-z)
      integer lda,n,m
      dimension abd(lda,n),b(n)
c
      do k = 1, n
         lm = min(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         t = ddot(lm,abd(la,k),1,b(lb),1)
         b(k) = (b(k) - t)/abd(m+1,k)
      end do
c
c.....solve  r*x = y
c
      do kb = 1, n
         k = n + 1 - kb
         lm = min(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         b(k) = b(k)/abd(m+1,k)
         t = -b(k)
         call daxpy(lm,t,abd(la,k),1,b(lb),1)
      end do
      return
      end
c---
      subroutine dppco(ap,n,rcond,z,info)
      integer n,info
      double precision ap(1),z(1)
      double precision rcond
c
c     dppco factors a double precision symmetric positive definite
c     matrix stored in packed form
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dppfa is slightly faster.
c     to solve  a*x = b , follow dppco by dppsl.
c     to compute  inverse(a)*c , follow dppco by dppsl.
c     to compute  determinant(a) , follow dppco by dppdi.
c     to compute  inverse(a) , follow dppco by dppdi.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      an upper triangular matrix  r , stored in packed
c                form, so that  a = trans(r)*r .
c                if  info .ne. 0 , the factorization is not complete.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.  if info .ne. 0 , rcond is unchanged.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is singular to working precision, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dppfa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dreal,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer i,ij,j,jm1,j1,k,kb,kj,kk,kp1
c
c
c     find norm of a
c
      j1 = 1
      do 30 j = 1, n
         z(j) = dasum(j,ap(j1),1)
         ij = j1
         j1 = j1 + j
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = z(i) + dabs(ap(ij))
            ij = ij + 1
   10    continue
   20    continue
   30 continue
      anorm = 0.0d0
      do 40 j = 1, n
         anorm = dmax1(anorm,z(j))
   40 continue
c
c     factor
c
      call dppfa(ap,n,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  trans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve trans(r)*w = e
c
         ek = 1.0d0
         do 50 j = 1, n
            z(j) = 0.0d0
   50    continue
         kk = 0
         do 110 k = 1, n
            kk = kk + k
            if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
            if (dabs(ek-z(k)) .le. ap(kk)) go to 60
               s = ap(kk)/dabs(ek-z(k))
               call dscal(n,s,z,1)
               ek = s*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = dabs(wk)
            sm = dabs(wkm)
            wk = wk/ap(kk)
            wkm = wkm/ap(kk)
            kp1 = k + 1
            kj = kk + k
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + dabs(z(j)+wkm*ap(kj))
                  z(j) = z(j) + wk*ap(kj)
                  s = s + dabs(z(j))
                  kj = kj + j
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  kj = kk + k
                  do 80 j = kp1, n
                     z(j) = z(j) + t*ap(kj)
                     kj = kj + j
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
c
c        solve r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (dabs(z(k)) .le. ap(kk)) go to 120
               s = ap(kk)/dabs(z(k))
               call dscal(n,s,z,1)
  120       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call daxpy(k-1,t,ap(kk+1),1,z(1),1)
  130    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
c
         ynorm = 1.0d0
c
c        solve trans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - ddot(k-1,ap(kk+1),1,z(1),1)
            kk = kk + k
            if (dabs(z(k)) .le. ap(kk)) go to 140
               s = ap(kk)/dabs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/ap(kk)
  150    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (dabs(z(k)) .le. ap(kk)) go to 160
               s = ap(kk)/dabs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call daxpy(k-1,t,ap(kk+1),1,z(1),1)
  170    continue
c        make znorm = 1.0
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0d0) rcond = ynorm/anorm
         if (anorm .eq. 0.0d0) rcond = 0.0d0
  180 continue
      return
      end
c*****
      subroutine dppdi(ap,n,det,job)
      integer n,job
      double precision ap(1)
      double precision det(2)
c
c     dppdi computes the determinant and inverse
c     of a double precision symmetric positive definite matrix
c     using the factors computed by dppco or dppfa .
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the output from dppco or dppfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        ap      the upper triangular half of the inverse .
c                the strict lower triangle is unaltered.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dpoco or dpofa has set info .eq. 0 .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal
c     fortran mod
c
c     internal variables
c
      double precision t
      double precision s
      integer i,ii,j,jj,jm1,j1,k,kj,kk,kp1,k1
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         s = 10.0d0
         ii = 0
         do 50 i = 1, n
            ii = ii + i
            det(1) = ap(ii)**2*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (det(1) .ge. 1.0d0) go to 20
               det(1) = s*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (det(1) .lt. s) go to 40
               det(1) = det(1)/s
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(r)
c
      if (mod(job,10) .eq. 0) go to 140
         kk = 0
         do 100 k = 1, n
            k1 = kk + 1
            kk = kk + k
            ap(kk) = 1.0d0/ap(kk)
            t = -ap(kk)
            call dscal(k-1,t,ap(k1),1)
            kp1 = k + 1
            j1 = kk + 1
            kj = kk + k
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = ap(kj)
               ap(kj) = 0.0d0
               call daxpy(k,t,ap(k1),1,ap(j1),1)
               j1 = j1 + j
               kj = kj + j
   80       continue
   90       continue
  100    continue
c
c        form  inverse(r) * trans(inverse(r))
c
         jj = 0
         do 130 j = 1, n
            j1 = jj + 1
            jj = jj + j
            jm1 = j - 1
            k1 = 1
            kj = j1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = ap(kj)
               call daxpy(k,t,ap(j1),1,ap(k1),1)
               k1 = k1 + k
               kj = kj + 1
  110       continue
  120       continue
            t = ap(jj)
            call dscal(j,t,ap(j1),1)
  130    continue
  140 continue
      return
      end
c*****
      subroutine dppfa(ap,n,info)
      integer n,info
      double precision ap(1)
c
c     dppfa factors a double precision symmetric positive definite
c     matrix stored in packed form.
c
c     dppfa is usually called by dppco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dppco) = (1 + 18/n)*(time for dppfa) .
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      an upper triangular matrix  r , stored in packed
c                form, so that  a = trans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas ddot
c     fortran dsqrt
c
c     internal variables
c
      double precision ddot,t
      double precision s
      integer j,jj,jm1,k,kj,kk
c     begin block with ...exits to 40
c
c
         jj = 0
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            kj = jj
            kk = 0
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               kj = kj + 1
               t = ap(kj) - ddot(k-1,ap(kk+1),1,ap(jj+1),1)
               kk = kk + k
               t = t/ap(kk)
               ap(kj) = t
               s = s + t*t
   10       continue
   20       continue
            jj = jj + j
            s = ap(jj) - s
c     ......exit
c            if (s .le. 0.0d0) go to 40
            if (s .le. 0.0d0) write(*,900) s
            if (s .le. 0.0d0) s = 0.00000001
            ap(jj) = dsqrt(s)
   30    continue
         info = 0
   40 continue
      return
  900 format (1x,'*** ERROR: FACTOR IN DPPFA < 0, set 0.0: s=',g12.5)
      end
      subroutine dppsl(ap,n,b)
      integer n
      double precision ap(1),b(1)
c
c     dppsl solves the double precision symmetric positive definite
c     system a * x = b
c     using the factors computed by dppco or dppfa.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the output from dppco or dppfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal.  technically this indicates
c        singularity but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dppco(ap,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dppsl(ap,n,c(1,j))
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,kk
c
      kk = 0
      do 10 k = 1, n
         t = ddot(k-1,ap(kk+1),1,b(1),1)
         kk = kk + k
         b(k) = (b(k) - t)/ap(kk)
   10 continue
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/ap(kk)
         kk = kk - k
         t = -b(k)
         call daxpy(k-1,t,ap(kk+1),1,b(1),1)
   20 continue
      return
      end
c@**********************************************************************
c                                                                      *
c     subroutine dptsl (n,d,e,b)                                       *
c                                                                      *
c     dptsl given a positive definite tridiagonal matrix and a right   *
c     hand side will find the solution of the associated set of        *
c     linear equations                                                 *
c     on entry                                                         *
c       n         integer                                              *
c                 is the order of the tridiagonal matrix               *
c                                                                      *
c       d         double(n)                                            *
c                 is the diagonal of the tridiagonal matrix            *
c                 on output d is destroyed                             *
c                                                                      *
c       e         double(n)                                            *
c                 is the offdiagonal of the tridiagonal matrix         *
c                 e(1) through e(n-1) should contain the offdiagonal   *
c                                                                      *
c       b         double(n) is the right hand side vektor              *
c                                                                      *
c     on return                                                        *
c                                                                      *
c       b         contains the solution                                *
c                                                                      *
c                                                                      *
c     linpack. this version dated 08/14/78                             *
c     jack dongarra, argonne national laboratory                       *
c                                                                      *
c                                                                      *
c***********************************************************************
      subroutine dptsl (n,d,e,b)
c
      integer n
      double precision d(n),e(n),b(n)
      integer k,kbm1,ke,kf,kp1,nm1,nm1d2
      double precision t1,t2
c
c.....check for 1x1 case
c
      if (n. eq. 1) then
         b(1)=b(1)/d(1)
      else
         nm1 = n-1
         nm1d2 = nm1/2
         if(n .ne. 2) then
            kbm1=n-1
c
c...........zero top half of subdiagonal and bottom half 
c...........of superdiagonal
c
            do k=1,nm1d2
               t1 = e(k)/d(k)
               d(k+1) = d(k+1) - t1*e(k)
               b(k+1) = b(k+1) - t1*b(k)
               t2 = e(kbm1)/d(kbm1+1)
               d(kbm1) = d(kbm1) - t2*e(kbm1)
               b(kbm1) = b(kbm1) - t2*b(kbm1+1)
               kbm1=kbm1-1
            end do
         end if
         kp1 = nm1d2 + 1 
c
c........clean up for possible 2x2 block at center
c
         if (mod(n,2) .eq. 0) then
            t1=e(kp1)/d(kp1)
            d(kp1+1) = d(kp1+1) - t1*e(kp1)
            b(kp1+1) = b(kp1+1) - t1*b(kp1)
            kp1 = kp1 + 1
         end if
c
c........back solve starting at the center, going 
c........towards the top and bottom
c
         b(kp1) = b(kp1)/d(kp1)
         if (n .ne. 2) then
            k=kp1-1
            ke=kp1+nm1d2-1
            do kf=kp1,ke
               b(k   ) = (b(k   )-e(k )*b(k+1)) / d(k   )
               b(kf+1) = (b(kf+1)-e(kf)*b(kf )) / d(kf+1)          
               k = k - 1
            end do
         end if
         if (mod(n,2) .eq. 0) b(1) = (b(1)-e(1)*b(2))/d(1)
      end if
      return
      end
