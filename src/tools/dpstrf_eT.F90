module dpstrf_eT
  contains


! ==================================================================
      SUBROUTINE dpstrf_e( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
!
!
!   .. Scalar Arguments ..
      DOUBLE PRECISION   TOL
      INTEGER(KIND=4)    INFO, RANK
      INTEGER            LDA, N
      CHARACTER          UPLO
!   ..
!   .. Array Arguments ..
      DOUBLE PRECISION   A( lda, * ), WORK( 2*n )
      INTEGER(KIND=4)    PIV( n )
!   ..
!
!=====================================================================
!
!   .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!   ..
!   .. Local Scalars ..
      DOUBLE PRECISION   AJJ, DSTOP, DTEMP
      INTEGER            I, ITEMP, J, JB, K, NB, PVT
      LOGICAL            UPPER
!   ..
!   .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            ILAENV
      LOGICAL            LSAME, DISNAN
      EXTERNAL           dlamch, ilaenv, lsame, disnan
!   ..
!   .. External Subroutines ..
    !  EXTERNAL           dgemv, dpstf2, dscal, dswap, dsyrk, xerbla
!   ..
!   .. Intrinsic Functions ..
    !  INTRINSIC          max, min, sqrt, maxloc
!   ..
!   .. Executable Statements ..
!
!   Test the input parameters.
!
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPSTRF', -info )
         RETURN
      END IF
!
!   Quick return if possible
!
      IF( n.EQ.0 )   RETURN
!
!   Get block size
!
      nb = ilaenv( 1, 'DPOTRF', uplo, n, -1, -1, -1 )
      IF( nb.LE.1 .OR. nb.GE.n ) THEN
!
!      Use unblocked code
!
         CALL dpstf2( uplo, n, a( 1, 1 ), lda, piv, rank, tol, work, &
                     info )
         GO TO 200
!
      ELSE
!
!   Initialize PIV
!
         DO 100 i = 1, n
            piv( i ) = i
  100    CONTINUE
!
!   Compute stopping value
!
         pvt = 1
         ajj = a( pvt, pvt )
         DO i = 2, n
            IF( a( i, i ).GT.ajj ) THEN
               pvt = i
               ajj = a( pvt, pvt )
            END IF
         END DO
         IF( ajj.LE.zero.OR.disnan( ajj ) ) THEN
            rank = 0
            info = 1
            GO TO 200
         END IF
!
!   Compute stopping value if not supplied
!
         IF( tol.LT.zero ) THEN
            dstop = n * dlamch( 'Epsilon' ) * ajj
         ELSE
            dstop = tol
         END IF
!
!
         IF( upper ) THEN
!
!         Compute the Cholesky factorization P**T * A * P = U**T * U
!
            DO 140 k = 1, n, nb
!
!            Account for last block not being NB wide
!
               jb = min( nb, n-k+1 )
!
!            Set relevant part of first half of WORK to zero,
!            holds dot products
!
               DO 110 i = k, n
                  work( i ) = 0
  110          CONTINUE
!
               DO 130 j = k, k + jb - 1
!
!            Find pivot, test for exit, else swap rows and columns
!            Update dot products, compute possible pivots which are
!            stored in the second half of WORK
!
                  DO 120 i = j, n
!
                     IF( j.GT.k ) THEN
                        work( i ) = work( i ) + a( j-1, i )**2
                     END IF
                     work( n+i ) = a( i, i ) - work( i )
!
  120             CONTINUE
!
                  IF( j.GT.1 ) THEN
                     itemp = maxloc( work( (n+j):(2*n) ), 1 )
                     pvt = itemp + j - 1
                     ajj = work( n+pvt )
                     IF( ajj.LE.dstop.OR.disnan( ajj ) ) THEN
                        a( j, j ) = ajj
                        GO TO 190
                     END IF
                  END IF
!
                  IF( j.NE.pvt ) THEN
!
!                  Pivot OK, so can now swap pivot rows and columns
!
                     a( pvt, pvt ) = a( j, j )
                     CALL dswap( j-1, a( 1, j ), 1, a( 1, pvt ), 1 )
                     IF( pvt.LT.n ) &
                        CALL dswap( n-pvt, a( j, pvt+1 ), lda, &
                                    a( pvt, pvt+1 ), lda )
                     CALL dswap( pvt-j-1, a( j, j+1 ), lda, &
                                 a( j+1, pvt ), 1 )
!
!                  Swap dot products and PIV
!
                     dtemp = work( j )
                     work( j ) = work( pvt )
                     work( pvt ) = dtemp
                     itemp = piv( pvt )
                     piv( pvt ) = piv( j )
                     piv( j ) = itemp
                  END IF
!
                  ajj = sqrt( ajj )
                  a( j, j ) = ajj
!
!               Compute elements J+1:N of row J.
!
                  IF( j.LT.n ) THEN
                     CALL dgemv( 'Trans', j-k, n-j, -one, a( k, j+1 ), &
                                 lda, a( k, j ), 1, one, a( j, j+1 ), &
                                 lda )
                     CALL dscal( n-j, one / ajj, a( j, j+1 ), lda )
                  END IF
!
  130          CONTINUE
!
!            Update trailing matrix, J already incremented
!
               IF( k+jb.LE.n ) THEN
                  CALL dsyrk( 'Upper', 'Trans', n-j+1, jb, -one, &
                              a( k, j ), lda, one, a( j, j ), lda )
               END IF
!
  140       CONTINUE
!
         ELSE
!
!      Compute the Cholesky factorization P**T * A * P = L * L**T
!
            DO 180 k = 1, n, nb
!
!            Account for last block not being NB wide
!
               jb = min( nb, n-k+1 )
!
!            Set relevant part of first half of WORK to zero,
!            holds dot products
!
               DO 150 i = k, n
                  work( i ) = 0
  150          CONTINUE
!
               DO 170 j = k, k + jb - 1
!
!            Find pivot, test for exit, else swap rows and columns
!            Update dot products, compute possible pivots which are
!            stored in the second half of WORK
!
                  DO 160 i = j, n
!
                     IF( j.GT.k ) THEN
                        work( i ) = work( i ) + a( i, j-1 )**2
                     END IF
                     work( n+i ) = a( i, i ) - work( i )
!
  160             CONTINUE
!
                  IF( j.GT.1 ) THEN
                     itemp = maxloc( work( (n+j):(2*n) ), 1 )
                     pvt = itemp + j - 1
                     ajj = work( n+pvt )
                     IF( ajj.LE.dstop.OR.disnan( ajj ) ) THEN
                        a( j, j ) = ajj
                        GO TO 190
                     END IF
                  END IF
!
                  IF( j.NE.pvt ) THEN
!
!                  Pivot OK, so can now swap pivot rows and columns
!
                     a( pvt, pvt ) = a( j, j )
                     CALL dswap( j-1, a( j, 1 ), lda, a( pvt, 1 ), lda )
                     IF( pvt.LT.n ) &
                        CALL dswap( n-pvt, a( pvt+1, j ), 1, &
                                    a( pvt+1, pvt ), 1 )
                     CALL dswap( pvt-j-1, a( j+1, j ), 1, a( pvt, j+1 ), &
                                 lda )
!
!                  Swap dot products and PIV
!
                     dtemp = work( j )
                     work( j ) = work( pvt )
                     work( pvt ) = dtemp
                     itemp = piv( pvt )
                     piv( pvt ) = piv( j )
                     piv( j ) = itemp
                  END IF
!
                  ajj = sqrt( ajj )
                  a( j, j ) = ajj
!
!               Compute elements J+1:N of column J.
!
                  IF( j.LT.n ) THEN
                     CALL dgemv( 'No Trans', n-j, j-k, -one, &
                                a( j+1, k ), lda, a( j, k ), lda, one, &
                                a( j+1, j ), 1 )
                     CALL dscal( n-j, one / ajj, a( j+1, j ), 1 )
                  END IF
!
  170          CONTINUE
!
!            Update trailing matrix, J already incremented
!
               IF( k+jb.LE.n ) THEN
                  CALL dsyrk( 'Lower', 'No Trans', n-j+1, jb, -one, &
                           a( j, k ), lda, one, a( j, j ), lda )
               END IF
!
  180       CONTINUE
!
         END IF
      END IF
!
!   Ran to completion, A has full rank
!
      rank = n
!
      GO TO 200
  190 CONTINUE
!
!   Rank is the number of steps completed.  Set INFO = 1 to signal
!   that the factorization cannot be used to solve a system.
!
      rank = j - 1
      info = 1
!
  200 CONTINUE
      RETURN
!
!   End of DPSTRF
!
      END SUBROUTINE
      end module dpstrf_eT
