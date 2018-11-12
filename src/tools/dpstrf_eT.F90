module dpstrf_et
  contains


! ==================================================================
      subroutine dpstrf_e( uplo, n, a, lda, piv, rank, tol, work, info )
!
!
!   .. scalar arguments ..
      double precision   tol
      integer            info, rank
      integer            lda, n
      character          uplo
!   ..
!   .. array arguments ..
      double precision   a( lda, * ), work( 2*n )
      integer            piv( n )
!   ..
!
!=====================================================================
!   .. integers kind=4 ..
!
      integer(kind=4)    info4, rank4
      integer(kind=4)    lda4, n4
      integer(kind=4)    piv4( n )
!
!=====================================================================
!
!   .. parameters ..
      double precision   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!   ..
!   .. local scalars ..
      double precision   ajj, dstop, dtemp
      integer            i, itemp, j, jb, k, nb, pvt
      logical            upper
!   ..
!   .. external functions ..
      double precision   dlamch
      integer            ilaenv
      logical            lsame, disnan
      external           dlamch, ilaenv, lsame, disnan
!   ..
!   .. external subroutines ..
    !  external           dgemv, dpstf2, dscal, dswap, dsyrk, xerbla
!   ..
!   .. intrinsic functions ..
    !  intrinsic          max, min, sqrt, maxloc
!   ..
!   .. executable statements ..
!
!   test the input parameters.
!
      info = 0
      upper = lsame( uplo, 'u' )
      if( .not.upper .and. .not.lsame( uplo, 'l' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, n ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'dpstrf', -info )
         return
      end if
!
!   quick return if possible
!
      if( n.eq.0 )   return
!
!   get block size
!
      nb = ilaenv( 1, 'dpotrf', uplo, n, -1, -1, -1 )
      if( nb.le.1 .or. nb.ge.n ) then
!
!      use unblocked code
!
         n4 = int(n,4)
         lda4 = int(lda,4)
         rank4 = int(rank,4)
         info4 = int(info,4)
!
         do i = 1,n
            piv4(i) = int(piv(i),4)
         enddo
!
         call dpstf2( uplo, n4, a( 1, 1 ), lda4, piv4, rank4, tol, work, &
                     info4 )
!
         n = int(n4,8)
         lda = int(lda4,8)
         rank = int(rank4,8)
         info = int(info4,8)
!
         do i = 1,n
            piv(i) = int(piv4(i),8)
         enddo
!
         go to 200
!
      else
!
!   initialize piv
!
         do 100 i = 1, n
            piv( i ) = i
  100    continue
!
!   compute stopping value
!
         pvt = 1
         ajj = a( pvt, pvt )
         do i = 2, n
            if( a( i, i ).gt.ajj ) then
               pvt = i
               ajj = a( pvt, pvt )
            end if
         end do
         if( ajj.le.zero.or.disnan( ajj ) ) then
            rank = 0
            info = 1
            go to 200
         end if
!
!   compute stopping value if not supplied
!
         if( tol.lt.zero ) then
            dstop = n * dlamch( 'epsilon' ) * ajj
         else
            dstop = tol
         end if
!
!
         if( upper ) then
!
!         compute the cholesky factorization p**t * a * p = u**t * u
!
            do 140 k = 1, n, nb
!
!            account for last block not being nb wide
!
               jb = min( nb, n-k+1 )
!
!            set relevant part of first half of work to zero,
!            holds dot products
!
               do 110 i = k, n
                  work( i ) = 0
  110          continue
!
               do 130 j = k, k + jb - 1
!
!            find pivot, test for exit, else swap rows and columns
!            update dot products, compute possible pivots which are
!            stored in the second half of work
!
                  do 120 i = j, n
!
                     if( j.gt.k ) then
                        work( i ) = work( i ) + a( j-1, i )**2
                     end if
                     work( n+i ) = a( i, i ) - work( i )
!
  120             continue
!
                  if( j.gt.1 ) then
                     itemp = maxloc( work( (n+j):(2*n) ), 1 )
                     pvt = itemp + j - 1
                     ajj = work( n+pvt )
                     if( ajj.le.dstop.or.disnan( ajj ) ) then
                        a( j, j ) = ajj
                        go to 190
                     end if
                  end if
!
                  if( j.ne.pvt ) then
!
!                  pivot ok, so can now swap pivot rows and columns
!
                     a( pvt, pvt ) = a( j, j )
                     call dswap( j-1, a( 1, j ), 1, a( 1, pvt ), 1 )
                     if( pvt.lt.n ) &
                        call dswap( n-pvt, a( j, pvt+1 ), lda, &
                                    a( pvt, pvt+1 ), lda )
                     call dswap( pvt-j-1, a( j, j+1 ), lda, &
                                 a( j+1, pvt ), 1 )
!
!                  swap dot products and piv
!
                     dtemp = work( j )
                     work( j ) = work( pvt )
                     work( pvt ) = dtemp
                     itemp = piv( pvt )
                     piv( pvt ) = piv( j )
                     piv( j ) = int(itemp,4)
                  end if
!
                  ajj = sqrt( ajj )
                  a( j, j ) = ajj
!
!               compute elements j+1:n of row j.
!
                  if( j.lt.n ) then
                     call dgemv( 'trans', j-k, n-j, -one, a( k, j+1 ), &
                                 lda, a( k, j ), 1, one, a( j, j+1 ), &
                                 lda )
                     call dscal( n-j, one / ajj, a( j, j+1 ), lda )
                  end if
!
  130          continue
!
!            update trailing matrix, j already incremented
!
               if( k+jb.le.n ) then
                  call dsyrk( 'upper', 'trans', n-j+1, jb, -one, &
                              a( k, j ), lda, one, a( j, j ), lda )
               end if
!
  140       continue
!
         else
!
!      compute the cholesky factorization p**t * a * p = l * l**t
!
            do 180 k = 1, n, nb
!
!            account for last block not being nb wide
!
               jb = min( nb, n-k+1 )
!
!            set relevant part of first half of work to zero,
!            holds dot products
!
               do 150 i = k, n
                  work( i ) = 0
  150          continue
!
               do 170 j = k, k + jb - 1
!
!            find pivot, test for exit, else swap rows and columns
!            update dot products, compute possible pivots which are
!            stored in the second half of work
!
                  do 160 i = j, n
!
                     if( j.gt.k ) then
                        work( i ) = work( i ) + a( i, j-1 )**2
                     end if
                     work( n+i ) = a( i, i ) - work( i )
!
  160             continue
!
                  if( j.gt.1 ) then
                     itemp = maxloc( work( (n+j):(2*n) ), 1 )
                     pvt = itemp + j - 1
                     ajj = work( n+pvt )
                     if( ajj.le.dstop.or.disnan( ajj ) ) then
                        a( j, j ) = ajj
                        go to 190
                     end if
                  end if
!
                  if( j.ne.pvt ) then
!
!                  pivot ok, so can now swap pivot rows and columns
!
                     a( pvt, pvt ) = a( j, j )
                     call dswap( j-1, a( j, 1 ), lda, a( pvt, 1 ), lda )
                     if( pvt.lt.n ) &
                        call dswap( n-pvt, a( pvt+1, j ), 1, &
                                    a( pvt+1, pvt ), 1 )
                     call dswap( pvt-j-1, a( j+1, j ), 1, a( pvt, j+1 ), &
                                 lda )
!
!                  swap dot products and piv
!
                     dtemp = work( j )
                     work( j ) = work( pvt )
                     work( pvt ) = dtemp
                     itemp = piv( pvt )
                     piv( pvt ) = piv( j )
                     piv( j ) = int(itemp,4)
                  end if
!
                  ajj = sqrt( ajj )
                  a( j, j ) = ajj
!
!               compute elements j+1:n of column j.
!
                  if( j.lt.n ) then
                     call dgemv( 'no trans', n-j, j-k, -one, &
                                a( j+1, k ), lda, a( j, k ), lda, one, &
                                a( j+1, j ), 1 )
                     call dscal( n-j, one / ajj, a( j+1, j ), 1 )
                  end if
!
  170          continue
!
!            update trailing matrix, j already incremented
!
               if( k+jb.le.n ) then
                  call dsyrk( 'lower', 'no trans', n-j+1, jb, -one, &
                           a( j, k ), lda, one, a( j, j ), lda )
               end if
!
  180       continue
!
         end if
      end if
!
!   ran to completion, a has full rank
!
      rank = n
!
      go to 200
  190 continue
!
!   rank is the number of steps completed.  set info = 1 to signal
!   that the factorization cannot be used to solve a system.
!
      rank = int(j - 1,4)
      info = 1
!
  200 continue
      return
!
!   end of dpstrf
!
   end subroutine
end module dpstrf_et
