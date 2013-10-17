      MODULE diagonalize
!
! diagonalizes real and complex matrices using lapack call
!
        USE nrtype
        IMPLICIT NONE

        INTERFACE diagh
           module procedure diagh_d, diagh_z
        END INTERFACE

      CONTAINS

        SUBROUTINE diagh_d( h, e, up )
          CHARACTER*1 :: up
          INTEGER :: n
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: h
          REAL(DP), DIMENSION(size(h,1)), INTENT(OUT) :: e
          INTEGER :: lwork, info
          REAL(DP), ALLOCATABLE :: work(:)
          n = size(h,1)
          if( size(h,1)/=size(h,2) ) then
             print*,'h is not a square'
             stop
          endif
          lwork = 3*n - 1
          allocate( work(lwork) )
          work = 0._dp
          call dsyev( 'V', up, n, h, n, e, work, lwork, info )
          if( info /= 0 ) then
             print*,'Something wrong in ddiagh',info
             stop
          endif
          deallocate( work )
        END SUBROUTINE diagh_d

        SUBROUTINE diagh_z( h, e, up )
          CHARACTER*1 :: up
          INTEGER :: n
          COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: h
          REAL(DP), DIMENSION(size(h,1)), INTENT(OUT) :: e
          INTEGER :: lwork, info
          REAL(DP), ALLOCATABLE :: work(:), rwork(:)
          n = size(h,1)
          if( size(h,1)/=size(h,2) ) then
             print*,'h is not a square'
             stop
          endif
          lwork = 2*n - 1
          allocate( work(2*lwork), rwork(3*n-2) )
          work = 0._dp
          rwork = 0._dp
          call zheev( 'V', up, n, h, n, e, work, lwork, rwork, info )
          if( info /= 0 ) then
             print*,'Something wrong in zdiagh',info
             stop
          endif
          deallocate( work, rwork )
        END SUBROUTINE diagh_z



      END MODULE diagonalize
