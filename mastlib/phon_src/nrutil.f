MODULE nrutil

  USE nrtype
  IMPLICIT NONE
  INTEGER, PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8

  INTERFACE inve
     MODULE PROCEDURE inve_s, inve_d, inve_z
  END INTERFACE

  INTERFACE swap
     MODULE PROCEDURE swap_i, swap_l, swap_d, swap_dv, swap_dm
  END INTERFACE

  INTERFACE vabs
     MODULE PROCEDURE vabs_dv, vabs_zv, vabs_dm, vabs_zm
  END INTERFACE

  INTERFACE vecprod
     MODULE PROCEDURE vecprod_d
  END INTERFACE

  INTERFACE outerprod
     MODULE PROCEDURE outerprod_d
  END INTERFACE

  INTERFACE arth
     MODULE PROCEDURE arth_d
  END INTERFACE

  INTERFACE assert
     MODULE PROCEDURE assert1, assert2
  END INTERFACE

  INTERFACE coth
     MODULE PROCEDURE coth_d
  END INTERFACE

CONTAINS

  !=====================================================================


  !============================================================
  SUBROUTINE inve_s(v,inv)
    !============================================================
    !
    ! inverts 3x3 matrices
    !
    USE nrtype

    IMPLICIT NONE

    REAL(SP) :: v(3,3), inv(3,3), d

    d = v(1,1) * ( v(2,2) * v(3,3) - v(2,3) * v(3,2) ) + &
         v(2,1) * ( v(3,2) * v(1,3) - v(1,2) * v(3,3) ) + &
         v(3,1) * ( v(1,2) * v(2,3) - v(1,3) * v(2,2) )

    inv(1,1) = ( v(2,2) * v(3,3) - v(2,3) * v(3,2) ) / d
    inv(1,2) = ( v(3,2) * v(1,3) - v(1,2) * v(3,3) ) / d
    inv(1,3) = ( v(1,2) * v(2,3) - v(1,3) * v(2,2) ) / d
    inv(2,1) = ( v(3,1) * v(2,3) - v(2,1) * v(3,3) ) / d
    inv(2,2) = - ( v(3,1) * v(1,3) - v(1,1) * v(3,3) ) / d
    inv(2,3) = ( v(2,1) * v(1,3) - v(1,1) * v(2,3) ) / d
    inv(3,1) = ( v(2,1) * v(3,2) - v(2,2) * v(3,1) ) / d
    inv(3,2) = ( v(3,1) * v(1,2) - v(1,1) * v(3,2) ) / d
    inv(3,3) = ( v(1,1) * v(2,2) - v(1,2) * v(2,1) ) / d

  END SUBROUTINE inve_s

  !============================================================
  SUBROUTINE inve_d(v,inv)
    !============================================================
    !
    ! inverts 3x3 matrices
    !
    USE nrtype

    IMPLICIT NONE

    REAL(DP) :: v(3,3), inv(3,3), d

    d = v(1,1) * ( v(2,2) * v(3,3) - v(2,3) * v(3,2) ) + &
         v(2,1) * ( v(3,2) * v(1,3) - v(1,2) * v(3,3) ) + &
         v(3,1) * ( v(1,2) * v(2,3) - v(1,3) * v(2,2) )

    inv(1,1) = ( v(2,2) * v(3,3) - v(2,3) * v(3,2) ) / d
    inv(1,2) = ( v(3,2) * v(1,3) - v(1,2) * v(3,3) ) / d
    inv(1,3) = ( v(1,2) * v(2,3) - v(1,3) * v(2,2) ) / d
    inv(2,1) = ( v(3,1) * v(2,3) - v(2,1) * v(3,3) ) / d
    inv(2,2) = - ( v(3,1) * v(1,3) - v(1,1) * v(3,3) ) / d
    inv(2,3) = ( v(2,1) * v(1,3) - v(1,1) * v(2,3) ) / d
    inv(3,1) = ( v(2,1) * v(3,2) - v(2,2) * v(3,1) ) / d
    inv(3,2) = ( v(3,1) * v(1,2) - v(1,1) * v(3,2) ) / d
    inv(3,3) = ( v(1,1) * v(2,2) - v(1,2) * v(2,1) ) / d

  END SUBROUTINE inve_d

  !============================================================
  SUBROUTINE inve_z(v,inv)
    !============================================================
    !
    ! inverts 3x3 matrices
    !
    USE nrtype

    IMPLICIT NONE

    COMPLEX(DP) :: v(3,3), inv(3,3), d

    d = v(1,1) * ( v(2,2) * v(3,3) - v(2,3) * v(3,2) ) + &
         v(2,1) * ( v(3,2) * v(1,3) - v(1,2) * v(3,3) ) + &
         v(3,1) * ( v(1,2) * v(2,3) - v(1,3) * v(2,2) )

    inv(1,1) = ( v(2,2) * v(3,3) - v(2,3) * v(3,2) ) / d
    inv(1,2) = ( v(3,2) * v(1,3) - v(1,2) * v(3,3) ) / d
    inv(1,3) = ( v(1,2) * v(2,3) - v(1,3) * v(2,2) ) / d
    inv(2,1) = ( v(3,1) * v(2,3) - v(2,1) * v(3,3) ) / d
    inv(2,2) = - ( v(3,1) * v(1,3) - v(1,1) * v(3,3) ) / d
    inv(2,3) = ( v(2,1) * v(1,3) - v(1,1) * v(2,3) ) / d
    inv(3,1) = ( v(2,1) * v(3,2) - v(2,2) * v(3,1) ) / d
    inv(3,2) = ( v(3,1) * v(1,2) - v(1,1) * v(3,2) ) / d
    inv(3,3) = ( v(1,1) * v(2,2) - v(1,2) * v(2,1) ) / d

  END SUBROUTINE inve_z

  SUBROUTINE swap_i( a, b )
    !         Swap the contents of a and b         
    INTEGER, INTENT(INOUT) :: a, b
    INTEGER :: dum
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_l( a, b )
    LOGICAL, INTENT(INOUT) :: a, b
    LOGICAL :: dum
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_l

  SUBROUTINE swap_d( a, b )
    REAL(DP), INTENT(INOUT) :: a, b
    REAL(DP) :: dum
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_d

  SUBROUTINE swap_dv( a, b )
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a, b
    REAL(DP), DIMENSION(size(a)) :: dum
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_dv

  SUBROUTINE swap_dm( a, b )
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a, b
    REAL(DP), DIMENSION(size(a,1),size(a,2)) :: dum
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_dm

  !===========================================================================

  FUNCTION vabs_dv( a )
    !       Compute |a|
    REAL(DP), DIMENSION(:), INTENT(IN) :: a
    REAL(DP) :: vabs_dv
    vabs_dv = sqrt( dot_product( a, a ) )
  END FUNCTION vabs_dv

  FUNCTION vabs_zv( a )
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
    REAL(DP) :: vabs_zv
    vabs_zv = sqrt( dot_product( a, conjg(a) ) )
  END FUNCTION vabs_zv

  FUNCTION vabs_dm( a )
    !       Compute |a|
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
    REAL(DP) :: vabs_dm
    INTEGER :: i, j
    vabs_dm = 0
    do i = 1, size(a,1)
       do j = 1, size(a,2)
          vabs_dm = vabs_dm + abs( a(i,j) )
       enddo
    enddo
  END FUNCTION vabs_dm

  FUNCTION vabs_zm( a )
    !       Compute |a|
    COMPLEX(DP), DIMENSION(:,:), INTENT(IN) :: a
    REAL(DP) :: vabs_zm
    INTEGER :: i, j
    vabs_zm = 0
    do i = 1, size(a,1)
       do j = 1, size(a,2)
          vabs_zm = vabs_zm + abs( a(i,j) )
       enddo
    enddo
  END FUNCTION vabs_zm


  !===========================================================================

  SUBROUTINE vecprod_d( a, b, c )
    REAL(DP), INTENT(IN) :: a(3), b(3)
    REAL(DP), INTENT(OUT) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  END SUBROUTINE vecprod_d

  !===========================================================================

  FUNCTION outerprod_d( a, b )
    REAL(DP), DIMENSION(:), INTENT(IN) :: a, b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a)) 
  END FUNCTION outerprod_d

  !===========================================================================

  FUNCTION arth_d( first, increment, n )
    REAL(DP), INTENT(IN) :: first, increment
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER :: k,k2
    REAL(DP) :: temp
    if( n > 0 ) arth_d(1) = first
    if( n <= NPAR_ARTH ) then
       do k=2,n
          arth_d(k) = arth_d(k-1) + increment
       enddo
    else
       do k=2,NPAR2_ARTH
          arth_d(k) = arth_d(k-1) + increment
       enddo
       temp = increment * NPAR2_ARTH
       k = NPAR2_ARTH
       do 
          if( k >= n ) exit
          k2 = k + k
          arth_d(k+1:min(k2,n)) = temp + arth_d(1:min(k,n-k))
          temp = temp + temp
          k = k2
       enddo
    endif
  END FUNCTION arth_d

  !===========================================================================

  SUBROUTINE assert1( n1, string )
    ! Report and die if any logical is false
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if( .not. n1 )then
       write(*,*) 'nrerror, an assertion failed with this tag:',&
            &        string
       STOP 'program terminated by assert1'
    endif
  END SUBROUTINE assert1

  SUBROUTINE assert2( n1, n2, string )
    ! Report and die if any logical is false
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1, n2
    if( .not. ( n1.and.n2 ) )then
       write(*,*) 'nrerror, an assertion failed with this tag: ',&
            &        string
       STOP ' program terminated by assert2'
    endif
  END SUBROUTINE assert2

  !==========================================================================

  FUNCTION coth_d( x )
    REAL(DP), INTENT(IN) :: x
    REAL(DP), PARAMETER :: eps = 1.d-6
    REAL(DP) :: coth_d

    if( x > eps ) then
       coth_d = ( exp(x) + exp(-x) ) / ( exp(x) - exp(-x) )
    else
       print*,'error in coth_d, argument too small ',x
    endif

  END FUNCTION coth_d

END MODULE nrutil






