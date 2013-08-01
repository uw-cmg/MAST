!--------------------------------------------------------------------------
! These two functions check the equivalence of two vectors in the lattice
!--------------------------------------------------------------------------
!==========================================================================
      FUNCTION equiva( a, b, eps )
!==========================================================================
      USE nrtype
      IMPLICIT NONE
      LOGICAL :: equiva
      REAL(DP) :: eps 
      REAL(DP) :: a(3), b(3)

      equiva = abs(  a(1) - b(1)  - nint( a(1) - b(1) ) ) < eps .and. &
     &     abs( a(2) - b(2)  - nint( a(2) - b(2) ) ) < eps .and. &
     &     abs( a(3) - b(3)  - nint( a(3) - b(3) ) ) < eps 
!      equiva = equiva .and. .not. ( any( abs ( nint( a - b) ) > 1 ) )

      END FUNCTION

!==========================================================================
    FUNCTION equiva1( a, b, eps )
!==========================================================================
      USE nrtype
      IMPLICIT NONE
      LOGICAL :: equiva1
      REAL(DP) :: eps
      REAL(DP) :: a(3), b(3)

      equiva1 = abs(  a(1) - b(1) ) < eps .and. &
           abs( a(2) - b(2) ) < eps .and. &
           abs( a(3) - b(3) ) < eps

    END FUNCTION equiva1
