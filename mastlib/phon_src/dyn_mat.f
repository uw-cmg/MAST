!----------------------------------------------------------------------------------
! Contructs the dynamical matrix at any q-vector q from the force constant matrix 
!----------------------------------------------------------------------------------
!===================================================================================
SUBROUTINE dyn_mat( dyn, latt, sh, q, ntot, list )
  !=================================================================================
  USE nrtype
  USE nrutil
  USE diagonalize
  USE data
  IMPLICIT NONE

  LOGICAL :: first=.true.
  SAVE :: first

  INTEGER :: na, nb, nr, nbranches, i, ntot, list(3,ntot), naa, nbb, natoms

  TYPE(dynmat), INTENT(INOUT) :: dyn
  TYPE(lattice), INTENT(IN) :: latt
  TYPE(shell), INTENT(IN)  :: sh

  REAL(DP), INTENT(IN) :: q(3)
  REAL(DP) :: temp(3), x(3), fac
  REAL(DP), PARAMETER :: eps = 1.d-6

   !=====================================================================
  ! On input, the forces are in eV/A and the distances in A.
  ! Conversion factor to THZ**2
  !=====================================================================

  fac = evtoj   &           !  D^2 U eV  --> Joule
       /1.d-20  &        !  D^2 R A^2 --> meter 
       /amtokg  &        !  proton mass --> Kg
       /1.d24/twopi**2   !  cicles in THZ^2

  nbranches = 3*latt%natoms

  if(first) then
     allocate( dyn%cmat(nbranches,nbranches), dyn%eig(nbranches) )
     first = .false.
  endif

  dyn%cmat = 0
  dyn%eig = 0
 
  do i = 1, ntot
     na = list(1,i)
     nb = list(2,i)
     nr = list(3,i)

!     if( dyn%weight(na,nb,nr) > 0 ) then

     x = sh%r(:,nr) + latt%x(:,nb) - latt%x(:,na)

        dyn%cmat( 1 + (na-1)*3 : 3 + (na-1)*3, 1 + (nb-1)*3 : 3 + (nb-1)*3 ) = &
             dyn%cmat( 1 + (na-1)*3 : 3 + (na-1)*3, 1 + (nb-1)*3 : 3 + (nb-1)*3 ) &
             + ( dyn%dmat(:,:,na,nb,nr) * exp(  -im * dot_product( q, x ) ) ) / &
             sqrt( latt%mass( latt%ityp(na) ) * latt%mass(latt%ityp(nb)) )

!     endif

  enddo

  dyn%cmat = dyn%cmat * fac

  !----------------------------
  ! Check Hermiticity
  !----------------------------
  temp = 0
  do na=1, nbranches
     do nb=1, nbranches
        if( nb /= na ) then
           temp(1) = temp(1) + abs( dyn%cmat(na,nb) - conjg(dyn%cmat(nb,na) ) )
        endif
        if( nb == na ) temp(2) = temp(2) + abs( dyn%cmat(na,na) - conjg(dyn%cmat(na,na) ) )
        temp(3) = temp(3) +  abs( dyn%cmat(na,nb) )
     enddo
  enddo
  if( abs(temp(1)) > eps .or. abs(temp(2)) > eps ) &
       print'(''Dynamical matrix not hermitean '',3e9.1,2x,2e9.1)', &
       temp(1),temp(2),temp(3),temp(1)/temp(3),temp(2)/temp(3)


END SUBROUTINE dyn_mat
