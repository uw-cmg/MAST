!=======================================================================
SUBROUTINE smooth( dos, sigma, ndos, emin, deltae )
  !=======================================================================
  !
  ! convolutes a Gaussian to the density of states
  !
  USE nrtype

  IMPLICIT NONE
  INTEGER :: i, j, ndos
  REAL(DP) :: dos(ndos), sigma, etempp, etempm, emin, deltae, peso, norma
  REAL(DP), ALLOCATABLE :: dostemp(:), table(:)

  allocate( dostemp(ndos), table(0:ndos-1) )

!  if( sigma < deltae/10 ) then
!     print*,'No smearing'
!     return
!  endif

  dostemp = 0
  do i = 0, ndos-1
     etempm = emin + deltae * i + deltae/2 
     table(i) = exp( - etempm**2 / sigma**2 )
  enddo
  do i = 1, ndos
     etempm = emin + deltae * ( i - 1 ) + deltae/2
     peso = 0
     do j = 1, ndos
        etempp = emin + deltae * ( j - 1 ) + deltae/2
        dostemp(i) = dostemp(i) + dos(j) * table(abs(i-j))
        peso = peso + table(abs(i-j))
     enddo
     dos(i) = dostemp(i) / peso
  enddo

!  dos = dostemp / sqrt( pi * sigma**2 ) * deltae


  deallocate( dostemp )

END SUBROUTINE smooth
