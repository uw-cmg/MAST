!============================================================================
SUBROUTINE rshell( sh, latt )
  !============================================================================
  !
  ! construct the lattice vectors (in A)
  !
  USE nrtype
  USE data
  IMPLICIT NONE

  INTEGER :: i1, i2, i3, imax1, imax2, imax3, nr, indsw, nsh
  REAL(DP) :: temp(3), eps = 1.d-6
  REAL(DP), POINTER :: rr(:)

  INTEGER, ALLOCATABLE :: nl(:), ishel(:)  

  TYPE(shell), INTENT(INOUT) :: sh
  TYPE(lattice), INTENT(IN)  :: latt

  imax1 = nint( sh%rmax * vabs( latt%bg(:,1) ) / latt%scale ) 
  imax2 = nint( sh%rmax * vabs( latt%bg(:,2) ) / latt%scale ) 
  imax3 = nint( sh%rmax * vabs( latt%bg(:,3) ) / latt%scale ) 

  sh%nrm = 0                    ! number of lattice points
  do i1 = -imax1, imax1
     do i2 = -imax2, imax2
        do i3 = -imax3, imax3
           temp = i1*latt%at(:,1) + i2*latt%at(:,2) + i3*latt%at(:,3)
           temp = temp*latt%scale
           if( vabs( temp ) <= sh%rmax )  sh%nrm = sh%nrm + 1
        enddo
     enddo
  enddo

  allocate( sh%r(3,sh%nrm), sh%rr(sh%nrm), nl(sh%nrm) )
  sh%r = 0
  sh%rr = 0

  nr = 0
  do i1 = -imax1, imax1
     do i2 = -imax2, imax2
        do i3 = -imax3, imax3
           temp = i1*latt%at(:,1) + i2*latt%at(:,2) + i3*latt%at(:,3)
           temp = temp * latt%scale
           if( vabs( temp ) <= sh%rmax ) then
              nr = nr + 1
              sh%r(:,nr) = temp
              sh%rr(nr) = dot_product( temp, temp )
           endif
        enddo
     enddo
  enddo

  rr => sh%rr

  nl=0
  print*,'sorting...', nr,' R-vectors'
  call hpsort( sh%nrm, rr, nl )
  !
  !   reorder also the R vectors 
  !
  do nr = 1, sh%nrm-1
20   indsw = nl(nr)
     if(indsw.ne.nr) then 
        call swap( sh%r(:,indsw), sh%r(:,nl(indsw)) )
        call swap( nl(nr), nl(indsw) )
        go to 20
     endif
  enddo
  !
  ! calculate nsh
  !
  nsh = 1
  do nr = 2, sh%nrm
     if( abs( sh%rr(nr) - sh%rr(nr-1) ) >= eps ) nsh = nsh + 1
  enddo

  print*,'number of shells: ',nsh
  allocate( ishel( nsh ), sh%rl( nsh ) )
  !
  !     calculate the number of R for each shel
  !
  ishel = 1
  nsh = 1
  sh%rl(1) = 0
  do nr = 2, sh%nrm
     if( abs( sh%rr(nr) - sh%rr(nr-1) ) < eps )then
        ishel(nsh) = ishel(nsh) + 1
     else
        nsh = nsh + 1
        sh%rl(nsh) = sqrt( sh%rr(nr) )
     endif
  enddo

  !      print*,ishel,sh%rl
  !      print*,sh%r
  deallocate( nl, ishel )

END SUBROUTINE rshell



