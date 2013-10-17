!==================================================================
SUBROUTINE set_lattice ( latt, symm, iprint, rmax )
  !==================================================================
  USE nrtype
  USE nrutil
  USE data

  IMPLICIT NONE

  INTEGER :: na, i, j, k, iprint

  REAL(DP) :: temp(3), temp1(3), rmax, eps

  TYPE (lattice) :: latt
  TYPE (symmetry):: symm

!  eps = symm%symprec
  eps = 0.11

  open(1,file='POSCAR',status='old') !  lattice vectors and positions
  read(1,*) latt%string2
  read(1,*) latt%scale                    ! from file POSCAR

  !======================================================================    
  ! lattice vectors are the columns of at
  !======================================================================

  read(1,*) latt%ats(:,1)
  read(1,*) latt%ats(:,2)
  read(1,*) latt%ats(:,3)

  !======================================================================
  ! V = latt%scale * |a1 * a2 x a3|
  !======================================================================

  call vecprod( latt%ats(:,2), latt%ats(:,3), temp )
  latt%volumes = abs( dot_product( latt%ats(:,1), temp ) )

  !======================================================================
  ! if latt%scale < 0 in the POSCAR file then V = latt%scale
  !======================================================================

  if( latt%scale < 0 ) latt%scale = exp(log(-latt%scale / latt%volumes)/3)
  latt%volumes = latt%volumes * latt%scale**3

  print'(/'' Cell Volume = '',f12.5/)',latt%volumes

  call inve(latt%ats,latt%bgs) 

  !=====================================================================
  ! the columns of latt%bgs contain the reciprocal vectors
  !=====================================================================      

  latt%bgs = transpose( latt%bgs )

  print*,'Lattice vectors:'
  print'(3f15.10)',latt%ats(:,1)
  print'(3f15.10)',latt%ats(:,2)
  print'(3f15.10)',latt%ats(:,3)
  print*,'Reciprocal vectors:'
  print'(3f15.10)',latt%bgs(:,1)
  print'(3f15.10)',latt%bgs(:,2)
  print'(3f15.10)',latt%bgs(:,3)
  read(1,*) ( latt%nions(i), i = 1, latt%ntypes )
  latt%natomss = sum( latt%nions )
  allocate ( latt%xs(3,latt%natomss), latt%ityps(latt%natomss) )
  na = 0
  do i = 1, latt%ntypes
     do j = 1, latt%nions(i)
        na = na + 1
        latt%ityps( na ) = i
     enddo
  enddo

  print'(/''latt%natomss = '',i10)', latt%natomss
  print'(/''latt%nions = '',i10)', latt%nions
  print'(/''latt%ntypes = '',i10)', latt%ntypes
  read(1,*) latt%string
  do na=1,latt%natomss
     read(1,*) latt%xs(:,na)
  enddo
  close(1)

  print'(/''Volume/atom = '',f10.5)', latt%volumes/latt%natomss

  !==============================================================================
  ! positions are needed in lattice coordinates to find the symmetry operations.
  !==============================================================================

  if( latt%string=='C' .or. latt%string=='c' ) &
       &     latt%xs = matmul( transpose(latt%bgs), latt%xs )


!--------------------
! Try to guess Rmax
! rmax = 0
  do i = -1, 1
     do j = -1, 1
        do k = -1, 1
           temp =  latt%scale*(i*latt%ats(:,1) + j*latt%ats(:,2) + k*latt%ats(:,3))
           if( rmax < vabs(temp) ) rmax = vabs(temp)
        enddo
     enddo
  enddo

  rmax = rmax + eps

  print'(''Rmax = '',f10.4)',rmax
!--------------------------------------------------------------------------------------

  call sgama ( symm, latt, iprint )


END SUBROUTINE set_lattice
