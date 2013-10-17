!-----------------------------------------------
! builds a supercell, dimensions set by NDIM
!-----------------------------------------------
!===========================================================================
SUBROUTINE build_supercell( latt )
!===========================================================================
  USE nrtype
  USE data
  IMPLICIT NONE

  REAL(DP) :: e1(3), e2(3), e3(3)
  REAL(DP), ALLOCATABLE :: x(:,:)

  INTEGER :: i1, i2, i3, na

  TYPE(lattice) :: latt

  allocate( x(3,latt%natoms) )

  e1 = 0; e2 = 0; e3 = 0; e1(1) = 1; e2(2) = 1; e3(3) = 1

  write(*,'(//''*********************************'')')
  write(*,'(  ''*     Building a super cell     *'')')
  write(*,'(''*********************************'')')
  open(1,file='SPOSCAR',status='unknown')
  write(1,'(''super cell'')')
  write(1,'(f20.10)') latt%scale 
  write(1,'(3f20.15)') latt%at(:,1) * latt%ndim(1)
  write(1,'(3f20.15)') latt%at(:,2) * latt%ndim(2)
  write(1,'(3f20.15)') latt%at(:,3) * latt%ndim(3)
  write(1,'(200i5)') ( latt%nions(na) * product( latt%ndim ), na = 1,latt%ntypes )
  write(1,*) 'Direct'

  write(*,'(/''New lattice vectors:'')')
  write(*,'(3f20.15)') latt%at(:,1) * latt%ndim(1)
  write(*,'(3f20.15)') latt%at(:,2) * latt%ndim(2)
  write(*,'(3f20.15)') latt%at(:,3) * latt%ndim(3)

  do na = 1, latt%natoms
     x(1,na) = latt%x(1,na) / latt%ndim(1)
     x(2,na) = latt%x(2,na) / latt%ndim(2)
     x(3,na) = latt%x(3,na) / latt%ndim(3)
     do i3 = 0, latt%ndim(3)-1
        do i2 = 0, latt%ndim(2)-1
           do i1 = 0, latt%ndim(1)-1
              write(1,'(3f20.15)') &
                   x(:,na) + &
                   i1 * e1 / latt%ndim(1)  +  &
                   i2 * e2 / latt%ndim(2)  +  &
                   i3 * e3 / latt%ndim(3)
           enddo
        enddo
     enddo
  enddo

  close(1)
  deallocate( x )


END SUBROUTINE build_supercell
