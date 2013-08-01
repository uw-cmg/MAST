!=========================================================================
SUBROUTINE sgama1 ( symm, latt, iprint )
!=========================================================================
!
! check if the symmetries of the bravais lattice are still  symmetries of the
! lattice with the basis, no fractionary traslations are allowed (only point 
! group )
!
  USE nrtype
  USE nrutil
  USE data
  
  IMPLICIT NONE
      
  LOGICAL :: found
  LOGICAL, ALLOCATABLE :: lrot(:)
  
  INTEGER :: nrot = 0, natoms, nsym, i, isym, iprint
  INTEGER, POINTER :: ityp(:)
  
  REAL(DP) :: ft(3), eps
  REAL(DP) :: at(3,3), bg(3,3)
  
  REAL(DP), POINTER :: x(:,:)
  
  TYPE (lattice) :: latt
  TYPE (symmetry):: symm
  

!=========================================================================
! now symm%is is the symmetry matrix in crystal coordinates
!=========================================================================
      
  eps = symm%symprec
  nrot = symm%nrot
  
  at = latt%at
  bg = latt%bg
  
  allocate( lrot(nrot) )
  nsym = 0
  lrot=.false.
  
  x => latt%x
  natoms = latt%natoms
  symm%is = symm%iscryst
  
  do i=1,nrot
     symm%is(:,:,i) =  matmul(transpose(latt%bgs),matmul(symm%is(:,:,i),latt%ats))
  enddo

!==============================================
! Find the symmetry operations of the lattice
!==============================================

  ft = 0
  do isym = 1, nrot

!---------------------------------------------------------------------
! check if every atom of cell is sent in some other atom under S + ft
!---------------------------------------------------------------------
     ityp => latt%ityp
     call checksym( symm%is(:,:,isym), ft, x, ityp, natoms, found, eps )

!-----------------------------------------------
! found the symmetry and print out 
!-----------------------------------------------
     if( found ) then
        lrot(isym) = .true.
        nsym = nsym + 1
        
        if( iprint > 0 ) then
           print*,' '

           print'(''symm'',i3)', isym
           if( iprint > 1 ) print'(3f15.10)', matmul( at, matmul( symm%is(:,:,isym), transpose(bg))) 
           if( iprint > 2 ) then
              print*,' '
              print'(3f15.10)', symm%is(:,:,isym)
           endif
           
        endif
     endif
     
  enddo
                  
  if(iprint>0) print'(''found    '',i5,'' symmetry operations'')',nsym
  symm%nsym = nsym

!=============================
! reorder simmetry operations
!=============================
  do isym = 1, nrot - 1
     if( .not.lrot(isym) )then
        i = 1
        do while ( .not.lrot(isym+i) .and. (isym+i+1)<=nrot )
           i = i + 1
        enddo
        call swap ( symm%is(:,:,isym), symm%is(:,:,isym+i) )
        call swap ( lrot(isym), lrot(isym+i) )
     endif
  enddo
  deallocate(lrot)

  do i=1,symm%nsym
     symm%is(:,:,i) =  matmul( latt%at, matmul(symm%is(:,:,i), transpose(latt%bg) ) )
  enddo

END SUBROUTINE sgama1
