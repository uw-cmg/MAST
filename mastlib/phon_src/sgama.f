!=========================================================================
SUBROUTINE sgama ( symm, latt, iprint )
  !=========================================================================
  !
  ! check if the symmetries of the bravais lattice are still  symmetries of the
  ! lattice with the basis
  !
  USE nrtype
  USE nrutil
  USE data

  IMPLICIT NONE

  LOGICAL :: found, equiva, first=.true.
  LOGICAL, ALLOCATABLE :: lrot(:), traslated(:)

  INTEGER :: nrot, natoms, nsym, trasl, i, j, isym, iprint, nat
  INTEGER :: nrot_cub = 0, nrot_hex = 0, ntau, nsuper(3)

  REAL(DP) :: rotx(3), ft(3)
  REAL(DP) :: at(3,3), bg(3,3)

  REAL(DP) :: eps

  TYPE (lattice) :: latt
  TYPE (symmetry):: symm

  REAL(DP), POINTER :: is(:,:,:), ftau(:,:), x(:,:)
  INTEGER, POINTER :: ityp(:)

  REAL(DP), ALLOCATABLE :: ftsuper(:,:)

! set symmetry precision
  eps = symm%symprec

  !==============================================================================
  ! find if the lattice is a cubic lattice or an hexagonal one
  !==============================================================================

  if(first) then 
     allocate( symm%is(3,3,48) )
  endif
  symm%is = 0
  nrot_cub = 48
  call cubicsym( latt%ats, symm%is, nrot_cub, eps ) 
  symm%is = 0
  nrot_hex = 24
  call hexsym( latt%ats, symm%is, nrot_hex, eps ) 

  if(iprint>0)print'(/''Rotations for cubic symmetry    '',i3,/ &
       &        ''Rotations for exagonal symmetry '',i3)', nrot_cub, nrot_hex

  nrot = max ( nrot_cub, nrot_hex )

  deallocate(symm%is)

  symm%nsym = nrot

  allocate( symm%is(3,3,nrot), symm%iscryst(3,3,nrot), symm%isstart(3,3,nrot), symm%ftau(3,nrot) )

  natoms = latt%natomss

  ityp => latt%ityps
  is => symm%is
  ftau => symm%ftau
  x => latt%xs
  at = latt%ats
  bg = latt%bgs

  !==============================================================================
  ! is contains the point group operations
  !==============================================================================
  if ( nrot == nrot_cub ) then
     if(iprint>0)print'(/''Found Cubic symmetry'')'
     call cubicsym( at, is, nrot, eps ) 
  else
     if(iprint>0)print'(/''Found Hexagonal symmetry'')'
     call hexsym( at, is, nrot, eps ) 
  endif

  !=========================================================================
  ! put symm%is is the symmetry matrix in crystal coordinates
  !=========================================================================

  !-----------------
  ! used in sgama1
  !-----------------
  symm%iscryst = symm%is


  do i=1,nrot
     symm%is(:,:,i) = matmul(transpose(latt%bgs),matmul(symm%is(:,:,i),latt%ats))
  enddo

  symm%nrot = nrot

  allocate( lrot(nrot), traslated(natoms) )
  nsym = 0; ntau = 0; trasl = 0; ftau = 0; lrot=.false.

  !==============================================================================
  ! search if the identity has fractionary traslations (the cell is a super cell)
  !==============================================================================
  trasl = 1
  do j = 2, natoms
     !--------------------------
     ! construct the traslation
     !--------------------------
     ft = x(:,1) - x(:,j) - nint( x(:,1) - x(:,j) )
     call checksym( is(:,:,1), ft, x, ityp, natoms, found, eps )
     if ( found ) then
        trasl = trasl + 1
     endif
  enddo

  if( trasl > 0 .and. iprint > 0 ) print'(/''Found'',i5,'' additional traslations'')', trasl

  !**********************************************************************
  !=====================================================================
  !             FIRST CALL ONLY
  !=====================================================================
  !
  if(first) then
     !----------------------------------------------------
     ! do it again, now fill ftsuper with the traslations
     !----------------------------------------------------
     allocate( ftsuper(3,trasl) )
     ftsuper = 0
     trasl = 1
     nsuper = 1
     do j = 2, natoms
        ft = x(:,1) - x(:,j) - nint( x(:,1) - x(:,j) )
        call checksym( is(:,:,1), ft, x, ityp, natoms, found, eps )
        if ( found ) then
           trasl = trasl + 1
           ftsuper(:,trasl) = ft
        endif
     enddo


     found = .false.
     do i = 1, 3
        do ntau = 2, trasl
           rotx = matmul( at, ftsuper(:,ntau) ) 
           call vecprod( rotx, latt%ats(:,i), ft )
           if( vabs( ft ) < eps ) then
              nsuper(i) = max( nsuper(i), nint( abs( 1._dp/ftsuper(i,ntau) ) ) )
              found = .true.
           endif
        enddo
     enddo
     !-------------------------------------------------------------------------
     ! if the traslation is not parallel to any lattice vector this is not a 
     ! supercell
     !-------------------------------------------------------------------------
     if(.not.found) trasl=1
     if( found ) print'(/''The cell is a supercell '',3i4)', nsuper

     !==============================
     ! Look for the primitive cell
     !==============================
     traslated = .false.
     do i = 1, natoms - 1
        do j = i+1, natoms
           do ntau = 2, trasl
              ft = x(:,j) - ftsuper(:,ntau)
              if( equiva( x(:,i), ft, eps ) ) traslated(j) = .true.
           enddo
        enddo
     enddo

     !--------------------------------------------------------------
     ! lattice vectors and reciprocal vectors of the primitive cell
     !--------------------------------------------------------------
     do i = 1, 3
        latt%at(:,i) = latt%ats(:,i) / nsuper(i)
        latt%bg(:,i) = latt%bgs(:,i) * nsuper(i)
     enddo

     !-----------------------------------------------------
     ! number of atoms and volume of the primitive cell
     !-----------------------------------------------------
     latt%natoms = latt%natomss / product(nsuper)
     latt%nions = latt%nions / product(nsuper)

     allocate( latt%ityp( latt%natoms ) )
     nat = 0
     do i = 1, latt%ntypes
        do j = 1, latt%nions(i)
           nat = nat + 1
           latt%ityp( nat ) = i
        enddo
     enddo

     latt%volume = latt%volumes / product(nsuper)

     if(iprint>0)print'(/''Primitive lattice vectors:'')'
     if(iprint>0)print'(3f20.15)',latt%at
     if(iprint>0)print'(/''Primitive reciprocal vectors:'')'
     if(iprint>0)print'(3f20.15)',latt%bg

     !-------------------------------
     ! atoms of the primitive cell
     !-------------------------------
     allocate( latt%x(3,latt%natoms) )
     nat = 0
     do i = 1, natoms
        if(.not.traslated(i)) then
           nat = nat + 1
           latt%x(:,nat) = x(:,i) * nsuper
        endif
     enddo

  endif
  !======================================================================
  !         END OF FIRST CALL 
  !======================================================================
  !**********************************************************************

  ntau = 0
  x => latt%x
  natoms = latt%natoms
  ityp => latt%ityp
  !==============================================
  ! Find the symmetry operations of the lattice
  !==============================================

  do isym = 1, nrot

     !----------------------------------------------------------------------
     ! search using all possible fractionary traslations (the first is zero)
     !----------------------------------------------------------------------
     traslation: do i=1,natoms

        !---------------
        ! rotate atom i
        !---------------
        rotx = matmul( symm%is(:,:,isym), x(:,i) ) 

        do j = 1, natoms

           !---------------------------------------------------------------------
           ! construct the traslation, mind, from the rotated to the non-rotated
           !---------------------------------------------------------------------
           ft = rotx - x(:,j) - nint( rotx - x(:,j) )

           !---------------------------------------------------------------------
           ! check if every atom of cell is sent in some other atom under S + ft
           !---------------------------------------------------------------------
           call checksym( symm%is(:,:,isym), ft, x, ityp, natoms, found, eps )

           !-----------------------------------------------------
           ! is it a symmetry operation? ok, go to the next rotation
           !-----------------------------------------------------
           if ( found ) exit traslation

        enddo
     enddo traslation

     !-----------------------------------------------
     ! found the symmetry, set ftau and print out 
     !-----------------------------------------------
     if( found ) then
        lrot(isym) = .true.
        ftau(:,isym) = ft
        if( vabs(ft) > eps ) ntau = ntau + 1
        nsym = nsym + 1
        if( iprint > 0 ) then
           print*,' '

           !-------------------------------------------------------------
           ! if this symmetry operation has no fractional traslation....
           !-------------------------------------------------------------
!           if( all ( abs(ft) < eps ) ) then
!              print'(''symm '',i3)', isym
!              if( iprint > 1 ) then
!                 print'(''Cartesian coordinates:'')'
!                 print'(3f15.10)', matmul( at, matmul( is(:,:,isym), transpose(bg))) 
!              endif
!              if( iprint > 2 ) then
!                 print*,' '
!                 print'(''Lattice coordinates:'')'
!                 print'(3f15.10)', is(:,:,isym)
!              endif

              !-------------------
              ! otherwise.......
              !-------------------
!           else
              print'(''symm '',i3)', isym
              if( iprint > 1 ) then
                 print'(''Cartesian coordinates:'')'
                 print'(''frac '',3f15.10)', latt%scale*matmul(at,ft)
                 print'(3f15.10)', matmul( at, matmul( is(:,:,isym), transpose(bg)))
              endif
              if( iprint > 2 ) then
                 print*,' '
                 print'(''Lattice coordinates'')'
                 print'(''frac '',3f15.10)', ft
                 print'(3f15.10)', is(:,:,isym)
              endif
!           endif
        endif
     endif

  enddo

  print'(/''found    '',i5,'' symmetry operations'')',nsym
  print'(''of which '',i5,'' have fractionary traslation''/)', ntau
  symm%nsym = nsym

  !-------------------------------
  ! Is there inversion symmetry ?
  !-------------------------------
  symm%invsym = .false.
  do isym = 2, nsym
     symm%invsym = ALL ( is(:,:,isym) == -is(:,:,1) ) .and. vabs( ftau(:,isym) ) < eps   
!     if( is(1,1,isym) == -1 .and. is(1,2,isym) == 0 .and. is(1,3,isym) == 0 .and. &
!          is(2,1,isym) == 0 .and. is(2,2,isym) == -1 .and. is(2,3,isym) == 0 .and. &
!          is(3,1,isym) == 0 .and. is(3,2,isym) == 0 .and. is(3,3,isym) == -1 .and. &
!          vabs( ftau(:,isym) ) < eps  ) then
!        symm%invsym = .true.
     if(symm%invsym) then 
       print'(''Found inversion symmetry ''/)'
       exit
     endif
  enddo
  



  !=============================
  ! reorder simmetry operations
  !=============================
  do isym = 1, nrot - 1
     if( .not.lrot(isym) )then
        i = 1
        do while ( .not.lrot(isym+i) .and. (isym+i+1)<=nrot )
           i = i + 1
        enddo
        call swap ( is(:,:,isym), is(:,:,isym+i) )
        call swap ( ftau(:,isym), ftau(:,isym+i) )
        call swap ( lrot(isym), lrot(isym+i) )
     endif
  enddo
  deallocate(lrot)

  !=============================
  ! set the map of the atoms
  !=============================
  if(first) then
     allocate( symm%irt(latt%natoms,nsym), symm%irts(latt%natomss,nsym) )

     symm%irt = -1
     symm%irts = -1
     
     do j = 1, latt%natoms
        do isym =  1, nsym
           !-----------------------------------------------------------
           ! rototraslate the atom (apply the symmetry operation isym)
           !------------------------------- ---------------------------
           rotx = matmul( is(:,:,isym), x(:,j) ) - ftau(:,isym)
           do i = 1, latt%natoms
              if( equiva( rotx, x(:,i), eps ) ) then
                 symm%irt(j,isym) = i
              endif
           enddo
           if( symm%irt(j,isym) < 0 ) stop 'Errore in sgama'
        enddo
     enddo
     
     do j = 1, latt%natomss
        do isym =  1, nsym
           !-----------------------------------------------------------
           ! rototraslate the atom (apply the symmetry operation isym)
           !------------------------------- ---------------------------
           rotx = matmul( is(:,:,isym), latt%xs(:,j) ) - ftau(:,isym) / nsuper
           do i = 1, latt%natomss
              if( equiva( rotx, latt%xs(:,i), eps ) ) then
                 symm%irts(j,isym) = i
              endif
           enddo
           if( symm%irts(j,isym) < 0 ) stop 'Errore in sgama'
        enddo
     enddo

  endif

  !--------------------------------------------------
  ! put the symmetry matrix in cartesian coordinates
  !--------------------------------------------------
  do i=1,symm%nsym
     symm%is(:,:,i) =  matmul( latt%at, matmul(symm%is(:,:,i), transpose(latt%bg) ) )
     symm%ftau(:,i) = latt%scale * matmul(latt%at,symm%ftau(:,i)) / latt%ndim
  enddo

  if( first ) then
     symm%isstart = 0
     symm%isstart = symm%is
     symm%nsymstart = symm%nsym
     first = .false.
  endif

  !---------------------
  ! for debug purpouses
  !---------------------
  if(iprint<0) write(*,*) 'exiting from sgama'

  RETURN
END SUBROUTINE sgama

!===================================================================
SUBROUTINE checksym ( is, ft, x, ityp, natoms, found, eps1 )
  !===================================================================
  !
  ! check if the roto-traslation is+ft is a symmetry operation, i.e. if the
  ! set of positions x is invariant
  !
  USE nrtype

  IMPLICIT NONE

  INTEGER :: i, j, natoms
  INTEGER :: ityp(natoms)

  REAL(DP) :: is(3,3), x(3,natoms), ft(3), rotx(3), eps(3), eps1

  LOGICAL :: found

  eps = eps1

  do i = 1, natoms
     rotx = matmul( is, x(:,i) ) - ft
     do j = 1, natoms
        found = all( abs ( rotx - x(:,j) - nint( rotx - x(:,j) ) ) < eps )
        found = found .and. ( ityp(j) == ityp(i) )
        if ( found ) exit      ! if found an equivalent atom proceed
     enddo
     if ( .not. found ) exit   ! otherwise this is not a symmetry operation
  enddo

END SUBROUTINE checksym
