!========================================================================================
SUBROUTINE numeric_force( latt, symm, sh, dyn, iprint, nti, lcentral, ldrift, lnosym )
  !======================================================================================
  
  USE nrtype
  USE nrutil
  USE data

  IMPLICIT NONE

  LOGICAL :: found, equiv, equiva, equiva1, esiste, lcentral, ldrift, lnosym

  INTEGER :: ndispl, i, ir, iir, na, nb, isym, naa, nbb, nbbb, nr, nrr, wsedge, wstot, &
       j, k, natoms, conta, iprint, startdisp, idisp, nti

  INTEGER, ALLOCATABLE :: adis(:), adis1(:), nions(:), ndis(:), ws(:)

!  REAL(DP), PARAMETER :: eps=1.d-8, large_disp = 0.2d0
  REAL(DP), PARAMETER :: large_disp = 0.2d0

  REAL(DP), ALLOCATABLE :: dis(:,:), force(:,:,:), force2(:,:,:), &
       force1(:,:,:), x(:,:), xtmp(:,:), xlatt_tmp(:,:), xlatt_tmp1(:,:), &
       tmpmat1(:,:,:), rmax(:), xws(:,:,:)

  REAL(DP) :: temp(3), temp1(3), tmp1(3), dx(3,3), invis(3,3), eps, eps1, &
       ftot(3,3), at(3,3), bg(3,3), scale, volume, tmpmat(3,3), spring_size, disp_size

  TYPE(lattice)  :: latt
  TYPE(shell)    :: sh
  TYPE(dynmat)   :: dyn
  TYPE(symmetry) :: symm

  allocate( nions(latt%ntypes), latt%super_atom(latt%natoms), &
       xlatt_tmp(3,latt%natoms), xlatt_tmp1(3,latt%natoms), rmax(latt%natoms) )

  eps = symm%symprec  / 100
  eps1 = symm%symprec 

  latt%super_atom = 0
  rmax = 0

  at = latt%ats
  call inve( at, bg )
  bg = transpose(bg)

  scale = latt%scale
  natoms = latt%natomss

  allocate ( x(3,natoms), xtmp(3,natoms), ws(natoms), xws(3,8,natoms) )

  x = latt%xs

  allocate( tmpmat1(3,3,natoms), dyn%tmpmat2(3,3,natoms,natoms) )

  call vecprod( at(:,2), at(:,3), tmp1 )
  volume = abs( dot_product( at(:,1), tmp1 ) )
  if( scale < 0 ) scale = exp(log(-scale / volume)/3 )


  !=================================================================
  ! open force file
  !=================================================================
  inquire(file='FORCES',exist=esiste)
  if(.not.esiste) stop 'cannot find file FORCES'
  open(1,file='FORCES',status='old')

  !-----------------------------------------------------------------
  ! How many displacements?
  !-----------------------------------------------------------------
  read(1,*) ndispl

  !------------------------------------
  ! allocate work space
  !------------------------------------
  allocate( adis(ndispl), adis1(ndispl), dis(3,ndispl), force(3,natoms,ndispl), &
       force1(3,natoms,3), force2(3,natoms,3) )

  !-----------------------------------------
  ! allocate space for the dynamical matrix
  !-----------------------------------------
  allocate( dyn%dmat(3,3,latt%natoms,latt%natoms,sh%nrm), &
       dyn%weight(latt%natoms,latt%natoms,sh%nrm) )

  !------------------------------------------------------------------
  ! readin from file FORCES the displacements and the induced forces
  !------------------------------------------------------------------
  force = 0
  !allocate( ndis(latt%natoms) ) !TTM read in more displacements
  allocate( ndis(latt%natomss) )


  !-------------------------------
  ! atoms in carthesian coordinates
  !-------------------------------
  x = scale * matmul( at, x )
  latt%x = latt%scale * matmul( latt%at, latt%x )

  do na = 1, latt%natoms
     do nb = 1, natoms
        if( equiva1( x(:,nb), latt%x(:,na), eps1 ) ) then
!-------------------------------------------------------------------
! set the corrispondence  cell ---> supercell
! the first atom of the supercell is the one we are interested in
!--------------------------------------------------------------------
           latt%super_atom(na) = nb
           exit
        endif
     enddo
  enddo

  ndis = 0
  idisp = 0
  do i=1,ndispl
     read(1,*) adis(i), dis(:,i)
     adis1(i) = adis(i)
     idisp = idisp + 1
     if( adis1(i) > adis1(i-1) .and. i > 1 ) idisp = 1

     !========================================================================
     ! Establish the corrispondence between the atom in the supercell
     ! and the one in the primitive cell
     !========================================================================
     do na = 1, latt%natoms
        if( equiva1( x(:,adis(i)), latt%x(:,na), eps1 ) ) then
           latt%super_atom(na) = adis(i)
           adis(i) = na
           exit
        endif
     enddo

     !-----------------------------------------------
     ! count how many displacements for each atom
     !-----------------------------------------------
     if(lcentral.and.i<=ndispl/2) ndis( adis(i) ) = ndis( adis(i) ) + 1
     if(.not.lcentral) ndis( adis(i) ) = ndis( adis(i) ) + 1
     disp_size = sqrt(dot_product(scale * matmul( at, dis(:,i) ),scale * matmul( at, dis(:,i) ))) 
     print'(''Atom '',i3,'', Displacement'',3f9.5,'', Size = '',f7.5'' (A)'')',adis(i),dis(:,i),&
          disp_size
     
     if( disp_size > large_disp ) then
        print'(/''  WARNING  size of displacement > '',f4.2,'' A'',/)', large_disp
     endif

     !------------------------------------------------
     ! put the displacements in carthesian axis
     !------------------------------------------------
     dis(:,i) = scale * matmul( at, dis(:,i) )

     !-------------------------------
     ! check the sum of the forces
     !-------------------------------
     temp = 0
     do na = 1, natoms
        read(1,*) force(:,na,i)
        temp = temp + force(:,na,i)
     enddo
     !-------------------------------
     ! added to check ldrift tag
     !-------------------------------
     if(ldrift) then
        if( iprint > 0 ) then
           print'(''Drift in forces: '',3f15.8)', temp
           print'(''Removing drift ...'')' 
        endif
        do na = 1, natoms
           force(:,na,i) = force(:,na,i) - temp / natoms 
        enddo
     endif
!
  enddo

!--------------------------------
! central differences
!--------------------------------
  if(lcentral) then
     ndispl = ndispl / 2
     do i = 1, ndispl
        force(:,:,i) = ( force(:,:,i) - force(:,:,ndispl + i) ) / 2
     enddo
  endif

  !---------------------------------------
  ! initialize the force constant matrix 
  !---------------------------------------  
  dyn%dmat = 0
  dyn%weight = 0

  !==================================================
  !==================================================
  !
  ! Now, for each atom in the primitive cell:
  !
  !==================================================
  !==================================================
  startdisp = 1

  !-----------------------------------------------------
  ! superlattice back to reciprocal space, and keep it
  !-----------------------------------------------------
  xlatt_tmp = matmul( transpose(latt%bg), latt%x ) / scale
  xlatt_tmp1 = matmul( transpose(latt%bg), latt%x ) / scale

  wstot = 0
  wsedge = 0

  do naa = 1, latt%natoms

     x = latt%xs
     latt%x = xlatt_tmp1
     xlatt_tmp = xlatt_tmp1
     !----------------------------------------------------------------------------
     ! set the symmetry of the crystal, fractionary traslations included,
     ! this is needed when for the non-moving atoms I look for the corresponding
     ! moving one
     !----------------------------------------------------------------------------
     symm%is = symm%isstart

     !-----------------------
     ! is the atom moving ?
     !-----------------------
     write(*,100)
     if( ndis( naa ) > 0 )then
        print'(''Atom #'',i4,'' MOVING'')', naa

        !--------------------------------------------------------
        ! xtmp contains the WS cell centred onto the moving atom
        !--------------------------------------------------------
        call get_wig( latt, sh, natoms, scale, at, bg, x, xtmp, &
             xlatt_tmp, naa, rmax(naa), ws, xws, eps1 ) 
        wstot = wstot + sum(ws) 
        do na = 1, natoms
           if( ws(na) > 1 ) wsedge = wsedge + ws(na)
        enddo

        !-------------------------------------------------------------------------
        ! find the point group symmetries of the crystal with naa in the origin
        !------------------------------------------------------------------------
        call sgama1 ( symm, latt, iprint )

        !------------------------------------------------
        ! construct the forces in carthesian coordinates
        !------------------------------------------------
        call set_forces( force, force1, force2, dis, &
             startdisp, symm, bg, scale, xtmp, ndis, naa, natoms, &
             ndispl, dx, iprint, lnosym )


        !---------------------------------------------------------------------
        ! construct the force constants for all the atoms in the WS supercell
        !---------------------------------------------------------------------
        dyn%tmpmat2(:,:,:,naa) = 0

        do na = 1, natoms

           !----------
           ! - F / dx
           !----------
           dyn%tmpmat2(:,1,na,naa) =  - force2(:,na,1) / dx(1,1)
           !----------
           ! - F / dy
           !----------
           dyn%tmpmat2(:,2,na,naa) =  - force2(:,na,2) / dx(2,2)
           !----------
           ! - F / dz
           !----------
           dyn%tmpmat2(:,3,na,naa) =  - force2(:,na,3) / dx(3,3)

        enddo

        !--------------------------------------------------------------------
        ! symmetrize the force constants: D(x) = sum_i S_i^-1 D(S_ix) S_i
        !--------------------------------------------------------------------
        tmpmat1 = 0
        do na = 1, natoms
           conta = 0
           do isym = 1, symm%nsym
              !------------------------------
              ! only point group operations
              !------------------------------
              call inve( symm%is(:,:,isym), invis )
              tmp1 = matmul( invis, xtmp(:,na) ) 
              tmp1 = matmul( transpose(bg), tmp1 ) / scale
              do nb = 1, natoms                       
                 temp = matmul( transpose(bg), xtmp(:,nb) ) / scale
                 !-----------------------------------------------------------------------------
                 ! isym sends (naa,na) in (naa,nb), so D(naa,na) = S_isym D(naa,nb) S_isym^-1
                 !-----------------------------------------------------------------------------
                 if( equiva( temp, tmp1, eps1 )  )then
                    tmpmat = dyn%tmpmat2(:,:,nb,naa)
                    tmpmat1(:,:,na) = tmpmat1(:,:,na) + &
                         matmul( matmul ( symm%is(:,:,isym), tmpmat ), invis )
                    conta = conta + 1
                 endif
              enddo
           enddo
           tmpmat1(:,:,na) = tmpmat1(:,:,na) / conta
           if( conta /= symm%nsym ) print*,'warning, not all symmetries ', &
                'have been found in symmetrization', conta, symm%nsym, na 
        enddo
        
        if( symm%lsymm ) dyn%tmpmat2(:,:,:,naa) = tmpmat1
        !===================================================================
        ! Fill the force constant matrix
        !===================================================================

        latt%x = matmul( latt%at, latt%x ) * latt%scale
        do nbb = 1, latt%natoms
           found = .false.
           do nrr = 1, sh%nrm
              !------------------------------------------------------
              ! generic vector of the crystal with naa in the origin
              !------------------------------------------------------
              temp = sh%r(:,nrr) + latt%x(:,nbb) - latt%x(:,naa)

              if( vabs(temp) <= sh%rmax )then
                 !-------------------------------------------------------
                 ! search for the corresponding vector in the WS cell
                 !-------------------------------------------------------
                 do j = 1, natoms
                    do k = 1, ws(j)                  
                       if( equiva1( temp, xws(:,k,j), eps1 ) ) then
                          found =.true.
                          dyn%dmat(:,:,naa,nbb,nrr) = dyn%tmpmat2(:,:,j,naa) / ws(j)
                          dyn%weight(naa,nbb,nrr) = ws(j)
                       endif
                    enddo
                 enddo
              endif
           enddo 
           if( .not. found ) stop 'Error, try to decrease SYMPREC or check your primitive cell'
        enddo

        ftot = 0
        do na = 1, natoms
           ftot = ftot + dyn%tmpmat2(:,:,na,naa) 
        enddo
        if( vabs(ftot) > 1.d-3 ) then
           print'(/''!!! WARNING The sum of the force constants is different from zero'')'
           print'(3f16.12)',ftot
        endif

     else

        print'(''Atom #'',i4,''  not moving '', 12x,3f10.6)', naa, xlatt_tmp(:,naa)

        !-----------------------------------------------------------------------------
        ! I am not moving this atom because there is a symmetry operation that sends
        ! it to an other atom of the basis: what is the symmetry operation and
        ! what is the target atom?
        !-----------------------------------------------------------------------------
        donb: do nbb = 1, naa-1

           found = .false.

           do isym = 1, symm%nsymstart
              if( naa == symm%irt(nbb,isym) ) then

                 !---------------------
                 ! found the target nb
                 !---------------------
                 found = .true.
                 write(*,'(''Because I already moved atom # '',i3, 2x,3f10.6)') nbb, xlatt_tmp(:,nbb)
                 write(*,'(''Symm #'',i3,'' translation: '',3f10.6)') &
                      isym, matmul( transpose(latt%bgs), symm%ftau(:,isym) ) / scale
                 !----------------------------------------------
                 ! I need the inverse of the symmetry operation
                 !----------------------------------------------
                 call inve( symm%isstart(:,:,isym), invis )

                 do na = 1, natoms
                    !-------------------------
                    ! S(nbb) = naa; S(na) = j
                    !-------------------------
                    j = symm%irts(na,isym)
                    dyn%tmpmat2(:,:,j,naa) = matmul( &
                         matmul( symm%isstart(:,:,isym), dyn%tmpmat2(:,:,na,nbb) ), invis )
                 enddo

                 !------------------------------------------------------
                 ! generate the WS cell (in xtmp) with nb in the center 
                 !------------------------------------------------------
                 call get_wig( latt, sh, natoms, scale, at, bg, x, &
                      xtmp, xlatt_tmp, nbb, rmax(naa), ws, xws, eps1 ) 
                 wstot = wstot + sum(ws) 
                 do na = 1, natoms
                    if( ws(na) > 1 ) wsedge = wsedge + ws(na)
                 enddo

                 !-------------------------------------------------------------
                 ! Now, choose one vector of the lattice with nbb in the origin
                 !-------------------------------------------------------------
                 latt%x = matmul( latt%at, latt%x ) * latt%scale
                 do nbbb = 1, latt%natoms
                    equiv = .false.
                    do nrr = 1, sh%nrm

                       temp = sh%r(:,nrr) + latt%x(:,nbbb) - latt%x(:,nbb)
                       
                       if( vabs(temp) <= sh%rmax ) then
                          temp = matmul(transpose(latt%bg), temp ) / latt%scale   ! to crystal
                          !---------------------------------------------------------
                          ! search the corresponding vector in the WS cell (of nbb)
                          !---------------------------------------------------------
                          do j = 1,natoms
                             do k = 1, ws(j)                  
                                temp1 = matmul(transpose(latt%bg), xws(:,k,j) ) / latt%scale   ! to crystal
                                if( equiva1( temp, temp1, eps1 ) ) then
!                                if( equiva1( temp, xws(:,k,j) ) ) then
                                   equiv = .true.
                                   !------------------------------------------------------------------------
                                   ! which is the vector in the non rotated crystal with naa in the center ?
                                   !------------------------------------------------------------------------
                                   donr: do nr = 1, sh%nrm
 
                                      do na = 1, latt%natoms
                                         tmp1 = sh%r(:,nr) + latt%scale * &
                                              matmul(latt%at,xlatt_tmp1(:,na) - xlatt_tmp1(:,naa))
                                         tmp1 = matmul( invis, tmp1 )
                                         tmp1 = matmul(transpose(latt%bg), tmp1 ) / latt%scale  ! to crystal

                                         if( equiva1( temp, tmp1, eps1 ) ) exit donr
                                      enddo
                                   enddo donr
                                   if( nr == sh%nrm + 1 ) &
                                        stop 'Error, nr not found 2, try to decrease SYMPREC or check your primitive cell'
                                   if( na == latt%natoms + 1 ) stop 'Error, na not found'
                                   dyn%dmat(:,:,naa,na,nr) =  matmul( symm%isstart(:,:,isym), matmul( &
                                        dyn%tmpmat2(:,:,j,nbb) / ws(j), invis ) ) 
                                   dyn%weight(naa,na,nr) = ws(j)
!                                   if(ws(j)==0) dyn%dmat(:,:,naa,na,nr) = 0
                                endif
                             enddo
                          enddo
                       endif
                    enddo
                 enddo
                 if( .not. equiv ) stop 'Error, ciccio vector not found'
              endif
              if(found) exit donb

           enddo
        enddo donb
        if( .not. found ) stop 'Error, equivalent atom not found'

        ftot = 0
        do na = 1, natoms
           ftot = ftot + dyn%tmpmat2(:,:,na,naa) 
        enddo

        if( vabs(ftot) > 1.d-3 ) then
           print'(/''!!! BAD NEWS !!! The sum of the force constants is different from zero'')'
           print'(3f16.12)',ftot
        endif

     endif

  enddo

100 format( '----------------------------------------------------------------------')

  latt%x = latt%scale * matmul( latt%at, xlatt_tmp1 )

  !-------------------------------
  ! set the cutoff in real space
  !-------------------------------
  sh%rmax = maxval( rmax ) 

!==========================================================================
! This has been implemented on the 28th of May 2002, on the plane from 
! London to Washington.
!
! When all the atoms are moved by the same amount, i.e. the crystal is 
! rigidly shifted, the force ON EACH ATOM must be zero. This is a stroger
! constraint than the usual one, in which it is the SUM of the forces on 
! each atom to be zero. The latter is expressed by:
!                \sum_(na,nb,nr) dmat(na,nb,nr) = 0                    (1)
! but the former is:
!                \sum_(na,nr) dmat(na,nb,nr) = 0;    nb = 1, natoms    (2)
! clearly, (2) ==> (1), but the opposite is not true in general. It is (2) 
! to imply that at q=(0,0,0) (gamma) the three acustic branches have
! identically zero frequencies.
!
! The constraint (2) is implemented in what follows. I also have to
! impose the symmetry D(-X) = Dt(X), and since at present I don't have 
! a better idea on how to do the two things together, I am doing it 
! iteratively one after the other. Anybody with a clever idea on how to
! avoid this iterative procedure is very wellcome to act.
!=========================================================================== 
  ! start iterative procedure
  if( nti > 1 ) then
     print'(/,''Imposing traslational invariance, iterative procedure...'')'
     print'(''If you do not want to do it (recover old behaviour) set NTI <=1'',/)'
  endif
  do i = 1,nti
     if( nti > 1 ) then
        ! this is the 'size of the spring' attached to the crystal, we want to zero it
        spring_size = 0._dp
        ! for each atom nbb in the primitive cell
        do nbb=1,latt%natoms
           conta = 0
           ftot = 0._dp
           ! calculate \sum_(naa,nr) dmat(na,nbb,nr) 
           do naa = 1, latt%natoms
              do nr = 1, sh%nrm
                 temp = sh%r(:,nr) + latt%x(:,nbb) - latt%x(:,naa) 
                 if( vabs(temp) <= sh%rmax ) then
                    if(dyn%weight(naa,nbb,nr)>0)ftot = ftot + dyn%dmat(:,:,naa,nbb,nr)
!                    ftot = ftot + dyn%dmat(:,:,naa,nbb,nr)
                    conta = conta + 1
                 endif
              enddo
           enddo
           spring_size = spring_size + vabs(ftot)
           ! impose the constraint (2)
           do naa = 1, latt%natoms
              do nr = 1, sh%nrm
                 temp = sh%r(:,nr) + latt%x(:,nbb) - latt%x(:,naa) 
                 if( vabs(temp) <= sh%rmax ) then
!                    dyn%dmat(:,:,naa,nbb,nr) = dyn%dmat(:,:,naa,nbb,nr) - ftot / conta
                    if(dyn%weight(naa,nbb,nr)>0)dyn%dmat(:,:,naa,nbb,nr) = dyn%dmat(:,:,naa,nbb,nr) - ftot / conta
                 endif
              enddo
           enddo
        enddo
     endif
     ! Now Impose the symmetry D(-X) = Dt(X)
     
     if( symm%lsymm )then
 
        do naa = 1, latt%natoms
           do nbb = 1, latt%natoms
              do nr = 1, sh%nrm

                 ! X = R + ( tau_b - tau_a ) 
                 temp = sh%r(:,nr) + latt%x(:,nbb) - latt%x(:,naa) 

                 if( vabs( temp ) <= sh%rmax ) then

                    ! I need this to find which is the nrr that gives  -X
                    do nrr = 1, sh%nrm
                  
                       ! -X = - [ R + ( tau_b - tau_a ) ] = - R + tau_a - tau_b   
                       tmp1 = sh%r(:,nrr) + latt%x(:,naa) - latt%x(:,nbb) 
                       tmp1 = -tmp1

                       if( vabs( tmp1 ) <= sh%rmax .and. equiva1( temp, tmp1, eps1 ) )then

                          if( vabs(dyn%dmat(:,:,naa,nbb,nr)) > eps .and. &
                               vabs(dyn%dmat(:,:,naa,nbb,nr))<vabs(dyn%dmat(:,:,1,1,1))*1e-30) then
                             print'(''Element '', 3i4,'' set to zero'')', naa,nbb,nr
                             dyn%dmat(:,:,naa,nbb,nr) = 0._dp
                             dyn%weight(naa,nbb,nr) = 0._dp
                          endif

                          if( vabs(dyn%dmat(:,:,nbb,naa,nrr)) > eps .and. &
                               vabs(dyn%dmat(:,:,nbb,naa,nrr))<vabs(dyn%dmat(:,:,1,1,1))*1e-30) then
                             print'(''Element '', 3i4,'' set to zero'')', naa,nbb,nrr
                             dyn%dmat(:,:,nbb,naa,nrr) = 0._dp
                             dyn%weight(nbb,naa,nrr) = 0._dp
                          endif

                          if( vabs(dyn%dmat(:,:,nbb,naa,nrr)) < eps .and. dyn%weight(naa,nbb,nr)>0 ) then
                             !                          dyn%dmat(:,:,nbb,naa,nrr) = transpose(dyn%dmat(:,:,naa,nbb,nr))
                             print*,'Zero matrix',dyn%weight(naa,nbb,nr), naa, nbb, nr
                          endif

                          ! General lattice with a basis: X =  R + ( tau_b - tau_a ) 
                          !                               D(-X) = Dt(X) 
                          tmpmat1(:,:,1) = 0.5_dp * ( dyn%dmat(:,:,nbb,naa,nrr) + &
                               transpose( dyn%dmat(:,:,naa,nbb,nr) ) )
                          dyn%dmat(:,:,naa,nbb,nr) = 0.5_dp * ( dyn%dmat(:,:,naa,nbb,nr) + &
                               transpose( dyn%dmat(:,:,nbb,naa,nrr) ) )
                          dyn%dmat(:,:,nbb,naa,nrr) = tmpmat1(:,:,1)
                          ! If there is inversion symmetry    D(-X) = D(X)     ==> D(X) = Dt(X)
                          if( symm%invsym ) then
                             dyn%dmat(:,:,naa,nbb,nr) = 0.5_dp * ( dyn%dmat(:,:,naa,nbb,nr) + &
                                  transpose( dyn%dmat(:,:,naa,nbb,nr) ) ) 
                             dyn%dmat(:,:,nbb,naa,nrr) = 0.5_dp * ( dyn%dmat(:,:,nbb,naa,nrr) + &
                                  transpose( dyn%dmat(:,:,nbb,naa,nrr) ) ) 
                          endif
                          if( dyn%weight(nbb,naa,nrr) /= dyn%weight(naa,nbb,nr) ) print*,'pesi diversi'
                            exit
                       endif
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
     if(iprint > 0 .and. nti > 1) &
          print'(''Iteration # '',i3,5x,''Spring size: '',e8.2)',i,spring_size 
     if( spring_size < eps ) exit
  enddo

  if(nti > 1 ) &
       print'(''I have been in the loop for '',i3,'' times; Size of the spring: '',e10.1,/)', &
       i, spring_size

  ftot = 0
  do naa = 1, latt%natoms
     do nbb = 1, latt%natoms
        do nr = 1, sh%nrm
           temp = sh%r(:,nr) + latt%x(:,nbb) - latt%x(:,naa) 
           if( vabs(temp) <= sh%rmax ) then
              ftot = ftot + dyn%dmat(:,:,naa,nbb,nr)
           endif
        enddo
     enddo
  enddo

  write(*,100)
  if( vabs(ftot) > 1.d-3 ) then
     print'(/''WARNING: Total sum of the force constants different from zero'')'
     print'(3f16.12)',ftot
  endif
  print'(/''Number of vectors on the edges of the WS cell = '',i5)', wsedge
  print'( ''Number of vectors inside the  WS cell         = '',i5)', wstot - wsedge 
  print'( ''Total number of vectors inside the WS cell    = '',i5)', wstot
  if(lcentral) then
     print'(/''Using central differences '')'
  else
     print'(/''Using forward differences '')'
  endif
  write(*,100)


END SUBROUTINE numeric_force
