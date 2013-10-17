!===========================================================================
SUBROUTINE set_forces( force, force1, force2, dis, startdisp, symm, &
     bg, scale, xtmp, ndis, naa, natoms, ndispl, dx, iprint )
  !===========================================================================
  !
  ! generate a set of three independent displacements and three independent forces,
  ! it is based on the point symmetry group of the crystal (sgama1 is called before
  ! the routine, fractional traslations are not allowed)
  !

  USE nrtype
  USE nrutil
  USE data

  IMPLICIT NONE

  LOGICAL :: found, equiva, found1

  INTEGER :: ndispl, i, na, nb, isym, naa, &
       j, natoms, iprint, startdisp, ndis(ndispl)

  REAL(DP) :: force(3,natoms,ndispl), force1(3,natoms,3), &
       force2(3,natoms,3), dx(3,3), dis(3,ndispl), tmp1(3), &
       invis(3,3), temp(3), xtmp(3,natoms), scale, bg(3,3), &
       invdx(3,3), dxnorm(3,3)


  REAL(DP) :: eps, eps1

  TYPE(symmetry) :: symm

!  eps = symm%symprec / 100
  eps = symm%symprec 
  eps1 = symm%symprec

  !---------------------
  ! first displacement
  !---------------------
  force1 = 0

  dx(:,1) = dis(:,startdisp)
  print'(/''Atom #'',i3,'' displacement #'',i3,3f10.5)',naa,startdisp,dx(:,1)

  force1(:,:,1) = force(:,:,startdisp)

  !---------------------------------------------------
  ! second displacement
  !---------------------------------------------------

  j = startdisp 
  dosym: do isym = 1, symm%nsym

     dx(:,2) = matmul( symm%is(:,:,isym), dis(:,j) )
     !----------------------------------------------------------------
     ! is dx(:,2) linearly independent from dx(:,1) ?
     !----------------------------------------------------------------
     call vecprod( dx(:,1)/vabs(dx(:,1)), dx(:,2)/vabs(dx(:,2)), tmp1 )
     if( vabs( tmp1 ) > eps ) then
        found = .true.
     else
        dx(:,2) = 0
        found = .false.
     endif
     if(found) then
        print'(/''Second displacement, found the symmetry #'',i2,/&
             &     ''The displacement is = '',3f8.4)', isym,dx(:,2)
        call inve( symm%is(:,:,isym), invis )
        !===========================================================================
        ! F_na ( S*u ) = S * F_{S^-1*na} ( u )                                 (0)
        !===========================================================================
        !----------------------
        ! So, for each na.....
        !----------------------
        do na=1,natoms
           !----------------------------
           ! construct temp = S^-1 * na
           !----------------------------
           temp = xtmp(:,na) 
           temp = matmul( invis, temp ) 
           temp = matmul( transpose(bg), temp ) / scale 
           !---------------------------------
           ! find nb so that nb = S^-1 * na
           !---------------------------------
           do nb=1,natoms
              found1 = .false. 
              !----------------------------------------------------------------------------
              ! it must be 'equiva' and not 'equiva1', i.e. nb = S^-1 * na + R and not
              ! nb = S^-1 * na, because not all the atoms on the edge of the Wigner-Seitz 
              ! cell are present in the xtmp. So, it could happen (and happens) 
              ! that S^-1 * na corresponds to the atom on the "wrong" side of the WS cell.
              ! This atom is not present in 'xtmp', but its relative on the other side of
              ! the WS cell is good, so take that one.
              !----------------------------------------------------------------------------  
              tmp1 = matmul( transpose(bg), xtmp(:,nb) ) / scale
              if( equiva( temp, tmp1, eps1 ) )then
                 found1 = .true.
                 !------------------------------------------
                 ! and calculate F_na ( S*u ) using Eq. (0)
                 !------------------------------------------
                 force1(:,na,2) = matmul( symm%is(:,:,isym), force(:,nb,j) )
                 exit  
              endif
           enddo
           if( .not. found1 )stop 'error in set_forces'
        enddo
        exit dosym
     endif
  enddo dosym
  if(.not.found) then
     j = j + 1
     dx(:,2) = dis(:,j)
     print'(/''Symmetry for second displacement not found''/''Using'', &
          &     '' next displacement #'',i4,3f8.4)',j,dx(:,2)
     force1(:,:,2) = force(:,:,j)
  endif

  !------------------------------------------------------------------
  ! Third displacement
  !------------------------------------------------------------------
  dosym1: do isym = 1, symm%nsym
     dx(:,3) = matmul( symm%is(:,:,isym), dis(:,j) )
     !----------------------------------------------------------------
     ! is dx(:,3) linearly independent from dx(:,1) and dx(:,2) ?
     !----------------------------------------------------------------
     call vecprod( dx(:,1), dx(:,2), tmp1 )
     temp = dot_product( tmp1/vabs(tmp1), dx(:,3)/vabs(dx(:,3)) )
     if( vabs( temp ) > eps ) then
        found=.true.
     else
        dx(:,3) = 0
        found=.false.
     endif
     if(found)then
        print'(/''Third displacement, found the symmetry #'',i2,/&
             &     ''The displacement is = '',3f8.4)', isym,dx(:,3)
        !------------------------------------------------
        ! see comments above for the second displacement
        !------------------------------------------------
        call inve( symm%is(:,:,isym), invis )
        do na = 1, natoms
           temp = xtmp(:,na) 
           temp = matmul( invis, temp ) 
           temp = matmul( transpose(bg), temp ) / scale
           found1 = .false.
           do nb=1,natoms
              tmp1 = matmul( transpose(bg), xtmp(:,nb) ) / scale
              if( equiva( temp, tmp1, eps1 ) )then
                 found1 = .true.
                 force1(:,na,3) = matmul( symm%is(:,:,isym), force(:,nb,j) )
                 exit 
              endif
           enddo 
           if( .not. found1 ) stop 'error in set_forces'
        enddo
        exit dosym1
     endif
  enddo dosym1
  if(.not.found) then
     j = j + 1
     dx(:,3) = dis(:,j)
     print'(/''Symmetry for third displacement not found,''/''Using'', &
          &     '' next displacement #'',i4,3f8.4)',j,dx(:,3)
     force1(:,:,3) = force(:,:,j)
  endif

  startdisp = startdisp + ndis( naa ) 

  !-------------------------------------------------------------------------
  ! Create a linear combination of the displacements and the forces so that 
  ! the three displacements are orientated along the three carthesian axis:
  ! Note that the matrix defining the linear combination is arbitrary. I use
  ! the inverse of the matrix whose columns are the displacements (and the
  ! linear combined displacements become e1, e2 and e3).
  !-------------------------------------------------------------------------

  invdx=0
  do i=1,3
     dxnorm(:,i) = dx(:,i) / vabs(dx(:,i))
  enddo
  call inve( dxnorm, invdx )

  do na = 1, natoms
     force2(:,na,1) = invdx(1,1) * force1(:,na,1) + &
          invdx(2,1) * force1(:,na,2) + invdx(3,1) * force1(:,na,3) 
     force2(:,na,2) = invdx(1,2) * force1(:,na,1) + &
          invdx(2,2) * force1(:,na,2) + invdx(3,2) * force1(:,na,3) 
     force2(:,na,3) = invdx(1,3) * force1(:,na,1) + &
          invdx(2,3) * force1(:,na,2) + invdx(3,3) * force1(:,na,3) 
  enddo

  if(iprint>1)then
     print*,'Displacement, and inverse'
     print'(3f10.5)',dx,invdx
     print*,'Force x'
     print'(3f10.5)',force2(:,:,1)
     print*,'Force y'
     print'(3f10.5)',force2(:,:,2)
     print*,'Force z'
     print'(3f10.5)',force2(:,:,3)
     print*,' '
  endif

  dx = matmul( invdx, dx )


END SUBROUTINE set_forces

