!==============================================================================
SUBROUTINE  get_wig( latt, sh, natoms, scale, at, bg,  x, xtmp, &
     &     xlatt_tmp, naa, rmax, ws, xws, eps ) 
  !=============================================================================
  !
  ! Given the set of positions xtmp the routine move the atoms so that on exit
  ! xtmp contains the vectors of the Wigner-Seitz cell
  !
  USE nrtype
  USE nrutil
  USE data

  IMPLICIT NONE

  TYPE(lattice)  :: latt
  TYPE(shell)    :: sh

  INTEGER :: natoms, na, naa, i, j, k,  ws(natoms)

  REAL(DP) :: at(3,3), bg(3,3), x(3,natoms), xtmp(3,natoms), &
       &     xlatt_tmp(3,natoms), temp(3), temp1(3), scale, rmax, xws(3,8,natoms)

  REAL(DP) :: e1(3), e2(3), e3(3)

!  REAL(DP), PARAMETER :: eps = 1.d-6
  REAL(DP) :: eps

  e1 = 0; e2 = 0; e3 = 0; e1(1) = 1; e2(2) = 1; e3(3) = 1

  !--------------------------------------------------
  ! set the origin on the moving atom
  !--------------------------------------------------
  do na = 1, natoms
     xtmp(:,na) = x(:,na) - x(:,latt%super_atom(naa))
  enddo
  do na = 1, latt%natoms
     latt%x(:,na) = xlatt_tmp(:,na) - xlatt_tmp(:,naa)
  enddo

  xtmp =  matmul( at, xtmp ) * scale

  !----------------------------------------------------------------------
  ! construct the WS cell (defined as the closest vectors to the origin)
  !----------------------------------------------------------------------
  do na = 1, natoms
     temp1 = xtmp(:,na)
     do i = -1, 1
        do j = -1, 1
           do k = -1, 1
              temp = xtmp(:,na) + scale*(i*at(:,1) + j*at(:,2) + k*at(:,3))
              if( vabs( temp ) - vabs( temp1 ) < - eps ) then
                 temp1 = temp
              endif
           enddo
        enddo
     enddo
     xtmp(:,na) = temp1

  enddo

  ws = 0
  xws = 0
  do na = 1, natoms
     temp1 = matmul( transpose(bg), xtmp(:,na) ) / scale
     do i = -1, 1
        do j = -1, 1
           do k = -1, 1
              temp = temp1 + i*e1 + j*e2 + k*e3
              temp = matmul( at, temp ) * scale
              if( abs( vabs(temp) - vabs(xtmp(:,na) ) ) < eps ) then
                 ws(na) = ws(na) + 1
                 xws(:,ws(na),na) = temp
              endif
           enddo
        enddo
     enddo
  enddo

  rmax = 0
  do na = 1, natoms
     if( rmax <= vabs(xtmp(:,na)) ) rmax = vabs(xtmp(:,na))
  enddo

  rmax = rmax + 0.05
  print'(/''Rmax = '',f12.6)', rmax
  if( 2*rmax >= sh%rmax ) then
     write(*,*) 'Rcut too small ', sh%rmax,' set it at least equal to ', 2*rmax
     stop
  endif


  RETURN
END SUBROUTINE get_wig
