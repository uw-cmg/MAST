! 
!-----------------------------------------------------------------------------
subroutine hexsym( at, is, nrot, eps )
  !-----------------------------------------------------------------------------
  !
  ! Provides symmetry operations for Hexagonal and Trigonal lattices.
  ! The c axis is assumed to be along the z axis
  !
  !
  !
  USE nrtype
  USE nrutil
  IMPLICIT NONE
  !
  !     first the input variables
  !
  integer :: &
       &      nrot                ! output: the number of symmetry matrices 
  !
  !    here the local parameters
  !  
  REAL(DP), PARAMETER :: sin3 = 0.866025403784438597_dp, cos3 = 0.5_dp, &
       &           msin3 =-0.866025403784438597_dp, mcos3 =-0.5_dp  
  !
  !   and the local variables
  !
  real(dp) ::  &
       &      s(3,3,12),   &       ! the s matrices in real variables
       &      overlap(3,3), &     ! overlap matrix between direct lattice vectors
       &      invover(3,3), &     ! overlap matrix between direct lattice vectors
       &      rat(3),     &        ! the rotated of a direct vector ( cartesian )
       &      rot(3,3),   &        ! the rotated of a direct vector ( in axis )
       &      value, eps               ! component of the s matrix in axis basis
  integer &
       &      jpol, &               ! counter over the polarizations
       &      kpol,  &             ! counter over the polarizations
       &      irot,   &            ! counter over the rotations
       &      mpol                ! counter over the polarizations 

  real(dp) :: &
       &      at(3,3), &             ! input: the direct lattice vectors
       &      is(3,3,nrot)        ! output: the symmetry matrices

  data s/ &
       &     1._dp, 0._dp, 0._dp,  0._dp, 1._dp, 0._dp,  0._dp, 0._dp, 1._dp, &
       &    -1._dp, 0._dp, 0._dp,  0._dp,-1._dp, 0._dp,  0._dp, 0._dp, 1._dp, &
       &    -1._dp, 0._dp, 0._dp,  0._dp, 1._dp, 0._dp,  0._dp, 0._dp,-1._dp, &
       &     1._dp, 0._dp, 0._dp,  0._dp,-1._dp, 0._dp,  0._dp, 0._dp,-1._dp, &
       &     cos3, sin3, 0._dp, msin3, cos3, 0._dp,  0._dp, 0._dp, 1._dp, &
       &     cos3,msin3, 0._dp,  sin3, cos3, 0._dp,  0._dp, 0._dp, 1._dp, &
       &    mcos3, sin3, 0._dp, msin3,mcos3, 0._dp,  0._dp, 0._dp, 1._dp, &
       &    mcos3,msin3, 0._dp,  sin3,mcos3, 0._dp,  0._dp, 0._dp, 1._dp, &
       &     cos3,msin3, 0._dp, msin3,mcos3, 0._dp,  0._dp, 0._dp,-1._dp, &
       &     cos3, sin3, 0._dp,  sin3,mcos3, 0._dp,  0._dp, 0._dp,-1._dp, &
       &    mcos3,msin3, 0._dp, msin3, cos3, 0._dp,  0._dp, 0._dp,-1._dp, &
       &    mcos3, sin3, 0._dp,  sin3, cos3, 0._dp,  0._dp, 0._dp,-1._dp /


  !
  !   first compute the overlap matrix between direct lattice vectors
  ! 

  do jpol = 1,3
     do kpol = 1,3
        overlap(kpol,jpol) = at(1,kpol)*at(1,jpol) + &
             &                           at(2,kpol)*at(2,jpol) + &
             &                           at(3,kpol)*at(3,jpol)
     enddo
  enddo
  !
  !    then its inverse
  ! 
  call inve(overlap,invover)

  nrot = 1
  do irot = 1,12
     !
     !   for each possible simmetry
     !
     do jpol = 1,3
        do mpol = 1,3
           !
           !   compute, in cartesian coordinates the rotated vector
           !
           rat(mpol) = s(mpol,1,irot)*at(1,jpol)  + &
                &                     s(mpol,2,irot)*at(2,jpol) + &
                &                     s(mpol,3,irot)*at(3,jpol)
        enddo

        do kpol = 1,3
           !
           !   the rotated vector is projected on the direct lattice
           !
           rot(kpol,jpol) = at(1,kpol)*rat(1) + &
                &                          at(2,kpol)*rat(2) + &
                &                          at(3,kpol)*rat(3)
        enddo
     enddo
     !
     !  and the inverse of the overlap matrix is applied
     !
     do jpol = 1,3
        do kpol = 1,3
           value = invover(jpol,1)*rot(1,kpol) + &
                &                 invover(jpol,2)*rot(2,kpol) + &
                &                 invover(jpol,3)*rot(3,kpol)
           if ( abs(float(nint(value))-value).gt.eps) then
              !
              ! if a noninteger is obtained, this implies that this operation
              ! is not a symmetry operation for the given lattice
              !
              go to 10
           end if
           !               is(kpol,jpol,nrot) = nint(value)
           is(kpol,jpol,nrot) = s(kpol,jpol,irot)
        enddo
     enddo
     nrot = nrot+1
10   continue
  enddo
  nrot = nrot-1
  !
  !     set the inversion symmetry ( Bravais lattices have always inversion )
  ! 
  do irot = 1, nrot
     do kpol = 1,3
        do jpol = 1,3
           is(kpol,jpol,irot+nrot) = -is(kpol,jpol,irot)
        end do
     end do
  end do

  nrot = 2*nrot

  return
end subroutine hexsym
