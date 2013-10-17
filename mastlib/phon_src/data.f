!------------------------------------------------------
! Definitions of the main data structures used by PHON
!------------------------------------------------------
MODULE data
  USE nrutil
  IMPLICIT NONE

  TYPE lattice  
     CHARACTER*1 :: string
     CHARACTER*3, POINTER :: name(:)
     CHARACTER*50 :: string2
     LOGICAL :: lsuper, lgamma
     INTEGER :: ntypes, natoms, natomss, ndim(3), qa, qb, qc, disp
     INTEGER, POINTER :: nions(:), ityp(:), ityps(:), super_atom(:)
     REAL(DP) :: at(3,3), bg(3,3), scale, volume
     REAL(DP) :: ats(3,3), bgs(3,3), volumes
     REAL(DP) :: dxstart(3)
     REAL(DP), POINTER :: x(:,:), xs(:,:), mass(:)
  END TYPE lattice

  TYPE shell
     INTEGER :: nrm, nsh
     REAL(DP) :: rmax, cutoff
     INTEGER, POINTER :: nl(:), ishel(:)
     REAL(DP), POINTER :: r(:,:), rr(:), rl(:)
  END TYPE shell

  TYPE dynmat
     LOGICAL :: lfree
     LOGICAL,      POINTER :: usethis(:), pdos_atom(:)
     INTEGER :: npoints
     INTEGER,      POINTER :: weight(:,:,:)
     REAL(DP),     POINTER :: tmpmat2(:,:,:,:), dmat(:,:,:,:,:), &
          qi(:,:), qf(:,:), eig(:)
     COMPLEX(DPC), POINTER :: cmat(:,:)
  END TYPE dynmat

  TYPE symmetry
     LOGICAL :: lsymm, invsym
     INTEGER :: nsym, nrot, nsymstart
     REAL(DP) :: symprec
     REAL(DP), POINTER :: is(:,:,:), isstart(:,:,:), iscryst(:,:,:), ftau(:,:)
     INTEGER, POINTER :: irt(:,:), irts(:,:)
  END TYPE symmetry


END MODULE data
