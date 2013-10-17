!-----------------------------------------------------------------------
SUBROUTINE simpson(mesh,dx,func,asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integrator 
  !
  USE nrtype
  IMPLICIT NONE

  INTEGER :: i, mesh
  REAL(DP) :: func(mesh), f1, f2, f3, r12, asum, dx

  asum = 0.0d0
  r12 = 1.0d0 / 12.0d0
  f3  = func(1) * r12 

  do i = 2,mesh-1,2
     f1 = f3
     f2 = func(i) * r12
     f3 = func(i+1) * r12
     asum = asum + 4.0d0*f1 + 16.0d0*f2 + 4.0d0*f3
  enddo
  asum = asum*dx

END SUBROUTINE simpson
