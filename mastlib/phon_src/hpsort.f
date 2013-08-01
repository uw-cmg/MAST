!---------------------------------------------------------------------
      subroutine  hpsort(n,ra,ind) 
!---------------------------------------------------------------------
! sort an array ra(1:n) into ascending order using heapsort algorithm.
! n is input, ra is replaced on output by its sorted rearrangement.
! create an index table (ind) by making an exchange in the index array 
! whenever an exchange is made on the sorted data array (ra).
! in case of equal values in the data array (ra) the values in the 
! index array (ind) are used to order the entries.
! if on input ind(1)  = 0 then indices are initialized in the routine,
! if on input ind(1) != 0 then indices are assumed to have been 
!                initialized before entering the routine and these
!                indices are carried around during the sorting process
!
! no work space needed !
! free us from machine-dependent sorting-routines !
!
! adapted from Numerical Recipes pg. 329 (new edition)
!
      USE nrtype
      IMPLICIT NONE
!-input/output variables
      integer n
      integer ind(n)
      real(dp) ra(n)
!-local variables
      integer i, ir, j, l,iind
      real(dp) rra
! initialize index array
      if (ind(1).eq.0) then
        do i=1,n
          ind(i)=i
        end do
      end if
      if (n.lt.2) return    ! nothing to order
! initialize indices for hiring and retirement-promotion phase
      l  = n/2 + 1
      ir = n
   10 continue
        if (l.gt.1) then    ! still in hiring phase
          l   = l-1
          rra = ra(l)
          iind= ind(l)
        else                ! in retirement-promotion phase.
          rra = ra(ir)      ! clear a space at the end of the array
          iind= ind(ir)     !  
          ra(ir) = ra(1)    ! retire the top of the heap into it
          ind(ir)= ind(1)   !
          ir     = ir -1    ! decrease the size of the corporation
          if (ir.eq.1) then ! done with the last promotion
             ra(1) = rra    ! the least competent worker at all !
             ind(1)= iind   !
             return
          end if 
        end if 
        i = l               ! wheter in hiring or promotion phase, we 
        j = l + l           ! set up to place rra in its proper level


        do while (j.le.ir) 
          if (j.lt.ir) then
            if (ra(j).lt.ra(j+1)) then ! compare to better underling
               j = j+1 
            else if (ra(j).eq.ra(j+1)) then
              if (ind(j).lt.ind(j+1)) j = j+1
            end if
          end if 
          if (rra.lt.ra(j)) then ! demote rra
            ra(i) = ra(j)
            ind(i)= ind(j)
            i = j
            j = j + j
          else if (rra.eq.ra(j)) then
            if (iind.lt.ind(j)) then ! demote rra
              ra(i) = ra(j)
              ind(i)= ind(j)
              i = j
              j = j + j
            else
              j = ir + 1         ! set j to terminate do-while loop
            end if 
          else                   ! this is the right place for rra
            j = ir + 1           ! set j to terminate do-while loop
          end if
        end do
        ra(i) = rra
        ind(i)= iind
      go to 10   

      end 


