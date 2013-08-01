!============================================================================
! Reads input variables
!============================================================================
SUBROUTINE reader( latt, dyn, temperature, ptemp, lforceout, lrecip, rmax, &
     dosin, dosend, dosstep, dossmear, nd, iprint, lsymm, nti, symprec, lcentral, &
     eigsize, ncycleseig, leigen, ldrift, lnosym )

  USE nrtype
  USE data

  IMPLICIT NONE

  LOGICAL :: LDUM, lopen, lforceout, lrecip, lsymm, lcentral, leigen, &
      ldrift, lnosym

  CHARACTER  :: CHARAC, c
  CHARACTER(5) :: name

  INTEGER :: IDUM, IERR, N, ntypes, nd, iprint, qa, qb, qc, nti, i, ncycleseig

  REAL(dp) :: RDUM, temperature, rmax, dosin, dosend, dosstep, dossmear, &
      dx, dy, dz, ptemp(3), symprec, eigsize

  COMPLEX(dpc) :: CDUM

  INTEGER, PARAMETER :: iu5 = 15, iu6 = 6

  TYPE(lattice) :: latt
  TYPE(dynmat)  :: dyn

  ! 'title'-string (defaults to 'unknown system'), keyword 'SYSTEM'
  lopen = .false.
  open(iu5,file='INPHON',status='old')

  !----------------------------------------------------------
  ! number of cycles for graphical representation of phonons
  !----------------------------------------------------------
  ncycleseig = 2
  call rdatab( lopen, 'INPHON', iu5, 'ncycleseig', '=', '#', ';', 'I', &
       ncycleseig, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''ncycleseig'' from file INPHON.'
     goto 150
  endif

  !----------------------------------------------------------------
  ! (inverse of) size of displacements in graphical representation
  !----------------------------------------------------------------
  eigsize = 1.0
  call rdatab( lopen, 'INPHON', iu5, 'eigsize', '=', '#', ';', 'F', &
       IDUM, eigsize, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''eigsize'' from file INPHON.'
     goto 150
  endif


  !------------------------------------------------------
  ! do you want to impose traslational invariance?
  ! How many loops in the procedure? (see numeric_force)
  !------------------------------------------------------
  nti = 1
  call rdatab( lopen, 'INPHON', iu5, 'nti', '=', '#', ';', 'I', &
       nti, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''nti'' from file INPHON.'
     goto 150
  endif
  if( nti < 1 ) nti = 1

  !--------------------------------------
  ! how many types of atoms? (no default)
  !--------------------------------------
  latt%ntypes = -1
  call rdatab( lopen, 'INPHON', iu5, 'ntypes', '=', '#', ';', 'I', &
       latt%ntypes, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''ntypes'' from file INPHON.'
     goto 150
  endif
  if(  latt%ntypes == -1 ) then
     write(iu6,*)'Error, there is no default for item ''ntypes'' '
     stop
  endif

  !----------------------------------------------
  ! allocate space from nions, mass, and usethis
  !----------------------------------------------
  allocate( latt%nions( latt%ntypes ), latt%mass( latt%ntypes), & 
            dyn%usethis(latt%ntypes), latt%name(latt%ntypes) )

  dyn%usethis = .true.
  call rdatab( lopen, 'INPHON', iu5, 'usethis', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, dyn%usethis, CHARAC, N, latt%ntypes, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''usethis'' from file INPHON.'
     goto 150
  endif

  !-----------------------------------------
  ! atom names
  !-----------------------------------------
  latt%name = 'H'
  do i = 1, latt%ntypes
     write(c,'(i1)')i
     name = 'name'//c
     call rdatab( lopen, 'INPHON', iu5, name, '=', '#', ';', 'S', &
          IDUM, RDUM, CDUM, LDUM, latt%name(i), N, 3*latt%ntypes, IERR )
     if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
          ((IERR==0).and.(N<1))) then
        write(iu6,*)'error reading item ''name'' from file INPHON.'
        goto 150
     endif
  enddo

  !-----------------------------------------
  ! readin symmetry precision (default 1e-5)
  !-----------------------------------------
  symprec = 1d-5
  call rdatab( lopen, 'INPHON', iu5, 'symprec', '=', '#', ';', 'F', &
       IDUM, symprec, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''symprec'' from file INPHON.'
     goto 150
  endif

  !----------------------------
  ! readin masses (default 1)
  !----------------------------
  latt%mass = 1
  call rdatab( lopen, 'INPHON', iu5, 'mass', '=', '#', ';', 'F', &
       IDUM, latt%mass(1), CDUM, LDUM, CHARAC, N, latt%ntypes, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''mass'' from file INPHON.'
     goto 150
  endif

  !---------------------------------
  ! Superlattice ? (default .true.)
  !---------------------------------
  latt%lsuper = .true.
  call rdatab( lopen, 'INPHON', iu5, 'lsuper', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, latt%lsuper, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lsuper'' from file INPHON.'
     goto 150
  endif

  !-------------------------------------------
  ! dimension of the supercell (default 1 1 1)
  !-------------------------------------------
  latt%ndim = 1  
  call rdatab( lopen, 'INPHON', iu5, 'ndim', '=', '#', ';', 'I', &
       latt%ndim(1), RDUM, CDUM, LDUM, CHARAC, N, 3, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''ndim'' from file INPHON.'
     goto 150
  endif

  !-----------------------------------------------------
  ! Thermodynamics (default .false.)
  !-----------------------------------------------------
  dyn%lfree = .false.
  call rdatab( lopen, 'INPHON', iu5, 'lfree', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, dyn%lfree, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lfree'' from file INPHON.'
     goto 150
  endif

  !----------------------------------------------------------------------------------------
  !  temperature ( default = -1, i.e. not set )
  !----------------------------------------------------------------------------------------
  temperature = -1
  call rdatab( lopen, 'INPHON', iu5, 'temperature', '=', '#', ';', 'F', &
       IDUM, temperature, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''temperature'' from file INPHON.'
     goto 150
  endif
  if( dyn%lfree .and. temperature < 0 ) stop 'TEMPERATURE not set, cannot proceed' 

  !-----------------------------------------------------------
  ! parameters for temperature ( INCREMENT NINCREMENT )
  !-----------------------------------------------------------
  ptemp = (/ 1.0, 0.0, 1.0 /)
  ptemp(1) = ptemp(1)*temperature
  call rdatab( lopen, 'INPHON', iu5, 'ptemp', '=', '#', ';', 'F', &
       IDUM, ptemp(2), CDUM, LDUM, CHARAC, N, 2, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''ptemp'' from file INPHON.'
     goto 150
  endif

  !------------------------------------------------
  ! Save force constant matrix ? (default .false.)
  !------------------------------------------------
  lforceout = .false.
  call rdatab( lopen, 'INPHON', iu5, 'lforceout', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, lforceout, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lforceout'' from file INPHON.'
     goto 150
  endif

  !----------------------------------------------------------------------------------
  ! Impose symmetry on the force constat matrix, D(-R) = D^T (R) ? (default .true. )
  !----------------------------------------------------------------------------------
  lsymm = .true.
  call rdatab( lopen, 'INPHON', iu5, 'lsymm', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, lsymm, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lsymm'' from file INPHON.'
     goto 150
  endif

  !----------------------------------------------------------------------------------
  ! Print out files MODE???.axsf and EIGEN.axsf for phonon visualisation with
  ! xcsysdens
  !----------------------------------------------------------------------------------
  leigen = .false.
  call rdatab( lopen, 'INPHON', iu5, 'leigen', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, leigen, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''leigen'' from file INPHON.'
     goto 150
  endif

  !----------------------------------------------------------------------------------
  ! Central differences (default .false.; forward differences )
  !----------------------------------------------------------------------------------
  lcentral = .false.
  call rdatab( lopen, 'INPHON', iu5, 'lcentral', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, lcentral, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lcentral'' from file INPHON.'
     goto 150
  endif

  !-------------------------------------------------------------------------------------------
  ! cutoff radius in real space (default 15 A). Setting this parameter should not be necessary
  !-------------------------------------------------------------------------------------------
  rmax = 0
  call rdatab( lopen, 'INPHON', iu5, 'rmax', '=', '#', ';', 'F', &
       IDUM, rmax, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''rmax'' from file INPHON.'
     goto 150
  endif

  !--------------------------------------------------------------------
  ! phonon q points in reciprocal lattice coordinates (default .true.)
  !--------------------------------------------------------------------
  lrecip = .true.
  call rdatab( lopen, 'INPHON', iu5, 'lrecip', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, lrecip, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lrecip'' from file INPHON.'
     goto 150
  endif


  !---------------------------------------------------------------------------------------
  ! do you want to suggest a displacement (default first displacement is along x; 1,0,0) ?
  !---------------------------------------------------------------------------------------
  latt%dxstart(1) = 1._dp; latt%dxstart(2) = 0._dp; latt%dxstart(3) = 0._dp
  call rdatab( lopen, 'INPHON', iu5, 'dxstart', '=', '#', ';', 'F', &
       IDUM, latt%dxstart, CDUM, LDUM, CHARAC, N, 3, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''dxstart'' from file INPHON.'
     goto 150
  endif

  !--------------------------------------------------------------
  ! how large must be the displacement? ( default 0.04 A )
  !--------------------------------------------------------------
  latt%disp = 25
  call rdatab( lopen, 'INPHON', iu5, 'disp', '=', '#', ';', 'I', &
       latt%disp, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''disp'' from file INPHON.'
     goto 150
  endif

  !--------------------------------------------------------------
  ! MP Q-points construction
  !--------------------------------------------------------------
  latt%qa = -1
  call rdatab( lopen, 'INPHON', iu5, 'qa', '=', '#', ';', 'I', &
       latt%qa, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''qa'' from file INPHON.'
     goto 150
  endif
  latt%qb = -1
  call rdatab( lopen, 'INPHON', iu5, 'qb', '=', '#', ';', 'I', &
       latt%qb, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''qb'' from file INPHON.'
     goto 150
  endif
  latt%qc = -1
  call rdatab( lopen, 'INPHON', iu5, 'qc', '=', '#', ';', 'I', &
       latt%qc, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''qc'' from file INPHON.'
     goto 150
  endif

  latt%lgamma = .false.
  call rdatab( lopen, 'INPHON', iu5, 'lgamma', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, latt%lgamma, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lgamma'' from file INPHON.'
     goto 150
  endif

  !------------------------------------------------------------------
  ! density of states, start, finish, step (default 0, 25, 0.1 THz)
  !------------------------------------------------------------------
  dosin = 0.0
  call rdatab( lopen, 'INPHON', iu5, 'dosin', '=', '#', ';', 'F', &
       IDUM, dosin, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''dosin'' from file INPHON.'
     goto 150
  endif
  dosend = 25.0
  call rdatab( lopen, 'INPHON', iu5, 'dosend', '=', '#', ';', 'F', &
       IDUM, dosend, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''dosend'' from file INPHON.'
     goto 150
  endif
  dosstep = 0.1
  call rdatab( lopen, 'INPHON', iu5, 'dosstep', '=', '#', ';', 'F', &
       IDUM, dosstep, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''dosstep'' from file INPHON.'
     goto 150
  endif
  dossmear = 0.02
  call rdatab( lopen, 'INPHON', iu5, 'dossmear', '=', '#', ';', 'F', &
       IDUM, dossmear, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''dossmear'' from file INPHON.'
     goto 150
  endif

  !----------------------------------
  ! number of q points (default = 0)
  !----------------------------------
  nd = 0
  call rdatab( lopen, 'INPHON', iu5, 'nd', '=', '#', ';', 'I', &
       nd, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''nd'' from file INPHON.'
     goto 150
  endif

  if( nd > 0 ) then 

     allocate( dyn%qi(3,nd), dyn%qf(3,nd) )

     !-----------------------------------------
     ! number of points between two q-points
     !-----------------------------------------
     dyn%npoints = 1
     call rdatab( lopen, 'INPHON', iu5, 'npoints', '=', '#', ';', 'I', &
          dyn%npoints, RDUM, CDUM, LDUM, CHARAC, N, 3, IERR )
     if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
          ((IERR==0).and.(N<1))) then
        write(iu6,*)'error reading item ''npoints'' from file INPHON.'
        goto 150
     endif

     !-------------
     ! starting qi
     !-------------
     dyn%qi = 0
     call rdatab( lopen, 'INPHON', iu5, 'qi', '=', '#', ';', 'F', &
          IDUM, dyn%qi(1,1), CDUM, LDUM, CHARAC, N, 3*nd, IERR )
     if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
          ((IERR==0).and.(N<1))) then
        write(iu6,*)'error reading item ''qi'' from file INPHON.'
        goto 150
     endif

     !---------------
     ! end qf
     !---------------
     dyn%qf = 0
     call rdatab( lopen, 'INPHON', iu5, 'qf', '=', '#', ';', 'F', &
          IDUM, dyn%qf(1,1), CDUM, LDUM, CHARAC, N, 3*nd, IERR )
     if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
          ((IERR==0).and.(N<1))) then
        write(iu6,*)'error reading item ''qf'' from file INPHON.'
        goto 150
     endif

  endif


  !---------------------------------------------------------
  ! output verbosity (default 0, keep verbosity at minimum)
  !---------------------------------------------------------
  iprint = 0
  call rdatab( lopen, 'INPHON', iu5, 'iprint', '=', '#', ';', 'I', &
       iprint, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''iprint'' from file INPHON.'
     goto 150
  endif

  !------------------------------------------------
  ! Remove drift forces on atoms ? (default .true.)
  !------------------------------------------------
  ldrift = .true.
  call rdatab( lopen, 'INPHON', iu5, 'ldrift', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, ldrift, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''ldrift'' from file INPHON.'
     goto 150
  endif

  !------------------------------------------------
  ! Use symmetry finder for forces ? (default .false.)
  !------------------------------------------------
  lnosym = .false.
  call rdatab( lopen, 'INPHON', iu5, 'lnosym', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, lnosym, CHARAC, N, 1, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(iu6,*)'error reading item ''lnosym'' from file INPHON.'
     goto 150
  endif


  close(iu5)

  RETURN
150 continue
  write(iu6,151) IERR, N

151 format(' Error code was IERR=',I1,' ... . Found N=',I5,' data.')
  STOP

END SUBROUTINE reader


