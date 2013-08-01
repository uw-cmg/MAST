!=========================================================================
!
! Program PHON: 
!
!
! Force constant matrix calculation program. Uses finite displacements.
!
!
! Dario Alfe`,      November  1998
!
! This program is freely distributed and comes with no warranty. Please
! send any comments to the author (d.alfe@ucl.ac.uk)
! 
! If you use this code to publish scientific results please include the following citation: 
! D. Alf\`e, "PHON: A program to calculate phonons using the small displacement method",
! Computer Physics Communications, Vol. 180, No. 12, pp. 2622-2633 (2009)
!
! 
! 26/2/2003  added partial density of states (USETHIS)   
! 21/8/2003  minor modification in 'numeric_force', which didn't compile
!            properly with the linux-pgf90 compiler
! 30/9/2003  added translational invariance symmetrization
! 16/12/2003 added makefile for linux ifc compiler
! 15/07/2004 removed generation of SPOSCART file (not needed anyway)
! 28/09/2005 corrected DXSTART; increased default of RMAX
! 03/03/2006 increased default of RMAX in the right routine!
! 12/03/2006 Drift in forces is not removed (for non periodic calculations)
! 18/03/2006 Corrected error in generate_punti
! 16/04/2006 Problem with symmetry in numeric_force hopefully fixed
!            (check performed in crystal coordinates rather than carthesian)
! 27/09/2006 Printing now the frequencies for iprint > 2 instead of the eigenvalues
! 06/10/2006 changed iu5=5 to iu5=15; added the example directory Al 
! 15/12/2006 Corrected bug in generate_punti
! 10/10/2007 Changed generation of displacements, now by default all 
!            have an amplitude of 0.04 A. This can be modified by changing
!            the value of the input variable DISP (default DISP=25)
! 06/12/2007 Changed warning message, now given when the sum of the force 
!            constants is larger than 1.d-3.
!            Corrected makefile.ibm 
! 06/05/2008 Corrected bug in get_displ.f, a confusion between carthesian and
!            crystal coordinates was suggesting wrong displacements in some cases
! 09/05/2008 Changed default amplitude of displacements to 0.02 A. This can be
!            modified by changing the value of the input variable DISP (default DISP=50)
! 20/05/2008 Added loop over temperature for thermodynamic properties (PTEMP).
!            File "THERMO" contains thermodynamic data as function of temperature
!            (suggestion by Peter Agoston).
! 23/05/2008 Corrected an inconstistency between get_displ and set_forces, which
!            in some special case would make the suggestion of displacements 
!            incorrect.
! 05/06/2008 Changed value of eps from 1.e-6 to 2.e-6 in get_wig; default value of 
!            DISP changed back to 25 (amplitude = 0.04 A)
!
! 03/03/2009 Polished for publication in Computer Physics Communications
! 22/10/2009 Removed printout of omegam2 and omegabar, left over from previous versions
! 24/11/2009 Introduced input variable SYMPREC (default 1d-6). Now eps is set equal to
!            SYMPREC in most parts of the code.
!            Reorganised FREQ? files, now they contain 48 branches each; FREQ and FREQ.cm
!            are not written if there are more than 48 branches (16 atoms in the primitive cell)
! 13/01/2010 Introduced possibility of doing central differences. Setting the variable 
!            LCENTRAL=.T. causes the code to produce + and - displacements. The FORCES used
!            by the code are then F = ( F(+) - F(-) ) / 2. 
!            Added interface with graphic program XCrysDen (http://www.xcrysden.org) for
!            graphical representations of vibrations at Gamma
! 20/01/2010 Corrected a bug for the creation of the EIGEN.axsf and MODE???.axsf files; added
!            creation of EIGEN.xyz if LEIGEN=.T., vibration can be visualised with jmol
! 27/05/2010 Corrected a bug in numeric_force.f appearing for central differences
! 4/06/2010  New way of calculating RMAX in set_lattice.f. Set RMAX=0 in reader.f.
!            Changed definition of rmax in get_wig.f. Removed LCENTRAL=T in
!            examples/Al/INPHON
! 19/06/2010 Changed back definition of rmax in get_wig.f. Drift of forces is
!            now removed in numeric_force.f
! 15/06/2010 .....
! 13/09/2010 Corrected definition of RMAX in set_lattice and get_wig. 
! 21/09/2010 Modified set_forces and get_wig which should now be completely compatible and
!            produce the same results independently from the size of the
!            displacement
! 24/07/2012 Reduced default value of SYMPREC to 1d-5
!            Added example of Graphene
! 15/12/2012 Removed extra symmetrisation with inversion symmetry in numeric_force.
! 31/01/2013 Improved (hopefully) installation procedure (see README)
! 10/04/2013 Corrected procedure to impose translational invariance (zero weight fcm are not used)
! 24/04/2013 Added capability of calculating partial density of states projected on single
!            atoms (pdos_atom) 
! 18/06/2013 Corrected definition of force constant matrix, which was the transpose of what it should 
!            have been. Bug found by David Olmsted. Note that the HARMONIC file still contains the 
!            force constant matrix in the old format.
!            Symmetrisation of force constant matrix Dt(R) = D(R) for systems with inversion symmetry
!            
!            
!            
! 22/07/2013 outside modification: add tag "ldrift" to control whether or not to remove drift forces
!            add tag "lnosym" to allow no symmetry finding for forces
!========================================================================= 

PROGRAM phon

  USE nrtype
  USE nrutil
  USE diagonalize
  USE data

  IMPLICIT NONE

  CHARACTER(1) :: fl
  CHARACTER(2) :: fl1
  CHARACTER(7) :: dum
  CHARACTER(3) :: nd_nmbr
  CHARACTER(80), PARAMETER :: version='  PHON, VERSION 1.36 (18/6/2013)'
  CHARACTER(80), PARAMETER :: citation='  Please cite: D. Alfe`, Comp. Phys. Comm. 180, 2622 (2009) '

  LOGICAL :: linverse, lforceout, lrecip, esiste, equiva1, lcentral, leigen, ldrift , lnosym

  INTEGER :: nbranches, nd, nq, iprint, ndos, j, k, kmax, ndiff, nresto, ncycleseig
  INTEGER :: i, i1, i2, na, nb, nrr, nr, ntot, idos, ntyp, conta, nti

  INTEGER :: itertemp
  REAL(DP), DIMENSION(3) :: ptemp

  REAL(DP), DIMENSION(:), ALLOCATABLE :: free, free_q, eint, x, cv

  INTEGER, ALLOCATABLE :: list(:,:)

  REAL(DP), ALLOCATABLE :: dos(:)

  REAL(DP) :: temp(3), peso, pesotot, pesopart, delta, slope(3), &
       temperature, q(3), dq(3), tmp1(3), q0, e, dosin, dosend, dosstep, &
       dossmear, tv, tc, tv0, tc0, omega, fakt, rcut, zeropoint, fac, freq, &
       dx, dy, dz, eigsize

  COMPLEX(DP) :: zdotc

  TYPE(lattice)  :: latt
  TYPE(shell)    :: sh
  TYPE(dynmat)   :: dyn
  TYPE(symmetry) :: symm

  LOGICAL :: LDUM
  CHARACTER  :: CHARAC
  INTEGER :: IDUM, IERR, N
  REAL(dp) :: RDUM
  COMPLEX(dpc) :: CDUM


  !==========================================================================
  ! read in run time variables
  !==========================================================================

  write(*,*)'----------------------------------------------------------------'
  write(*,*) version
  write(*,*) citation
  write(*,*)'----------------------------------------------------------------'
  call reader( latt, dyn, temperature, ptemp, lforceout, lrecip, sh%rmax, &
       dosin, dosend, dosstep, dossmear, nd, iprint, symm%lsymm, nti, symm%symprec, &
       lcentral, eigsize, ncycleseig, leigen, ldrift, lnosym )


  ALLOCATE( free(nint(ptemp(3))), free_q(nint(ptemp(3))), eint(nint(ptemp(3))), &
       x(nint(ptemp(3))), cv(nint(ptemp(3))) ) 

  fac = evtoj   &        !  D^2 U eV  --> Joule
       /1.d-20  &        !  D^2 R A^2 --> meter 
       /amtokg  &        !  proton mass --> Kg
       /1.d24/twopi**2   !  cicles in THZ^2

  ndos = (dosend - dosin)/dosstep + 1
  allocate( dos(ndos) )
  dos = 0

  !======================================================================
  ! set lattice
  !=====================================================================
  call set_lattice ( latt, symm, iprint, sh%rmax )
  call generate_punti( symm, latt )

  !---------------------------
  ! Number of phonon branches
  !---------------------------
  nbranches = latt%natoms * 3

  !===========================
  ! Build a supercell ?
  !===========================
  if(latt%lsuper)then
     call build_supercell( latt )
     call get_displ( latt, symm, iprint, lcentral )
     stop
  endif

  !=============================================================================
  ! construct the shell of vectors in real space
  !=============================================================================
  call rshell( sh, latt ) 

  call numeric_force( latt, symm, sh, dyn, iprint, nti, lcentral, ldrift, lnosym )

! partial density of states
  open(15,file='INPHON',status='old')
  allocate (dyn%pdos_atom(latt%natoms))
  dyn%pdos_atom = .true.
  call rdatab( .false., 'INPHON', 15, 'pdos_atom', '=', '#', ';', 'L', &
       IDUM, RDUM, CDUM, dyn%pdos_atom, CHARAC, N, latt%natoms, IERR )
  if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
       ((IERR==0).and.(N<1))) then
     write(*,*)'error reading item ''pdos_atom'' from file INPHON.'
     stop
  endif
  if(any(dyn%pdos_atom)) dyn%usethis = .false.
  close(15)


  !--------------------------
  ! set the list of vectors
  !--------------------------
  ntot = 0
  do nrr=1,sh%nrm
     do na=1,latt%natoms
        do nb=1,latt%natoms
           tmp1 = sh%r(:,nrr) + latt%x(:,nb) - latt%x(:,na)
           if( vabs(tmp1) < sh%rmax ) ntot = ntot + 1
        enddo
     enddo
  enddo
  allocate( list(3,ntot) )
  ntot = 0
  do nrr=1,sh%nrm
     do na=1,latt%natoms
        do nb=1,latt%natoms
           tmp1 = sh%r(:,nrr) + latt%x(:,nb) - latt%x(:,na)
           if( vabs(tmp1) < sh%rmax ) then
              ntot = ntot + 1
              list(1,ntot) = na
              list(2,ntot) = nb
              list(3,ntot) = nrr
           endif
        enddo
     enddo
  enddo

  !==========================================================
  ! Write the force constant matrix in the file HARMONIC ?
  !==========================================================
  if(lforceout)then
     open(1,file='HARMONIC',status='unknown')
     write(1,*) sh%rmax,'   cutoff radius'
     write(1,*) ntot,'   number of vectors'
     do i = 1, ntot
        na = list(1,i)
        nb = list(2,i)
        nrr = list(3,i)
        write(1,'(''vector:'',4i20)') na, nb, nrr, dyn%weight(na,nb,nrr)
        write(1,'(3f20.12)') sh%r(:,nrr) + latt%x(:,nb) - latt%x(:,na)
        write(1,'(''force constant matrix:'')')
        write(1,'(3f20.12)') transpose(dyn%dmat(:,:,na,nb,nrr))
     enddo
     write(*,'(/''Force constant matrix file HARMONIC written''/)')
     close(1)
  endif

  !==========================================================================
  ! calculate the free energy as the sum of frequencies, you must provide
  ! the file 'QPOINTS' which contains the special points of the BZ 
  ! IN RECIPROCAL LATTICE COORDINATES
  !=========================================================================

  if( dyn%lfree ) then

     inquire(file='QPOINTS',exist=esiste)
     if(.not.esiste) stop 'cannot find file QPOINTS, use QA, QB and QC to generate one'
     open(1,file='QPOINTS',status='old')

     read(1,*) nd
     print'(/''Integrating frequencies...'')'
     print'(/''Using '',i6,'' q-points from file QPOINTS (in reciprocal space coordinates)''/)',nd
     pesotot = 0; 
     omega = 0; zeropoint = 0
     do itertemp = 1,nint(ptemp(3))
        eint(itertemp) = 0
        free_q(itertemp) = 0
        free(itertemp) = 0
        cv(itertemp) = 0
     end do
     conta = 0
     do nq = 1, nd
        read(1,*) q, peso
        pesotot = pesotot + peso
     enddo
     rewind(1)
     read(1,*) nd
     do nq = 1, nd
        read(1,*) q, peso
        q = twopi/latt%scale*matmul(latt%bg,q) ! in cartesian coordinates

        !==========================================================================
        ! construct the dynamical matrix
        !==========================================================================
        call dyn_mat( dyn, latt, sh, q, ntot, list )

        if(iprint>2)then
           print'(''Dynamical matrix:'')'
           do j=1,latt%natoms
              do k=1,latt%natoms
                 print'(2i3)',j,k
                 do i1=1,3
                    print'(6f13.8)', ( dyn%cmat( i1 + (j-1)*3, i2 + (k-1)*3 ), i2=1,3) 
                 enddo
              enddo
           enddo
        endif

        !==========================================================================
        ! diagonalize the dynamical matrix
        !==========================================================================
        call diagh( dyn%cmat, dyn%eig, 'L' )

        if(iprint>2)then
           do j = 1, nbranches
              if( dyn%eig(j) < 0.d0 ) then
                 freq = - sqrt(-dyn%eig(j))
              else
                 freq = sqrt(dyn%eig(j))
              endif
              print'(''Eigenvalue'',i3,''   Frequency (Thz) '',f20.12)',j, freq
              print'(''Eigenvector '',i3)',j
              do k=1,latt%natoms
                 print'(''atom '',i3)',k
                 do i1=1,3
                    print'(2f20.12)', dyn%cmat(i1+(k-1)*3,j)
                 enddo
              enddo
           enddo
        endif

        !=====================
        ! density of states
        !=====================
        do i1 = 1, nbranches
           pesopart = 0
           do k = 1, latt%natoms
              if(dyn%usethis(latt%ityp(k)).or.dyn%pdos_atom(k)) pesopart = pesopart + &
                   abs(zdotc(3,dyn%cmat(1+(k-1)*3,i1),1, &
                   dyn%cmat(1+(k-1)*3,i1),1))
           enddo
           idos = 0
           kmax = (dosend - dosin)/dosstep
           do k = 1, kmax + 1
              e = dosin + (k-1)*dosstep + dosstep/2
              idos = idos + 1
              if( dyn%eig(i1) >= 0 )then
                 if( sqrt(dyn%eig(i1)) > e - dosstep/2 .and. sqrt(dyn%eig(i1))< e + dosstep/2 )then
                    dos(idos) = dos(idos) + peso * pesopart
                 endif
              endif
           enddo

           !==========================================================================
           ! the phonon frequencies are the squares of the matrix eigenvalues
           !==========================================================================
           !-----------------------------------------------------------------------
           ! exclude the acustic mode at gamma and any imaginary frequency 
           !-----------------------------------------------------------------------
           if( dyn%eig(i1) < -zero ) then
              print*,'negative eigenvalues',dyn%eig(i1),i1
              !                  stop 'negative eigenvalues'
           elseif( dyn%eig(i1) > zero ) then
              conta = conta + 1
              fakt = 1.d4 * twopi**2 / evtoj 
              !----------------------------------------------------------
              ! omega is the average logaritmic frequency.
              ! nu = sqrt( dyn%eig(i) ) in THZ
              ! omega = 2pi * nu
              !--------------------------------------------------------------
              omega = omega + peso * pesopart * log( dyn%eig(i1) * fakt )

              zeropoint = zeropoint + peso * pesopart * &
                   0.5_dp * hplank * sqrt( dyn%eig(i1) )*1.d12 / evtoj

              do itertemp = 1,nint(ptemp(3))
                 temperature = (itertemp-1)*ptemp(2) + ptemp(1)

                 free(itertemp) = free(itertemp) + peso * pesopart * bolkev * temperature * &
                      log( hplank * sqrt( dyn%eig(i1) )*1.d12 / temperature / bolk )

                 free_q(itertemp) = free_q(itertemp)  + peso * pesopart * ( &
                      0.5_dp * hplank * sqrt( dyn%eig(i1) )*1.d12 / evtoj + &
                      bolkev * temperature * &
                      log( 1 - exp( -hplank * sqrt( dyn%eig(i1) )*1.d12 / temperature / bolk) ) )

                 eint(itertemp) = eint(itertemp) + &
                      peso * pesopart * (hplank/2)*sqrt( dyn%eig(i1) )*1.d12 * &
                      coth( hplank*sqrt( dyn%eig(i1) )*1.d12 / ( 2*bolk*Temperature ) ) 
              
                 x(itertemp) = hplank*sqrt( dyn%eig(i1) )*1.d12 / ( bolk*Temperature ) 

                 cv(itertemp) = cv(itertemp) + &
                      peso * pesopart * x(itertemp)**2*exp(x(itertemp))/(exp(x(itertemp)) - 1)**2
              end do

           endif
        enddo
     enddo

     omega = omega / pesotot / nbranches
     zeropoint = zeropoint / pesotot

     do itertemp = 1,nint(ptemp(3))
        free(itertemp) = free(itertemp) / pesotot
        free_q(itertemp) = free_q(itertemp) / pesotot
        cv(itertemp) = cv(itertemp) / pesotot
        eint(itertemp) = eint(itertemp) / pesotot / evtoj
     end do

     if( conta /= nbranches*nd ) &
          write(*,'(/''WARNING, Found '',i8, &
          & '' frequencies different from zero, out of a total of '',i8/)') conta, nbranches*nd

     !------------------------------------------------------------
     ! print out thermodynamic quantities and writes file THERMO 
     !------------------------------------------------------------
     write(6,'(/''Your primitive cell contains N ='',i3,'' atoms''/)') latt%natoms

     write(6,'(''Zero point energy           '',f15.5,''  (eV/cell)'')') zeropoint

     open(23,FILE='THERMO',STATUS='replace')
     write(23,'("#  T(K)     E(eV/cell)    F(eV/cell)    Fc(eV/cell)   S(kB/cell)    Cv(kB/cell)")')
     do itertemp = 1,nint(ptemp(3))
      temperature = (itertemp-1)*ptemp(2) + ptemp(1)

      write(6,'(''Temperature                   '',f10.2,'' K'')') temperature
      write(6,'(''Free energy                 '',f15.5,''  (eV/cell)'', &
           & f12.5,''  (eV/atom)'')') free_q(itertemp), free_q(itertemp)/latt%natoms 
      write(6,'(''Free energy(classical_limit)'',f15.5,''  (eV/cell)'', &
           & f12.5,''  (eV/atom)'')') free(itertemp), free(itertemp)/latt%natoms
      write(6,'(''Internal Energy             '',f15.5,''  (eV/cell)'')') eint(itertemp)
      write(6,'(''3 N kB T                    '',f15.5,''  (eV/cell)'')') &
           3*temperature*bolkev*latt%natoms
      write(6,'(''Cv                          '',f15.5,''  (kB/cell)'')') cv(itertemp)
      write(6,'(''S                           '',f15.5,''  (kB/cell)''/)') &
           (eint(itertemp) - free_q(itertemp)) / (temperature*bolkev)
 
      write (23,'(f8.2,5f14.8)') temperature,eint(itertemp),free_q(itertemp),free(itertemp), &
       (eint(itertemp)-free_q(itertemp))/(temperature*bolkev),cv(itertemp)
     end do
     close(23)
     close(1)

     !---------------------------------
     ! writes density of states files
     !---------------------------------
     open(1,file='DOS',status='unknown')
     open(2,file='DOS.cm',status='unknown')
     open(3,file='DOS.meV',status='unknown')
     idos = 0
     dos = dos / pesotot / nbranches / dosstep
     call smooth( dos, dossmear, ndos, dosin, dosstep )
     write(6,'(/''Density of states integral = '',f10.5)') sum(dos)*dosstep
     if( abs(sum(dos)*dosstep-1.0_dp) > 1.d-6 ) then
        print*,'renormalizing DOS..'
        dos = dos / ( sum(dos)*dosstep )
     endif        

     kmax = (dosend - dosin)/dosstep
     do k = 1, kmax + 1
        e = dosin + (k-1)*dosstep + dosstep/2
        idos = idos + 1
        write(1,'(2f20.5)') e, dos(idos)
        write(2,'(2f20.5)') e*convert_thz_to_cm, dos(idos)/convert_thz_to_cm
        write(3,'(2f20.5)') e*convert_thz_to_mev, dos(idos)/convert_thz_to_mev
     enddo

     close(1); close(2); close(3)
     stop

  endif

  !============================================================================
  ! calculate phonon dispersions, written in files FREQ and FREQ?
  !===========================================================================

  if(nbranches < 48)then
    open(2,file='FREQ',status='unknown')
    open(3,file='FREQ.cm',status='unknown')
  else
    do i=1,nbranches/48
       if( i <=9 ) then
          write(fl,'(i1)')i
          open(10+i,file='FREQ'//fl,status='unknown')
       else
          write(fl1,'(i2)')i
          open(10+i,file='FREQ'//fl1,status='unknown')
       endif
    enddo
    nresto = mod(nbranches,48)
    if(nresto > 0) then
       if( i <=9 ) then
          write(fl,'(i1)')i
          open(10+i,file='FREQ'//fl,status='unknown')
       else
          write(fl1,'(i2)')i
          open(10+i,file='FREQ'//fl1,status='unknown')
       endif
    endif
  endif

  q=0
  if(leigen)then
     call dyn_mat( dyn, latt, sh, q, ntot, list )
     call diagh( dyn%cmat, dyn%eig, 'L' )
     open(4,file='EIGEN.axsf',status='unknown')
     write(4,'(''ANIMSTEPS '',i4)') nbranches
     write(4,'(''CRYSTAL'')')
     write(4,'(''PRIMVEC'')')
     write(4,'(3f10.6)')latt%at(1,:)*latt%scale
     write(4,'(3f10.6)')latt%at(2,:)*latt%scale
     write(4,'(3f10.6)')latt%at(3,:)*latt%scale
     do j=1,nbranches
        write(4,'(''PRIMCOORD'',i4)')j
        write(4,'(i4,i2)')latt%natoms, 1
        i2 = 0
        do k = 1, latt%ntypes
           do i1 = 1, latt%nions(k)
              i2 = i2 + 1
              write(4,'(a3,3f9.5,3f12.7)')latt%name(k),latt%x(:,i2),&
                   dble(dyn%cmat(1+(i2-1)*3:3+(i2-1)*3,j))*nbranches*&
                   eigsize/20.0/sqrt(latt%mass(k))
           enddo
        enddo
        if(j<10)then
           nd_nmbr='00'
           write(nd_nmbr(3:3),'(i1)')j
        elseif(j<100)then
           nd_nmbr='0'
           write(nd_nmbr(2:3),'(i2)')j
        else
           write(nd_nmbr,'(i3)')j
        endif
        open(7,file='MODE'//nd_nmbr//'.axsf',status='unknown')
        write(7,'(''ANIMSTEPS '',i4)') 20*ncycleseig
        write(7,'(''CRYSTAL'')')
        write(7,'(''PRIMVEC'')')
        write(7,'(3f10.6)')latt%at(1,:)*latt%scale
        write(7,'(3f10.6)')latt%at(2,:)*latt%scale
        write(7,'(3f10.6)')latt%at(3,:)*latt%scale
        do na = 1, 20*ncycleseig
           write(7,'(''PRIMCOORD'',i4)')na
           write(7,'(i4,i2)')latt%natoms, 1
           i2 = 0
           do k = 1, latt%ntypes
              do i1 = 1, latt%nions(k)
                 i2 = i2 + 1
                 write(7,'(a3,3f9.5,3f12.7)')latt%name(k),latt%x(:,i2)+ &
                   dble(dyn%cmat(1+(i2-1)*3:3+(i2-1)*3,j))*nbranches*&
                   eigsize/10.0*sin((na-1)/20.0*2.*3.14159265)/sqrt(latt%mass(k))
              enddo
           enddo
        enddo
        close(7)
     enddo
     close(4)
     open(4,file='EIGEN.xyz',status='unknown')
     do j=1,nbranches
        write(4,'(i4,i2)')latt%natoms
        write(4,'(''MODE'',i4)')j
        i2 = 0
        do k = 1, latt%ntypes
           do i1 = 1, latt%nions(k)
              i2 = i2 + 1
              write(4,'(a3,3f9.5,3f12.7)')latt%name(k),latt%x(:,i2),&
                   dble(dyn%cmat(1+(i2-1)*3:3+(i2-1)*3,j))*nbranches*&
                   eigsize/20.0/sqrt(latt%mass(k))
           enddo
        enddo
     enddo
     close(4)
  endif

  write(*,'(/''Distances:'')')
  do nq = 1, nd
     print'(''point '',i4,f10.6,10x,3f10.6)', nq, q0, dyn%qi(:,nq)
     if( dyn%npoints > 1 ) then
        dq = ( dyn%qf(:,nq) - dyn%qi(:,nq) ) / ( dyn%npoints - 1 )      
     else
        dq = 0
     endif
     if( lrecip ) dq = matmul(latt%bg,dq)
     q0 = q0 - vabs( dq )
     do i = 1, dyn%npoints
        q = dyn%qi(:,nq) 
        if( lrecip ) then
           q = twopi/latt%scale*matmul(latt%bg,q) + &
                (i-1)*dq*twopi/latt%scale
        else
! this only works if latt%scale is the lattice parameter
           q = twopi/latt%scale*q + (i-1)*dq*twopi/latt%scale
        endif
        if(iprint > 1)then
           write(*,'(/''Point #'',i3,'' q-vector: '',3f10.5)') i, q*latt%scale/twopi
        endif
        q0 = q0 + vabs(dq)
        !==========================================================================
        ! construct the dynamical matrix
        !==========================================================================
        call dyn_mat( dyn, latt, sh, q, ntot, list )

        if(iprint>2)then
           print'(''Dynamical matrix:'')'
           do j=1,latt%natoms
              do k=1,latt%natoms
                 print'(2i3)',j,k
                 do i1=1,3
                    print'(6f13.8)', &
                         ( dyn%cmat( i1 + (j-1)*3, i2 + (k-1)*3 )/fac, i2=1,3) 
                 enddo
              enddo
           enddo
        endif
        !==========================================================================
        ! diagonalize the dynamical matrix
        !==========================================================================
        call diagh( dyn%cmat, dyn%eig, 'L' )

        do i1=1,nbranches

           !==========================================================================
           ! the phonon frequencies are the square roots of the matrix eigenvalues
           !==========================================================================
           if( dyn%eig(i1) < 0 ) then
              dyn%eig(i1) = -sqrt(-dyn%eig(i1))
           else
              dyn%eig(i1) = sqrt(dyn%eig(i1))
           endif
        enddo
        if(nbranches < 48)then
          write(2,'(f12.8,48f12.5)')q0,dyn%eig
          write(3,'(f12.8,48f12.5)')q0,dyn%eig*convert_thz_to_cm
        else
          do i1=1,nbranches/48
             write(10+i1,'(f12.8,48f12.5)')q0,dyn%eig( 1 + (i1 - 1)*48 : i1*48 )
          enddo
          if(nresto>0)write(10+i1,'(f12.8,48f12.5)')q0,dyn%eig( 1 + (i1 - 1)*48 : (i1-1)*48+nresto )
        endif

        !---------------------------------------------------------
        ! estimates the slopes of the first three branches
        ! this only works if latt%scale is the lattice parameter
        !---------------------------------------------------------
        if( iprint > 0 ) then
           if( i == 1 .or. i == dyn%npoints - 1 ) then
              slope = dyn%eig(1:3)
           elseif( i == 2 .or. i == dyn%npoints ) then
              slope = ( dyn%eig(1:3) - slope ) / vabs( dq/latt%scale ) * 100
              write(*,'(''Slopes (m/s) : '',3f10.1)') slope
           endif
        endif

        !==============================================================================
        ! density of states is calculated also here, but this is not the right way
        ! of calculating the DOS
        !===============================================================================
        do i1 = 1, nbranches
           idos = 0
           do k = 1, ndos
              e = dosin + (k-1)*dosstep
              idos = idos + 1
              if( dyn%eig(i1) >= 0 )then
                 if( dyn%eig(i1) > e - dosstep/2 .and. dyn%eig(i1) < e + dosstep/2 )then
                    dos(idos) = dos(idos) + 1
                 endif
              endif
           enddo
        enddo

        if(iprint>2)then
           do j=1,nbranches
              print'(''Eigenvalue '',i3,'' : '',f12.5)',j,dyn%eig(j) 
              print'(''Eigenvector '',i3,'' (divided by sqrt of atomic masses) '')',j
              do k=1,latt%natoms
                 print'(''atom '',i3)',k
                 do i1=1,3
                    print'(2f14.8)', dyn%cmat(i1+(k-1)*3,j)/sqrt(latt%mass(latt%ityp(k)))
                 enddo
              enddo
           enddo
        endif

     enddo
  enddo

  print'(''end point '',f10.6,10x,3f10.6)',  q0, dyn%qf(:,nd)

  open(4,file='DOS',status='unknown')
  open(7,file='DOS.cm',status='unknown')
  open(8,file='DOS.meV',status='unknown')
  idos = 0
  dos = dos / nd / dyn%npoints / nbranches / dosstep
  call smooth( dos, dossmear, ndos, dosin, dosstep )
  write(6,'(/''Density of states integral = '',f10.5)') sum(dos)*dosstep
  if( abs(sum(dos)*dosstep-1.0_dp) > 1.d-6 ) then
     print*,'renormalizing DOS..'
     dos = dos / ( sum(dos)*dosstep )
  endif
  do k = 1, ndos
     e = dosin + (k-1)*dosstep
     idos = idos + 1
     write(4,'(2f20.5)') e, dos(idos)
     write(7,'(2f20.5)') e*convert_thz_to_cm, dos(idos)/convert_thz_to_cm
     write(8,'(2f20.5)') e*convert_thz_to_meV, dos(idos)/convert_thz_to_mev
  enddo

END PROGRAM phon
