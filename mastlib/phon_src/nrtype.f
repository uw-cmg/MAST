MODULE nrtype
!-----------------------------------------------------
  ! Data type definitions and constant settings
  !---------------------------------------------------

  INTEGER, PARAMETER :: I4B = selected_int_kind(9) 
  INTEGER, PARAMETER :: I2B = selected_int_kind(4) 
  INTEGER, PARAMETER :: I1B = selected_int_kind(2) 

  INTEGER, PARAMETER :: SP = kind(1.0) 
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(10) 

  INTEGER, PARAMETER :: SPC = kind( (1.0,1.0) ) 
  INTEGER, PARAMETER :: DPC = kind( (1.0_dp,1.0_dp) ) 

  INTEGER, PARAMETER :: LGT = kind( .true. ) 

  REAL(DP), PARAMETER :: PI = 3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: TWOPI = 2._dp*PI

  COMPLEX(DPC), PARAMETER :: im = (0._dp,1._dp)

  !  Some important Parameters, to convert to a.u.
  !  - AUTOA  = 1. a.u. in Angstroem
  !  - RYTOEV = 1 Ry in Ev
  !  - EVTOJ  = 1 eV in Joule
  !  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
  !  - BOLKEV = Boltzmanns constant in eV/K
  !  - BOLK   = Boltzmanns constant in Joule/K

  REAL(DP), PARAMETER :: EVTOJ=1.60217733E-19_dp,AMTOKG=1.6605402E-27_dp, &
       BOLKEV=8.6173857E-5_dp,BOLK=BOLKEV*EVTOJ, HPLANK=6.6262E-34_dp, &
       zero = 1.d-6, convert_thz_to_cm = 33.357, &
       convert_thz_to_meV = 4.1357

END MODULE nrtype

