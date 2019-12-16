! Author: Simon Mason
Module data_cpt_constants
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Parameters
!
! Integer parameters
  Integer, Parameter, Public :: ia_cca = 1             ! - CCA identifier -
  Integer, Parameter, Public :: ia_pcr = 2             ! - PCR identifier -
  Integer, Parameter, Public :: ia_mlr = 3             ! - MLR identifier -
  Integer, Parameter, Public :: ia_gcm = 4             ! - GCM identifier -
  Integer, Parameter, Public :: ia_pfv = 5             ! - PFV identifier -
  Integer, Parameter, Public :: ia_for = 11            ! - forecast identifier -
  Integer, Parameter, Public :: ia_dia = 12            ! - diagnostics identifier -
  Integer, Parameter, Public :: ipm = 100              ! - maximum percentage of missing values -
  Integer, Parameter, Public :: mnt = 4                ! - minimum length of training period -
  Integer, Parameter, Public :: mre = 8                ! - maximum reasonable number of EOF modes -
  Integer, Parameter, Public :: na = 5                 ! - number of analytical methods -
  Integer, Parameter, Public :: na10 = 3               ! - number of analytical methods up to CPT version 10 -
  Integer, Parameter, Public :: nbs = 2                ! - number of probability bin options -
  Integer, Parameter, Public :: ncorr = 3              ! - number of correlations -
  Integer, Parameter, Public :: ncos = 4               ! - number of correction options -
  Integer, Parameter, Public :: ndcon = 8              ! - number of continuous deterministic scores -
  Integer, Parameter, Public :: ndcat = 8              ! - number of categorised deterministic scores -
  Integer, Parameter, Public :: nel = 2                ! - number of expertise levels -
  Integer, Parameter, Public :: nensc = 3              ! - number of ensemble forecast options -
  Integer, Parameter, Public :: nep = 1000             ! - number of points on exceedance probability curve -
  Integer, Parameter, Public :: nev = 3                ! - number of error-variance options -
  Integer, Parameter, Public :: n_file = 3             ! - number of input files -
  Integer, Parameter, Public :: ng = 3                 ! - number of categories -
  Integer, Parameter, Public :: nmcs = 6               ! - number of model combination options -
  Integer, Parameter, Public :: nmo = 4                ! - number of missing value replacement options -
  Integer, Parameter, Public :: npall = 10             ! - number of all-category probabilistic scores -
  Integer, Parameter, Public :: npcat = 9              ! - number of category-specific probabilistic scores -
  Integer, Parameter, Public :: npo = 3                ! - number of principal component calculation options -
  Integer, Parameter, Public :: nstruct = 3            ! - number of data structures -
  Integer, Parameter, Public :: nstd = 5               ! - number of standardization options -
  Integer, Parameter, Public :: nthr = 3               ! - number of threshold options -
  Integer, Parameter, Public :: ntl = 5                ! - number of time options -
  Integer, Parameter, Public :: nto = 2                ! - number of tranformation options -
#if GUI == 0
  Integer, Parameter, Public :: nueo = 3               ! - number of user-error options -
#endif
  Integer, Parameter, Public :: nwr = 3                ! - number of weather roulette options -
!
  Integer, Parameter, Public :: m_file = n_file + 3    ! - total number of input files -
  Integer, Parameter, Public :: nscore = ndcon + ndcat ! - number of scores -
  Integer, Parameter, Public :: nps = npall + npcat    ! - number of probabilistic scores -
  Integer, Parameter, Public :: mps = npall + npcat*ng ! - total number of probabilistic scores -
  Integer, Parameter, Public :: nts = ng - 1           ! - number of thresholds -
  Integer, Parameter, Public :: nskill = 9 + nts       ! - number of skill scores -
!
  Integer, Dimension(   nps), Parameter, Public :: ipsi = &  ! - probabilistic score identifiers -
     (/101, 111, 121, 122, 125, 131, 141, 151, 152, 153, &
     & 201, 202, 203, 221, 222, 223, 224, 225, 231/)
  Integer, Dimension(   nbs), Parameter, Public :: nbss = &  ! - number of probability bins -
     (/11, 21/)
!
! Character parameters
  Character(Len=*), Dimension(  ncos), Parameter, Public :: ccorrects = & ! - correction options -
     (/'No correction                   ', 'Correct mean biases             ', 'Correct mean and variance biases', &
       'Correct for skill               '/)
  Character(Len=*), Dimension( nensc), Parameter, Public :: censc = &     ! - ensemble forecasting options -
     (/'None               ', 'Ensemble counting  ', 'Fitted distribution'/)
  Character(Len=*), Dimension(m_file), Parameter, Public :: c_file = &    ! - input file types -
     (/'X', 'Y', 'Z', '1', '2', '3'/)
  Character(Len=*), Dimension( nmcs), Parameter, Public :: cmcs = &       ! - model combination options -
     (/'Average uncalibrated  ', 'Average calibrated    ', 'Average recalibrated  ', 'Skill-weighted average', &
       'Best model by location', 'Best overall model    '/)
  Character(Len=*), Dimension(  nstd), Parameter, Public :: cstds = &     ! - standardization options -
     (/'No standardization      ', 'Anomalies               ', 'Standardized anomalies  ', 'Standardized index (SPI)', &
       '% of average            '/)
  Character(Len=*), Dimension(  nthr), Parameter, Public :: cthrs = &     ! - threshold options -
     (/'Climatological probabilities', 'Absolute thresholds         ', 'Reference dates             '/)
  Character(Len=*), Dimension(   ntl), Parameter, Public :: ctls = &      ! - time options -
     (/'Gregorian calendar (British) ', &
       'Julian calendar              ', &
       'Ignore leap years            ', &
       'All months 30 days           ', &
       'Gregorian calendar (Catholic)'/)
  Character(Len=*), Dimension(   nto), Parameter, Public :: ctos = &      ! - transformation options -
     (/'Empirical distribution', 'Gamma distribution    '/)
!
! Logical parameters
  Logical, Dimension(nscore), Parameter, Public :: lposcore = & ! - positively oriented score? -
     (/  .true.,  .true.,  .true.,  .true.,  .true., .false., .false., .false., &
         .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true./)
  Logical, Dimension(nscore), Parameter, Public :: lbmscore = & ! - beginner mode score? -
     (/  .true., .false.,  .true., .false., .false., .false., .false.,  .true., &
        .false., .false., .false., .false., .false.,  .true.,  .true.,  .true./)
!
End Module data_cpt_constants
