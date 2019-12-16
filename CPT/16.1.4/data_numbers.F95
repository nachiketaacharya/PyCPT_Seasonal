! Author: Simon Mason
Module data_numbers
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Parameters
!
! Data kind numbers
  Integer, Parameter, Public :: dp = Kind(1.0d0) ! - double precision -
  Integer, Parameter, Public :: sp = Kind(1.0e0) ! - single precision -
#if DP == 1
  Integer, Parameter, Public :: rp = dp          ! - default precision -
  Integer, Parameter, Public :: ap = sp          ! - alternative precision -
#else
  Integer, Parameter, Public :: rp = sp          ! - default precision -
  Integer, Parameter, Public :: ap = dp          ! - alternative precision -
#endif
!
! Whole numbers
  Real(Kind=rp), Parameter, Public :: zero =          0.0_rp ! - zero -
  Real(Kind=rp), Parameter, Public :: one =           1.0_rp ! - one -
  Real(Kind=rp), Parameter, Public :: two =           2.0_rp ! - two -
  Real(Kind=rp), Parameter, Public :: three =         3.0_rp ! - three -
  Real(Kind=rp), Parameter, Public :: four =          4.0_rp ! - four -
  Real(Kind=rp), Parameter, Public :: five =          5.0_rp ! - five -
  Real(Kind=rp), Parameter, Public :: six =           6.0_rp ! - six -
  Real(Kind=rp), Parameter, Public :: seven =         7.0_rp ! - seven -
  Real(Kind=rp), Parameter, Public :: eight =         8.0_rp ! - eight -
  Real(Kind=rp), Parameter, Public :: nine =          9.0_rp ! - nine -
  Real(Kind=rp), Parameter, Public :: ten =          10.0_rp ! - ten -
  Real(Kind=rp), Parameter, Public :: eleven =       11.0_rp ! - eleven -
  Real(Kind=rp), Parameter, Public :: twelve =       12.0_rp ! - twelve -
  Real(Kind=rp), Parameter, Public :: sixten =       16.0_rp ! - sixteen -
  Real(Kind=rp), Parameter, Public :: twenty =       20.0_rp ! - twenty -
  Real(Kind=rp), Parameter, Public :: thirty =       30.0_rp ! - thirty -
  Real(Kind=rp), Parameter, Public :: fifty =        50.0_rp ! - fifty -
  Real(Kind=rp), Parameter, Public :: oneh =        100.0_rp ! - one hundred -
  Real(Kind=rp), Parameter, Public :: onet =       1000.0_rp ! - one thousand -
  Real(Kind=rp), Parameter, Public :: tent =      10000.0_rp ! - ten thousand -
  Real(Kind=rp), Parameter, Public :: onem =    1000000.0_rp ! - one million -
  Real(Kind=rp), Parameter, Public :: onetm =  10000000.0_rp ! - ten million -
  Real(Kind=rp), Parameter, Public :: onehm = 100000000.0_rp ! - one hundred million -
  Real(Kind=rp), Parameter, Public :: oneb = 1000000000.0_rp ! - one billion -
!
! - single precision - 
  Real(Kind=sp), Parameter, Public :: zero_sp =   0.0_sp ! - zero -
  Real(Kind=sp), Parameter, Public :: one_sp =    1.0_sp ! - one -
  Real(Kind=sp), Parameter, Public :: ten_sp =   10.0_sp ! - ten -
  Real(Kind=sp), Parameter, Public :: oneh_sp = 100.0_sp ! - one hundred -
!
! - double precision - 
  Real(Kind=dp), Parameter, Public :: zero_dp = 0.0_dp ! - zero -
  Real(Kind=dp), Parameter, Public :: one_dp =  1.0_dp ! - one -
  Real(Kind=dp), Parameter, Public :: ten_dp = 10.0_dp ! - ten -
!
! Fractions
  Real(Kind=rp), Parameter, Public :: onettth = 0.0001_rp ! - one ten thousandth -
  Real(Kind=rp), Parameter, Public :: onetth =   0.001_rp ! - one thousandth -
  Real(Kind=rp), Parameter, Public :: onehth =    0.01_rp ! - one hundredth -
  Real(Kind=rp), Parameter, Public :: tenth =      0.1_rp ! - one tenth -
  Real(Kind=rp), Parameter, Public :: third = one/three   ! - one third -
  Real(Kind=rp), Parameter, Public :: half =       0.5_rp ! - half -
  Real(Kind=rp), Parameter, Public :: onep5 =      1.5_rp ! - one and a half -
!
! Angles
  Real(Kind=rp), Parameter, Public :: r90 =   90.0_rp ! - 90 degrees -
  Real(Kind=rp), Parameter, Public :: r180 = 180.0_rp ! - 180 degrees -
  Real(Kind=rp), Parameter, Public :: r270 = 270.0_rp ! - 270 degrees -
  Real(Kind=rp), Parameter, Public :: r360 = 360.0_rp ! - 360 degrees -
!
! Irrational numbers
  Real(Kind=rp), Parameter, Public :: pi = 3.141592653589793_rp ! - pi -
  Real(Kind=rp), Parameter, Public :: hpi = pi/two              ! - half pi -
!
! Ordinal numbers
  Character(Len = 2), Parameter, Public :: cordns(0:3) = & ! - ordinal number suffixes -
     (/'th', 'st', 'nd', 'rd'/)
!
! Geographical parameters
  Real(Kind=rp), Parameter, Public :: re = 6372795.0_rp ! - radius of Earth (in metres) -
!
! Physical parameters
  Real(Kind=rp), Parameter, Public :: h2okgm3 = 998.0_rp ! - mass of water per cubic metre - 
!
! Temperature parameters
  Real(Kind=rp), Parameter, Public :: zero_K = 273.15_rp ! - zero degrees in Kelvin -
  Real(Kind=rp), Parameter, Public :: zero_F =  32.00_rp ! - zero degrees in Farenheit -
!
! Rounding errors
  Real(Kind=rp), Parameter, Public :: tolp = 0.0201_rp ! - probability rounding tolerance -
!
  Real(Kind=sp), Parameter, Public :: tolv = 0.0001_sp ! - version-number rounding tolerance -
!
! Digits
  Character(Len=*), Parameter, Public :: cdigits = '0123456789' ! - digits -
!
! Scalars
!
! Integer scalars
  Integer, Public :: ihuge  ! - huge integer -
!
! Real scalars
  Real(Kind=rp), Public :: ln2      ! - log of 2 -
  Real(Kind=rp), Public :: ln2pi    ! - log of 2*pi -
  Real(Kind=rp), Public :: lnsqr2pi ! - log of square root of 2*pi -
  Real(Kind=rp), Public :: rhuge    ! - huge real -
  Real(Kind=rp), Public :: sqrt2    ! - square root of two -
  Real(Kind=rp), Public :: sqrt3    ! - square root of three -
  Real(Kind=rp), Public :: sqrpi    ! - square root of pi -
  Real(Kind=rp), Public :: sqr2pi   ! - square root of 2*pi -
  Real(Kind=rp), Public :: base     ! - base -
  Real(Kind=rp), Public :: elim     ! - maximum exponent -
  Real(Kind=rp), Public :: eps      ! - machine precision -
  Real(Kind=rp), Public :: bt       ! - machine precision -
  Real(Kind=rp), Public :: sfmin    ! - safe minimum -
  Real(Kind=rp), Public :: sfmax    ! - safe maximum -
  Real(Kind=rp), Public :: smlnum   ! - small number -
  Real(Kind=rp), Public :: bignum   ! - big number -
  Real(Kind=rp), Public :: tol      ! - tolerance for identifying probable floating errors -
!
End Module data_numbers
