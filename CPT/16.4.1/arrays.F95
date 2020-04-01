! Author: Simon Mason
Module arrays
!
! Declarations
!
! Modules
  Use data_numbers, Only: ap, rp
  Use maths,        Only: interval
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: convert_units, get_cvsample, get_flag, insertion_sort, moments, proj_read_results, proj_write_results, &
            rank_data, set_flag
!
! Arrays
!
! Integer arrays
  Integer, Dimension(:), Allocatable, Public :: iusex ! - indices of used X variables -
  Integer, Dimension(:), Allocatable, Public :: iusey ! - indices of used Y variables -
!
! Real arrays
  Real(Kind=rp), Dimension(:), Allocatable, Public :: c   ! - climatological probabilities of exceedance -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: e   ! - empirical probabilities of exceedance -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: f   ! - probabilities of exceedance given forecast -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: tex ! - probability of exceedance thresholds -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: xm  ! - means of explanatory variables -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: xsd ! - standard deviations of explanatory variables -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: ym  ! - means of response variables -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: ysd ! - standard deviations of response variables -
!
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: b      ! - regression coefficients -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: b2     ! - regression coefficients -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: b0     ! - regression constants -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: fcastt ! - transformed forecast values -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: pev    ! - prediction error variance -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: xvp    ! - predictors during forecast period -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yc     ! - cross-validated training period response data -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yfit   ! - fitted values of response variables -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yhat   ! - cross-validated hindcasts -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yhatt  ! - transformed cross-validated hindcasts -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yret   ! - retroactive predictions -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yrett  ! - transformed retroactive predictions -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: yt     ! - transformed response data -
!
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: fcast ! - forecast values -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: vcopy ! - copy of data -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: x     ! - explanatory variables -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: xc    ! - cross-validated training period explanatory data -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: xiny  ! - predictors interpolated to predictand locations -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: xk    ! - copy of explanatory data -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: y     ! - response variables -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: yopt  ! - cross-validated hindcasts for optimization -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: ypc   ! - persisted component response variables -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: z     ! - updated explanatory variables -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: ziny  ! - updated predictors interpolated to predictand locations -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: zyp   ! - persistence component for forecasts -
!
! Derived-type arrays
  Type(interval), Dimension(:,:), Allocatable, Public :: fpls ! - prediction limits -
  Type(interval), Dimension(:,:), Allocatable, Public :: yrpls ! - retroactive forecast prediction limits -
!
! Workspace
!
! Workspace dimensions
  Integer, Public :: liwk ! - integer workspace dimensions -
  Integer, Public :: lrwk ! - real workspace dimensions -
!
! Integer workspace
  Integer, Dimension(:), Allocatable, Public :: iwk ! - integer workspace -
!
! Real workspace
! - single precision -
  Real(Kind=ap), Dimension(:), Allocatable, Public :: awk ! - alternative precision workspace -
! - double precision -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: rwk ! - default precision workspace -
!
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: rwk2 ! - default precision workspace -
!
! Logical workspace
  Logical, Dimension(:), Allocatable, Public :: lwk ! - logical workspace -
!
! Interfaces
!
! Generic interfaces
  Interface insertion_sort
   Module Procedure insertion_sort_1
   Module Procedure insertion_sort_2
  End Interface insertion_sort
!
Contains
!
!
 Pure Subroutine moments (m, n, v, vm, vsd)
!
! Calculates arithmetic column means and standard deviations of a 2-D array
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: m ! - number of variables -
  Integer, Intent(In) :: n ! - number of cases -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v ! - variables; shape(m,n) -
!
! Output arrays
  Real(Kind=rp), Dimension(:), Intent(Out) :: vm  ! - means; shape(m) -
  Real(Kind=rp), Dimension(:), Intent(Out) :: vsd ! - standard deviations; shape(m) -
!
! Locals
!
! Local scalars
  Integer :: i ! - indices -
!
  Real(Kind=rp) :: df  ! - number of cases -
  Real(Kind=rp) :: dof ! - degrees of freedom -
!
! Functions and subroutines
!
! Intrinsic functions
  Intrinsic Real, Sqrt, Sum
!
! Executable Statements
!
! Calculate degrees of freedom
  df = Real(n, Kind=rp)
  dof = Real(n-1, Kind=rp)
!
! Calculate means and standard deviations
  ForAll (i=1:m)
     vm(i) = Sum(v(i,1:n))/df
     vsd(i) = Sqrt(Sum((v(i,1:n) - vm(i))**2)/dof)
  End ForAll
!
  Return
 End Subroutine moments
!
!
!
 Subroutine get_cvsample (hcw, icv, v1, v1cv, &
            v2, v2cv, v3, v3cv)
!
! Generates a cross-validated sample
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: hcw ! - half cross-validation width -
  Integer, Intent(In) :: icv ! - cross-validation index -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v1 ! - sample data -
!
! - optional input arrays -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: v2 ! - additional 2-D sample data -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In), Optional :: v3 ! - additional 3-D sample data -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: v1cv ! - cross-validated sample data -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: v2cv ! - additional 2-D cross-validated sample data -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(Out), Optional :: v3cv ! - additional 3-D cross-validated sample data -
!
! Locals
!
! Local scalars
  Integer :: i1 ! - start of cross-validation sample -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Cshift, Present
!
! Executable Statements
!
! Identify start of cross-validation sample
  i1 = hcw + icv - 1
!
! Cross-validate
  v1cv(:,:) = Cshift(v1(:,:), i1, 2)
  If (Present(v2)) v2cv(:,:) = Cshift(v2(:,:), i1, 2)
  If (Present(v3)) v3cv(:,:,:) = Cshift(v3(:,:,:), i1, 2)
!
  Return
 End Subroutine get_cvsample
!
!
!
 Subroutine insertion_sort_1 (n, v, order)
!
! Sorts data in ascending or descending order using insertion sort algorithm
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of cases -
!
  Character(Len=1), Intent(In) :: order ! - order -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:), Intent(InOut) :: v ! - data to be sorted -
!
! Locals
!
! Local scalars
  Integer :: i,j ! - indices -
!
  Real(Kind=rp) :: swp ! - swapped value -
!
! Executable Statements
!
! Sort data
! - ascending -
  Select Case (order)
   Case ('A', 'a')
     Do i = 2, n
        swp = v(i)
        Do j = 1, i-1
           If (v(i) < v(j)) Exit
        End Do
        If (j < i) Then
           v(j+1:i) = v(j:i-1)
           v(j) = swp
        End If
     End Do
! - descending -
   Case ('D', 'd')
     Do i = 2, n
        swp = v(i)
        Do j = 1, i-1
           If (v(i) > v(j)) Exit
        End Do
        If (j < i) Then
           v(j+1:i) = v(j:i-1)
           v(j) = swp
        End If
     End Do
  End Select
!
  Return
 End Subroutine insertion_sort_1
!
!
!
 Subroutine insertion_sort_2 (m, n, v, order)
!
! Sorts 2-D data in ascending or descending order of second dimension using insertion sort algorithm
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: m ! - number of locations -
  Integer, Intent(In) :: n ! - number of cases -
!
  Character(Len=1), Intent(In) :: order ! - order -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: v ! - data to be sorted -
!
! Locals
!
! Local scalars
  Integer :: i, j ! - indices -
  Integer :: i1   ! - indices -
!
  Real(Kind=rp) :: swp ! - swapped value -
!
! Executable Statements
!
! Sort data
! - ascending -
  Select Case (order)
   Case ('A', 'a')
     Do i1 = 1, m
        Do i = 2, n
           swp = v(i1,i)
           Do j = 1, i-1
              If (v(i1,i) < v(i1,j)) Exit
           End Do
           If (j < i) Then
              v(i1,j+1:i) = v(i1,j:i-1)
              v(i1,j) = swp
           End If
        End Do
     End Do
! - descending -
   Case ('D', 'd')
     Do i1 = 1, m
        Do i = 2, n
           swp = v(i1,i)
           Do j = 1, i-1
              If (v(i1,i) > v(i1,j)) Exit
           End Do
           If (j < i) Then
              v(i1,j+1:i) = v(i1,j:i-1)
              v(i1,j) = swp
           End If
        End Do
     End Do
  End Select
!
  Return
 End Subroutine insertion_sort_2
!
!
!
 Subroutine rank_data (v, n, order, r)
!
! Ranks data in ascending or descending order
!
! Modules
  Use data_numbers, Only: half, one
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of cases -
!
  Character(Len=1), Intent(In) :: order ! - order -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: v ! - data to be ranked -
!
! Output arrays
  Real(Kind=rp), Dimension(:), Intent(Out) :: r ! - ranks -
!
! Locals
!
! Local scalars
  Integer :: i, j ! - indices -
!
! Executable Statements
!
! Initialise ranks
  r(1:n) = one
!
! Rank data
! - ascending -
  Select Case (order)
   Case ('A', 'a')
     Do i = 1, n-1
        Do j = i+1, n
           If (v(j) > v(i)) Then
              r(j) = r(j) + one
           Else If (v(j) < v(i)) Then
              r(i) = r(i) + one
           Else
              r(i) = r(i) + half
              r(j) = r(j) + half
           End If
        End Do
     End Do
! - descending -
   Case ('D', 'd')
     Do i = 1, n-1
        Do j = i+1, n
           If (v(j) < v(i)) Then
              r(j) = r(j) + one
           Else If (v(j) > v(i)) Then
              r(i) = r(i) + one
           Else
              r(i) = r(i) + half
              r(j) = r(j) + half
           End If
        End Do
     End Do
  End Select
!
  Return
 End Subroutine rank_data
!
!
!
 Function get_flag(iflags)
!
! Identifies the index of the first flagged item in an integer array
!
! Function type
  Integer :: get_flag
!
! Arguments
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iflags ! - flags -
!
! Executable Statements
!
! Identify flagged index
  get_flag = 1
  Do
     If (iflags(get_flag) == 1) Exit
     get_flag = get_flag + 1
  End Do
!
  Return
 End Function get_flag
!
!
!
 Subroutine set_flag (iflag,iflags)
!
! Sets the index in a flagged integer array
!
! Arguments
!
! Input scalar
  Integer, Intent(In) :: iflag ! - flag -
!
! Output arrays
  Integer, Dimension(:), Intent(Out) :: iflags ! - flags -
!
! Executable Statements
!
! Set flagged index
  iflags(:) = 0
  iflags(iflag) = 1
!
  Return
 End Subroutine set_flag
!
!
!
 Subroutine convert_units (aunits, bunits, v, ifail)
!
! Converts units of measurement
!
! Modules
  Use data_numbers,        Only: five, nine, ten, thirty, oneh, onet, zero_F, zero_K, h2okgm3
  Use data_time_constants, Only: hrpdy, minphr, secpmin
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: bunits ! - new units -
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: aunits ! - original units -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: v ! - data to be converted -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Temperature conversions
  Select Case (Trim(aunits))
! - from Kelvin -
   Case ('K', 'Kelvin', 'Kelvin_scale')
     Select Case (Trim(bunits))
      Case ('K', 'Kelvin', 'Kelvin_scale') ! - alternative labelling -
        Continue
      Case ('C', 'degree_Celsius', 'Celsius_scale') ! - to Celsius -
        v(:,:) = v(:,:) - zero_K
      Case ('F', 'Farenheit', 'Farenheit_scale') ! - to Farenheit -
        v(:,:) = zero_F + (v(:,:) - zero_K)*nine/five
      Case Default
        GoTo 1
     End Select
! - from Celsius -
   Case ('C', 'degree_Celsius', 'Celsius_scale')
     Select Case (Trim(bunits))
      Case ('C', 'degree_Celsius', 'Celsius_scale') ! - alternative labelling -
        Continue
      Case ('K', 'Kelvin', 'Kelvin_scale') ! - to Kelvin -
        v(:,:) = zero_K + v(:,:)
      Case ('F', 'Farenheit', 'Farenheit_scale') ! - to Farenheit -
        v(:,:) = zero_F + v(:,:)*nine/five
      Case Default
        GoTo 1
     End Select
! - from Farenheit -
   Case ('F', 'Farenheit', 'Farenheit_scale')
     Select Case (Trim(bunits))
      Case ('F', 'Farenheit', 'Farenheit_scale') ! - alternative labelling -
        Continue
      Case ('K', 'Kelvin', 'Kelvin_scale') ! - to Kelvin -
        v(:,:) = zero_K + (v(:,:) - zero_F)*five/nine
      Case ('C', 'degree_Celsius', 'Celsius_scale') ! - to Celsius -
        v(:,:) = (v(:,:) - zero_F)*nine/five
      Case Default
        GoTo 1
     End Select
!
! Rainfall conversions
! - from kilogrammes per metre squared per second -
   Case ('kg m-2 s-1')
     Select Case (Trim(bunits))
      Case ('m/s', 'meters/second', 'm s-1') ! - to metres per second -
        v(:,:) = v(:,:)/h2okgm3
      Case ('mm/day') ! - to mm per day -
        v(:,:) = v(:,:)*onet*secpmin*minphr*hrpdy/h2okgm3
      Case ('mm/month') ! - to mm per month -
        v(:,:) = v(:,:)*thirty*onet*secpmin*minphr*hrpdy/h2okgm3
      Case ('cm/day') ! - to cm per day -
        v(:,:) = v(:,:)*oneh*secpmin*minphr*hrpdy/h2okgm3
      Case Default
        GoTo 1
     End Select
! - from metres per second -
   Case ('m/s', 'meters/second', 'm s-1')
     Select Case (Trim(bunits))
      Case ('m/s', 'meters/second', 'm s-1') ! - alternative labelling -
        Continue
      Case ('kg m-2 s-1') ! - to kilogrammes per metre squared per second -
        v(:,:) = v(:,:)*h2okgm3
      Case ('mm/day') ! - to mm per day -
        v(:,:) = v(:,:)*onet*secpmin*minphr*hrpdy
      Case ('mm/month') ! - to mm per month -
        v(:,:) = v(:,:)*thirty*onet*secpmin*minphr*hrpdy
      Case ('cm/day') ! - to cm per day -
        v(:,:) = v(:,:)*oneh*secpmin*minphr*hrpdy
      Case Default
        GoTo 1
     End Select
! - from mm per day -
   Case ('mm/day')
     Select Case (Trim(bunits))
      Case ('kg m-2 s-1') ! - to kilogrammes per metre squared per second -
        v(:,:) = v(:,:)*h2okgm3/(onet*hrpdy*minphr*secpmin)
      Case ('m/s', 'meters/second', 'm s-1') ! - to metres per second -
        v(:,:) = v(:,:)/(onet*hrpdy*minphr*secpmin)
      Case ('mm/month') ! - to mm per month -
        v(:,:) = v(:,:)*thirty
      Case ('cm/day') ! - to cm per day -
        v(:,:) = v(:,:)/ten
      Case Default
        GoTo 1
     End Select
! - from mm per month -
   Case ('mm/month')
     v(:,:) = v(:,:)/thirty
     Select Case (Trim(bunits))
      Case ('kg m-2 s-1') ! - to kilogrammes per metre squared per second -
        v(:,:) = v(:,:)*h2okgm3/(onet*hrpdy*minphr*secpmin)
      Case ('m/s', 'meters/second', 'm s-1') ! - to metres per second -
        v(:,:) = v(:,:)/(onet*hrpdy*minphr*secpmin)
      Case ('cm/day') ! - to cm per day -
        v(:,:) = v(:,:)/ten
      Case Default
        GoTo 1
     End Select
! - from cm per day -
   Case ('cm/day')
     Select Case (Trim(bunits))
      Case ('kg m-2 s-1') ! - to kilogrammes per metre squared per second -
        v(:,:) = v(:,:)*h2okgm3/(oneh*hrpdy*minphr*secpmin)
      Case ('m/s', 'meters/second', 'm s-1') ! - to metres per second -
        v(:,:) = v(:,:)/(oneh*hrpdy*minphr*secpmin)
      Case ('cm/day') ! - to mm per day -
        v(:,:) = v(:,:)*ten
      Case ('cm/month') ! - to mm per month -
        v(:,:) = v(:,:)*ten*thirty
      Case Default
        GoTo 1
     End Select
   Case Default
     GoTo 1
  End Select
!
! Set new units
  aunits = bunits
  ifail = 0
  Return
!
! Unable to convert
1 ifail = 1
!
  Return
 End Subroutine convert_units
!
!
!
 Subroutine proj_read_results (iin, ianal, igauss, ihind, iretro, ifail)
!
! Reads results from project file
!
! Modules
  Use data_cpt_constants, Only: ia_pcr, ia_mlr, ia_gcm
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin    ! - input file unit number -
  Integer, Intent(In) :: ianal  ! - analysis indicator -
  Integer, Intent(In) :: igauss ! - transform data flag -
  Integer, Intent(In) :: ihind  ! - hindcasts flag -
  Integer, Intent(In) :: iretro ! - retroactive forecasts flag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Read results
  If (ihind == 1) Then
     Read (Unit=iin, Err=1, End=2) yhat
     Read (Unit=iin, Err=1, End=2) xm
     Read (Unit=iin, Err=1, End=2) xsd
     Read (Unit=iin, Err=1, End=2) ym
     Read (Unit=iin, Err=1, End=2) ysd
     If (igauss > 0) Then
        Read (Unit=iin, Err=1, End=2) yt
        Read (Unit=iin, Err=1, End=2) yhatt
        If (iretro == 1) Read (Unit=iin, Err=1, End=2) yrett
     End If
     If (iretro == 1) Then
        Read (Unit=iin, Err=1, End=2) yret
        Read (Unit=iin, Err=1, End=2) yrpls
     End If
  End If
  Select Case (ianal)
   Case (ia_pcr) ! - PCR -
     Read (Unit=iin, Err=1, End=2) b
   Case (ia_mlr) ! - MLR -
     Read (Unit=iin, Err=1, End=2) b
     Read (Unit=iin, Err=1, End=2) b0
   Case (ia_gcm) ! - GCM -
     Read (Unit=iin, Err=1, End=2) b
     Read (Unit=iin, Err=1, End=2) b0
     Read (Unit=iin, Err=1, End=2) xiny
  End Select
!
  ifail = 0
  Return
!
! Errors in project file
! - problem reading file -
1 ifail = 1
  Return
!
! - end of file -
2 ifail = 2
  Return
!
 End Subroutine proj_read_results
!
!
!
 Subroutine proj_write_results (iout,ianal,igauss,ihind,iretro,ifail)
!
! Writes results to project file
!
! Modules
  Use data_cpt_constants, Only: ia_pcr, ia_mlr, ia_gcm
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout   ! - output file unit number -
  Integer, Intent(In) :: ianal  ! - analysis indicator -
  Integer, Intent(In) :: igauss ! - transform data flag -
  Integer, Intent(In) :: ihind  ! - hindcasts flag -
  Integer, Intent(In) :: iretro ! - retroactive forecasts flag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write results
  If (ihind == 1) Then
     Write (Unit=iout, Err=1) yhat
     Write (Unit=iout, Err=1) xm
     Write (Unit=iout, Err=1) xsd
     Write (Unit=iout, Err=1) ym
     Write (Unit=iout, Err=1) ysd
     If (igauss > 0) Then
        Write (Unit=iout, Err=1) yt
        Write (Unit=iout, Err=1) yhatt
        If (iretro == 1) Write (Unit=iout, Err=1) yrett
     End If
     If (iretro == 1) Then
        Write (Unit=iout, Err=1) yret
        Write (Unit=iout, Err=1) yrpls
     End If
  End If
  Select Case (ianal)
   Case (ia_pcr) ! - PCR -
     Write (Unit=iout, Err=1) b
   Case (ia_mlr) ! - MLR -
     Write (Unit=iout, Err=1) b
     Write (Unit=iout, Err=1) b0
   Case (ia_gcm) ! - GCM -
     Write (Unit=iout, Err=1) b
     Write (Unit=iout, Err=1) b0
     Write (Unit=iout, Err=1) xiny
  End Select
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
!
 End Subroutine proj_write_results
End Module arrays
