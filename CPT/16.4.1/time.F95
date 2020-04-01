! Author: Simon Mason
Module time
!
! Modules
  Use data_time_constants, Only: ndw, nmn
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: add_to_month, adjust_pdate, adjust_pprd, check_lead, check_pprd, date_diff, day_of_week, day_of_year, get_cprd,       &
            get_month, get_pdate, get_sequence, idate_to_pdate, init_time, ndays, nyears, pdate_to_idate, reset_iseq, reset_time, &
            set_bdate, set_fdate, set_refdate, valid_pdate
!
! Derived-type definitions
!
! - parsed date -
  Public :: pdate
  Type pdate
     Sequence
     Integer :: iyr ! - year -
     Integer :: imn ! - month -
     Integer :: idy ! - day -
  End Type pdate
!
! - parsed period -
  Public :: pprd
  Type pprd
     Sequence
     Type(pdate) :: sdate ! - start date -
     Type(pdate) :: edate ! - end date -
  End Type pprd
!
! - relative date -
  Private :: rdate
  Type rdate
     Sequence
     Integer :: idate ! - time since reference -
     Integer :: iunit ! - units of time (1=year; 2=month; 3=day; 4=seconds) -
  End Type rdate
!
! - relative period -
  Private :: rprd
  Type rprd
     Sequence
     Type(rdate) :: sdate ! - start date -
     Integer :: lprd      ! - length of period -
  End Type rprd
!
! Arrays
!
! Integer arrays
  Integer, Dimension(ndw), Public :: idywk   ! - used days of the week (for bi-weekly data -
  Integer, Dimension(nmn), Public :: ims     ! - initial months of seasons -
  Integer, Dimension(nmn), Public :: ims_bk  ! - backup initial months of seasons -
  Integer, Dimension(nmn), Public :: mpsn    ! - number of months in seasons -
  Integer, Dimension(nmn), Public :: mpsn_bk ! - backup number of months in seasons -
!
! Scalars
!
! Integer scalars
  Integer, Public :: iseq        ! - time sequence identifier -
  Integer, Public :: iafter = -2 ! - A's date should be before (-1) equal to (0) or after (1) B's date? -
  Integer, Public :: islim       ! - start-of-season limits indicator -
!                                    0 = no season limits
!                                    1 = same calendar date
!                                    2 = same day of year
!                                    3 = same day of week
!                                    4 = evenly spaced pentads (or other interval)
  Integer, Public :: itime       ! - time option for setting length of months -
  Integer, Public :: lead        ! - lead-time -
  Integer, Public :: nse         ! - number of sub-seasons -
  Integer, Public :: nss         ! - number of seasons -
!
  Integer, Private :: isn         ! - current season -
  Integer, Private :: iyr0 = 2000 ! - reference year -
  Integer, Private :: imn0 = 1    ! - reference month -
  Integer, Private :: idy0 = 1    ! - reference day -
!
  Logical, Public :: lmsq  ! - permit monthly sequencing? -
  Logical, Public :: lssq  ! - permit seasonal sequencing? -
  Logical, Public :: lslim ! - use season limits? -
!
! Derived types
  Type(pdate), Private :: date0 ! - reference date -
!
  Type(pdate), Public :: dnow ! - current date -
!
  Type(pprd), Public :: ssnlims    ! - season limits (for sub-seasonal data) -
  Type(pprd), Public :: ssnlims_bk ! - backup season limits (for sub-seasonal data) -
!
! Interfaces
!
! Interface assignments
  Public :: Assignment(=)
  Interface Assignment(=)
     Module Procedure init_pdate
     Module Procedure init_pdates
     Module Procedure init_pprd
     Module Procedure init_pprds
  End Interface
!
! Interface operators
  Public :: Operator(+)
  Interface Operator(+)
     Module Procedure add_pdate
     Module Procedure add_pprd
  End Interface
!
  Public :: Operator(-)
  Interface Operator(-)
     Module Procedure minus_pdate
     Module Procedure minus_pprd
  End Interface
!
  Public :: Operator(==)
  Interface Operator(==)
     Module Procedure same_pdates
     Module Procedure same_pdate
     Module Procedure eq_pdate
     Module Procedure same_pprd
  End Interface
!
  Public :: Operator(/=)
  Interface Operator(/=)
     Module Procedure diff_pdate
     Module Procedure diff_pprd
     Module Procedure neq_pdate
     Module Procedure neq_pprd
  End Interface
!
  Public :: Operator(<)
  Interface Operator(<)
     Module Procedure lt_pdate
  End Interface
!
  Public :: Operator(<=)
  Interface Operator(<=)
     Module Procedure le_pdate
  End Interface
!
  Public :: Operator(>)
  Interface Operator(>)
     Module Procedure gt_pdate
  End Interface
!
  Public :: Operator(>=)
  Interface Operator(>=)
     Module Procedure ge_pdate
  End Interface
!
! Generic interfaces
  Public same_ssn
  Interface same_ssn
   Module Procedure same_ssn_date
   Module Procedure same_ssn_period
  End Interface same_ssn
!
  Public get_cdate
  Interface get_cdate
   Module Procedure get_cdate_date
   Module Procedure get_cdate_period
   Module Procedure get_cdate_lsn
   Module Procedure get_cdate_lsns
  End Interface get_cdate
!
Contains
!
!
 Subroutine init_pdate (d, i)
!
! Initialises a parsed date
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - assignment value -
!
! Output scalars
  Type(pdate), Intent(Out) :: d ! - date -
!
! Executable Statements
!
! Assign dates
  d%iyr = i
  d%imn = i
  d%idy = i
!
  Return
 End Subroutine init_pdate
!
!
!
 Subroutine init_pdates (d, i)
!
! Initialises parsed dates
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - assignment value -
!
! Output arrays
  Type(pdate), Dimension(:), Intent(Out) :: d ! - dates -
!
! Locals
!
! Local scalars
  Integer :: j ! - index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Size
!
! Executable Statements
!
! Assign dates
  Do j = 1, Size(d)
     d(j) = i
  End Do
!
  Return
 End Subroutine init_pdates
!
!
!
 Subroutine init_pprd (p, i)
!
! Initialises a parsed period
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - assignment value -
!
! Output scalars
  Type(pprd), Intent(Out) :: p ! - period -
!
! Executable Statements
!
! Assign period
  p%sdate = i
  p%edate = i
!
  Return
 End Subroutine init_pprd
!
!
!
 Subroutine init_pprds (p, i)
!
! Initialises parsed periods
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - assignment value -
!
! Output arrays
  Type(pprd), Dimension(:), Intent(Out) :: p ! - periods -
!
! Locals
!
! Local scalars
  Integer :: j ! - index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Size
!
! Executable Statements
!
! Assign periods
  Do j = 1, Size(p)
     p(j) = i
  End Do
!
  Return
 End Subroutine init_pprds
!
!
!
 Recursive Function add_pdate(d, i) &
           Result (dnew)
!
! Increments a parsed date by a constant
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr
!
! Function type
  Type(pdate) :: dnew
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - increment -
!
  Type(pdate), Intent(In) :: d ! - date -
!
! Locals
!
! Local scalars
  Integer :: iid    ! - number of dekadal increments -
  Integer :: iim    ! - number of monthly increments -
  Integer :: iiy    ! - number of annual increments -
  Integer :: ia     ! - absolute increment -
  Integer :: id     ! - direction -
  Integer :: ima    ! - previous season index -
  Integer :: imb    ! - current season index -
  Integer :: isq_bk ! - backup time sequencing -
  Integer :: k      ! - index -
  Integer :: l      ! - index -
  Integer :: nrs    ! - number of remaining seasons -
!
  Type(pprd) :: cssnlim ! - current season limits -
!
! Functions and Subroutines
!
! Instrinsic functions
  Intrinsic Abs, Any, Mod, Sign
!
! Executable Statements
!
! Increment year
  dnew = d
  If (i == 0) Return
  Select Case (iseq)
   Case (isq_yr)
     dnew%iyr = d%iyr + i
!
! Increment season
   Case (isq_sn)
     dnew%iyr = d%iyr + i/nss ! - add difference in whole years -
     nrs = Sign(Mod(i, nss), i) ! - identify partial year differences -
     ima = get_ssn(d%imn)
     If (ima == 0) Return
     imb = ima
     Select Case (nrs)
      Case (0) ! - same season -
        Continue
      Case (1:) ! - increment seasons -
        Do l = 1, nrs
           imb = imb + 1
           If (imb > nss) imb = imb - nss
           If (ims(imb) < ims(ima)) dnew%iyr = dnew%iyr + 1
           ima = imb
        End Do
      Case (:-1) ! - decrement seasons -
        Do l = 1, -nrs
           imb = imb - 1
           If (imb < 1) imb = imb + nss
           If (ims(imb) > ims(ima)) dnew%iyr = dnew%iyr - 1
           ima = imb
        End Do
     End Select
     dnew%imn = ims(imb)
     isn = ima ! - set current season -
!
! Increment month
   Case (isq_mn, isq_so)
     dnew%imn = add_to_month(d%imn, i, &
                             iyr=dnew%iyr)
!
! Increment day, pentad, week
   Case (1, 5, 7)
! - no seasonal limits -
     dnew = d
     If (.not.lslim .or. ssnlims%sdate == ssnlims%edate) Then
        Call add_days (i)
! - seasonal limits -
     Else
        cssnlim = ssnlims
        If (d%imn > ssnlims%sdate%imn .or. (d%imn == ssnlims%sdate%imn .and. d%idy >= ssnlims%sdate%idy)) Then
           cssnlim%sdate%iyr = d%iyr
        Else
           cssnlim%sdate%iyr = d%iyr - 1
        End If
        If (ssnlims%sdate < ssnlims%edate) Then ! - season does not span year-end -
           cssnlim%edate%iyr = cssnlim%sdate%iyr
        Else ! - season spans year-end -
           cssnlim%edate%iyr = cssnlim%sdate%iyr + 1
        End If
        isq_bk = iseq
        iseq = 1
        If (i > 0) Then ! - increment dates -
           Do k = 1, i*isq_bk
              Call add_days (1)
              If (dnew > cssnlim%edate) Then
                 cssnlim%sdate%iyr = cssnlim%sdate%iyr + 1
                 cssnlim%edate%iyr = cssnlim%edate%iyr + 1
                 dnew = cssnlim%sdate
              End If
           End Do
        Else ! - subtract dates -
           Do k = 1, -i*isq_bk
              Call add_days (-1)
              If (dnew < cssnlim%sdate) Then
                 cssnlim%sdate%iyr = cssnlim%sdate%iyr - 1
                 cssnlim%edate%iyr = cssnlim%edate%iyr - 1
                 dnew = cssnlim%edate
              End If
           End Do
        End If
        iseq = isq_bk
     End If
!
! Increment half-weeks
   Case (3, 4)
     isq_bk = iseq
     iseq = 1
     dnew = d + (i/2)*7
     If (Mod(i,2) /= 0) Then ! - adjust to nearest used day-of-week if a partial week is incremented -
        If (i > 0) Then
           dnew = dnew + 3
        Else
           dnew = dnew - 4
        End If
        If (.not.Any(idywk(1:2) == day_of_week(dnew))) dnew = dnew + 1 ! - adjust if day of the week is incorrect -
     End If
     iseq = isq_bk
!
! Increment dekads
   Case (10)
     dnew = d
     ia = Abs(i)
     id = Sign(1,i)
     If (ia < 3) Then
        iid = i
     Else
        iim = ia/3
        If (iim < nmn) Then
           iiy = 0
        Else
           iiy = iim/nmn
           dnew%iyr = d%iyr + iiy*id
           iim = iim - iiy*nmn
        End If
        dnew%imn = add_to_month(d%imn, iim*id, &
                                iyr=dnew%iyr)
        iid = (ia - (iiy*nmn + iim)*3)*id
        If (iid == 0) Return
     End If
     Do k = id, iid, id
        dnew%idy = dnew%idy + id*10
        Select Case (dnew%idy)
         Case (1, 11, 21)
           Continue
         Case (31)
           dnew%idy = 1
           dnew%imn = add_to_month(dnew%imn, 1, &
                                   iyr=dnew%iyr)
         Case (-9)
           dnew%idy = 21
           dnew%imn = add_to_month(dnew%imn, -1, &
                                   iyr=dnew%iyr)
        End Select
     End Do
  End Select
!
  Return
!
 Contains
!
!
  Subroutine add_days (i)
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - increment -
!
! Increment
  dnew%idy = dnew%idy + i*iseq
  If (dnew%idy < 1) Then ! - adjust if now earlier month -
     Do
        dnew%imn = add_to_month(dnew%imn, -1, &
                                iyr=dnew%iyr)
        dnew%idy = dnew%idy + ndays(dnew%iyr, dnew%imn)
        If (dnew%idy >= 1) Exit
     End Do
  Else If (dnew%idy > ndays(dnew%iyr, dnew%imn)) Then ! - adjust if now later month -
     Do
        dnew%idy = dnew%idy-ndays(dnew%iyr, dnew%imn)
        dnew%imn = add_to_month(dnew%imn, 1, &
                                iyr=dnew%iyr)
        If (dnew%idy <= ndays(dnew%iyr, dnew%imn)) Exit
     End Do
  End If
!
  Return
  End Subroutine add_days
 End Function add_pdate
!
!
!
 Function minus_pdate(d, i)
!
! Decreases a parsed date by a constant
!
! Function type
  Type(pdate) :: minus_pdate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - decrease -
!
  Type(pdate), Intent(In) :: d ! - date -
!
! Executable Statements
!
! Decrease date
  minus_pdate = add_pdate(d, -i)
!
  Return
 End Function minus_pdate
!
!
!
 Function add_to_month(i0, i, &
          iyr) &
          Result(imn)
!
! Adds (or subtracts) a specified number of months
!
! Function type
  Integer :: imn
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i0 ! - initial month -
  Integer, Intent(In) :: i  ! - increment -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: iyr ! - year -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Add increment to initial month
  imn = i0 + i
  Select Case (imn)
   Case (1:nmn)
     Continue
   Case (:0) ! - adjust if now earlier year -
     Do
        If (Present(iyr)) iyr = iyr - 1
        imn = imn + nmn
        If (imn > 0) Exit
     End Do
   Case (nmn+1:) ! - adjust if now later year -
     Do
        If (Present(iyr)) iyr = iyr + 1
        imn = imn - nmn
        If (imn <= nmn) Exit
     End Do
  End Select
!
  Return
 End Function add_to_month
!
!
!
 Function add_pprd(p, i)
!
! Increments a parsed period by a constant
!
! Modules
  Use data_time_constants, Only: isq_dy, isq_mn, isq_sn, isq_so, isq_yr
!
! Function type
  Type(pprd) :: add_pprd
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - increment -
!
  Type(pprd), Intent(In) :: p ! - period -
!
! Locals
!
! Local scalars
  Integer :: isq_bk ! - back-up time sequence -
!
  Logical :: ls_bk ! - back-up use seasonal limits? -
!
! Executable Statements
!
! Increment period
  add_pprd%sdate = p%sdate + i
  Select Case (iseq)
   Case (isq_mn, isq_so, isq_yr)
     add_pprd%edate = p%edate + i
   Case (isq_sn)
     If (date_diff(p%sdate, p%edate, isq_mn) <= nmn) Then ! - check if "season" is multi-year -
        add_pprd%edate = add_pprd%sdate
     Else
        add_pprd%edate = p%edate
        add_pprd%edate%iyr = add_pprd%edate%iyr + i/nss
        add_pprd%edate%imn = add_pprd%sdate%imn
        If (add_pprd%edate < add_pprd%sdate) add_pprd%edate%iyr = add_pprd%edate%iyr + 1
     End If
     isn = get_ssn(add_pprd%edate%imn)
     If (isn > 0) Then
        iseq = isq_mn
        add_pprd%edate = add_pprd%edate + (mpsn(isn) - 1)
        iseq = isq_sn
     End If
   Case (10)
     Select Case (add_pprd%sdate%idy)
      Case (1)
        add_pprd%edate%iyr = add_pprd%sdate%iyr
        add_pprd%edate%imn = add_pprd%sdate%imn
        add_pprd%edate%idy = 10
      Case (11)
        add_pprd%edate%iyr = add_pprd%sdate%iyr
        add_pprd%edate%imn = add_pprd%sdate%imn
        add_pprd%edate%idy = 20
      Case (21)
        add_pprd%edate%idy = ndays(add_pprd%sdate%iyr, add_pprd%sdate%imn)
        add_pprd%edate%iyr = add_pprd%sdate%iyr
        add_pprd%edate%imn = add_pprd%sdate%imn
      Case Default
        add_pprd%edate = add_pprd%sdate + 10
     End Select
   Case Default ! - add period length in days if data are sub-seasonal -
     isq_bk = iseq
     ls_bk = lslim
     iseq = isq_dy
     lslim = .false.
     add_pprd%edate = add_pprd%sdate + date_diff(p%sdate, p%edate, isq_dy)
     iseq = isq_bk
     lslim = ls_bk
  End Select
!
  Return
 End Function add_pprd
!
!
!
 Function minus_pprd(p, i)
!
! Decreases a parsed period by a constant
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn
!
! Function type
  Type(pprd) :: minus_pprd
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - increment -
!
  Type(pprd), Intent(In) :: p ! - period -
!
! Executable Statements
!
! Increment period
  minus_pprd%sdate = p%sdate - i
  Select Case (iseq)
   Case Default
     minus_pprd%edate = p%edate - i
   Case (isq_sn)
     iseq = isq_mn
     minus_pprd%edate = minus_pprd%sdate + (mpsn(isn) - 1)
     iseq = isq_sn
   Case (10)
     Select Case (minus_pprd%sdate%idy)
      Case (1)
        minus_pprd%edate%iyr = minus_pprd%sdate%iyr
        minus_pprd%edate%imn = minus_pprd%sdate%imn
        minus_pprd%edate%idy = 10
      Case (11)
        minus_pprd%edate%iyr = minus_pprd%sdate%iyr
        minus_pprd%edate%imn = minus_pprd%sdate%imn
        minus_pprd%edate%idy = 20
      Case (21)
        minus_pprd%edate%idy = ndays(minus_pprd%sdate%iyr, minus_pprd%sdate%imn)
        minus_pprd%edate%iyr = minus_pprd%sdate%iyr
        minus_pprd%edate%imn = minus_pprd%sdate%imn
      Case Default
        minus_pprd%edate = minus_pprd%sdate + 10
     End Select
  End Select
!
  Return
 End Function minus_pprd
!
!
!
 Function same_pdate(d1, d2) &
          Result (same)
!
! Identical parsed dates?
!
! Function type
  Logical :: same
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  same = .false.
  If (d1%iyr /= d2%iyr) Return
  If (d1%imn /= d2%imn) Return
  If (d1%idy /= d2%idy) Return
  same = .true.
!
  Return
 End Function same_pdate
!
!
!
 Function same_pdates(d1, d2) &
          Result (same)
!
! Array of identical parsed dates?
!
! Function type
  Logical :: same
!
! Arguments
!
! Input arrays
  Type(pdate), Dimension(:), Intent(In) :: d1 ! - first dates -
  Type(pdate), Dimension(:), Intent(In) :: d2 ! - second dates -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any
!
! Executable Statements
!
! Compare dates
  same = .false.
  If (Any(d1(:)%iyr /= d2(:)%iyr)) Return
  If (Any(d1(:)%imn /= d2(:)%imn)) Return
  If (Any(d1(:)%idy /= d2(:)%idy)) Return
  same = .true.
!
  Return
 End Function same_pdates
!
!
!
 Function same_ssn_date(d1, d2) &
          Result (same)
!
! Identical seasons?
!
! Function type
  Logical :: same
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  same = .false.
  If (d1%imn /= d2%imn) Return
  If (d1%idy /= d2%idy) Return
  same = .true.
!
  Return
 End Function same_ssn_date
!
!
!
 Function same_ssn_period(p1, p2) &
          Result (same)
!
! Identical seasons?
!
! Function type
  Logical :: same
!
! Arguments
!
! Input scalars
  Type(pprd), Intent(In) :: p1 ! - first period -
  Type(pprd), Intent(In) :: p2 ! - second period -
!
! Executable Statements
!
! Compare dates
  same = .false.
  If (p1%sdate%imn /= p2%sdate%imn) Return
  If (p1%sdate%idy /= p2%sdate%idy) Return
  If (p1%edate%imn /= p2%edate%imn) Return
  If (p1%edate%idy /= p2%edate%idy) Return
  same = .true.
!
  Return
 End Function same_ssn_period
!
!
!
 Function eq_pdate(d, i) &
          Result (eq)
!
! Initialised parsed date?
!
! Function type
  Logical :: eq
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - constant -
!
  Type(pdate), Intent(In) :: d ! - date -
!
! Executable Statements
!
! Compare dates
  eq = .false.
  If (d%iyr /= i) Return
  If (d%imn /= i) Return
  If (d%idy /= i) Return
  eq = .true.
!
  Return
 End Function eq_pdate
!
!
!
 Function neq_pdate(d, i) &
          Result (neq)
!
! Uninitialised parsed date?
!
! Function type
  Logical :: neq
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - constant -
!
  Type(pdate), Intent(In) :: d ! - date -
!
! Executable Statements
!
! Compare dates
  neq = .true.
  If (d%iyr /= i) Return
  If (d%imn /= i) Return
  If (d%idy /= i) Return
  neq = .false.
!
  Return
 End Function neq_pdate
!
!
!
 Function neq_pprd(p, i) &
          Result (neq)
!
! Uninitialised parsed period?
!
! Function type
  Logical :: neq
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - constant -
!
  Type(pprd), Intent(In) :: p ! - period -
!
! Executable Statements
!
! Compare dates
  neq = (p%sdate /= i) .or. (p%edate /= i)
!
  Return
 End Function neq_pprd
!
!
!
 Function same_pprd(p1, p2) &
          Result (same)
!
! Identical parsed periods?
!
! Function type
  Logical :: same
!
! Arguments
!
! Input scalars
  Type(pprd), Intent(In) :: p1 ! - first period -
  Type(pprd), Intent(In) :: p2 ! - second period -
!
! Executable Statements
!
! Compare dates
  same = .false.
  If (p1%sdate%iyr /= p2%sdate%iyr) Return
  If (p1%sdate%imn /= p2%sdate%imn) Return
  If (p1%sdate%idy /= p2%sdate%idy) Return
  If (p1%edate%iyr /= p2%edate%iyr) Return
  If (p1%edate%imn /= p2%edate%imn) Return
  If (p1%edate%idy /= p2%edate%idy) Return
  same = .true.
!
  Return
 End Function same_pprd
!
!
!
 Function diff_pdate(d1, d2) &
          Result (differ)
!
! Different parsed dates?
!
! Function type
  Logical :: differ
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  differ = .true.
  If (d1%iyr /= d2%iyr) Return
  If (d1%imn /= d2%imn) Return
  If (d1%idy /= d2%idy) Return
  differ = .false.
!
  Return
 End Function diff_pdate
!
!
!
 Function diff_pprd(p1, p2) &
          Result (differ)
!
! Different parsed periods?
!
! Function type
  Logical :: differ
!
! Arguments
!
! Input scalars
  Type(pprd), Intent(In) :: p1 ! - first period -
  Type(pprd), Intent(In) :: p2 ! - second period -
!
! Executable Statements
!
! Compare periods
  differ = .true.
  If (p1%sdate%iyr /= p2%sdate%iyr) Return
  If (p1%sdate%imn /= p2%sdate%imn) Return
  If (p1%sdate%idy /= p2%sdate%idy) Return
  If (p1%edate%iyr /= p2%edate%iyr) Return
  If (p1%edate%imn /= p2%edate%imn) Return
  If (p1%edate%idy /= p2%edate%idy) Return
  differ = .false.
!
  Return
 End Function diff_pprd
!
!
!
 Function lt_pdate(d1,d2)
!
! Is first parsed date less than second parsed date?
!
! Function type
  Logical :: lt_pdate
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  lt_pdate = .true.
  If (d1%iyr < d2%iyr) Then
     Return
  Else If (d1%iyr == d2%iyr) Then
     If (d1%imn < d2%imn) Then
        Return
     Else If (d1%imn == d2%imn) Then
        If (d1%idy < d2%idy) Return
     End If
  End If
  lt_pdate = .false.
!
  Return
 End Function lt_pdate
!
!
!
 Function le_pdate(d1, d2)
!
! Is first parsed date less than or equal to second parsed date?
!
! Function type
  Logical :: le_pdate
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  If (d1 == d2 .or. d1 < d2) Then
     le_pdate = .true.
  Else
     le_pdate = .false.
  End If
!
  Return
 End Function le_pdate
!
!
!
 Function gt_pdate(d1, d2)
!
! Is first parsed date greater than second parsed date?
!
! Function type
  Logical :: gt_pdate
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  gt_pdate = .true.
  If (d1%iyr > d2%iyr) Then
     Return
  Else If (d1%iyr == d2%iyr) Then
     If (d1%imn > d2%imn) Then
        Return
     Else If (d1%imn == d2%imn) Then
        If (d1%idy > d2%idy) Return
     End If
  End If
  gt_pdate = .false.
!
  Return
 End Function gt_pdate
!
!
!
 Function ge_pdate(d1, d2)
!
! Is first parsed date greater than or equal to second parsed date?
!
! Function type
  Logical :: ge_pdate
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  If (d1 == d2 .or. d1 > d2) Then
     ge_pdate = .true.
  Else
     ge_pdate = .false.
  End If
!
  Return
 End Function ge_pdate
!
!
!
 Subroutine init_time ()
!
! Initialises time constants
!
! Modules
  Use data_version, Only: cyr1, cyr2, date_this, date_ver1
!
! Locals
!
! Local scalars
  Integer :: iyr1 ! - index -
!
! Local arrays
  Integer, Dimension(8) :: idattime ! - date and time -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Date_and_Time, Index
!
! Executable Statements
!
! Identify year of first version
  iyr1 = Index(date_ver1, ', ') + 2
  cyr1 = date_ver1(iyr1:iyr1+3)
!
! Identify year of current version
  iyr1 = Index(date_this, ', ') + 2
  cyr2 = date_this(iyr1:iyr1+3)
!
! Initialise season settings
  Call reset_time (.true.)
!
! Use British Gregorian calendar
  itime = 1
!
! Set current date and time
  Call Date_and_Time (values=idattime)
  dnow%iyr = idattime(1)
  dnow%imn = idattime(2)
  dnow%idy = idattime(3)
!
  Return
 End Subroutine init_time
!
!
!
 Subroutine reset_time (lhard)
!
! Resets definitions of season start dates and lengths
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lhard ! - hard reset? -
!
! Locals
!
! Local scalars
  Integer :: l ! - index -
!
! Executable Statements
!
! Set used days of the week
  idywk(:) = 0
!
! Set initial months of seasons
  nss = 1
  ForAll (l=1:nmn)
     ims(l)=l
  End ForAll
  mpsn(:) = 1
  isn = 1
!
! Reset all time settings
  If (lhard) Then
     iseq = 0
     ssnlims = 0
     islim = 0
     lslim = .false.
     Call set_refdate ()
  End If
!
  Return
 End Subroutine reset_time
!
!
!
 Subroutine set_refdate ( &
            iy0)
!
! Sets reference date
!
! Arguments
!
! Input scalars
! - optional inpnut scalars -
  Integer, Optional :: iy0
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Set reference date
  If (Present(iy0)) Then
     date0%iyr = iy0
  Else
     date0%iyr = iyr0
     date0%imn = imn0
     date0%idy = idy0
  End If
!
  Return
 End Subroutine set_refdate
!
!
!
 Subroutine reset_iseq (isq_new)
!
! Resets time-sequencing. Restores backup if no argument is present
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Integer, Intent(In), Optional :: isq_new ! - new time seuencing -
!
! Locals
!
! Local scalars
  Integer, Save :: isq_bk ! - back-up time-sequencing -
!
  Logical, Save :: lseq_bk = .false. ! - backed-up? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Reset time sequencing
! - set new time sequencing -
  If (Present(isq_new)) Then
     If (.not.lseq_bk) isq_bk = iseq
     iseq = isq_new
     lseq_bk = .true.
! - restore previous time sequencing -
  Else
     iseq = isq_bk
     lseq_bk = .false.
  End If
!
  Return
 End Subroutine reset_iseq
!
!
!
 Function get_month(cmn)
!
! Identifies month number
!
! Modules
  Use data_text,           Only: nlang
  Use data_time_constants, Only: cmon, lcmon, umon
!
! Function type
  Integer :: get_month
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cmn ! - month -
!
! Locals
!
! Local scalars
  Integer :: i ! - language index -
  Integer :: l ! - month index -
!
! Executable Statements
!
! Identify month
  Do l = 1, nmn
     Do i = 1, nlang
        If (cmn(1:lcmon) == cmon(l,i) .or. cmn(1:lcmon) == umon(l,i)) Then
           get_month = l
           Return
        End If
     End Do
  End Do
  get_month = 0
!
  Return
 End Function get_month
!
!
!
 Function get_ssn(imn)
!
! Identifies current season
!
! Function type
  Integer :: get_ssn
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imn ! - start month of season -
!
! Locals
!
! Local scalars
  Integer :: l ! - season index -
!
! Executable Statements
!
! Identify season
  get_ssn = 0
  Do l = 1, nss
     If (imn == ims(l)) Then
        get_ssn = l
        Exit
     End If
  End Do
!
  Return
 End Function get_ssn
!
!
!
 Function ndays(iyr, imn)
!
! Calculates number of days in the month
! NB - assumes the Gregorian calendar as implemented by Britain and the British Empire
! For other calendars, set itime to:
! 1 British Gregorian calendar
! 2 Julian calendar
! 3 no leap years
! 4 30-day months
! 5 Catholic Gregorian calendar
!
! Function type
  Integer :: ndays
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iyr ! - year -
  Integer, Intent(In) :: imn ! - month -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod
!
! Executable Statements
!
! Define number of days
  If (itime == 4) Then ! - 30-day months -
     ndays = 30
     Return
  End If
  Select Case (imn)
   Case (1)  ! January
     ndays = 31
   Case (2)  ! February
     Select Case (itime)
      Case (1, 5) ! - include Gregorian leap years -
        If (Mod(iyr, 4) == 0) Then
           If (Mod(iyr, 100) == 0) Then
              If (Mod(iyr, 400) == 0) Then
                 ndays = 29
              Else
                 ndays = 28
              End If
           Else
              ndays = 29
           End If
        Else
           ndays = 28
        End If
      Case (2) ! - include Julian leap years -
        If (Mod(iyr, 4) == 0) Then
           ndays = 29
        Else
           ndays = 28
        End If
      Case (3) ! - exclude leap years -
        ndays = 28
     End Select
   Case (3)  ! March
     ndays = 31
   Case (4)  ! April
     ndays = 30
   Case (5)  ! May
     ndays = 31
   Case (6)  ! June
     ndays = 30
   Case (7)  ! July
     ndays = 31
   Case (8)  ! August
     ndays = 31
   Case (9)  ! September
     Select Case (itime)
      Case (1)
        If (iyr /= 1752) Then
           ndays = 30
        Else
           ndays = 19
        End If
     Case Default
        ndays = 30
     End Select
   Case (10) ! October
     Select Case (itime)
      Case (5)
        If (iyr /= 1582) Then
           ndays = 31
        Else
           ndays = 21
        End If
     Case Default
        ndays = 31
     End Select
   Case (11) ! November
     ndays = 30
   Case (12) ! December
     ndays = 31
   Case Default
     ndays = 0
  End Select
!
  Return
 End Function ndays
!
!
!
 Function nyears(n, lssy)
!
! Identifies number of years
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr, nmn
!
! Function type
  Integer :: nyears
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of cases -
!
  Logical, Intent(In) :: lssy ! - convert sub-seasonal data to years? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Ceiling, Real
!
! Executable Statements
!
! Identify number of years
  Select Case (iseq)
   Case (isq_yr)
     nyears = n
   Case (isq_sn)
     nyears = Ceiling(Real(n)/Real(nss))
   Case (isq_mn, isq_so)
     nyears = Ceiling(Real(n)/Real(nmn))
   Case (1:)
     If (lssy) Then
        nyears = Ceiling(Real(n)/Real(nse))
     Else
        nyears = n
     End If
  End Select
!
  Return
 End Function nyears
!
!
!
 Function day_of_year(d1)
!
! Calculates day of the year
!
! Function type
  Integer :: day_of_year
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - date -
!
! Locals
!
! Local scalars
  Logical :: lslim_bk ! - backup use season limits? -
!
  Type(pdate) :: d0 ! - first day of the year -
!
! Executable Statements
!
! Set first day of the year
  d0 = d1
  d0%imn = 1
  d0%idy = 1
!
! Calculate day of the year
  lslim_bk = lslim
  lslim = .false.
  day_of_year = 1 + date_diff(d0, d1, 1)
  lslim = lslim_bk
!
  Return
 End Function day_of_year
!
!
!
 Function day_of_week(d1)
!
! Calculates day of the week, using Zeller's Rule
! NB - only works for the following calendars
! 1 British Gregorian calendar
! 2 Julian calendar
! 5 Catholic Gregorian calendar
! Returns:
!  1 Sunday
!  2 Monday
!  3 Tuesday
!  4 Wednesday
!  5 Thursday
!  6 Friday
!  7 Saturday
! -1 Invalid calendar
! -2 Invalid day of month
! -3 Invalid month
!
! Modules
  Use data_time_constants, Only: nmn
!
! Function type
  Integer :: day_of_week
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - date -
!
! Locals
!
! Local scalars
  Integer :: iy ! - adjusted year -
  Integer :: im ! - adjusted month -
!
! Executable Statements
!
! Return if calendard is invalid
  Select Case (itime)
   Case (1, 2, 5)
     Continue
   Case Default
     day_of_week = -1
     Return
  End Select
!
! Check for valid date
! - invalid month -
  If (d1%imn < 1 .or. d1%imn > nmn) Then
     day_of_week = -3
     Return
  End If
! - invalid day -
  If (d1%idy < 1 .or. d1%idy > ndays(d1%iyr, d1%imn)) Then
     day_of_week = -2
     Return
  End If
!
! Identify day using Zeller's Rule
! - adjust to a March start-of-year -
  iy = d1%iyr
  im = d1%imn
  If (im < 3) Then
     im = im + nmn
     iy = iy - 1
  End If
! - adjust if calendar is not Gregorian -
  Select Case (itime)
   Case (1) ! - British Gregorian -
     Select Case (d1%iyr)
      Case (:1751)
        day_of_week = day_of_wk(.false.)
      Case (1752)
        Select Case (d1%imn)
         Case (1:8)
           day_of_week = day_of_wk(.false.)
         Case (9)
           Select Case (d1%idy)
            Case (1:2)
              day_of_week = day_of_wk(.false.)
            Case (3:13)
              day_of_week = -2
            Case (14:30)
              day_of_week = day_of_wk(.true.)
           End Select
         Case (10:12)
           day_of_week = day_of_wk(.true.)
        End Select
      Case (1753:)
        day_of_week = day_of_wk(.true.)
     End Select
   Case (2) ! - Julian -
     day_of_week = day_of_wk(.false.)
   Case (5) ! - Catholic Gregorian -
     Select Case (d1%iyr)
      Case (:1581)
        day_of_week = day_of_wk(.false.)
      Case (1582)
        Select Case (d1%imn)
         Case (1:9)
           day_of_week = day_of_wk(.false.)
         Case (10)
           Select Case (d1%idy)
            Case (1:4)
              day_of_week = day_of_wk(.false.)
            Case (5:14)
              day_of_week = -2
            Case (15:31)
              day_of_week = day_of_wk(.true.)
           End Select
         Case (11:12)
           day_of_week = day_of_wk(.true.)
        End Select
      Case (1583:)
        day_of_week = day_of_wk(.true.)
     End Select
  End Select
!
 Contains
!
!
 Function day_of_wk(lgc)
!
! Calculates day of week
!
! Function type
  Integer :: day_of_wk
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lgc ! - Gregorian Calendar? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod
!
! Executable Statements
!
! Calculate day-of-week
  If (lgc) Then ! - Gregorian -
     day_of_wk = d1%idy + 13*(im + 1)/5 + iy + iy/4 - iy/100 + iy/400
  Else ! - Julian -
     day_of_wk = d1%idy + 13*(im + 1)/5 + iy + iy/4 + 5
  End If
  day_of_wk = Mod(day_of_wk, 7)
  If (day_of_wk == 0) day_of_wk = 7
!
  Return
 End Function day_of_wk
 End Function day_of_week
!
!
!
 Subroutine get_pdate (cdate, ctag, sdate, ifail, edate)
!
! Identifies start and end dates of a parsed period from a character string
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Requested CPT date tag not present
!    ifail =  2 Problem with date format
!    ifail =  3 Invalid month
!    ifail =  4 Invalid day of month
!    ifail =  5 Invalid month and day of month
!
! Modules
  Use data_numbers, Only: cdigits
!
! Arguments
!
! Input scalars
  Character(Len=1), Intent(In) :: cdate ! - date type -
  Character(Len=*), Intent(In) :: ctag  ! - tags -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
  Type(pdate), Intent(Out) :: sdate ! - start date -
! - optional output scalars -
  Type(pdate), Intent(Out), Optional :: edate ! - end date -
!
! Locals
!
! Local scalars
  Integer :: i1 ! - start of date tag -
  Integer :: i2 ! - locator -
  Integer :: i3 ! - end of date tag -
  Integer :: i4 ! - locator -
  Integer :: it ! - date -
!
  Logical :: lday ! - day present? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Len_Trim, Present, Scan, Verify
!
! Executable Statements
!
! Initialise dates
  sdate=0
!
! Locate date tag
  Select Case (cdate)
   Case (' ')
     i1 = 1
     i3 = Verify(ctag, cdigits//'/-T:') - 1
     If (i3 == -1) i3 = Len_Trim(ctag)
   Case ('F', 'S', 'T')
     i1 = Index(ctag(1:), 'cpt:'//cdate//'=')
     If (i1 == 0) Then
        ifail = 1
        Return
     End If
     i1 = i1 + Len('cpt:'//cdate//'=')
     i3 = Index(ctag(i1:), ', ')
     If (i3 > 0) Then
        i3 = i1 + i3 - 2
     Else
        i3 = Len_Trim(ctag)
     End If
  End Select
  If (i3 == 0) GoTo 3
!
! Identify start year
  i2 = Index(ctag(i1:i3), '-')
  If (i2 == 0) Then
     Read (Unit=ctag(i1:i3), Fmt=*, Err=3) sdate%iyr
     If (Present(edate)) edate = sdate
     GoTo 2
  End If
  i2 = i1 + i2 - 2
  Read (Unit=ctag(i1:i2), Fmt=*, Err=3) sdate%iyr
!
! Read month
  i1 = i2 + 2
  i2 = Scan(ctag(i1:i3), '/-T:')
  If (i2 == 0) Then
     Read (Unit=ctag(i1:i3), Fmt=*, Err=3) sdate%imn
     If (Present(edate)) edate = sdate
     GoTo 2
  End If
  i2 = i1 + i2 - 2
  Read (Unit=ctag(i1:i2), Fmt=*, Err=3) sdate%imn
!
! Identify whether day is present
  lday = .false.
  i1 = i2 + 1
  If (ctag(i1:i1) == '-') Then ! - day -
     i1 = i1 + 1
     i2 = Index(ctag(i1:i3), '/')
     If (i2 == 0) Then
        i2 = Index(ctag(i1:i3), 'T')
        If (i2 == 0) Then
           i2 = i3
        Else
           i2 = i1 + i2 - 2
        End If
        Read (Unit=ctag(i1:i2), Fmt=*, Err=3) sdate%idy
        If (Present(edate)) edate = sdate
        GoTo 2
     End If
     i2 = i1 + i2 - 2
     i4 = Index(ctag(i1:i3), 'T')
     If (i4 == 0) Then
        i4 = i2
     Else
        i4 = i1 + i4 - 2
     End If
     Read (Unit=ctag(i1:i4), Fmt=*, Err=3) sdate%idy
     lday = .true.
  End If
!
! Identify end date
  If (Present(edate)) Then
     edate = sdate
     i1 = i2 + 1
     If (ctag(i1:i1) == '/') Then
        i1 = i1 + 1
1       i2 = Index(ctag(i1:i3), '-')
        If (i2 == 0) Then
           i2 = Index(ctag(i1:i3), 'T')
           If (i2 == 0) Then
              i2 = i3
           Else
              i2 = i1 + i2 - 2
           End If
           Read (Unit=ctag(i1:i2), Fmt=*, Err=3) it
           If (lday) Then
              edate%idy = it
           Else
              edate%imn = it
           End If
        Else
           i2 = i1 + i2 - 2
           Read (Unit=ctag(i1:i2), Fmt=*, Err=3) it
           i1 = i2 + 2
           i2 = Index(ctag(i1:i3), '-')
           If (i2 == 0) Then
              If (lday) Then
                 edate%imn = it
              Else
                 edate%iyr = it
              End If
           Else
              edate%iyr = it
           End If
           GoTo 1
        End If
     Else
        GoTo 3
     End If
  End If
!
! Check validity of dates
2 ifail = valid_pdate(sdate)
  If (ifail /= 0) Then
     ifail = ifail + 2
     Return
  End If
  If (Present(edate)) Then
     ifail = valid_pdate(edate)
     If (ifail/=0) Then
        ifail = ifail + 2
        Return
     End If
  End If
  ifail = 0
  Return
!
! Errors
3 ifail = 2
  Return
 End Subroutine get_pdate
!
!
!
 Function valid_pdate(d)
!
! Checks that parsed date is valid
!
! On exit:
!    valid_pdate =  0 Valid
!    valid_pdate =  1 Invalid month
!    valid_pdate =  2 Invalid day of month
!    valid_pdate =  3 Invalid month and day of month
!
! Modules
  Use data_time_constants, Only: mdm
!
! Function type
  Integer :: valid_pdate
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d ! - date -
!
! Executable Statements
!
! Check validity of date
  valid_pdate = 0
  If (d%imn >= 1 .and. d%imn <= nmn) Then
     If (d%idy < 0 .or. d%idy > ndays(d%iyr, d%imn)) valid_pdate = valid_pdate + 2
  Else If (d%imn /= 0) Then
     valid_pdate = valid_pdate + 1
     If (d%idy < 0 .or. d%idy > mdm) valid_pdate = valid_pdate + 2
  Else
     If (d%idy /= 0) Then
        valid_pdate = valid_pdate + 1
        If (d%idy < 1 .or. d%idy > ndays(d%iyr, d%imn)) valid_pdate = valid_pdate + 2
     End If
  End If
!
  Return
 End Function valid_pdate
!
!
!
 Function valid_pprd(p, isq)
!
! Checks that parsed period is valid
!
! On exit:
!    valid_pprd =  0 Valid
!    valid_pprd =  1 Invalid length of period
!    valid_pprd =  2 Invalid start date for dekad
!    valid_pprd =  3 Invalid end date for dekad
!
! Function type
  Integer :: valid_pprd
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - sequencing -
!
  Type(pprd), Intent(In) :: p ! - period -
!
! Locals
!
! Local scalars
  Integer :: i ! - date difference -
!
  Logical :: lslim_bk ! - backup use season limits? -
!
! Executable Statements
!
! Check validity of pentads and weeks
  valid_pprd = 0
  Select Case (isq)
   Case (5, 7)
     lslim_bk = lslim
     lslim = .false.
     i = 1 + date_diff(p%sdate, p%edate, 1)
     lslim = lslim_bk
     If (i /= isq .and. i /= 1) valid_pprd = 1
!
! Check validity of dekads
   Case (10)
     If (p%sdate%iyr /= p%edate%iyr .or. p%sdate%imn /= p%edate%imn) Then
        valid_pprd = 1
        Return
     End If
     Select Case (p%sdate%idy)
      Case (1)
        If (p%edate%idy /= 10) valid_pprd = 3
      Case (11)
        If (p%edate%idy /= 20) valid_pprd = 3
      Case (21)
        If (p%edate%idy /= ndays(p%edate%iyr, p%edate%imn)) valid_pprd = 3
      Case Default
        valid_pprd = 2
     End Select
  End Select
!
  Return
 End Function valid_pprd
!
!
!
 Function check_pprd(cs1, cs2, p1, p2, isq) &
          Result (check)
!
! Checks that periods are consistent
!
! On exit:
!    check = 0 Successful
!    check = 1 Inconsistent lengths of seasons
!    check = 2 Invalid length of second period
!    check = 3 Invalid start date for second dekad
!    check = 4 Invalid end date for second dekad
!    check = 5 Invalid spacing for subseasonal data
!
! Modules
  Use data_time_constants, Only: isq_dy, isq_yr
!
! Function type
  Integer :: check
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - sequence -
!
  Character(Len=*), Intent(In) :: cs1 ! - first season -
  Character(Len=*), Intent(In) :: cs2 ! - second season -
!
  Type(pprd), Intent(In) :: p1 ! - first period -
  Type(pprd), Intent(In) :: p2 ! - second period -
!
! Locals
!
! Local scalars
  Logical :: lslim_bk ! - backup use seasonal limits? -
!
! Executable Statements
!
! Check that date differences are identical
  check = 0
  Select Case (isq)
   Case (isq_yr)
     If (cs1 /= cs2) check = 1
   Case (1:)
     check = valid_pprd(p2, isq)
     Select Case (check)
      Case (0, 1)
        check = 0
      Case Default
        check = check + 1
        Return
     End Select
     If (isq /= 10) Then
        lslim_bk = lslim
        lslim = .false.
        If (date_diff(p1%sdate, p1%edate, isq_dy) /= date_diff(p2%sdate, p2%edate, isq_dy)) check = 1
        lslim = lslim_bk
     End If
  End Select
!
  Return
 End Function check_pprd
!
!
!
 Subroutine check_lead (xfdate, yfdate, xmdate, it1, isq, isynch, ifail)
!
! Warns of unreasonable lead-times
!
! On exit:
!    ifail = 0 All checks passed
!    ifail = 1 negative lead
!    ifail = 2 synchronous predictors and predictands with synchronous predictors switched off
!    ifail = 3 large lead-time
!    ifail = 4 forecast and target periods do not match
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_yr, nmn
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: it1    ! - index of first X date of interest -
  Integer, Intent(In) :: isq    ! - time sequencing -
  Integer, Intent(In) :: isynch ! - synchronous predictors flag -
!
  Type(pdate), Intent(In) :: xfdate ! - first X date of interest -
  Type(pdate), Intent(In) :: yfdate ! - first Y date of interest -
  Type(pdate), Intent(In) :: xmdate ! - first X model date -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: isq_u  ! - used time-sequencing -
  Integer :: isyn_u ! - used synchronisation flag -
  Integer :: lead1  ! - lead-time -
  Integer :: lead2  ! - target period offset -
!
  Logical :: lslim_bk ! - backup use seasonal limts? -
!
! Executable Statements
!
! Check for yearly sequencing
  Select Case (isq)
   Case (isq_yr)
     isq_u = isq_mn
   Case Default
     isq_u = isq
  End Select
!
! Calculate lead
  lslim_bk = lslim
  lslim = .true.
  If (xmdate == 0) Then
     lead1 = date_diff(xfdate, yfdate, isq_u)
     lead2 = 0
     lead = lead1
     isyn_u = isynch
  Else
     lead1 = date_diff(xmdate+(it1-1), yfdate, isq_u)
     lead2 = date_diff(xfdate, yfdate, isq_u)
     lead = lead2
     isyn_u = 1
  End If
  lslim = lslim_bk
!
! Check for negative leads
  If (lead1 < 0) Then
     ifail = 1
!
! Check for synchronous predictors
  Else If (lead == 0) Then
     ifail = 2*(1 - isyn_u)
!
! Check for large leads
  Else
     ifail = 0
     Select Case (isq)
      Case (isq_yr) ! - yearly -
        If (lead1/nmn > 0) ifail = 3
        If (xmdate == 0) Then
           lead = yfdate%iyr - xfdate%iyr
        Else
           lead = yfdate%iyr - (xmdate%iyr + it1 - 1)
        End If
      Case Default
        If (lead1 > 1) ifail = 3
     End Select
  End If
!
! Check for offset
  If (ifail == 0) Then
     If (lead2 /= 0) ifail = 4
  End If
!
  Return
 End Subroutine check_lead
!
!
!
 Recursive Function date_diff(d1, d2, isq) &
           Result (d)
!
! Calculates difference in dates
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr
!
! Function type
  Integer :: d
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - date sequence -
!
  Type(pdate), Intent(In) :: d1 ! - first date -
  Type(pdate), Intent(In) :: d2 ! - second date -
!
! Locals
!
! Local scalars
  Integer :: idadj  ! - date adjustment -
  Integer :: incd   ! - date increment -
  Integer :: iya    ! - current year -
  Integer :: ima    ! - current month -
  Integer :: imb    ! - next month -
  Integer :: iendy  ! - end year -
  Integer :: iendm  ! - end month -
  Integer :: irds   ! - remaining days -
  Integer :: isq_bk ! - backup time sequencing -
  Integer :: k      ! - time index -
  Integer :: nms    ! - number of months -
  Integer :: nyr    ! - number of years -
!
  Logical :: lsl ! - season limits? -
!
  Type(pdate) :: d0 ! - first date -
!
  Type(pprd) :: cssnlim ! - current season limits -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod, Sign
!
! Executable Statements
!
! Return if dates are identical
  If (d2 == d1) Then
     d = 0
     Return
  End If
!
! Compare yearly dates
  Select Case (isq)
   Case (isq_yr)
     d = d2%iyr - d1%iyr
     If (d1%imn > 0 .and. d2%imn > 0) Then
        If (d2%imn < d1%imn) Then
           d = d - 1
        Else If (d2%imn == d1%imn) Then
           If (d1%idy > 0 .and. d2%idy > 0) Then
              If (d2%idy < d1%idy) d = d - 1
           End If
        End If
     End If
!
! Compare seasonal dates
   Case (isq_sn)
     d = date_diff(d1, d2, isq_mn) ! - calculate difference in months -
     If (d == 0) Return
     nyr = d/nmn ! - calculate difference in whole years -
     nms = Mod(d, nmn) ! - identify partial year differences -
     d = nyr*nss ! - multiply whole years by number of seasons per year -
     If (nms == 0) Return
! - calculate number of seasons in partial year difference -
     ima = 0
     Do k = 1, nss
        If (d1%imn == ims(k)) Then
           ima = k
           Exit
        End If
     End Do
     If (ima == 0) Then
        d = 0
        Return
     End If
     imb = 0
     Do k = 1, nss
        If (d2%imn == ims(k)) Then
           imb = k
           Exit
        End If
     End Do
     If (imb == 0) Then
        d = 0
        Return
     End If
     If (nms > 0) Then
        nms = imb - ima
        If (nms < 0) nms = nms + nss
        d = d + nms
     Else
        nms = ima - imb
        If (nms < 0) nms = nms + nss
        d = d - nms
     End If
!
! Compare monthly dates
   Case (isq_mn, isq_so)
     d = (d2%iyr - d1%iyr)*nmn + d2%imn - d1%imn
!
! Compare daily, bi-weekly, pentadal, and weekly dates, and other sub-seasonal other than dekadal
! - determine direction of difference -
   Case Default
     If (d2 > d1) Then
        incd = isq
     Else
        incd = -isq
     End If
     lsl = lslim
     If (lsl) lsl = (ssnlims /= 0) .and. (ssnlims%sdate /= ssnlims%edate)
! - no seasonal limits -
     If (.not.lsl) Then
        d = d2%idy - d1%idy
        If (d2%iyr /= d1%iyr .or. d2%imn /= d1%imn) Then
           If (incd > 0) Then
              iya = d1%iyr
              ima = d1%imn
              iendy = d2%iyr
              iendm = d2%imn
           Else
              iya = d2%iyr
              ima = d2%imn
              iendy = d1%iyr
              iendm = d1%imn
              d = -d
           End If
           Do
              d = d + ndays(iya, ima)
              ima = add_to_month(ima, 1, &
                    iyr=iya)
              If (iya == iendy .and. ima == iendm) Exit
           End Do
        End If
        Select Case (isq)
         Case Default
           d = d/incd
         Case (3, 4)
           irds = Mod(d, 7)
           d = 2*(d/7)
           Select Case (irds)
            Case (0:2)
              Continue
              d = d/7
            Case Default
              d = d + Sign(irds, incd)
           End Select
        End Select
! - season limits -
     Else
        isq_bk = iseq
        iseq = isq
        cssnlim = ssnlims
        d0 = d1
        If (d1%imn > ssnlims%sdate%imn .or. (d1%imn == ssnlims%sdate%imn .and. d1%idy >= ssnlims%sdate%idy)) Then
           cssnlim%sdate%iyr = d0%iyr
        Else
           cssnlim%sdate%iyr = d0%iyr - 1
        End If
        If (ssnlims%sdate < ssnlims%edate) Then ! - season does not span year-end -
           cssnlim%edate%iyr = cssnlim%sdate%iyr
        Else ! - season spans year-end -
           cssnlim%edate%iyr = cssnlim%sdate%iyr + 1
        End If
        If (incd > 0) Then ! - second date is later than first -
           d = 1
           If (d0 < cssnlim%sdate) Then ! - adjust if date is outside season limits -
              d0 = cssnlim%sdate
              If (d0 < d2) d = d + 1
           End If
           If (d0 < d2) Then
              Do
                 d0 = d0 + 1
                 If (d0 > cssnlim%edate) Then
                    cssnlim%sdate%iyr = cssnlim%sdate%iyr + 1
                    cssnlim%edate%iyr = cssnlim%edate%iyr + 1
                 End If
                 If (d0 >= d2) Exit
                 d = d + 1
              End Do
           End If
        Else ! - second date is earlier than first -
           d = -1
           If (d0 > cssnlim%edate) Then ! - adjust if date is outside season limits -
              idadj = date_diff(cssnlim%sdate, cssnlim%edate, isq)
              If (isq > 1) idadj = idadj - 1
              d0 = cssnlim%sdate + idadj
              If (d0 > d2) d = d - 1
           End If
           If (d0 > d2) Then
              Do
                 d0 = d0 - 1
                 If (d0 < cssnlim%sdate) Then
                    cssnlim%sdate%iyr = cssnlim%sdate%iyr - 1                                                        
                    cssnlim%edate%iyr = cssnlim%edate%iyr - 1
                 End If
                 If (d0 <= d2) Exit
                 d = d - 1
              End Do
           End If
        End If
        iseq = isq_bk
     End If
! - dekadal -
   Case (10)
     d = date_diff(d1, d2, -1)*3 + (d2%idy - d1%idy)/10
  End Select
!
  Return
 End Function date_diff
!
!
!
 Subroutine get_sequence (p1, p2, isq, ifail)
!
! Identifies period sequencing. On exit:
!    isq = isq_yr if years are consecutive
!    isq = isq_sn if seasons are consecutive
!    isq = isq_so if seasons are overlapping and are monthly consecutive
!    isq = isq_mn if months are consecutive
!    isq =      1 if days are consecutive
!    isq =      5 if pentads are consecutive
!    isq =      7 if weeks are consecutive
!    isq =     10 if dekads are consecutive
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Dates are identical
!    ifail =  2 Second date is before first
!    ifail =  3 Unrecognized sequencing
!    ifail =  4 Inconsistent period length
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr, mdy, ndw
!
! Arguments
!
! Input scalars
  Type(pprd), Intent(In) :: p1 ! - period 1 -
  Type(pprd), Intent(In) :: p2 ! - period 2 -
!
! Output scalars
  Integer, Intent(Out) :: isq   ! - date sequence indicator -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: isq1 ! - length of first period -
  Integer :: isq2 ! - length of second period -
!
! Executable Statements
!
! Check for start date sequences
  Call get_sequence_dates (p1%sdate, p2%sdate, isq, ifail)
  If (ifail /= 0) Return
!
! Check for consistent periods
  If (isq < 1) Then
     If (p1%sdate == p1%edate .and. p2%sdate == p2%edate) Then
        If (isq == isq_sn) isq = isq_mn
        ifail = 0
     Else
        Call get_sequence_dates (p1%sdate, p1%edate, isq1, ifail) ! - check length of period -
        Call get_sequence_dates (p2%sdate, p2%edate, isq2, ifail)
        If (isq1 == isq_mn .and. isq2 == isq_mn) Then
           ifail = 0
        Else If (isq1 < 0 .and. isq2 < 0) Then ! - set overlapping seasons as monthly sequencing -
           If (isq == isq_mn) isq = isq_so
           ifail = 0
        Else If (p1%sdate%imn == p2%sdate%imn .and. &
                 p1%sdate%idy == p2%sdate%idy) Then ! - allow annual sequencing of non-calendar-monthly periods -
           If (p1%edate%imn == p2%edate%imn .and. &
               p1%edate%idy == p2%edate%idy) Then ! - same end dates -
              ifail = 0
           Else If (date_diff(p1%sdate, p1%edate, 1) == date_diff(p2%sdate, p2%edate, 1)) Then ! - same length of period -
              ifail = 0
           Else
              ifail = 3
           End If
        Else
           ifail = 3
        End If
     End If
!
! Check for standard periods of days
  Else If (isq > 1) Then
     Select Case (isq)
      Case (1, 3, 4, 5, 7) ! - check for half-weeks, pentads and weeks -
        Call get_sequence_dates (p1%sdate, p1%edate, isq1, ifail) ! - check length of period -
        isq1 = isq1 + 1
        Call get_sequence_dates (p2%sdate, p2%edate, isq2, ifail)
        isq2 = isq2 + 1
        If (isq1 == isq2) Then
           ifail = 0
        Else
           isq = 0
           ifail = 4
        End If
      Case (8:13) ! - check for dekads -
        If (((p1%sdate%idy ==  1 .and. p1%edate%idy == 10) .or.  &
             (p1%sdate%idy == 11 .and. p1%edate%idy == 20) .or.  &
             (p1%sdate%idy == 21 .and. p1%edate%idy == ndays(p1%sdate%iyr, p1%sdate%imn))) .and.  &
            ((p2%sdate%idy ==  1 .and. p2%edate%idy == 10) .or.  &
             (p2%sdate%idy == 11 .and. p2%edate%idy == 20) .or.  &
             (p2%sdate%idy == 21 .and. p2%edate%idy == ndays(p2%sdate%iyr, p2%sdate%imn)))) Then
           isq = 10
           ifail = 0
        Else
           isq = 0
           ifail = 3
        End If
      Case (mdy-ndw:mdy+ndw-1) ! - check for approximate annual sequencing -
        If (day_of_year(p1%sdate) == day_of_year(p2%sdate) .and. &
            day_of_year(p1%edate) == day_of_year(p2%edate)) Then ! - allow for identical day-of-year -
           isq = isq_yr
           ifail = 0
        Else If (day_of_week(p1%sdate) == day_of_week(p2%sdate) .and. &
                 day_of_week(p1%edate) == day_of_week(p2%edate)) Then ! - allow for identical day-of-week -
           isq = isq_yr
           ifail = 0
        Else
           isq = 0
           ifail = 3
        End If
      Case Default
        isq = 0
        ifail = 3
     End Select
  End If
!
  Return
!
 Contains
!
!
  Subroutine get_sequence_dates (d1, d2, isq, ifail)
!
! Identifies date sequencing.
! On exit:
!    isq = isq_yr if years are consecutive
!    isq = isq_mn if seasons are consecutive
!    isq = isq_sn if months are consecutive
!    isq >      0 if days are n days apart
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Dates are identical
!    ifail =  2 Second date is before first
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - date 1 -
  Type(pdate), Intent(In) :: d2 ! - date 2 -
!
! Output scalars
  Integer, Intent(Out) :: isq   ! - date sequence indicator -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Check for date sequences
  isq = 0
  If (d2%iyr < d1%iyr) Then
     ifail = 2
     Return
  End If
! - day is constant -
  If (d2%idy == d1%idy) Then
     If (d2%imn == d1%imn) Then ! - same month -
        If (d2%iyr > d1%iyr) Then
           isq = isq_yr
           ifail = 0
        Else
           ifail = 1
        End If
     Else If (d2%imn == d1%imn+1) Then ! - next month -
        If (d2%iyr == d1%iyr) Then
           isq = isq_mn
           ifail = 0
        Else
           isq = isq_sn
           ifail = 0
        End If
     Else If (d1%imn == nmn .and. d2%imn == 1) Then ! - next month, next year -
        If (d2%iyr == d1%iyr+1) Then
           isq = isq_mn
           ifail = 0
        Else
           isq = isq_sn
           ifail = 0
        End If
     Else ! - multiple month gap -
        If (date_diff(d1, d2, -1) > 0) Then
           isq = isq_sn
           ifail = 0
        Else
           ifail = 2
        End If
     End If
! - check for consecutive days, pentads, weeks, and dekads -
  Else
     isq = date_diff(d1, d2, 1)
     If (isq > 0) Then
        ifail = 0
     Else
        ifail = 2
     End If
  End If
!
  Return
  End Subroutine get_sequence_dates
 End Subroutine get_sequence
!
!
!
 Function get_cdate_date(d, ifmt) &
          Result (cdate)
!
! Formats date as a character string
!
! If ifmt = 1, ISO format
! If ifmt = 2, common format
! If ifmt = 3, GrADS format
! If ifmt = 4, International text format
!
! Modules
  Use data_text,           Only: ilang
  Use data_io_constants,   Only: ldat
  Use data_time_constants, Only: cmon
!
! Function type
  Character(Len=ldat) :: cdate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ifmt ! - format indicator -
!
  Type(pdate), Intent(In) :: d ! - date -
!
! Executable Statements
!
! Create date
  Select Case (ifmt)
! - ISO format -
   Case (1)
     If (d%imn > 0) Then
        If (d%idy > 0) Then
           Write (Unit=cdate, Fmt='(I4,A,I2.2,A,I2.2)') d%iyr, '-', d%imn, '-', d%idy
        Else
           Write (Unit=cdate, Fmt='(I4,A,I2.2)') d%iyr, '-', d%imn
        End If
     Else
        Write (Unit=cdate, Fmt='(I4)') d%iyr
     End If
! - common format -
   Case (2)
     If (d%imn > 0) Then
        If (d%idy > 0) Then
           Write (Unit=cdate, Fmt='(I2.2,3A,I4)') d%idy, ' ', cmon(d%imn,ilang), ' ', d%iyr
        Else
           Write (Unit=cdate, Fmt='(2A,I4)') cmon(d%imn,ilang), ' ', d%iyr
        End If
     Else
        Write (Unit=cdate, Fmt='(I4)') d%iyr
     End If
! - GrADS format -
   Case (3)
     If (d%imn > 0) Then
        If (d%idy > 0) Then
           Write (Unit=cdate, Fmt='(I2.2,A,I4)') d%idy, cmon(d%imn,1), d%iyr
        Else
           Write (Unit=cdate, Fmt='(2A,I4)') '1', cmon(d%imn,1), d%iyr
        End If
     Else
        Write (Unit=cdate, Fmt='(2A,I4)') '1', cmon(d%imn,1), d%iyr
     End If
! - International text format -
   Case (4)
     If (d%imn > 0) Then
        If (d%idy > 0) Then
           Write (Unit=cdate, Fmt='(I4,3A,I2.2)') d%iyr, ' ', cmon(d%imn,ilang), ' ', d%idy
        Else
           Write (Unit=cdate, Fmt='(I4,2A)') d%iyr, ' ', cmon(d%imn,ilang)
        End If
     Else
        Write (Unit=cdate, Fmt='(I4)') d%iyr
     End If
  End Select
!
  Return
 End Function get_cdate_date
!
!
!
 Function get_cdate_period(p, lsn, ifmt) &
          Result (cdate)
!
! Formats period as a character string
!
! If ifmt = 1, ISO format
! If ifmt = 2, common format
! If ifmt = 3, GrADS format
! If ifmt = 4, International text format
!
! Modules
  Use data_text,           Only: ilang
  Use data_io_constants,   Only: lprd
  Use data_time_constants, Only: cma, cmon
!
! Function type
  Character(Len=lprd) :: cdate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lsn  ! - length of season (for seasonally sequenced data) -
  Integer, Intent(In) :: ifmt ! - format indicator -
!
  Type(pprd), Intent(In) :: p ! - period -
!
! Executable Statements
!
! Create date
  If (p%edate == p%sdate) Then
     cdate = get_cdate_date(p%sdate, ifmt)
  Else If (p%sdate == 0) Then
     cdate = '0'
  Else
     Select Case (ifmt)
! - ISO format -
      Case (1)
        If (p%edate%iyr > p%sdate%iyr) Then
           cdate = Trim(get_cdate_date(p%sdate, ifmt))//'/'//get_cdate_date(p%edate, ifmt)
        Else If (p%edate%imn /= p%sdate%imn) Then
           If (p%sdate%idy > 0) Then
              Write (Unit=cdate, Fmt='(I4,4(A,I2.2))') &
                 p%sdate%iyr, '-', p%sdate%imn, '-', p%sdate%idy, '/', p%edate%imn, '-', p%edate%idy
           Else
              Write (Unit=cdate, Fmt='(I4,2(A,I2.2))') p%sdate%iyr, '-', p%sdate%imn, '/', p%edate%imn
           End If
        Else
           Write (Unit=cdate, Fmt='(I4,3(A,I2.2))') &
              p%sdate%iyr, '-', p%sdate%imn, '-', p%sdate%idy, '/', p%edate%idy
        End If
! - common format -
      Case (2)
        If (iseq < 0) Then
           Select Case (lsn)
            Case (:1) ! - single month -
              If (p%edate%iyr == p%sdate%iyr) Then
                 If (p%sdate%idy == 0) Then
                    If (p%sdate%imn > 0) Then
                       Write (Unit=cdate, Fmt='(2A,I4)') cma(ilang)(p%sdate%imn:p%edate%imn), ' ', p%sdate%iyr
                    Else
                       Write (Unit=cdate, Fmt='(I4)') p%sdate%iyr
                    End If
                 Else
                    If (p%edate%imn == p%sdate%imn) Then
                       Write (Unit=cdate, Fmt='(2(I2.2,A),2A,I4)') &
                          p%sdate%idy, '-', p%edate%idy, ' ', cmon(p%sdate%imn,ilang), ' ', p%sdate%iyr
                    Else
                       Write (Unit=cdate, Fmt='(2(I2.2,3A),I4)') &
                          p%sdate%idy, ' ', cmon(p%sdate%imn,ilang), ' - ', &
                          p%edate%idy, ' ', cmon(p%edate%imn,ilang), ' ', p%sdate%iyr
                    End If
                 End If
              Else
                 If (p%sdate%idy == 0) Then
                    Write (Unit=cdate, Fmt='(2A,I4,A,I2.2)') &
                       cma(ilang)(p%sdate%imn:p%edate%imn+nmn), ' ', p%sdate%iyr, '/', Mod(p%edate%iyr, 100)
                 Else
                    cdate = Trim(get_cdate_date(p%sdate, ifmt))//' - '//get_cdate_date(p%edate, ifmt)
                 End If
              End If
            Case (2:4) ! - standard season -
              If (p%edate%iyr == p%sdate%iyr) Then
                 Write (Unit=cdate, Fmt='(2A,I4)') &
                    cma(ilang)(p%sdate%imn:p%edate%imn), ' ', p%sdate%iyr
              Else
                 Write (Unit=cdate, Fmt='(2A,I4,A,I2.2)') &
                    cma(ilang)(p%sdate%imn:p%edate%imn+nmn), ' ', p%sdate%iyr, '/', Mod(p%edate%iyr, 100)
              End If
            Case (5:) ! - long season -
              If (p%edate%iyr == p%sdate%iyr) Then
                 Write (Unit=cdate, Fmt='(4A,I4)') &
                    cmon(p%sdate%imn,ilang), '-', cmon(p%edate%imn,ilang), ' ', p%sdate%iyr
              Else
                 Write (Unit=cdate, Fmt='(4A,I4,A,I2.2)') &
                    cmon(p%sdate%imn,ilang), '-', cmon(p%edate%imn,ilang), ' ', p%sdate%iyr, '/', Mod(p%edate%iyr, 100)
              End If
           End Select
        Else
           If (p%sdate == p%edate) Then ! - single day -
              cdate = Trim(get_cdate_date(p%sdate, ifmt))
           Else ! - multiple days -
              If (p%edate%iyr == p%sdate%iyr) Then
                 If (p%sdate%imn == p%edate%imn) Then
                    Write (Unit=cdate, Fmt='(2(I2.2,A),I4)') &
                         p%sdate%idy, ' - ', p%edate%idy, ' '//cmon(p%sdate%imn,ilang)//' ', p%sdate%iyr
                 Else
                    Write (Unit=cdate, Fmt='(2(I2.2,A),I4)') &
                         p%sdate%idy, ' '//cmon(p%sdate%imn,ilang)//' - ', &
                         p%edate%idy, ' '//cmon(p%edate%imn,ilang)//' ', p%sdate%iyr
                 End If
              Else
                 Write (Unit=cdate, Fmt='(I2.2,A,I4,A,I2.2,A,I4)') &
                      p%sdate%idy, ' '//cmon(p%sdate%imn,ilang)//' ', p%sdate%iyr, ' - ', &
                      p%edate%idy, ' '//cmon(p%edate%imn,ilang)//' ', p%edate%iyr
              End If
           End If
        End If
! - GrADS format -
      Case (3)
        cdate = Trim(get_cdate_date(p%sdate, ifmt))
! - International text format -
      Case (4)
        If (iseq < 0) Then
           Select Case (lsn)
            Case (:1) ! - single month -
              If (p%edate%iyr == p%sdate%iyr) Then
                 If (p%sdate%idy == 0) Then
                    If (p%sdate%imn > 0) Then
                       Write (Unit=cdate, Fmt='(I4,2A)') p%sdate%iyr, ' ', cmon(p%sdate%imn,ilang)
                    Else
                       Write (Unit=cdate, Fmt='(I4)') p%sdate%iyr
                    End If
                 Else
                    If (p%edate%imn == p%sdate%imn) Then
                       Write (Unit=cdate, Fmt='(I4,2A,2(A,I2.2))') &
                          p%sdate%iyr, ' ', cmon(p%sdate%imn,ilang), ' ', p%sdate%idy, '/', p%edate%idy
                    Else
                       Write (Unit=cdate, Fmt='(I4,2(3A,I2.2))') &
                          p%sdate%iyr, ' ', cmon(p%sdate%imn,ilang), ' ', p%sdate%idy, ' / ', &
                                            cmon(p%edate%imn,ilang), ' ', p%edate%idy
                    End If
                 End If
              Else
                 If (p%sdate%idy == 0) Then
                    Write (Unit=cdate, Fmt='(I4,A,I2.2,2A)') &
                       p%sdate%iyr, '/', Mod(p%edate%iyr, 100), ' ', cmon(p%sdate%imn,ilang)
                 Else
                    cdate = Trim(get_cdate_date(p%sdate, ifmt))//' / '//get_cdate_date(p%edate, ifmt)
                 End If
              End If
            Case (2:4) ! - standard season -
              If (p%edate%iyr == p%sdate%iyr) Then
                 Write (Unit=cdate, Fmt='(I4,2A)') &
                    p%sdate%iyr, ' ', cma(ilang)(p%sdate%imn:p%edate%imn)
              Else
                 Write (Unit=cdate, Fmt='(I4,A,I2.2,2A)') &
                    p%sdate%iyr, '/', Mod(p%edate%iyr, 100), ' ', cma(ilang)(p%sdate%imn:p%edate%imn+nmn)
              End If
            Case (5:) ! - long season -
              If (p%edate%iyr == p%sdate%iyr) Then
                 Write (Unit=cdate, Fmt='(I4,4A)') &
                    p%sdate%iyr, ' ', cmon(p%sdate%imn,ilang), '-', cmon(p%edate%imn,ilang)
              Else
                 Write (Unit=cdate, Fmt='(I4,A,I2.2,4A)') &
                    p%sdate%iyr, '/', Mod(p%edate%iyr, 100), ' ', cmon(p%sdate%imn,ilang), '-', cmon(p%edate%imn,ilang)
              End If
           End Select
        Else
           If (p%sdate == p%edate) Then ! - single day -
              cdate = Trim(get_cdate_date(p%sdate, ifmt))
           Else ! - multiple days -
              If (p%edate%iyr == p%sdate%iyr) Then
                 If (p%sdate%imn == p%edate%imn) Then
                    Write (Unit=cdate, Fmt='(I4,2(A,I2.2))') &
                         p%sdate%iyr, ' '//cmon(p%sdate%imn,ilang), p%sdate%idy, '/', p%edate%idy
                 Else         
                    Write (Unit=cdate, Fmt='(I4,2(A,I2.2))') &
                         p%sdate%iyr, ' '//cmon(p%sdate%imn,ilang)//' ', p%sdate%idy, &
                                      '/'//cmon(p%edate%imn,ilang)//' ', p%edate%idy
                 End If
              Else
                 Write (Unit=cdate, Fmt='(I4,A,I2.2,A,I4,A,I2.2)') &
                         p%sdate%iyr, ' '//cmon(p%sdate%imn,ilang)//' ', p%sdate%idy, '/', &
                         p%edate%iyr, ' '//cmon(p%edate%imn,ilang)//' ', p%edate%iyr
              End If
           End If
        End If
     End Select
  End If
!
  Return
 End Function get_cdate_period
!
!
!
 Function get_cdate_lsn(d1, lper, ifmt) &
          Result (cdate)
!
! Formats period as a character string given start date and length of period
!
! If ifmt = 1, ISO format
! If ifmt = 2, common format
! If ifmt = 3, GrADS format
! If ifmt = 4, International text format
!
! Modules
  Use data_io_constants,   Only: lprd
  Use data_time_constants, Only: isq_dy, isq_mn, isq_sn, isq_so, isq_yr
!
! Function type
  Character(Len=lprd) :: cdate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lper ! - length of period -
  Integer, Intent(In) :: ifmt ! - format indicator -
!
  Type(pdate), Intent(In) :: d1 ! - start date -
!
! Locals
!
! Local scalars
  Integer :: iseq_bk ! - backup time sequencing -
  Integer :: lp      ! - length of period -
!
  Logical :: lreset   ! - reset time-sequencing? -
  Logical :: lsl      ! - season limits? -
  Logical :: lslim_bk ! - backup use season limits? -
!
  Type(pdate) :: d2 ! - end date -
!
  Type(pprd) :: p ! - date -
!
! Executable Statements
!
! Calculate end date
  Select Case (iseq)
   Case (isq_sn, isq_so, isq_yr)
     lp = lper
     iseq_bk = iseq
     iseq = isq_mn
     lreset = .true.
   Case (isq_mn)
     lp = lper
     lreset = .false.
   Case Default
     lp = lper*iseq
     iseq_bk = iseq
     iseq = isq_dy
     lreset = .true.
  End Select
  lsl = .not.(ssnlims%sdate == ssnlims%edate)
  If (lsl) Then
     lslim_bk = lslim
     lslim = .false.
  End If
  d2 = d1 + (lp - 1)
  If (lsl) lslim = lslim_bk
  p%sdate = d1
  p%edate = d2
  If (lreset) iseq = iseq_bk
!
! Construct date
  cdate = get_cdate_period(p, 1, ifmt)
!
  Return
 End Function get_cdate_lsn
!
!
!
 Function get_cdate_lsns(d1, lb, lsns, ifmt) &
          Result (cdate)
!
! Formats period as a character string given persistence component
!
! If ifmt = 1, ISO format
! If ifmt = 2, common format
! If ifmt = 3, GrADS format
! If ifmt = 4, International text format
!
! Modules
  Use data_text,           Only: ilang
  Use data_io_constants,   Only: lprd
  Use data_time_constants, Only: cma, isq_mn, isq_sn, isq_yr
!
! Function type
  Character(Len=lprd) :: cdate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lb   ! - lower bound on lsns -
  Integer, Intent(In) :: ifmt ! - format indicator -
!
  Type(pdate), Intent(In) :: d1 ! - start date (of persistence component, not forecast component) -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsns ! - length of seasons -
!
! Locals
!
! Local scalars
  Integer :: iseq_bk ! - backup time sequencing -
!
  Logical :: lsqc ! - time sequencing changed? -
!
  Type(pprd) :: p1 ! - persistence component dates -
  Type(pprd) :: p2 ! - forecast component dates -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Sum
!
! Executable Statements
!
! Construct date
  Select Case (ifmt)
! - ISO format -
   Case (1)
     cdate = get_cdate_lsn(d1, Sum(lsns(:)), ifmt)
! - common format -
   Case (2)
     If (lb == 1 .or. lsns(lb) == 0) Then
        cdate = get_cdate_lsn(d1, Sum(lsns(:)), ifmt)
     Else If (iseq < 0) Then
        Select Case (iseq)
         Case (isq_sn, isq_yr)
           iseq_bk = iseq
           iseq = isq_mn
           lsqc = .true.
         Case Default
           lsqc = .false.
        End Select
        p1%sdate = d1
        p1%edate = p1%sdate + (lsns(lb) - 1)
        p2%sdate = p1%edate + 1
        p2%edate = p2%sdate + (lsns(lb+1) - 1)
        If (lsqc) iseq = iseq_bk
        If (p2%edate%iyr == p1%sdate%iyr) Then
           Write (Unit=cdate, Fmt='(5A,I4)') &
              '(', cma(ilang)(p1%sdate%imn:p1%edate%imn), ')', &
                   cma(ilang)(p2%sdate%imn:p2%edate%imn),' ', p1%sdate%iyr
        Else
           If (p1%edate%iyr == p1%sdate%iyr) Then
              If (p2%edate%iyr == p2%sdate%iyr) Then
                 Write (Unit=cdate, Fmt='(5A,I4,A,I2.2)') &
                    '(', cma(ilang)(p1%sdate%imn:p1%edate%imn), ')', &
                         cma(ilang)(p2%sdate%imn:p2%edate%imn),' ', p1%sdate%iyr, '/', Mod(p2%edate%iyr, 100)
              Else
                 Write (Unit=cdate, Fmt='(5A,I4,A,I2.2)') &
                    '(', cma(ilang)(p1%sdate%imn:p1%edate%imn), ')', &
                         cma(ilang)(p2%sdate%imn:p2%edate%imn+nmn),' ', p1%sdate%iyr, '/', Mod(p2%edate%iyr, 100)
              End If
           Else
              If (p2%edate%iyr == p2%sdate%iyr) Then
                 Write (Unit=cdate, Fmt='(5A,I4,A,I2.2)') &
                    '(', cma(ilang)(p1%sdate%imn:p1%edate%imn+nmn), ')', &
                         cma(ilang)(p2%sdate%imn:p2%edate%imn),' ', p1%sdate%iyr, '/', Mod(p2%edate%iyr, 100)
              Else
                 Write (Unit=cdate, Fmt='(5A,I4,A,I2.2)') &
                    '(', cma(ilang)(p1%sdate%imn:p1%edate%imn+nmn), ')', &
                         cma(ilang)(p2%sdate%imn:p2%edate%imn+nmn),' ', p1%sdate%iyr, '/', Mod(p2%edate%iyr, 100)
              End If
           End If
        End If
     Else ! - persistence not yet fully implemented for daily to sub-seasonal scales -
        cdate = get_cdate_lsn(d1, Sum(lsns(:)), ifmt)
     End If
! - GrADS format -
   Case (3)
     cdate = get_cdate_lsn(d1, Sum(lsns(:)), ifmt)
! - International text format -
   Case (4)
     cdate = get_cdate_lsn(d1, Sum(lsns(:)), ifmt)
  End Select
!
  Return
 End Function get_cdate_lsns
!
!
!
 Function get_cprd(d, d1, isq, lse, lb, lsn, nsn) &
          Result (cprd)
!
! Creates selected season as a character string
!
! Modules
  Use data_io_constants,   Only: lprd
  Use data_time_constants, Only: isq_dy, isq_sn, isq_so
!
! Function type
  Character(Len=lprd) :: cprd
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - time-sequence -
  Integer, Intent(In) :: lb  ! - lower-bound -
  Integer, Intent(In) :: lse ! - length of sub-season -
  Integer, Intent(In) :: nsn ! - number of seasons -
!
  Type(pdate), Intent(In) :: d  ! - beginning date of period -
  Type(pdate), Intent(In) :: d1 ! - beginning date of first period -
!
! Input scalars
  Integer, Dimension(lb:), Intent(In) :: lsn ! - length of seasons -
!
! Locals
!
! Local scalars
  Integer :: isnc   ! - current season -
  Integer :: isq_bk ! - time sequencing backup -
!
  Type(pprd) :: p ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod
!
! Executable Statements
!
! Identify and create season
  Select Case (isq)
   Case Default
     cprd = get_cdate(d, lb, lsn(lb:1), 2)
   Case (1:)
     p%sdate = d
     p%edate = d
     isq_bk = iseq
     iseq = isq_dy
     p%edate = p%edate + (lse - 1)
     iseq = isq_bk
     cprd = get_cdate(p, 0, 2)
   Case (isq_sn, isq_so)
     isnc = date_diff(d1, d, isq) + 1
     isnc = Mod(isnc, nsn)
     If (isnc == 0) isnc = nsn
     Select Case (lb)
      Case (1)
        cprd = get_cdate(d, lsn(isnc), 2)
      Case (0)
        cprd = get_cdate(d, lb, (/lsn(0), lsn(isnc)/), 2)
     End Select
  End Select
!
  Return
 End Function get_cprd
!
!
!
 Subroutine set_fdate (iafter, isq, bfdate, adate1, abdate, afdate, ifail)
!
! Sets appropriate start date
!
! Modules
  Use data_time_constants, Only: isq_yr
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iafter ! - A's date should be before (-1) equal to (0) or after (1) B's date? -
  Integer, Intent(In) :: isq    ! - date sequencing -
!
  Type(pdate), Intent(In) :: bfdate ! - B's start date -
  Type(pdate), Intent(In) :: adate1 ! - A's first date -
!
! Input/output scalars
  Type(pdate), Intent(InOut) :: abdate ! - A's beginning date -
  Type(pdate), Intent(InOut) :: afdate ! - A's start date -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: lag ! - lag -
!
! Executable Statements
!
! Determine current difference
  Select Case (iafter)
   Case (1)
     lag = date_diff(bfdate, afdate, isq)
     If (afdate%imn == bfdate%imn) lag = lag - 1
   Case (0)
     If (isq /= isq_yr .or. afdate%imn == bfdate%imn) Then
        lag = date_diff(bfdate, afdate, isq)
     Else
        lag = 1 + date_diff(bfdate, afdate, isq)
     End If
   Case (-1)
     lag = -date_diff(afdate, bfdate, isq)
     If (afdate%imn == bfdate%imn) lag = lag + 1
  End Select
!
! Adjust if A start date is not the same as or immediately after B start date
  Call reset_iseq (isq_new=isq)
  If (lag /= 0) Then
     afdate = afdate - lag
     abdate = abdate - lag
  End If
!
! Adjust if A start date is out of range
  ifail = 0
  Do
     If (abdate >= adate1) Exit
     abdate = abdate + 1
     afdate = afdate + 1
     ifail = 1
  End Do
  Call reset_iseq ()
!
  Return
 End Subroutine set_fdate
!
!
!
 Subroutine set_bdate (lb, lsn, nsn, sdate, fdate, bdate, &
       d1, d2, cprdb)
!
! Sets beginning date from first date
!
! Modules
  Use data_time_constants, Only: isq_mn
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lb  ! - lower-bound -
  Integer, Intent(In) :: nsn ! - number of seasons -
!
  Type(pdate), Intent(In) :: sdate ! - first date in file -
  Type(pdate), Intent(In) :: fdate ! - first date -
!
! - optional input scalars -
  Type(pdate), Intent(In), Optional :: d1 ! - start date -
!
! Input/output scalars
! - optional input/output scalars -
  Type(pdate), Intent(InOut), Optional :: d2 ! - end date -
!
! Output scalars
  Type(pdate), Intent(Out) :: bdate ! - beginning date -
!
! - optional output scalars -
  Character(Len=*), Intent(Out), Optional :: cprdb ! - period -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsn ! - length of seasons -
!
! Locals
!
! Local scalars
  Integer :: lse ! - length of subseason -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Sum
!
! Execuable Statements
!
! Set beginning date
  If (lb == 0) Then
     If (iseq < 0) Then
        Call reset_iseq (isq_new=isq_mn)
     Else
        Call reset_iseq (isq_new=iseq)
     End If
     bdate = fdate - lsn(0)
     If (Present(d2) .and. Present(d1)) d2 = d1 + (Sum(lsn(0:)) - 1)
     Call reset_iseq ()
  Else
     bdate = fdate
  End If
  If (Present(cprdb)) Then
     Select Case (iseq)
      Case Default
        lse = 0
      Case (1:2, 5:)
        lse = iseq
      Case (3:4)
        lse = 7
     End Select
     cprdb = get_cprd(bdate, sdate, iseq, lse, lb, lsn(lb:), nsn)
  End If
!
  Return
 End Subroutine set_bdate
!
!
!
 Subroutine adjust_pdate (d1, td, isq, lafter)
!
! Adjusts the date to be the first after or last before a target date
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - date sequence -
!
  Logical, Intent(In) :: lafter ! - date should be immediately after first date available? (before last date if false) - 
!
  Type(pdate), Intent(In) :: td ! - target date (first date available if lafter is true; last date available if lafter is false) -
!
! Input/output scalars
  Type(pdate), Intent(InOut) :: d1 ! - date to adjust -
!
! Locals
!
! Local scalars
  Integer :: iseq_bk ! - backup time sequencing -
  Integer :: lag     ! - lag -
!
! Executable Statements
!
! Adjust if A start date is not the earliest possible
  If (lafter) Then
     lag = date_diff(td, d1, isq)
     If (lag /= 0) Then
        iseq_bk = iseq
        iseq = isq
        d1 = d1 - lag
        iseq = iseq_bk
     End If
!
! Adjust if A start date is not the latest possible
  Else
     lag = date_diff(d1, td, isq)
     If (lag /= 0) Then
        iseq_bk = iseq
        iseq = isq
        d1 = d1 + lag
        iseq = iseq_bk
     End If
  End If
!
  Return
 End Subroutine adjust_pdate
!
!
!
 Subroutine adjust_pprd (p1, td, isq, lafter)
!
! Adjusts a period to be the first after or last before a target date
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - date sequence -
!
  Logical, Intent(In) :: lafter ! - period should be immediately after first date available? (before last date if false) - 
!
  Type(pdate), Intent(In) :: td ! - target date (first date available if lafter is true; last date available if lafter is false) -
!
! Input/output scalars
  Type(pprd), Intent(InOut) :: p1 ! - period to adjust -
!
! Locals
!
! Local scalars
  Integer :: iseq_bk ! - backup time sequencing -
  Integer :: lag     ! - lag -
!
! Executable Statements
!
! Adjust if A start date is not the earliest possible
  If (lafter) Then
     lag = date_diff(td, p1%sdate, isq)
     If (lag /= 0) Then
        iseq_bk = iseq
        iseq = isq
        p1%sdate = p1%sdate - lag
        p1%edate = p1%edate - lag
        iseq = iseq_bk
     End If
!
! Adjust if A start date is not the latest possible
  Else
     lag = date_diff(p1%edate, td, isq)
     If (lag /= 0) Then
        iseq_bk = iseq
        iseq = isq
        p1%sdate = p1%sdate + lag
        p1%edate = p1%edate + lag
        iseq = iseq_bk
     End If
  End If
!
  Return
 End Subroutine adjust_pprd
!
!
!
 Function pdate_to_idate(iunit, apdate) &
          Result (idate)
!
! Converts parsed date to relative date
!
! Function type
  Integer :: idate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iunit ! - units of date -
!
  Type(pdate), Intent(In) :: apdate ! - parsed date -
!
! Locals
!
! Local scalars
  Type(pdate) :: d0 ! - reference date -
!
! Executable Statements
!
! Convert to relative date
  d0 = date0
  If (iunit == 10) Then
     d0%idy = (d0%idy/10)*10 + 1
     If (d0%idy > 21) d0%idy = 21
  End If
  idate = date_diff(d0, apdate, iunit)
!
  Return
 End Function pdate_to_idate
!
!
!
 Function idate_to_pdate(idate) &
          Result (adate)
!
! Converts relative date to parsed date
!
! Function type
  Type(pdate) :: adate
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: idate ! - relative date -
!
! Executable Statements
!
! Convert to parsed date
  adate = date0 + idate
!
  Return
 End Function idate_to_pdate
End Module time
