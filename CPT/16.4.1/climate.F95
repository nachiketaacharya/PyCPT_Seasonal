! Author: Simon Mason
Module climate
!
! Modules
  Use data_numbers,      Only: rp
  Use data_io_constants, Only: lprd, lprd_v15
  Use time,              Only: pdate
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: calc_climates, check_climates, climatology, close_climate, get_clim_limits, init_climate, init_climatep, &
            label_climate, print_climatology, proj_read_climate, proj_write_climate, reset_climate, set_climate, valid_climate
!
! Arrays
!
! Integer arrays
  Integer, Dimension(:), Allocatable, Public :: ncu ! - number of used climatological cases -
!
! Real arrays
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: ave   ! - climatological averages -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: climx ! - model climatological (sorted) data -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: sdev  ! - climatological standard deviations -
!
  Real(Kind=rp), Dimension(:,:,:), Pointer, Public :: clim => Null() ! - climatological (sorted) response data -
!
! Scalars
!
! Integer scalars
  Integer, Public :: icex = 0  ! - climatological period lies outside the range of the training period flag -
!                                  icex = :
!                                     0 climatological period lies within the range of the training period -
!                                     1 climatological period starts before the training period -
!                                     2 climatological period ends after the training period -
!                                     3 climatological period spans the range of the training period -
!                                    -1 climatological period starts and ends before the training period with a gap -
!                                    -2 climatological period starts and ends after the training period with a gap -
  Integer, Public :: iclim     ! - setable climatological period flag -
  Integer, Public :: ieclim    ! - climatological period can extend beyond training period flag -
  Integer, Public :: ieclim_bk ! - climatological period can extend beyond training period backup flag -
  Integer, Public :: itrain1   ! - index of first date of training period (relative to combined climate and training dataset) -
  Integer, Public :: itrain2   ! - index of last date of training period (relative to combined climate and training dataset) -
  Integer, Public :: nc        ! - length of climatological period -
  Integer, Public :: ncy       ! - length of climatological period in years -
  Integer, Public :: ntc       ! - length of climatological period with length of training period -
!
  Integer, Private :: iclim1 ! - index of first date of climatological period -
  Integer, Private :: iclim2 ! - index of last date of climatological period -
  Integer, Private :: lcsn   ! - length of climatological season -
!
! Derived types
!
! Derived type definitions
! - climatological period -
  Public clim_per
  Type clim_per
     Sequence
     Type(pdate) :: d1        ! - first date of climatology -
     Type(pdate) :: d2        ! - last date of climatology -
!
     Character(lprd) :: clim1 ! - first date of climatology -
     Character(lprd) :: clim2 ! - last date of climatology -
  End Type clim_per
!
! - version 15 climatological period -
  Type clim_per_v15
     Sequence
     Type(pdate) :: d1        ! - first date of climatology -
     Type(pdate) :: d2        ! - last date of climatology -
!
     Integer :: ic1           ! - index of first date of climatological period -
     Integer :: ic2           ! - index of last date of climatological period -
     Integer :: nc            ! - length of climatological period -
     Integer :: il1           ! - index of first date of climatological period relative to independent training period -
     Integer :: il2           ! - index of last date of climatological period relative to independent training period -
     Integer :: lsn           ! - length of season -
!
     Character(lprd_v15) :: clim1 ! - first date of climatology -
     Character(lprd_v15) :: clim2 ! - last date of climatology -
  End Type clim_per_v15
!
! Interface operators
  Private :: Operator(==)
!
  Public :: Operator(/=)
!
  Interface Operator(==)
     Module Procedure same_clim_per
  End Interface
  Interface Operator(/=)
     Module Procedure diff_clim_per
  End Interface
!
! Derived type scalars
  Type(clim_per), Public :: climate_per ! - climatological period -
  Type(clim_per), Public :: climate_bak ! - backup climatological period -
!
Contains
!
!
 Function same_clim_per(cp1, cp2) &
          Result(same)
!
! Are climate periods identical?
!
! Modules
  Use time, Only: Operator(==)
!
! Function type
  Logical same
!
! Arguments
!
! Input scalars
  Type(clim_per), Intent(In) :: cp1 ! - first climatological period -
  Type(clim_per), Intent(In) :: cp2 ! - second climatological period -
!
! Executable Statements
!
! Test whether climate periods are identical
  same = (cp1%d1 == cp2%d1 .and. cp1%d2 == cp2%d2)
!
  Return
 End Function same_clim_per
!
!
!
 Function diff_clim_per(cp1, cp2) &
          Result(differ)
!
! Are climate periods different?
!
! Modules
  Use time, Only: Operator(==)
!
! Function type
  Logical differ
!
! Arguments
!
! Input scalars
  Type(clim_per), Intent(In) :: cp1 ! - first climatological period -
  Type(clim_per), Intent(In) :: cp2 ! - second climatological period -
!
! Executable Statements
!
! Test whether climate periods are different
  differ = .not.(cp1%d1 == cp2%d1 .and. cp1%d2 == cp2%d2)
!
  Return
 End Function diff_clim_per
!
!
!
 Subroutine reset_climate ()
!
! Resets climatological period
!
! Executable Statements
!
! Switch off climatological period
  iclim = 0
  iclim1 = 0
  iclim2 = 0
  itrain1 = 0
  itrain2 = 0
  lcsn = 0
  nc = 0
  ncy = 0
  ntc = 0
  Call init_climatep ()
!
  Return
 End Subroutine reset_climate
!
!
!
 Subroutine init_climatep ()
!
! Initialises climatological period
!
! Modules
  Use time, Only: Assignment(=)
!
! Executable Statements
!
! Initialise climatological period
  climate_per%d1 = 0      ! - first date of climatology -
  climate_per%d2 = 0      ! - last date of climatology -
  climate_per%clim1 = ' ' ! - first period of climatology -
  climate_per%clim2 = ' ' ! - last period of climatology -
!
  Return
 End Subroutine init_climatep
!
!
!
 Subroutine init_climate (i0, nsn, ifail)
!
! Initialises memory for climatologies
!
! Modules
  Use arrays,   Only: rwk
  Use iofiles,  Only: mya
  Use settings, Only: igauss
  Use distribs, Only: alpha, beta, pcen
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i0  ! - lower bound -
  Integer, Intent(In) :: nsn ! - number of seasons -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Allocated
!
! Executable Statements
!
! Free existing memory
  If (Allocated(rwk)) Deallocate (rwk)
!
! Allocate memory
! - climatological data -
  Allocate (clim(mya,nc,i0:nsn), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - number of used climatological cases -
  Allocate (ncu(i0:nsn), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - climatological averages -
  Allocate (ave(mya,i0:nsn), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - climatological standard deviations -
  Allocate (sdev(mya,i0:nsn), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - gamma parameters -
  If (igauss == 2) Then
     Allocate (alpha(mya,-1:nsn), Stat=ifail)
     If (ifail /= 0) GoTo 1
     Allocate (beta(mya,-1:nsn), Stat=ifail)
     If (ifail /= 0) GoTo 1
     Allocate (pcen(mya,-1:nsn), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
!
  Return
!
! Errors
1 ifail = 1
!
  Return
 End Subroutine init_climate
!
!
!
 Function climatology(icalc, isn0, lsn0, nsn, lprst)
!
! Redefines climatological period
!
! Modules
  Use data_numbers,   Only: one
  Use progress_meter, Only: end_progress_meter, set_progress_increment, start_progress_meter
  Use arrays,         Only: y
  Use iofiles,        Only: mya
  Use settings,       Only: igauss, nu, &
                            record_change
  Use categories,     Only: ithr, tobs, tobst, &
                            set_percentiles, set_thresholds
!
! Function type
  Integer :: climatology
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icalc ! - calculation completed flag -
  Integer, Intent(In) :: isn0  ! - first season -
  Integer, Intent(In) :: lsn0  ! - length of persistence seasons -
  Integer, Intent(In) :: nsn   ! - number of seasons -
!
  Logical, Intent(In) :: lprst ! - include persistence component? -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
!
! Executable Statements
!
! Recalculate climatology
  Call set_climate (icalc, lsn0)
  If (iclim1 < 1 .or. iclim2 > nu) Then
     climatology = -1
     Return
  End If
  Select Case (ithr)
   Case (1, 3)
     Select Case (igauss)
      Case (0)
        Call set_progress_increment (one/Real(3, Kind=rp), .false.)
      Case Default
        Call set_progress_increment (one/Real(4, Kind=rp), .false.)
     End Select
   Case (2)
     Select Case (igauss)
      Case (0)
        Call set_progress_increment (one/Real(2, Kind=rp), .false.)
      Case Default
        Call set_progress_increment (one/Real(3, Kind=rp), .false.)
     End Select
  End Select
  Call start_progress_meter ('Calculating climatologies', .false.)
  Call calc_climates (.true., isn0, nsn, mya, y(1:mya,:,:), ave(:,:), sdev(:,:), ncu(:), clim(1:mya,:,:), ifail)
  If (lprst) ncu(0) = ncu(1)
!
! Calculate thresholds
  If (ifail == 0) Then
     Select Case (ithr)
      Case (1)
        Call set_thresholds (igauss, mya, nsn, ncu(:), clim(:,:,1:), tobs, tobst, ifail)
        climatology = ifail
      Case (2)
        climatology = set_percentiles(nsn, ncu(:), clim(:,:,1))
        Call set_thresholds (igauss, mya, nsn, ncu(:), clim(:,:,1:), tobs, tobst, ifail)
        climatology = ifail
      Case (3)
        climatology = set_percentiles(nsn, ncu(:), clim(:,:,1))
     End Select
     If (climatology /= 0) Then
        climate_per = climate_bak
        Call end_progress_meter (.false., .false.)
        climatology = 1
     End If
     climatology = record_change() ! - record change for project file -
     Call end_progress_meter (.false., .true.)
     climatology = 0
  Else
     climate_per = climate_bak
     Call end_progress_meter (.false., .false.)
     climatology = 1
  End If
!
  Return
 End Function climatology
!
!
!
 Function valid_climate()
!
! Checks that end of climatological period is later than beginning
!
! Modules
  Use time,     Only: Operator(+), &
                      date_diff, nyears
  Use screen,   Only: window_update
  Use errors,   Only: cpt_error
  Use iofiles,  Only: yfile
  Use settings, Only: nt
!
! Function type
  Integer :: valid_climate
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Check for valid settings
  nc = date_diff(climate_per%d1, climate_per%d2, yfile%iseq) + 1
  If (nc <= 0) Then
     ncy = nyears(nc, .true.)
     climate_per%d2 = yfile%prd1%sdate + (yfile%it1 + nt - 2)
     Call window_update (climate_per%d2%iyr)
     ifail = 1
     Call cpt_error ('valid_climate', ifail, .false., &
          i_arg1=yfile%iseq)
     valid_climate = 2
  Else
     valid_climate = -1
  End If
  Call window_update (ncy)
!
  Return
 End Function valid_climate
!
!
!
 Subroutine check_climates (cfile, mfile, lsn0, lpfv)
!
! Checks that climate period is within data range
!
! Modules
  Use data_time_constants, Only: isq_sn, isq_so
  Use time,                Only: Operator(+), Operator(-), Operator(>), Operator(<), Operator(>=), Operator(<=), Operator(==), &
                                 iseq, nss, &
                                 add_to_month, adjust_pdate, date_diff, nyears
  Use screen,              Only: window_update
  Use iofiles,             Only: ifile
  Use settings,            Only: nt
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lsn0 ! - length of persistence season -
!
  Logical, Intent(In) :: lpfv ! - PFV? -
!
  Type(ifile), Intent(In) :: cfile ! - file containing climatological data -
  Type(ifile), Intent(In) :: mfile ! - file containing monthly climatological data -
!
! Locals
!
! Local scalars
  Integer :: idc ! - difference between first and last seasons of climatology -
!
  Type(pdate) :: dbdate ! - default beginning date -
  Type(pdate) :: dldate ! - default last date -
  Type(pdate) :: d1     ! - first date -
  Type(pdate) :: d2     ! - second date -
!
! Executable Statements
!
! Set default climatological date limits
  Select Case (ieclim)
! - set default climatological period to the training period -
   Case (0)
     dbdate = cfile%bdate
     dldate = cfile%bdate + (nt - 1)
   Case (1)
! - set default climatological period immediately before training period, if possible -
     dbdate = cfile%prd1%sdate
     If (lpfv) Then
        dldate = cfile%bdate - 1
        If (dldate <= dbdate) dldate = cfile%prdn%sdate ! - set to the end of the full period if end of period is before start -
! - set default climatological period to the full period -
     Else
        dldate = cfile%prdn%sdate
     End If
  End Select
!
! Set to default if climatological period has not been previously set
  If (nc <= 1) Then
     climate_per%d1 = dbdate
     climate_per%d2 = dldate
!
! Apply season settings to climatology
  Else
     climate_per%d1%imn = cfile%prd1%sdate%imn
     If (lsn0 > 0) climate_per%d1%imn = add_to_month(climate_per%d1%imn, -lsn0, &
                                        iyr=climate_per%d1%iyr) ! - adjust climatology start date for persistence component -
     Select Case (ieclim)
      Case (0)
        If (climate_per%d1 < mfile%bdate) Call adjust_pdate (climate_per%d1, mfile%bdate, cfile%iseq, .true.)
      Case (1)
        If (climate_per%d1 < mfile%prd1%sdate) Call adjust_pdate (climate_per%d1, mfile%prd1%sdate, cfile%iseq, .true.)
     End Select
     climate_per%d2%imn = cfile%prdn%sdate%imn
     If (climate_per%d2 > mfile%prdn%edate) Call adjust_pdate (climate_per%d2, mfile%prdn%edate, cfile%iseq, .false.)
! - check that previously set limits are within range of data -
     If (climate_per%d1 < dbdate) climate_per%d1 = dbdate
     If (climate_per%d2 > dldate) climate_per%d2 = dldate
     If (climate_per%d1 > climate_per%d2) Then
        climate_per%d1 = dbdate
        climate_per%d2 = dldate
     End If
  End If
! - adjust end date to ensure equal number of each season -
  Select Case (iseq)
   Case (isq_sn, isq_so)
     idc = date_diff(climate_per%d1, climate_per%d2, iseq) + 1
     idc = Mod(idc, nss)
     If (idc /= 0) climate_per%d2 = climate_per%d2 - idc
  End Select
!
! Check whether climatological period extends outside of training period
  icex = get_icex(cfile%bdate, cfile%bdate+(nt-1))
  Select Case (icex)
   Case (0)
     ntc = nt
   Case (-1) ! - climatological period ends before start of training period -
     iclim1 = 1
   Case (1) ! - climatological period starts before training period -
     iclim1 = 1
     d1 = climate_per%d1
     d2 = cfile%bdate + (nt - 1)
   Case (3) ! - climatological period contains training period -
     iclim1 = 1
     d1 = climate_per%d1
     d2 = climate_per%d2
   Case (-2) ! - climatological period ends after training period -
     iclim1 = nt + 1
   Case (2) ! - climatological period ends after training period -
     iclim1 = date_diff(cfile%bdate, climate_per%d1, cfile%iseq) + 1
     d1 = cfile%bdate
     d2 = climate_per%d2
  End Select
!
! Determine length of climatological and combined periods
  nc =  date_diff(climate_per%d1, climate_per%d2, iseq) + 1
  Select Case (icex)
   Case (1:)
     ntc = date_diff(d1, d2, iseq) + 1
   Case (:-1)
     ntc = nt + nc
  End Select
  ncy = nyears(nc, .true.)
  Call window_update (ncy)
!
  Return
 End Subroutine check_climates
!
!
!
 Subroutine set_climate (icalc, lsn0)
!
! Sets climatological parameters given date ranges
!
! Modules
  Use data_time_constants, Only: isq_dy, isq_sn, isq_so
  Use time,                Only: Operator(+), Operator(-), &
                                 iseq, lslim, &
                                 date_diff, nyears, reset_iseq
  Use iofiles,             Only: yfile
  Use settings,            Only: nt, nu
  Use fields,              Only: yfield
  Use season,              Only: len_ssn
  Use missing,             Only: kcuse
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icalc ! - calculation completed flag -
  Integer, Intent(In) :: lsn0  ! - length of persistence seasons -
!
! Locals
!
! Local scalars
  Logical :: ls_bk ! - backup use season limits? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count
!
! Executable Statements
!
! Identify climatological data
  icex = get_icex(yfile%bdate, yfile%bdate+(nt-1))
  Call get_clim_limits (icex, iseq, &
       iclim1=iclim1, iclim2=iclim2, itrain1=itrain1, itrain2=itrain2)
! - adjust for compression if analysis is running -
  If (icalc == 1) Then
     If (iseq /= isq_sn .and. iseq /= isq_so) Then
        If (iclim1 >= 1) Then
           iclim1 = iclim1 - Count(.not.kcuse(1:iclim1))
           If (iclim1 == 0) iclim1 = 1
        End If
        iclim2 = iclim2 - Count(.not.kcuse(1:iclim2))
        nc = iclim2 + 1 - iclim1
        ncu(:) = nc
        If (itrain1 >= 1) Then
           itrain1 = itrain1 - Count(.not.kcuse(1:itrain1))
           If (itrain1 == 0) itrain1 = 1
        End If
        itrain2 = itrain1 + nu - 1
     Else
! - calculate separately for each season (data are not yet compressed) -
        ncu(:) = nyears(nc, .false.)
     End If
  End If
!
! Calculate length of climatological season
  Select Case (yfile%nse)
   Case Default
     lcsn = len_ssn(yfile%prd1%sdate, yfile%prd1%edate, yfile%iseq)
   Case (2:)
     Call reset_iseq (isq_new=isq_dy)
     ls_bk = lslim
     lslim = .false.
     lcsn = date_diff(yfield(1)%tprd%sdate, yfield(yfile%nse)%tprd%edate, isq_dy) + 1
     lslim = ls_bk
     Call reset_iseq ()
  End Select
!
! Construct climatological period labels
  Call label_climate (lsn0)
!
  Return
 End Subroutine set_climate
!
!
!
 Function get_icex(d1, d2, &
          lsn0) &
          Result (icex)
!
! Identifies whether climatological period lies outside training period
!
! Modules
  Use data_time_constants, Only: isq_mn
  Use time,                Only: Operator(+), Operator(-), Operator(>), Operator(<), &
                                 date_diff
!
! Function type
  Integer :: icex
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d1 ! - first date of training period -
  Type(pdate), Intent(In) :: d2 ! - end date training period -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: lsn0 ! - length of persistence season -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Identify whether climatological period lies outside training period
  icex = 0
  If (climate_per%d1 < d1) Then
     icex = 1
  Else If (Present(lsn0)) Then
     If (lsn0 > 0) Then ! - adjust for persistence component -
        If (date_diff(climate_per%d1, d1, isq_mn) > lsn0) icex = 1
     End If
  End If
  If (climate_per%d2 > d2) icex = icex + 2
!
! Identify whether there is a gap between the climatological and training periods 
  Select Case (icex)
   Case (1) ! - climatological period is before training period -
     If (climate_per%d2 < (d1 - 1)) icex = -icex
   Case (2) ! - climatological period is after training period -
     If (climate_per%d1 > (d2 + 1)) icex = -icex
  End Select
!
  Return
 End Function get_icex
!
!
!
 Subroutine get_clim_limits (icex, isq, &
            iclim1, iclim2, itrain1, itrain2)
!
! Modules
  Use data_time_constants, Only: isq_sn, isq_so, isq_yr
  Use time,                Only: date_diff, nyears
  Use iofiles,             Only: yfile
  Use settings,            Only: nt
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icex ! - climatological period extends beyond training period flag -
  Integer, Intent(In) :: isq  ! - time-sequencing -
!
! Output scalars
! - optional output scalars -
  Integer, Intent(Out), Optional :: iclim1  ! - index of first date of climatological period -
  Integer, Intent(Out), Optional :: iclim2  ! - index of last date of climatological period -
  Integer, Intent(Out), Optional :: itrain1 ! - index of first date of training period (relative to combined climate and training dataset) -
  Integer, Intent(Out), Optional :: itrain2 ! - index of last date of training period (relative to combined climate and training dataset) -
!
! Local scalars
  Integer :: i1 ! - index of first date -
!
! Executable Statements
!
! Identify limits of the climatological and training periods
  Select Case (isq)
   Case Default
     Select Case (icex)
      Case (0)
        If (Present(iclim1)) iclim1 = date_diff(yfile%bdate, climate_per%d1, yfile%iseq) + 1
        If (Present(iclim2)) iclim2 = date_diff(yfile%bdate, climate_per%d2, yfile%iseq) + 1
        If (Present(itrain1)) itrain1 = 1
        If (Present(itrain2)) itrain2 = nt
      Case (1) ! - climatological period is before training period (no gap) -
        If (Present(iclim1)) iclim1 = 1
        If (Present(iclim2)) iclim2 = nc
        If (Present(itrain1) .or. Present(itrain2)) i1 = date_diff(climate_per%d1, yfile%bdate, yfile%iseq) + 1
        If (Present(itrain1)) itrain1 = i1
        If (Present(itrain2)) itrain2 = i1 + nt - 1
      Case (-1) ! - climatological period is before training period (gap) -
        If (Present(iclim1)) iclim1 = 1
        If (Present(iclim2)) iclim2 = nc
        If (Present(itrain1)) itrain1 = nc + 1
        If (Present(itrain2)) itrain2 = ntc
      Case (2) ! - climatological period is after training period (no gap) -
        If (Present(iclim1) .or. Present(iclim2)) i1 = date_diff(yfile%bdate, climate_per%d1, yfile%iseq) + 1
        If (Present(iclim1)) iclim1 = i1
        If (Present(iclim2)) iclim2 = i1 + nc - 1
        If (Present(itrain1)) itrain1 = 1
        If (Present(itrain2)) itrain2 = nt
      Case (-2) ! - climatological period is after training period (gap) -
        If (Present(iclim1)) iclim1 = nt + 1
        If (Present(iclim2)) iclim2 = ntc
        If (Present(itrain1)) itrain1 = 1
        If (Present(itrain2)) itrain2 = nt
      Case (3) ! - training period is within climatological period -
        If (Present(iclim1)) iclim1 = 1
        If (Present(iclim2)) iclim2 = ntc
        If (Present(itrain1)) itrain1 = date_diff(climate_per%d1, yfile%bdate, yfile%iseq) + 1
        If (Present(itrain2)) itrain2 = itrain1 + nt - 1
     End Select
! - correct for seasonally sequenced data -
   Case (isq_sn, isq_so)
     Select Case (icex)
      Case (0)
        If (Present(iclim1)) iclim1 = date_diff(yfile%bdate, climate_per%d1, isq_yr) + 1
        If (Present(iclim2)) iclim2 = date_diff(yfile%bdate, climate_per%d2, isq_yr) + 1
        If (Present(itrain1)) itrain1 = 1
        If (Present(itrain2)) itrain2 = nyears(nt, .false.)
      Case (1) ! - climatological period is before training period (no gap) -
        If (Present(iclim1)) iclim1 = 1
        If (Present(iclim2)) iclim2 = iclim1 + nyears(nc, .false.) - 1
        If (Present(itrain1)) itrain1 = date_diff(climate_per%d1, yfile%bdate, isq_yr) + 1
        If (Present(itrain2)) itrain2 = itrain1 + nyears(nt, .false.) - 1
      Case (-1) ! - climatological period is before training period (gap) -
        If (Present(iclim1)) iclim1 = 1
        If (Present(iclim2)) iclim2 = nyears(nc, .false.)
        If (Present(itrain1)) itrain1 = nyears(nc, .false.) + 1
        If (Present(itrain2)) itrain2 = nyears(ntc, .false.)
      Case (2) ! - climatological period is after training period (no gap) -
        If (Present(iclim1) .or. Present(iclim2)) i1 = date_diff(yfile%bdate, climate_per%d1, isq_yr) + 1
        If (Present(iclim1)) iclim1 = i1
        If (Present(iclim2)) iclim2 = i1 + nyears(nc, .false.) - 1
        If (Present(itrain1)) itrain1 = 1
        If (Present(itrain2)) itrain2 = nyears(nt, .false.)
      Case (-2) ! - climatological period is after training period (gap) -
        If (Present(iclim1)) iclim1 = nyears(nt, .false.) + 1
        If (Present(iclim2)) iclim2 = nyears(ntc, .false.)
        If (Present(itrain1)) itrain1 = 1
        If (Present(itrain2)) itrain2 = nyears(nt, .false.)
      Case (3) ! - training period is within climatological period -
        If (Present(iclim1)) iclim1 = 1
        If (Present(iclim2)) iclim2 = nyears(ntc, .false.)
        If (Present(itrain1) .or. Present(itrain2)) i1 = date_diff(climate_per%d1, yfile%bdate, isq_yr) + 1
        If (Present(itrain1)) itrain1 = i1
        If (Present(itrain2)) itrain2 = i1 + nyears(nt, .false.) - 1
     End Select
  End Select
!
  Return
 End Subroutine get_clim_limits
!
!
!
 Subroutine label_climate (lsn0)
!
! Constructs climatological period labels
!
! Modules
  Use data_time_constants, Only: isq_dy, isq_mn, isq_sn, nmn
  Use time,                Only: pdate, &
                                 iseq, lslim, nse, nss, &
                                 date_diff, get_cdate, reset_iseq, same_ssn
  Use iofiles,             Only: yfile
  Use fields,              Only: yfield
  Use season,              Only: which_season
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lsn0 ! - length of persistence seasons -
!
! Locals
!
! Local scalars
  Integer :: isn ! - season index -
  Integer :: l   ! - season counter -
  Integer :: ll  ! - season counter -
  Integer :: ls1 ! - length of season -
  Integer :: ls2 ! - length of season -
!
  Logical :: ls_bk ! - backup use season limits? -
  Logical :: lseq  ! - reset time sequencing -
!
  Type(pdate) :: d2 ! - end date of climatology -
!
! Executable Statements
!
! Determine length of climatological season
  lseq = (iseq > 0)
  If (lseq) Then
     Call reset_iseq (isq_new=isq_dy)
     If (nse > 0) Then
        ls_bk = lslim
        lslim = .false.
        ls1 = date_diff(yfile%prd1%sdate, yfield(yfile%nse)%tprd%edate, isq_dy) + 1
        ls1 = lcsn + lsn0
     Else
        ls1 = iseq
     End If
  Else If (lcsn <= nmn) Then
     ls1 = lcsn + lsn0
  Else
     ls1 = lcsn + lsn0 + 1 - nmn
  End If
  Select Case (iseq)
   Case Default
     ls2 = ls1
     d2 = climate_per%d2
   Case (isq_sn)
     isn = which_season(yfile%fdate, nss, yfield(:)%tprd%sdate) ! - identify first season used -
     isn = isn - 1                                              ! - extend climatological year to the previous period -
     If (isn == 0) isn = nss
     lcsn = date_diff(yfield(isn)%tprd%sdate, yfield(isn)%tprd%edate, isq_mn) + 1
     If (lcsn <= nmn) Then ! - there is no persistence component with seasonal data -
        ls2 = lcsn
     Else
        ls2 = lcsn + 1 - nmn
     End If
     d2 = climate_per%d2
   Case (1:)
     If (nse > 1) Then
        l = 0 
        Do ll = 1, nse
           l = l + 1
           If (same_ssn(climate_per%d2, yfield(ll)%tprd%sdate)) Exit
        End Do
        If (l == 0) l = nse
        ls2 = date_diff(yfield(1)%tprd%sdate, yfield(l)%tprd%sdate, isq_dy) + &
              date_diff(yfield(1)%tprd%sdate, yfield(1)%tprd%edate, isq_dy) + 1
        d2 = yfield(1)%tprd%sdate
        d2%iyr = climate_per%d2%iyr - (yfield(l)%tprd%sdate%iyr - yfield(1)%tprd%sdate%iyr)
     Else
        d2 = climate_per%d2
        ls2 = iseq
     End If
  End Select
!
! Construct climatological period labels
  climate_per%clim1 = get_cdate(climate_per%d1, ls1, 2)
  climate_per%clim2 = get_cdate(d2, ls2, 2)
  If (lseq) Then
     Call reset_iseq ()
     If (nse > 0) lslim = ls_bk
  End If
!
  Return
 End Subroutine label_climate
!
!
!
 Subroutine calc_climates (lstd, isn0, nsn, m, y, ybar, ystd, ncu, clim, ifail)
!
! Calculates mean and variance for a climate period
!
! Modules
  Use data_numbers,   Only: zero, one
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: yhat, yret, yrpls, &
                            insertion_sort
  Use settings,       Only: igauss, iretro, istd, nu, nur
  Use distribs,       Only: alpha, beta, pcen, &
                            stdize, thom_gamma, ustdize
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isn0 ! - first season -
  Integer, Intent(In) :: nsn  ! - number of seasons -
  Integer, Intent(In) :: m    ! - number of variables -
!
  Logical, Intent(In) :: lstd ! - recalculate climatology? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(isn0:), Intent(In) :: ncu ! - number of used climatological cases -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:,isn0:), Intent(InOut) :: y ! - data -
!
! Output arrays
  Real(Kind=rp), Dimension(:,isn0:), Intent(Out) :: ybar ! - climatological mean -
  Real(Kind=rp), Dimension(:,isn0:), Intent(Out) :: ystd ! - climatological standard deviation -
!
  Real(Kind=rp), Dimension(:,:,isn0:), Intent(Out) :: clim ! - climatological (sorted) data -
!
! Locals
!
! Local scalars
  Integer :: i   ! - variable index -
  Integer :: isn ! - season index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Sqrt, Sum
!
! Executable Statements
!
! Unstandardise, if requested
  ifail = -1
  Do isn = isn0, nsn
     If (lstd .and. istd /= 0) Then
        Call ustdize (istd, igauss, isn, m, nu, ybar(:,isn), ystd(:,isn), y(:,:,isn))
        If (isn == nsn) Then
           Call ustdize (istd, igauss, isn, m,  nu, ybar(:,isn), ystd(:,isn), yhat(:,:))
           If (iretro == 1) Then
              Call ustdize (istd, igauss, isn, m, nur, ybar(:,isn), ystd(:,isn), yret(:,:))
              Call ustdize (istd, igauss, isn, m, nur, ybar(:,isn), ystd(:,isn), yrpls(:,:)%lower)
              Call ustdize (istd, igauss, isn, m, nur, ybar(:,isn), ystd(:,isn), yrpls(:,:)%upper)
           End If
        End If
     End If
!
! Extract climatological data
     Do i = 1, m
        clim(i,1:ncu(isn),isn) = y(i,iclim1:iclim2,isn)
!
! Calculate climatological moments
! - means -
        ybar(i,isn) = Sum(clim(i,1:ncu(isn),isn))/Real(ncu(isn), Kind=rp)
! - standard deviations -
        If (nc > 1) Then
           ystd(i,isn) = Sum((clim(i,1:ncu(isn),isn) - ybar(i,isn))**2)/Real(ncu(isn)-1, Kind=rp)
           If (ystd(i,isn) > zero) Then
              ystd(i,isn) = Sqrt(ystd(i,isn))
           Else
              ystd(i,isn) = one
           End If
        Else
           ystd(i,isn) = one
        End If
!
! Create climatology
! - sort data -
        Call insertion_sort (ncu(isn), clim(i,:,isn), 'a')
! - calculate gamma parameters -
        If (igauss == 2) Call thom_gamma (ncu(isn), clim(i,:,isn), alpha(i,isn), beta(i,isn), pcen(i,isn))
!
! Update progress meter
        If (update_progress_meter(.false.) /= 0) Return
     End Do
!
! Restandardise, if requested
     If (lstd .and. istd /= 0) Then
        Call stdize (istd, igauss, isn, m,       nu, ybar(:,isn), ystd(:,isn), y(:,:,isn))
        Call stdize (istd, igauss, isn, m, ncu(isn), ybar(:,isn), ystd(:,isn), clim(:,:,isn))
        If (isn == nsn) Then
           Call stdize (istd, igauss, isn, m,  nu, ybar(:,isn), ystd(:,isn), yhat(:,:))
           If (iretro == 1) Then
              Call stdize (istd, igauss, isn, m, nur, ybar(:,isn), ystd(:,isn), yret(:,:))
              Call stdize (istd, igauss, isn, m, nur, ybar(:,isn), ystd(:,isn), yrpls(:,:)%lower)
              Call stdize (istd, igauss, isn, m, nur, ybar(:,isn), ystd(:,isn), yrpls(:,:)%upper)
           End If
        End If
     End If
!
! Update progress meter
     If (update_progress_meter(.false.) /= 0) Return
  End Do
  ifail = 0
!
  Return
 End Subroutine calc_climates
!
!
!
 Subroutine print_climatology (lprmpt, &
            iout)
!
! Prints climatology
!
! Modules
  Use gui, Only: print_text
!
! Arguments
!
! Input scalars
  Logical, Intent(In), Optional :: lprmpt ! - print prompt? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
!
! Functions nad Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Print climatology
  If (Present(iout)) Then
     If (lprmpt) Then
        Call print_text ('Climatology:', &
             iout=iout)
        Call print_text ('    '//Trim(climate_per%clim1)//' to '//Trim(climate_per%clim2), &
             iout=iout)
     Else
        Call print_text (Trim(climate_per%clim1)//' to '//Trim(climate_per%clim2), &
             iout=iout)
     End If
  Else
     If (lprmpt) Then
        Call print_text ('Climatology:', &
             lcolour=.true.)
        Call print_text ('    '//Trim(climate_per%clim1)//' to '//Trim(climate_per%clim2))
     Else
        Call print_text (Trim(climate_per%clim1)//' to '//Trim(climate_per%clim2))
     End If
  End If
!
  Return
 End Subroutine print_climatology
!
!
!
 Subroutine proj_read_climate (iin, pver, lprst, nt, nsn, lsn0, lres, ifail)
!
! Reads climatology settings from project file
!
! Modules
  Use data_numbers, Only: sp, tolv
  Use time,         Only: Operator(+)
  Use iofiles,      Only: yfile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin  ! - input file unit number -
  Integer, Intent(In) :: lsn0 ! - length of persistence season -
  Integer, Intent(In) :: nsn  ! - number of seasons -
  Integer, Intent(In) :: nt   ! - number of cases in training period -
!
  Real(Kind=sp), Intent(In) :: pver ! - project file version -
!
  Logical, Intent(In) :: lprst ! - persistence component? -
  Logical, Intent(In) :: lres  ! - read results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Type(clim_per_v15) :: climate_per_v15 ! - version 15 climate -
!
! Executable Statements
!
! Read climatology settings
  If (.not.lres) Then
     If (pver > 16.0_sp) Then
        Read (Unit=iin, Err=1, End=2) climate_per
        Read (Unit=iin, Err=1, End=2) ieclim
     Else
        Read (Unit=iin, Err=1, End=2) climate_per_v15, ieclim
        Call clim_per_v15_to_v16 (climate_per_v15, climate_per)
     End If
!
! Check whether climatological period extends beyond training period
     If (yfile%cfile%lset) Then
        icex = get_icex(yfile%fdate, yfile%fdate+(nt-1), &
               lsn0=lsn0)
     Else
        icex = 0
     End If
!
! Read results
  Else
     If (pver > 16.02_sp-tolv) Then
        Read (Unit=iin, Err=1, End=2) nc
        If (.not.lprst) Then
           Call init_climate (1, nsn, ifail)
        Else
           Call init_climate (0, nsn, ifail)
        End If
        If (ifail /= 0) GoTo 3
        Read (Unit=iin, Err=1, End=2) ave
        Read (Unit=iin, Err=1, End=2) clim
        Read (Unit=iin, Err=1, End=2) sdev
        Read (Unit=iin, Err=1, End=2) ncu
     Else
        ifail = 1
        Return
     End If
  End If
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
! - memory allocation error -
3 ifail = 3
  Return
!
 End Subroutine proj_read_climate
!
!
!
 Subroutine proj_write_climate (iout, lres, ifail)
!
! Writes climatology settings to project file
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output file unit number -
!
  Logical, Intent(In) :: lres ! - write results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write climatology settings
  If (.not.lres) Then
     Write (Unit=iout, Err=1) climate_per
     Write (Unit=iout, Err=1) ieclim
!
! Write results
  Else
     Write (Unit=iout, Err=1) nc
     Write (Unit=iout, Err=1) ave
     Write (Unit=iout, Err=1) clim
     Write (Unit=iout, Err=1) sdev
     Write (Unit=iout, Err=1) ncu
  End If
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
!
 End Subroutine proj_write_climate
!
!
!
 Subroutine clim_per_v15_to_v16 (climate_v15, climate_v16)
!
! Converts version 15 climate to version 16
!
! Arguments
!
! Input scalars
  Type(clim_per_v15), Intent(In) :: climate_v15 ! - version 15 climate -
!
! Output scalars
  Type(clim_per), Intent(Out) :: climate_v16 ! - version 16 climate -
!
! Executable Statements
!
! Convert climatological period
  climate_v16%d1 = climate_v15%d1       ! - first date of climatology -
  climate_v16%d2 = climate_v15%d2       ! - last date of climatology -
  climate_v16%clim1 = climate_v15%clim1 ! - first date of climatology -
  climate_v16%clim2 = climate_v15%clim2 ! - last date of climatology -
!
  Return
 End Subroutine clim_per_v15_to_v16
!
!
!
 Function close_climate()
!
! Frees memory allocated for climate
!
! Modules
  Use distribs, Only: alpha, beta, pcen
!
! Function type
  Integer :: close_climate
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated, Associated
!
! Executable Statements
!
! Free memory allocated for categories
  If (Allocated(pcen))  Deallocate (pcen)
  If (Allocated(beta))  Deallocate (beta)
  If (Allocated(alpha)) Deallocate (alpha)
  If (Allocated(sdev))  Deallocate (sdev)
  If (Allocated(ave))   Deallocate (ave)
  If (Allocated(ncu))   Deallocate (ncu)
  If (Associated(clim)) Then
     Deallocate (clim)
     NullIfY (clim)
  End If
  close_climate = 0
!
  Return
 End Function close_climate
End Module climate
