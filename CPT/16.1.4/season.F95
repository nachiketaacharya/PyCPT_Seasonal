! Author: Simon Mason
Module season
!
! Modules
  Use data_numbers,      Only: rp
  Use data_io_constants, Only: lprd
  Use time,              Only: pprd
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: available_ssn, close_season, get_cssel, get_cssn, get_seasons, len_ssn, proj_read_season, proj_write_season, &
            prompt_season, reset_season, seasonal, seasonal_zyp, set_icsv, set_seasons, set_ssnlimits, which_season
!
! Parameters
!
! Local parameters
  Integer, Parameter, Public :: lssel = 40 ! - length of current season selection -
!
! Arrays
!
! Integer arrays
  Integer, Dimension(:), Allocatable, Public :: lsn ! - length of season -
!
! Real arrays
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Private :: vs ! - seasonal aggregates -
!
! Logical arrays
  Logical, Dimension(:), Allocatable, Public :: lssr ! - seasons have shortened records? -
!
  Logical, Dimension(:), Allocatable, Private :: ldsn ! - distinct season? -
!
! Scalars
!
! Integer scalars
  Integer, Public :: icsv  ! - calculate seasonal values flag -
  Integer, Public :: ispi  ! - SPI period -
  Integer, Public :: lensn ! - total length of season -
  Integer, Public :: lb    ! - lower-bound of lsn -
  Integer, Public :: nsn   ! - number of seasons -
!
  Integer, Private :: l0  ! - first month of persistence component -
  Integer, Private :: nfd ! - total number of fields -
!
! Character scalars
  Character(Len=lssel), Public :: cssel
!
! Logical scalars
  Logical, Public :: lprst ! - include persistence component? -
!
! Derived-type scalars
  Type(pprd), Pointer, Private :: pperiod => Null()  ! - pointer period -
!
Contains
!
!
 Function len_ssn(sdate, edate, isq_in)
!
! Calculates season length
!
! Modules
  Use data_time_constants, Only: isq_mn
  Use time,                Only: pdate, &
                                 date_diff
!
! Function type
  Integer :: len_ssn
!
! Arguments
!
! Input scalars
  Integer :: isq_in ! - input sequence -
!
  Type(pdate), Intent(In) :: sdate ! - start date -
  Type(pdate), Intent(In) :: edate ! - end date -
!
! Locals
!
! Local scalars
  Integer :: isq
!
! Executable Statements
!
! Return appropriate sequencing
  Select Case (isq_in)
   Case (:0)
     isq = isq_mn
   Case Default
     isq = isq_in
  End Select
! 
! Calculate season length
  len_ssn = date_diff(sdate, edate, isq) + 1
!
  Return
 End Function len_ssn
!
!
!
 Subroutine reset_season (lresetlims)
!
! Resets season definitions
!
! Modules
  Use time, Only: Assignment(=), &
                  ssnlims
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lresetlims
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Reset limits
  If (lresetlims) ssnlims = 0
!
! Reset season settings
  icsv = 0
  ispi = 0
  lb = 1
  lensn = 0
  nsn = 1
  lprst = .false.
!
! Reset memory
  If (Allocated(lsn)) Deallocate (lsn)
  Allocate (lsn(0:nsn))
  lsn(0) = 0
  lsn(1) = 1
!
  Return
 End Subroutine reset_season
!
!
!
 Function which_season(ssn, nsn, ssns) &
          Result (issn)
!
! Identifies which season the current one is in a sequence
!
! Modules
  Use time, Only: pdate, &
                  same_ssn
!
! Function type
  Integer :: issn 
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nsn ! - number of seasons -
!
  Type(pdate), Intent(In) :: ssn ! - current season -
!
! Input arrays
  Type(pdate), Dimension(:), Intent(In) :: ssns ! - seasons -
!
! Locals
!
! Local scalars
  Integer :: isn ! - season index -
!
! Exeutable Statements
!
! Compare seasons
  issn = 0
  Do isn = 1, nsn
     If (same_ssn(ssn, ssns(isn))) Then
        issn = isn
        Exit
     End If
  End Do
!
  Return
 End Function which_season
!
!
!
 Subroutine get_seasons (ianal, xfield, yfield, ifail)
!
! Identifies (PFV and model predictors) or prompts (otherwise) for target season
!
! On exit:
!    ifail =  0 Successful
!    ifail = -1 Cancelled
!    ifail =  1 Problem allocating memory
!
! Modules
  Use data_cpt_constants,  Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr, nmn
  Use time,                Only: Operator(==), Operator(+), &
                                 iafter, iseq, &
                                 reset_iseq, same_ssn
  Use errors,              Only: cpt_error
  Use iofiles,             Only: mfile, xfile, yfile
  Use settings,            Only: isynch
  Use fields,              Only: field, &
                                 mfield, nmfd, &
                                 init_field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis indicator -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Type(field), Dimension(:), Pointer :: xfield ! - X field settings -
  Type(field), Dimension(:), Pointer :: yfield ! - Y field settings -
!
! Locals
!
! Local scalars
  Integer :: i      ! - field index -
  Integer :: isn    ! - season/field index -
  Integer :: lsn_bk ! - back-up length of season -
  Integer :: l      ! - month index -
  Integer :: il     ! - field/month index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Max, Size, Ubound
!
! Exeutable Statements
!
! Compare seasons
  If (ianal/=ia_pfv .and. isynch == 0 .and. xfile%cfile%lset) Then
     If (xfield(1)%mdate == 0) Then
        iafter = 1
     Else If (.not.same_ssn(xfield(1)%tprd, yfield(1)%tprd)) Then ! - force a lag if Y season is incompatible with X -
        If (xfile%iseq == mfile%iseq) Then ! - seasons are necessarily incompatible if the sequencing is identical -
           iafter = 1
        Else If (mfile%iseq == isq_mn .and. xfile%iseq < 0) Then ! - no incompatibility if Y-file is monthly -
           iafter = 0
        Else
           iafter = 0
        End If
     Else
        iafter = 0
     End If
  Else
     iafter = 0
  End If
!
! Count number of different seasons
  If (xfile%cfile%lset) Then
     Call count_seasons (xfield)
  Else
     Call count_seasons (yfield)
  End If
  Select Case (ianal)
    Case (ia_cca, ia_pcr, ia_mlr)
     If (xfield(1)%mdate == 0) nsn = 1
  End Select
!
! Allocate workspace
  If (Ubound(lsn, Dim=1) /= nsn) Then
     lsn_bk = lsn(1)
     Deallocate (lsn)
     Allocate (lsn(0:nsn))
     lsn(1) = lsn_bk
  End If
!
! Set X seasons if Y file is unopened
  If (.not.yfile%cfile%lset) Then
     Select Case (iseq)
      Case (isq_sn)
        isn = 0
        Do l = 1, nfd
           If (.not.ldsn(l)) Cycle
           isn = isn + 1
           lsn(isn) = len_ssn(xfield(l)%tprd%sdate, xfield(l)%tprd%edate, xfile%iseq)
        End Do
      Case (isq_so)
        lsn(:) = len_ssn(xfield(1)%tprd%sdate, xfield(1)%tprd%edate, xfile%iseq)
     End Select
  End If
!
! Set seasons and start-at dates
  If (icsv == 0) Then
     xfile%bdate = xfile%fdate
     lb = 1
     lensn = 0
     lsn(0) = 0
     If (xfile%cfile%lset) Then
        Select Case (mfile%iseq)
         Case (isq_sn)
           isn = 0
           Do l = 1, Size(yfield)
              isn = isn + 1
              lsn(isn) = len_ssn(yfield(l)%tprd%sdate, yfield(l)%tprd%edate, isq_mn)
           End Do
         Case (isq_so)
           lsn(:) = len_ssn(yfield(1)%tprd%sdate, yfield(1)%tprd%edate, isq_mn)
           lensn = nmn + lsn(1) - 1
        End Select
     End If
     lprst = .false. ! - switch off persistence component if seasonal calculations are not required -
     If (.not.xfile%cfile%lset) yfile%bdate = yfile%fdate
     Call label_seasons (lsn(1:))
     ifail = 0
     Return
  End If
!
! Backup Y fields as monthly fields
  Nullify (mfield)
  Select Case (mfile%iseq)
   Case (isq_yr)
     nmfd = Size(yfield)
   Case (isq_mn)
     nmfd = Size(yfield)*Max(1,lensn)
   Case (1:)
     nmfd = mfile%nse
   Case Default
     nmfd = nsn
  End Select
  Call init_field (mfield, nmfd, yfield(:)%rmiss, ifail)
  If (ifail /= 0) Return
  Select Case (iseq)
   Case (isq_yr)
     Select Case (mfile%iseq)
      Case Default
        mfield(:) = yfield(:)
      Case (isq_mn)
        Call reset_iseq (isq_mn)
        Do i = 1, Size(yfield)
           il = (i-1)*Max(lensn, 1) + 1
           mfield(il) = yfield(i)
           Do l = 2, lensn
              il = il + 1
              mfield(il) = yfield(i)
              mfield(il)%tprd = mfield(il-1)%tprd + 1
           End Do
        End Do
        Call reset_iseq ()
     End Select
   Case Default
     mfield(1) = yfield(1)
     mfield(2:) = mfield(1)
  End Select
!
! Identify season
  Select Case (ianal)
   Case (ia_cca, ia_pcr, ia_mlr, ia_gcm)
     If (xfield(1)%mdate == 0 .or. lb == 0) Then
        nsn = 1
        Call prompt_season (nsn, 1, yfield, ifail)
     Else
        Call init_field (yfield, nsn, mfield(1:1)%rmiss, ifail) ! - re-initialise Y fields -
        If (ifail == 0) Then
           yfield(:) = mfield(1)
           Call identify_seasons (nsn, xfield(:)%tprd, yfield(:)%tprd)
        End If
     End If
   Case (ia_pfv)
     Call init_field (yfield, nsn, mfield(1:1)%rmiss, ifail) ! - re-initialise Y fields -
     If (ifail == 0) Then
        yfield(:) = mfield(1)
        Call identify_seasons (nsn, xfield(:)%tprd, yfield(:)%tprd)
     End If
  End Select
!
! Set seasons
  Select Case (ifail)
   Case (0)
     Call set_seasons (nsn, yfield(:)%cssn, yfield(:)%tprd)
   Case (-1) ! - cancelled -
     Return
   Case Default ! - problem allocating memory -
     ifail = 1
     Call cpt_error ('get_seasons', ifail, .false.)
  End Select
!
  Return
!
 Contains
!
!
 Subroutine count_seasons (afield)
!
! Counts number of seasons
!
! Modules
  Use data_time_constants, Only: isq_sn, isq_so, isq_yr, nmn
  Use time,                Only: iseq, nss
  Use fields,              Only: field
!
! Arguments
!
! Pointer arrays
  Type(field), Dimension(:), Pointer :: afield ! - field settings -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated, Any, Size
!
! Exectuable Statements
!
! Count number of seasons
  nfd = Size(afield)
  If (Allocated(ldsn)) Deallocate (ldsn)
  Allocate (ldsn(nfd))
  ldsn(1) = .true.
  Select Case (iseq)
   Case (isq_yr, isq_sn)
     nsn = 1
     If (nfd > 1) Then
        Do isn = 2, nfd
           If (Any(afield(1:isn-1)%cssn == afield(isn)%cssn)) Then
              ldsn(isn) = .false.
           Else
              ldsn(isn) = .true.
              nsn = nsn + 1
           End If
        End Do
     End If
     If (iseq == isq_sn) nss = nsn
   Case (isq_so)
     nsn = nmn
     nss = 1
  End Select
!
  Return
  End Subroutine count_seasons
 End Subroutine get_seasons
!
!
!
 Subroutine prompt_season (nsn, lby, yfield, ifail)
!
! Prompts for season definitions
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_yr, nmn
  Use labels,              Only: cg_cselection, &
                                 l_tseason
  Use time,                Only: Operator(+), Operator(-), &
                                 adjust_pprd, reset_iseq
  Use gui,                 Only: iw, &
                                 box_close, box_open, init_win, print_cvalue, print_text, prompt_integer, win_prompt
#if GUI == 1
  Use gui,                 Only: print_advisory, prompt_boolean
#endif
  Use errors,              Only: cpt_error
  Use iofiles,             Only: mfile, xfile, yfile
#if GUI == 1
  Use settings,            Only: istd
#endif
  Use settings,            Only: record_change, set_options
  Use fields,              Only: field, &
                                 init_field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nsn ! - number of distinct seasons -
  Integer, Intent(In) :: lby ! - lower bound for yfield -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Type(field), Dimension(:), Pointer :: yfield ! - Y field settings -
!
! Locals
!
! Local scalars
  Integer :: iprst ! - use persistence component -
  Integer :: istdu ! - used standardization option -
!
  Character(Len=64) :: cse ! - error seasons -
!
! Local arrays
  Type(field), Dimension(:), Allocatable :: tfield ! - temporary field information -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Kind, Lbound, Trim, Ubound
!
! Executable Statements
!
! Backup season settings
  Call reset_iseq (isq_new=isq_mn)
  If (ispi == 0) ispi = 1
  ifail = set_options(i1=yfield(1)%tprd%sdate%imn, i2=lsn(1), l1=lprst)
  If (lprst) Then
     iprst = 1
  Else
     iprst = 0
  End If
!
! Prompt for season
  cssel = get_cssel(0, yfield(1)%tprd%sdate%imn, lsn(1), 1, 0)
  pperiod => yfield(1)%tprd
  Call init_win (title=Trim(l_tseason%c))
  Call box_open (.false., .false., &
       title='Set season')
  Call prompt_integer ('First month of season to forecast', yfield(1)%tprd%sdate%imn, &
       itab=60, ilow=1, ihgh=nmn, fcb1=record_change, fcb2=update_ctsel, lfocus = .false.)
!
! Prompt for season
  Call prompt_integer ('Length of season to forecast', lsn(1), &
       itab=60, ilow=1, ihgh=nmn, fcb1=record_change, fcb2=update_ctsel, lfocus = .false.)
!
! Prompt for persistence component
#if GUI == 1
  Call prompt_boolean ('Add persistence component?', iprst)
#endif
  Call print_text (' ')
  Call print_cvalue (Trim(cg_cselection), cssel, &
       iwid=lssel)
  Call print_text (' ')
  Call box_close (.false.)
  iw = win_prompt(.true.)
!
! Set option
  lprst = (iprst == 1)
  If (set_options(iw=iw, i1=yfield(1)%tprd%sdate%imn, i2=lsn(1), l1=lprst) >= 3) Then
     ifail = -1
     GoTo 2
  End If
!
! Prompt for persistence component
1 Continue
#if GUI == 1
  If (lprst) Then
     If (lsn(1) < nmn) Then
        istdu = istd
#else
        istdu = 3
#endif
        l0 = yfield(lby)%tprd%sdate%imn
        If (istdu /= 3) Then
           ispi = pperiod%sdate%imn - l0
           If (ispi < 0) ispi = ispi + nmn
           ispi = ispi + lsn(1)
        End If
        cssel = get_cssel(0, pperiod%sdate%imn, lsn(1), 1, ispi)
        ifail = set_options(i1=ispi, i2=l0)
        Call init_win (title=Trim(l_tseason%c))
        Call box_open (.false., .false., &
             title='Persistence component')
        Select Case (istdu)
         Case Default
           Call prompt_integer ('First month of persistence component', l0, &
                ilow=1, ihgh=nmn, fcb1=record_change, fcb2=update_cl0, lfocus = .false.)
         Case (3) ! - SPI standardization -
           Call prompt_integer ('Length of SPI', ispi, &
                ilow=lsn(1), ihgh=nmn, fcb1=record_change, fcb2=update_cspi, lfocus = .false.)
        End Select
        Call print_text (' ')
        Call print_cvalue (Trim(cg_cselection), cssel, &
             iwid=lssel)
        Call box_close (.false.)
        iw = win_prompt(.true.)
        If (istdu /= 3) Then
           ispi = pperiod%sdate%imn - l0
           If (ispi < 0) ispi = ispi + nmn
           ispi = ispi + lsn(1)
        End If
        If (set_options(iw=iw, i1=ispi, i2=l0) >= 3) ispi = lsn(1)
#if GUI == 1
     Else
        Call print_advisory ('It is not possible to add a persistence component when the length of the season is one year')
        ispi = lsn(1)
     End If
  Else
     ispi = lsn(1)
  End If
#endif
! - check for invalid SPI period -
  Allocate (tfield(0:Ubound(yfield, Dim=1)))
  tfield(0) = yfield(1)
  If (ispi > nmn) Then
     ifail = 1
     tfield(1)%tprd = yfield(1)%tprd
     tfield(1)%tprd%edate = tfield(1)%tprd%sdate + (lsn(1) - 1)
     tfield(0)%tprd%sdate = tfield(1)%tprd%edate - (ispi - 1)
     tfield(0)%tprd%edate = tfield(1)%tprd%sdate - 1
     lsn(0) = ispi - lsn(1)
     cse='Persistence season: '//Trim(get_cssn(tfield(0)%tprd, isq_yr, 1))//&
        &'; Target season: '//Trim(get_cssn(tfield(1)%tprd, isq_yr, 1))
     Call cpt_error ('prompt_season', ifail, .false., &
          c_arg1=Trim(cse))
     ispi = lsn(1)
     Deallocate (tfield)
     GoTo 1
! - check for persistence component -
  Else If (ispi > lsn(1)) Then
     lsn(0) = ispi - lsn(1)
     lprst = .true.
     lb = 0
  Else
     lsn(0) = 0
     lprst = .false.
     lb = 1
  End If
!
! Adjust start year if month is unavailable
  yfield(1)%tprd%edate = yfield(1)%tprd%sdate + (lsn(1) - 1)
  Call adjust_pprd (yfield(1)%tprd, mfile%prd1%sdate, xfile%iseq, .true.)
!
! Reset Y fields
  If (Lbound(yfield, Dim=1) /= lb .or. Ubound(yfield, Dim=1)/yfile%nls /= nsn) Then
     tfield(1:) = yfield(1:)
     Call init_field (yfield, Ubound(tfield, Dim=1, Kind=Kind(0)), (/tfield(1)%rmiss/), ifail, &
          lb=lb)
     If (ifail /= 0) GoTo 2
     yfield(lb:) = tfield(lb:)
  End If
  If (lb == 0) yfield(0)%tprd = yfield(1)%tprd - lsn(0)
  ifail = 0
!
! Restore time sequencing
2 Call reset_iseq ()
!
  Return
 End Subroutine prompt_season
!
!
!
 Subroutine set_ssnlimits (date1, nses, isq, ssnlims)
!
! Sets season limits (for sub-seasonal data)
!
! Modules
  Use data_time_constants, Only: isq_dy
  Use time,                Only: pdate, pprd, &
                                 Operator(+), &
                                 lslim, nse, &
                                 reset_iseq
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq  ! - time sequencing -
  Integer, Intent(In) :: nses ! - number of sub-seasons -
!
  Type(pdate), Intent(In) :: date1 ! - first date -
!
! Output scalars
  Type(pprd), Intent(Out) :: ssnlims ! - season limits -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod
!
! Exeutable Statements
!
! Set season limits
  Call reset_iseq (isq_new=isq_dy)
  ssnlims%sdate%iyr = 0
  ssnlims%sdate%imn = date1%imn
  ssnlims%sdate%idy = date1%idy
  Select Case (isq)
   Case Default
     ssnlims%edate = ssnlims%sdate + (nses*isq - 1)
   Case (3, 4)
     If (Mod(nses,2) == 0) Then
        ssnlims%edate = ssnlims%sdate + ((nses/2)*7 - 2)
     Else
        ssnlims%edate = ssnlims%sdate + ((nses/2)*7 + 2)
     End If
  End Select
  ssnlims%edate%iyr = 0
  lslim = .true.
  nse = nses
  Call reset_iseq ()
!
  Return
 End Subroutine set_ssnlimits
!
!
!
 Subroutine identify_seasons (nsn, xprd, yprd)
!
! Identifies season definitions from X field settings
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_so
  Use time,                Only: pprd, &
                                 Operator(+), Operator(>), &
                                 iseq, &
                                 adjust_pprd, reset_iseq
  Use iofiles,             Only: mfile, xfile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nsn ! - number of distinct seasons -
!
! Input/output arrays
  Type(pprd), Dimension(:), Intent(InOut) :: xprd ! - X field target periods -
  Type(pprd), Dimension(:), Intent(InOut) :: yprd ! - Y field target periods -
!
! Locals
!
! Local scalars
  Integer :: ifd    ! - field index -
  Integer :: isn    ! - season index -
  Integer :: isn2   ! - season index -
  Integer :: iyd    ! - year difference -
!
  Logical :: lso ! - seasonally-overlapping sequencing? -
!
! Executable Statements
!
! Identify season
  lso = (iseq == isq_so)
  Call reset_iseq (isq_new=isq_mn)
  isn = 0
  Do ifd = 1, nfd
     If (.not.ldsn(ifd)) Cycle
     isn = isn + 1
     lsn(isn) = len_ssn(xprd(ifd)%sdate, xprd(ifd)%edate, xfile%iseq)
     yprd(isn)%sdate%imn = xprd(isn)%sdate%imn
     yprd(isn)%edate = yprd(isn)%sdate + (lsn(isn) - 1)
     If (ifd == 1) Then
        If (mfile%prd1%sdate > xprd(ifd)%sdate) Then
           Call adjust_pprd (yprd(isn), mfile%prd1%sdate, xfile%iseq, .true.)
        End If
     Else
        iyd = xprd(isn)%sdate%iyr - xprd(1)%sdate%iyr
        yprd(isn)%sdate%iyr = yprd(1)%sdate%iyr + iyd
        yprd(isn)%edate = yprd(isn)%sdate + (lsn(isn) - 1)
        If (iyd < 0) Then ! - adjust if season is earlier than first -
           Call adjust_pprd (yprd(isn), mfile%prd1%sdate, xfile%iseq, .true.)
           iyd = xprd(1)%sdate%iyr - xprd(isn)%sdate%iyr
           Do isn2 = 1, isn-1
              yprd(isn2)%sdate%iyr = yprd(isn2)%sdate%iyr + iyd
              yprd(isn2)%edate = yprd(isn2)%sdate + (lsn(isn2) - 1)
           End Do
        End If
     End If
  End Do
!
! Set overlapping seasons
  If (lso) Then
     lsn(2:) = lsn(1)
     Do isn = 2, nsn
        yprd(isn) = yprd(isn-1) + 1
     End Do
  End If
  Call reset_iseq ()
!
! Switch off persistence component
  ispi = lsn(1)
  lsn(0) = 0
  lprst = .false.
!
  Return
 End Subroutine identify_seasons
!
!
!
 Subroutine set_seasons (nsn, cssny, tperiod)
!
! Sets season definitions
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr
  Use time,                Only: Operator(==), Operator(<), Operator(>), Operator(-), Operator(+), &
                                 iafter, iseq, &
                                 add_to_month, adjust_pprd, date_diff, set_bdate, set_fdate
  Use iofiles,             Only: mfile, xfile, yfile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nsn ! - number of seasons -
!
! Input/output scalars
  Type(pprd), Dimension(lb:), Intent(InOut) :: tperiod ! - Y field target periods -
!
! Output arrays
  Character(Len=*), Dimension(lb:), Intent(Out) :: cssny ! - Y field seasons -
!
! Locals
!
! Local scalars
  Integer :: ilsn  ! - season length -
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Ubound
!
! Executable Statements
!
! Set Y-file season
  yfile%prd1%sdate = tperiod(1)%sdate
  yfile%prd1%edate = tperiod(1)%edate
  If (lprst) Then
     lensn = ispi
  Else
     If (nsn == 1) Then
        lensn = lsn(1)
     Else
        lensn = date_diff(tperiod(1)%sdate, tperiod(nsn)%edate, isq_mn) + 1
     End If
  End If
  If (nsn == 1 .and. Ubound(tperiod, Dim=1) > 1) tperiod(1:) = tperiod(1) ! - apply to any additional fields -
!
! Determine last season available
  Select Case (iseq)
   Case Default
     yfile%prdn%edate%imn = yfile%prd1%edate%imn
     yfile%prdn%sdate = yfile%prdn%edate
     yfile%prdn%sdate%imn = add_to_month(yfile%prdn%sdate%imn, 1-lensn, &
                            iyr=yfile%prdn%sdate%iyr)
     Call adjust_pprd (yfile%prdn, mfile%prdn%edate, xfile%iseq, .false.)
   Case (isq_sn, isq_so)
     yfile%prdn = yfile%prd1
     Do
        If (yfile%prdn%edate < mfile%prdn%edate) Then
           yfile%prdn = yfile%prdn + 1
        Else
           Exit
        End If
     End Do
  End Select
!
! Set appropriate start date
  yfile%fdate%imn = yfile%prd1%sdate%imn
  yfile%bdate = yfile%fdate
  If (xfile%fdate%imn == 0 .and. iseq == isq_yr) Then
     xfile%fdate%imn = xfile%prd1%sdate%imn
     xfile%fdate%idy = xfile%prd1%sdate%idy
     xfile%bdate%imn = xfile%prd1%sdate%imn
     xfile%bdate%idy = xfile%prd1%sdate%idy
  End If
! - set appropriate lag if necessary -
  Call set_bdate (lb, lsn(lb:), nsn, yfile%prd1%sdate, yfile%fdate, yfile%bdate, &
       d1=tperiod(lb)%sdate, d2=tperiod(lb)%edate)
  Call set_fdate (iafter, xfile%iseq, xfile%fdate, mfile%prd1%sdate, yfile%bdate, yfile%fdate, ifail)
  If (ifail /= 0) Call set_fdate (-iafter, xfile%iseq, yfile%fdate, xfile%prd1%sdate, xfile%bdate, xfile%fdate, ifail)
  mfile%bdate = yfile%bdate
  mfile%fdate = yfile%fdate
!
! Set season labels
  Select Case (iseq)
   Case Default
     ilsn = len_ssn(xfile%prd1%sdate, xfile%prd1%edate, xfile%iseq)
     Call label_seasons ((/ilsn/), &
          ilsn=lsn(1), tperiod=tperiod, cssny=cssny)
   Case (isq_sn, isq_so)
     Call label_seasons (lsn(1:), &
          ilsn=1, tperiod=tperiod, cssny=cssny)
  End Select
!
  Return
 End Subroutine set_seasons
!
!
!
 Subroutine available_ssn (isq, nt)
!
! Identifies months / sesons with fewer cases than implied by number of years
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, nmn
  Use time,                Only: iseq, nss, &
                                 date_diff
  Use fields,              Only: yfield
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - time sequencing -
  Integer, Intent(In) :: nt  ! - number of cases -
!
! Locals
!
! Local scalars
  Integer :: i1 ! - lower bound -
  Integer :: ns ! - number of seasons -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Min, Mod
!
! Executable Statements
!
! Identify any seasons with fewer cases
! - identify number of seasons -
  Select Case (isq)
   Case (isq_mn)
     ns = Min(nmn, lensn)
     i1 = lb
   Case (isq_sn)
     ns = nss
     i1 = 1
   Case (isq_so)
     ns = nmn
     i1 = 1
   Case Default
     ns = 1
     i1 = 1
  End Select
  Allocate (lssr(i1:ns))
  lssr(:) = .false.
! - single seasons have the specified number of years by design -
  If (ns == 1) Then
     Return
  Else
! - some seasons may have fewer years if there are multiple seasons -
     Select Case (iseq)
      Case (isq_sn)
        i1 = Mod(nt, nss)
        If (i1 == 0) i1 = nsn
        Select Case (isq)
         Case (isq_mn)
           i1 = date_diff(yfield(1)%tprd%sdate, yfield(i1)%tprd%edate, isq_mn) + 1
           i1 = Mod(i1, lensn)
        End Select
      Case (isq_so)
        Select Case (isq)
         Case (isq_mn)
           i1 = Mod(nt+lensn-nmn, nmn)
         Case (isq_so)
           i1 = Mod(nt, nmn)
        End Select
      Case (isq_mn)
        i1 = Mod(nt, nmn)
     End Select
     If (i1 > 0) lssr(i1+1:) = .true.
  End If
!
  Return
 End Subroutine available_ssn
!
!
!
 Subroutine label_seasons (lsns, &
            ilsn, tperiod, cssny)
!
! Labels seasons
!
! Modules
  Use time,    Only: Operator(+), &
                     get_cdate, get_cprd
  Use screen,  Only: window_update
  Use iofiles, Only: mfile, xfile, yfile
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Integer, Intent(In), Optional :: ilsn  ! - season length -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: lsns ! - lengths of seasons -
!
! Input arrays
! - optional input arrays -
  Type(pprd), Dimension(lb:), Intent(In), Optional :: tperiod ! - Y field target periods -
!
! Output arrays
! - optional output arrays -
  Character(Len=*), Dimension(lb:), Intent(Out), Optional :: cssny ! - Y field seasons -
!
! Locals
!
! Local scalars
  Integer :: isn ! - season index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!

  If (nsn == 1) Then
     xfile%cssn = get_cssn(xfile%prd1, xfile%iseq, nsn)
  Else
     xfile%cssn = get_cssn(xfile%prd1, xfile%iseq, 1)
     Do isn = 2, nsn
        xfile%cssn = Trim(xfile%cssn)//' '//get_cssn(xfile%prd1+(isn-1), xfile%iseq, 1)
     End Do
  End If
  If (xfile%cfile%lset) xfile%cprdb = get_cprd(xfile%bdate, xfile%prd1%sdate, xfile%iseq,  1,  lsns(:), nsn)
  If (yfile%cfile%lset) yfile%cprdb = get_cprd(yfile%bdate, yfile%prd1%sdate, yfile%iseq, lb, lsn(lb:), nsn)
  Call window_update (xfile%cprdb)
  Call window_update (yfile%cprdb)
  If (Present(cssny)) Then
     yfile%cprd1 = get_cdate(yfile%prd1, ilsn, 2)
     yfile%cprdn = get_cdate(yfile%prdn, ilsn, 2)
     yfile%cssn = get_cssn(yfile%prd1, yfile%iseq, nsn)
     mfile%cssn = get_cssn(mfile%prd1, mfile%iseq, 1)
     Do isn = lb, nsn
        cssny(isn) = get_cssn(tperiod(isn), yfile%iseq, 1)
     End Do
  End If
!
  Return
 End Subroutine label_seasons
!
!
!
 Function get_cssn(p, isq, nsn) &
          Result (cssn)
!
! Creates season as a character string
!
! Modules
  Use data_text,           Only: ilang
  Use data_time_constants, Only: cma, cmon, isq_mn, nmn
  Use time,                Only: pprd, &
                                 Operator(==), &
                                 date_diff
!
! Function type
  Character(Len=lprd) :: cssn
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - sequence -
  Integer, Intent(In) :: nsn ! - number of seasons -
!
  Type(pprd), Intent(In) :: p ! - period -
!
! Locals
!
! Local scalars
  Integer :: ls ! - length of season -
!
! Executable Statements
!
! Create season
  If (isq <= 0 .and. p%sdate%imn /= 0 .and. p%edate%imn /= 0) Then
     If (p%edate == p%sdate) Then
        cssn = cmon(p%sdate%imn, ilang)
     Else If (nsn <= 1) Then
        If (p%edate%iyr == p%sdate%iyr) Then
           If (p%sdate%imn == p%edate%imn) Then
              cssn = cmon(p%sdate%imn,ilang)
           Else
              If (p%edate%idy == p%sdate%idy) Then
                 Write (Unit=cssn, Fmt='(A)') cma(ilang)(p%sdate%imn:p%edate%imn)
              Else
                 Write (Unit=cssn, Fmt='(I2.2,3A,I2.2,2A)') &
                    p%sdate%idy, ' ', cmon(p%sdate%imn,ilang), ' / ', p%edate%idy, ' ', cmon(p%edate%imn,ilang)
              End If
           End If
        Else
           If (p%sdate%idy == 0) Then
              ls = date_diff(p%sdate, p%edate, isq_mn)
              Select Case (ls)
               Case (:4)
                 Write (Unit=cssn, Fmt='(A)') cma(ilang)(p%sdate%imn:p%edate%imn+nmn)
               Case Default
                 Write (Unit=cssn, Fmt='(3A)') cmon(p%sdate%imn,ilang), ' - ', cmon(p%edate%imn,ilang)
              End Select
           Else
              Write (Unit=cssn, Fmt='(I2.2,3A,I2.2,2A)') &
                 p%sdate%idy, ' ', cmon(p%sdate%imn,ilang), ' / ', p%edate%idy, ' ', cmon(p%edate%imn,ilang)
           End If
        End If
     Else
        Write (Unit=cssn, Fmt='(3A)') &
           cma(ilang)(p%sdate%imn:p%sdate%imn+lsn(1)-1), ' - ', cma(ilang)(p%edate%imn+nmn-(lsn(nsn)-1):p%edate%imn+nmn)
     End If
  Else
     cssn = ' '
  End If
!
  Return
 End Function get_cssn
!
!
!
 Function get_cssel(imd0, imd1, ls, nld, ispi) &
          Result (cssel)
!
! Returns current season selection
!
! Modules
  Use data_text,           Only: ilang
  Use data_time_constants, Only: cma, cmon, nmn
  Use labels,              Only: cg_initin_l
  Use time,                Only: add_to_month
!
! Function type
  Character(Len=lssel) :: cssel
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imd0 ! - month of forecast initialization -
  Integer, Intent(In) :: imd1 ! - first month of season -
  Integer, Intent(In) :: ls   ! - length of season -
  Integer, Intent(In) :: nld  ! - number of seasons -
  Integer, Intent(In) :: ispi ! - length of SPI -
!
! Locals
!
! Local scalars
  Integer :: l  ! - season index -
  Integer :: lp ! - length of persistence season -
  Integer :: m0 ! - zeroth month -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Current season(s)
  If (ls == 1) Then
     cssel = cmon(imd1, ilang)
  Else
     cssel = cma(ilang)(imd1:imd1+ls-1)
  End If
  Select Case (nld)
   Case (1)
     Continue
   Case (2:4)
     Do l = 2, nld
        If (ls == 1) Then
           If (imd1+l-1 <= nmn) Then
              cssel = Trim(cssel)//', '//cmon(imd1+l-1,ilang)
           Else
              cssel = Trim(cssel)//', '//cmon(imd1+l-(nmn+1),ilang)
           End If
        Else
           m0 = add_to_month(imd1, (l-1)*ls-1)
           cssel = Trim(cssel)//', '//cma(ilang)(m0+1:m0+ls)
        End If
     End Do
   Case Default
     If (ls == 1) Then
        If (imd1+nld-1 <= nmn) Then
           cssel = Trim(cssel)//' - '//cmon(imd1+nld-1,ilang)
        Else
           cssel = Trim(cssel)//' - '//cmon(imd1+nld-(nmn+1),ilang)
        End If
     Else
        m0 = add_to_month(imd1, (nld-1)*ls-1)
        cssel = Trim(cssel)//' - '//cma(ilang)(m0+1:m0+ls)
     End If
  End Select
!
! Add start date
  If (imd0 > 0) cssel = Trim(cssel)//' '//Trim(cg_initin_l)//' '//cmon(imd0,ilang)
  If (ispi > 0) Then
     cssel = 'Forecast: '//Trim(cssel)
     lp = ispi - ls
     m0 = add_to_month(imd1, -(lp+1))
     If (m0 ==12) m0 = 0
     Select Case (lp)
      Case (:0)
        cssel = 'Persistence: N/A; '//Trim(cssel)
      Case (1)
        cssel = 'Persistence: '//cmon(m0+1,ilang)//'; '//Trim(cssel)
      Case Default
        cssel = 'Persistence: '//cma(ilang)(m0+1:m0+lp)//'; '//Trim(cssel)
     End Select
  End If
!
  Return
 End Function get_cssel
!
!
!
 Function update_ctsel() &
          Result (update)
!
! Updates target season selection indicator
!
! Modules
  Use screen, Only: window_update
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Update season selection indicator
  cssel = get_cssel(0, pperiod%sdate%imn, lsn(1), 1, 0)
  Call window_update (cssel)
  update = 2
!
  Return
 End Function update_ctsel
!
!
!
 Function update_cl0() &
          Result (update)
!
! Updates target season selection indicator given first month of persistence
!
! Modules
  Use data_time_constants, Only: nmn
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Update season selection indicator
  ispi = pperiod%sdate%imn - l0
  If (ispi < 0) ispi = ispi + nmn
  ispi = ispi + lsn(1)
  update = update_cspi()
!
  Return
 End Function update_cl0
!
!
!
 Function update_cspi() &
          Result (update)
!
! Updates target season selection indicator given SPI length
!
! Modules
  Use screen, Only: window_update
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Update season selection indicator
  cssel = get_cssel(0, pperiod%sdate%imn, lsn(1), 1, ispi)
  Call window_update (cssel)
  update = 2
!
  Return
 End Function update_cspi
!
!
!
 Subroutine set_icsv (icsv, units)
!
! Determines whether to calculate seasonal averages or seasonal aggregations
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: units ! - field units -
!
! Input/output scalars
  Integer, Intent(InOut) :: icsv ! - calculate seasonal values flag -
!
! Executable Statements
!
! Determine appropriate means of seasonal aggregation
  Select Case (Trim(units))
   Case ('K', 'Kelvin_scale', 'C', 'degree_Celsius', 'Celsius_scale', 'F', 'Farenheit') ! - temperatures -
     icsv = 1
   Case ('m/s', 'meters/second', 'm s-1', 'mm/day', 'cm/day') ! - rates -
     icsv = 1
   Case ('m', 'meters', 'cm', 'mm', 'mm/month') ! - depths and monthly rates -
     icsv = 2
   Case Default
     If (icsv == 0) icsv = 1
  End Select
!
  Return
 End Subroutine set_icsv
!
!
!
 Subroutine seasonal (bdate, n, lb, imn1, rmiss, &
            iyrs)
!
! Converts monthly to seasonal data
!
! Modules
  Use data_time_constants, Only: isq_sn, isq_so, isq_yr
  Use time,                Only: pdate, &
                                 Operator(+), Operator(-), &
                                 iseq
  Use arrays,              Only: y
  Use iofiles,             Only: ny
  Use settings,            Only: ntm
  Use fields,              Only: yfield
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lb ! - lower bound -
  Integer, Intent(In) :: n  ! - number of cases -
!
  Type(pdate), Intent(In) :: bdate ! - first date of data -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: imn1 ! - first month of season -
!
  Real(Kind=rp), Dimension(lb:), Intent(In) :: rmiss ! - missing values -
!
! - optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: iyrs ! - years of data -
!
! Locals
!
! Local scalars
  Integer :: i     ! - locator -
  Integer :: isn   ! - season index -
  Integer :: ltrim ! - trimmed length -
  Integer :: ns    ! - number of seasons -
  Integer :: ifail ! - error indicator -
!
  Character(Len=32) :: cprd ! - period -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any, Index, Len_Trim, Trim
!
! Executable Statements
!
! Determine parameters
  If (lprst) Then
     ns = 1
  Else
     ns = nsn
  End If
!
! Calculate seasonal aggregates
  Allocate (vs(ny,n,lb:nsn), Stat=ifail)
  If (ifail /= 0) GoTo 1
  Call calc_seasonal (icsv, lprst, ny, n, ns, lssr(lb:), lsn(lb:), lb, imn1(lb:), bdate, rmiss(lb:), y(:,:,:), vs(:,:,lb:), &
       iyrs=iyrs)
!
! Reallocate memory
  Deallocate (y)
  Allocate (y(ny,n,lb:nsn), Stat=ifail)
  If (ifail /= 0) GoTo 1
  y(:,:,:) = vs(:,:,:)
  Deallocate (vs)
!
! Adjust units if monthly are accumulated
  If (icsv == 2) Then
     Do isn = lb, nsn-(1-lb)
        i = Index(yfield(isn)%unit, '/')
        If (i > 0) Then
           ltrim = Len_Trim(yfield(isn)%unit)
           cprd = yfield(isn)%unit(i+1:ltrim)
           Write (yfield(isn)%unit(i+1:), Fmt='(I0,2A)') lsn(isn), '-', Trim(cprd)
        End If
        i = Index(yfield(isn)%unitc, '/')
        If (i > 0) Then
           ltrim = Len_Trim(yfield(isn)%unitc)
           cprd = yfield(isn)%unitc(i+1:ltrim)
           Write (yfield(isn)%unitc(i+1:), Fmt='(I0,2A)') lsn(isn), '-', Trim(cprd)
        End If
     End Do
  End If
!
! Reset numbers of cases once seasonal totals are calculated
  Select Case (iseq)
   Case (isq_yr)
     If (.not.Any(lssr(:))) ntm = n
   Case (isq_sn, isq_so)
     ntm = n
  End Select
!
  Return
!
1 ifail = 1
  Return
 End Subroutine seasonal
!
!
!
 Subroutine seasonal_zyp (imn1, rmiss)
!
! Converts monthly to seasonal data for persistence components of forecasts
!
! Modules
  Use data_time_constants, Only: isq_yr
  Use time,                Only: Operator(+), Operator(-), &
                                 iseq, &
                                 nyears
  Use arrays,              Only: zyp
  Use iofiles,             Only: ny, pfile
  Use settings,            Only: nmf, nf
!
! Arguments
!
! Input arrays
  Integer, Dimension(0:), Intent(In) :: imn1 ! - first month of season -
!
  Real(Kind=rp), Dimension(0:), Intent(In) :: rmiss ! - missing values -
!
! Locals
!
! Local scalars
  Integer :: nyr   ! - number of years -
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Calculate seasonal aggregates from training data
  nyr = nyears(nf, .false.)
  Allocate (vs(ny,nyr,0:0), Stat=ifail)
  If (ifail /= 0) GoTo 1
  Call calc_seasonal (icsv, lprst, ny, nyr, 1, (/.false./), lsn(0:), 0, imn1(0:), pfile%bdate, rmiss(0:), zyp(:,:,:), vs(:,:,0:))
! - reassign memory -
  Deallocate (zyp)
  Allocate (zyp(ny,nyr,0:0), Stat=ifail)
  If (ifail /= 0) GoTo 1
  zyp(:,:,:) = vs(:,:,:)
  nf = nyr
  Deallocate (vs)
  If (iseq == isq_yr) nmf = nyr
!
  Return
!
1 ifail = 1
  Return
 End Subroutine seasonal_zyp
!
!
!
 Subroutine calc_seasonal (icsv, lcs, m, nyr, nsn, lssr, lsn, lb, imn1, d1, rmiss, vm, vs, &
            iyrs)
!
! Calculates seasonal averages/totals
!
! Modules
  Use data_numbers,        Only: zero
  Use data_time_constants, Only: nmn
  Use maths,               Only: approx
  Use time,                Only: pdate, &
                                 add_to_month, ndays
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icsv ! - seasonal aggregation method -
  Integer, Intent(In) :: m    ! - number of locations -
  Integer, Intent(In) :: nyr  ! - number of years -
  Integer, Intent(In) :: nsn  ! - number of seasons -
  Integer, Intent(In) :: lb   ! - lower bound -
!
  Logical, Intent(In) :: lcs ! - consecutive seasons? -
!
  Type(pdate), Intent(In) :: d1 ! - first date available -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsn  ! - length of season -
  Integer, Dimension(lb:), Intent(In) :: imn1 ! - first month of season -
!
  Real(Kind=rp), Dimension(lb:), Intent(In) :: rmiss ! - missing values -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: vm ! - monthly values (dimensions: m,nyr,lensn) -
!
  Logical, Dimension(lb:), Intent(In) :: lssr ! - season has shortened record? -
!
! - optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: iyrs ! - years of data -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:,lb:), Intent(Out) :: vs ! - seasonal values (dimensions: m,nyr,nsn) -
!
! Locals
!
! Local scalars
  Integer :: i   ! - location index -
  Integer :: k   ! - year index -
  Integer :: kk  ! - year index -
  Integer :: l   ! - month index -
  Integer :: ll  ! - month index -
  Integer :: isn ! - season index -
  Integer :: is1 ! - index of first month of season -
  Integer :: isl ! - index of last month of season -
  Integer :: imn ! - current month -
  Integer :: iyr ! - current year -
  Integer :: lmn ! - length of month -
  Integer :: lsd ! - length of season in days -
  Integer :: nya ! - number of years adjustment -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Real
!
! Executable Statements
!
! Initialise
  If (icsv == 2) Then
     lmn = 1
     lsd = 1
  End If
  Do isn = lb, nsn
     vs(:,:,isn) = rmiss(isn)
     If (lcs .and. isn > lb) Then
        is1 = isl + 1
     Else
        is1 = add_to_month(imn1(isn), 1-d1%imn)
     End If
     isl = is1 + lsn(isn) - 1
!
! Calculate seasonal aggregate
     If (lssr(isn)) Then
        nya = 1
     Else
        nya = 0
     End If
     Do k = 1, nyr-nya
        If (.not.Present(iyrs)) Then
           iyr = d1%iyr + k - 1
        Else
           iyr = iyrs(k)
        End If
        imn = imn1(isn) - 1
        If (icsv == 1) lsd = 0
        vs(:,k,isn) = zero
        Do l = is1, isl
           imn = add_to_month(imn, 1, &
                 iyr=iyr)
           If (icsv == 1) Then ! - seasonal average -
              lmn = ndays(iyr, imn)
              lsd = lsd + lmn
           End If
           Select Case (l)
            Case (1:nmn)
              ll = l
              kk = k
            Case (nmn+1:)
              ll = l - nmn
              kk = k + 1
           End Select
           Do i = 1, m
              If (approx(vm(i,kk,ll), rmiss(isn))) Then
                 vs(i,k,isn) = rmiss(isn)
              Else If (vs(i,k,isn) /= rmiss(isn)) Then
                 vs(i,k,isn) = vs(i,k,isn) + vm(i,kk,ll)*Real(lmn, Kind=rp)
              End If
           End Do
        End Do
        Where (vs(:,k,isn) /= rmiss(isn)) vs(:,k,isn) = vs(:,k,isn)/Real(lsd, Kind=rp)
     End Do
  End Do
!
  Return
 End Subroutine calc_seasonal
!
!
!
 Subroutine proj_read_season (iin, pver, ifail)
!
! Reads program settings from project file
!
! Modules
  Use data_numbers,        Only: sp, tolv
  Use data_time_constants, Only: isq_so
  Use time,                Only: iseq, lslim, nse, nss, ssnlims
  Use iofiles,             Only: xfile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin ! - input file unit number -
!
  Real(Kind=sp), Intent(In) :: pver ! - project file version -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Clear existing settings
  Call close_season ()
!
! Read program settings
  If (pver > 16.0_sp) Then
     Read (Unit=iin, Err=1, End=2) icsv, ispi, lb, lensn, lprst, lensn, nsn, nss, nse, iseq, lslim, ssnlims
  Else If (pver > 14.07_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) icsv, ispi, lb, lensn, lprst, lensn, nsn, nss, iseq
     nse = 0
     lslim = .false.
  Else
     Read (Unit=iin, Err=1, End=2) icsv, ispi, lb, lensn, lprst, lensn, nsn
     nse = 0
     lslim = .false.
  End If
!
! Read season settings if saved
  If (icsv > 0) Then
     If (pver < 14.07_sp-tolv) Then
        Select Case (iseq)
         Case Default
           nss = nsn
         Case (isq_so)
           nss = 1
        End Select
     End If
     nfd = xfile%nfs*nss
     Allocate (ldsn(nfd), Stat=ifail)
     If (ifail /= 0) GoTo 3
     Deallocate (lsn)
     Allocate (lsn(0:nsn), Stat=ifail)
     If (ifail /= 0) GoTo 3
     Read (Unit=iin, Err=1, End=2) ldsn, lsn
  Else
     If (pver < 14.07_sp-tolv) nss = 1
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
! - problem initializing memory -
3 ifail = 3
  Return
!
 End Subroutine proj_read_season
!
!
!
 Subroutine proj_write_season (iout, ifail)
!
! Writes program settings to project file
!
! Modules
  Use time, Only: iseq, lslim, nse, nss, ssnlims
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output file unit number -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write program settings
  Write (Unit=iout, Err=1) icsv, ispi, lb, lensn, lprst, lensn, nsn, nss, nse, iseq, lslim, ssnlims
!
! Write season settings if saved
  If (icsv > 0) Write (Unit=iout, Err=1) ldsn, lsn
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
!
 End Subroutine proj_write_season
!
!
!
 Subroutine close_season ()
!
! Frees memory allocated for calculating seasonal aggregates
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory allocated for calculating seasonal averages
  Call reset_season (.false.)
  If (Allocated(ldsn)) Deallocate (ldsn)
  If (Allocated(lssr)) Deallocate (lssr)
!
  Return
 End Subroutine close_season
End Module season
