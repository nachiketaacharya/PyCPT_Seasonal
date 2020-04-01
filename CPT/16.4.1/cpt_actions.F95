! Author: Simon Mason
Module cpt_actions
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: cpt_reset, reset, reset_mem, run_analysis
#if GUI == 1
  Public :: check_reset
#endif
!
Contains
!
!
 Function run_analysis()
!
! Reads data and fits selected model
!
! Modules
#if GUI == 1
  Use clrwin$,             Only: clear_window$
#endif
  Use data_numbers,        Only: rp, zero, one
  Use data_cpt_constants,  Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv, mnt, ng, nts
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, nmn
  Use maths,               Only: zero_bound
  Use labels,              Only: ca_atypes_a, &
                                 l_climper
  Use time,                Only: pprd, &
                                 Operator(+), Operator(>), Operator(<), &
                                 iseq, nss, &
                                 date_diff, get_cdate, nyears, reset_iseq
  Use screen,              Only: window_update
  Use gui,                 Only: set_cursor_waiting
#if GUI == 1
  Use gui,                 Only: ih_con
#endif
  Use menus,               Only: im_ov
  Use progress_meter,      Only: lstop, prog, &
                                 end_progress_meter, set_progress_increment, start_progress_meter, stop_progress_meter, &
                                 update_progress_meter
  Use errors,              Only: cpt_error
  Use arrays,              Only: vcopy, x, y, yhat, ypc, yret
  Use iofiles,             Only: mfile, mxa, mya, xfile, yfile, zfile
  Use settings,            Only: ieav, igauss, igto, iretro, istd, istdo, izero, lcw, lzero, &
                                 nt, ntr, nu, nur, nyr, &
                                 check_n, get_cv, get_retro, record_results, set_nused
  Use space,               Only: cstndy
#if GUI == 1
  Use fields,              Only: iffx, iffy, ilimx0, ilimx1, ilimx2, ilimy0, ilimy1, ilimy2, &
                                 set_ilims
#endif
  Use fields,              Only: xfield, yfield
  Use distribs,            Only: set_gamma, stdize, ustdize
  Use season,              Only: icsv, lb, lensn, lprst, lsn, lssr, nsn, &
                                 available_ssn, get_cssn, seasonal
  Use pcs,                 Only: lcde, &
                                 check_pcs, get_cca_settings, get_eofx_settings, get_eofy_settings, set_scree
  Use missing,             Only: kcuse, ktuse, &
                                 compress, replace_missing_xy
  Use categories,          Only: irobs, pthr, rfps, tobs, tobst, &
                                 calc_thresholds, init_thresholds, set_cv_categories, set_ra_categories
  Use climate,             Only: ave, clim, climate_per, icex, itrain1, itrain2, nc, ncu, ncy, ntc, sdev, &
                                 calc_climates, check_climates, get_clim_limits, init_climate, set_climate
  Use ensemble,            Only: ensemble_ave
  Use analysis,            Only: ianal, icalc, &
                                 analysis_flags, close_analysis, init_results, perform_analysis
  Use pfv,                 Only: perform_pfv
  Use cpt_input,           Only: get_data
  Use forecast,            Only: ifc, ifcast
  Use verif,               Only: irv
!
! Function type
  Integer :: run_analysis
!
! Locals
!
! Local scalars
  Integer :: j     ! - location index -
  Integer :: k     ! - time index -
  Integer :: kk    ! - time index -
  Integer :: l     ! - field index -
  Integer :: iarg1 ! - first integer error argument -
  Integer :: ierr  ! - error parameter -
  Integer :: isn   ! - season index -
  Integer :: is0   ! - zeroth season -
  Integer :: it1   ! - index of first date of training period (relative to combined climate and training dataset) -
  Integer :: it2   ! - index of last date of training period (relative to combined climate and training dataset) -
  Integer :: iz    ! - zero-bound option -
  Integer :: nact  ! - number of actions -
  Integer :: nems  ! - total number of ensemble members -
  Integer :: ntcy  ! - number of years in climatological and training periods -
  Integer :: nuy   ! - number of available years in training period -
  Integer :: ifail ! - error indicator -
!
  Logical :: lpfv ! - PFV? -
!
  Character(Len=32) :: cproc ! - error procedure -
  Character(Len=64) :: carg1 ! - first character error argument -
  Character(Len=64) :: carg2 ! - second character error argument -
  Character(Len=64) :: carg3 ! - third character error argument -
!
  Type(pprd) :: error_prd ! - period of erroneous data -
!
! Local arrays
  Integer, Dimension(3) :: iloc ! - location of first negative value -
!
  Integer, Dimension(:), Allocatable :: iyrs ! - years of data -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count, Max, MaxVal, Real, Sum, Trim
!
! Executable Statements
!
! Initialise errors
  iarg1 = 0
  carg1 = ' '
  carg2 = ' '
  carg3 = ' '
!
! Check for invalid length of training period
! - check for too few -
  lpfv = (ianal == ia_pfv)
  ifail = check_n(nt, lcw, .true., ianal, iretro, ierr)
  If (ifail /= 0) Then
     If (ifail /= 2 .and. .not.lpfv) Then
        Call window_update (lcw)
        Call cpt_error ('check_n', ifail, .false., &
             i_arg1=ierr)
        run_analysis = 1
        Return
     End If
  End If
!
! Check for cancelled retroactive calculations
  Select Case (iretro)
   Case (0) ! - cross-validation -
     run_analysis = get_cv()
   Case (1) ! - retroactive -
     If (.not.lpfv) Then
        run_analysis = get_retro()
        If (iretro == -1) Then
           iretro = 0
           run_analysis = 1
           Return
        End If
     End If
  End Select
  run_analysis = 1
!
! Initialise analysis
! - switch on hour glass -
  Call set_cursor_waiting (1)
! - initialise calculation flags -
  Call analysis_flags ('on')
  ifc = 0
!
! Check climatologies
! - set climatological probabilities -
  If (lpfv) Then
     If (Sum(xfield(1)%prob(:)) /= zero) Then
        pthr(1) = xfield(1)%prob(1)
        Do j = 2, nts
           pthr(j) = pthr(j-1) + xfield(1)%prob(j)
        End Do
     Else
        xfield(1)%prob(1) = pthr(1)
        Do j = 2, nts
           xfield(1)%prob(j) = pthr(j) - pthr(j-1)
        End Do
        xfield(1)%prob(ng) = one - pthr(nts)
     End If
  End If
! - set initial climatological period settings -
  Call check_climates (yfile, mfile, lsn(0), lpfv)
!
! Read data
  Call start_progress_meter ('Reading data', .false.)
  ifail = get_data(icsv, icex, .true.)
  If (ifail /= 0 .or. lstop) GoTo 2
  Call end_progress_meter (.false., .false.)
!
! Adjust number of cases for seasonal sequencing
  nyr = nyears(nt, .false.)
  ntcy = nyears(ntc, .false.)
!
! Check for zero-bound
  If (igto == 2 .and. igauss > 0) Then
     iz = 2
  Else If (izero == 1) Then
     iz = 3
  Else
     iz = 0
  End If
! - search for negative values if zero bound is on -
  If (iz /= 0) Then
     Call zero_bound (yfile%nfl, ntcy, yfield(:)%nv, yfield(:)%rmiss, y(:,:,:), iloc(:))
 ! - negative values found -
     If (iloc(2) /= 0) Then
        If (icsv == 0) Then
           error_prd = yfile%prd1 + (yfile%it1 + iloc(2) - 2)
        Else
           Call reset_iseq (isq_new=isq_mn)
           error_prd = mfile%prd1 + (mfile%it1 + iloc(2) - 2)
           Call reset_iseq ()
        End If
        ifail = iz
        cproc = 'run_analysis'
        carg1 = Trim(yfile%cfile%ffile)
        carg2 = get_cdate(error_prd,lsn(1),2)
        If (yfile%idstr /= 1) Then
           l = 1
           Do
              If (iloc(1) <= yfield(l)%nv) Exit
              iloc(1) = iloc(1) - yfield(l)%nv
           End Do
           j = iloc(1)
           carg3 = cstndy(j,l)
        End If
        GoTo 1
     End If
  End If
!
! Identify training period
  Select Case (icex)
   Case (0, 2, -2)
     itrain1 = 1
   Case (1, 3) ! - climatological period starts before training period -
     itrain1 = date_diff(climate_per%d1, yfile%bdate, yfile%iseq) + 1
   Case (-1) ! - climatological period is before training period (gap) -
     itrain1 = nc + 1
  End Select
!
! Replace missing values and check for missing variables
  Call start_progress_meter ('Checking for missing values', .false.)
  Call available_ssn (mfile%iseq, ntc)
  ifail = replace_missing_xy(ianal, icex, nt, ntc, ntcy, itrain1, nsn, lensn,  1, lssr(:),  .true., (icsv > 0 .and. lensn > 1), &
                             .false., nu)
  DeAllocate (lssr)
  If (ifail /= 0 .or. lstop) GoTo 2
!
! Calculate seasonal values
  If (icsv > 0 .and. lensn > 1) Then
     Print *, 'Calculating seasonal values ...'
     Call available_ssn (yfile%iseq, ntc)
     Select Case (icex) ! - (use nt or ntc instead of nu because data have not yet been compressed -
      Case (0)
        Call seasonal (mfile%bdate, nyr, lb, yfield(lb:)%tprd%sdate%imn, yfield(lb:)%rmiss)
      Case (1, 3)
        Call seasonal (climate_per%d1, ntcy, lb, yfield(lb:)%tprd%sdate%imn, yfield(lb:)%rmiss)
      Case (2)
        Call seasonal (mfile%bdate, ntcy, lb, yfield(lb:)%tprd%sdate%imn, yfield(lb:)%rmiss)
      Case (-1)
        Allocate (iyrs(ntcy))
        iyrs(1) = climate_per%d1%iyr
        Do k = 2, ncy
           iyrs(k) = iyrs(k-1) + 1
        End Do
        iyrs(ncy+1) = mfile%bdate%iyr
        Do k = ncy+2, ntcy
           iyrs(k) = iyrs(k-1) + 1
        End Do
        Call seasonal (climate_per%d1, ntcy, lb, yfield(lb:)%tprd%sdate%imn, yfield(lb:)%rmiss, &
             iyrs=iyrs)
        Deallocate (iyrs)
      Case (-2)
        Allocate (iyrs(ntcy))
        iyrs(1) = mfile%bdate%iyr
        Do k = 2, nyr
           iyrs(k) = iyrs(k-1) + 1
        End Do
        iyrs(nyr+1) = climate_per%d1%iyr
        Do k = nyr+2, ntcy
           iyrs(k) = iyrs(k-1) + 1
        End Do
        Call seasonal (climate_per%d1, ntcy, lb, yfield(lb:)%tprd%sdate%imn, yfield(lb:)%rmiss, &
             iyrs=iyrs)
        Deallocate (iyrs)
     End Select
! - repeat missing value estimation for seasonal values -
     Call start_progress_meter ('Re-checking for missing values', .false.)
     ifail = replace_missing_xy(ianal, icex, nt, ntc, ntcy, itrain1, nsn,     1, lb, lssr(1:nsn), .false.,    .false., lprst, nu)
     DeAllocate (lssr)
     If (ifail /= 0 .or. lstop) GoTo 2
  End If
!
! Identify number of non-missing years
  Select Case (iseq)
   Case Default
     nyr = Count(kcuse(itrain1:itrain1+nt-1))
   Case (isq_sn, isq_so) ! - seasonal data have not yet been compressed -
     nyr = nyears(nt, .true.)
     nur = nt
     nu = nt
  End Select
!
! Double check for invalid length of training period given missing values
  If (nu < nt) Then
     ifail = check_n(nu, lcw, .false., ianal, iretro, ierr)
     Select Case (ifail)
      Case (0)
        Continue
      Case Default
        If (lpfv .and. ifail == 2) Then
           Continue
        Else
           ifail = 5
           cproc = 'check_n'
           iarg1 = ierr
           GoTo 1
        End If
     End Select
  End If
  ifail = set_nused((iseq /= isq_sn .and. iseq /= isq_so), kcuse(itrain1:))
  If (ifail /= 0) Then
     cproc = 'set_nused'
     iarg1 = mnt - 1
     GoTo 1
  End If
  Call end_progress_meter (.false., .false.)
!
! Calculate ensemble average
  If (ieav == 1 .and. .not.lpfv) Then
     If (xfile%nem > 1 .or. xfile%nms > 1) Then ! - average ensemble members -
        nems = Max(1,xfile%nem)*Max(1,xfile%nms)
        Call ensemble_ave (xfile%nfs/nems, xfile%nls, nu, xfield(:)%nv, nems, x(:,:,1))
     End If
  End If
!
! Initialise progress meter
  Call start_progress_meter ('Calculating climatologies and thresholds', .false.)
  nact = 1              ! - initial checks -
  nact = nact + mya + 1 ! - calc_climate -
  nact = nact + 1       ! - standardize -
  nact = nact + mya     ! - set_percentiles -
  Select Case (igauss)  ! - set_thresholds -
   Case (0)
     nact = nact + mya
   Case Default
     nact = nact + mya + 1
  End Select
  Call set_progress_increment (one/Real(nact+1, Kind=rp), .false.)
!
! Set climatology
  If (.not.lprst) Then
     Call init_climate (1, nsn, ifail)
  Else
     Call init_climate (0, nsn, ifail)
  End If
  If (lstop) GoTo 2
  If (ifail /= 0) Then
     cproc = 'init_climate'
     GoTo 1
  End If
  If (init_thresholds(nsn) /= 0) Then
     cproc = 'init_thresholds'
     GoTo 1
  End If
  Call set_climate (1, lsn(0))
  Print *, Trim(l_climper%c)//': '//Trim(climate_per%clim1)//' to '//Trim(climate_per%clim2)
  If (update_progress_meter(.false.) /= 0) GoTo 2
!
! Calculate climatologies
  nuy = nyears(nu, .false.)
! - add persistence and predicted components of season -
  If (lprst) Then
     Allocate (ypc(mya,nuy,lb:nsn), Stat=ifail)
     If (ifail /= 0) Then
        ifail = 4
        cproc = 'run_analysis'
        carg1 = 'calculating persistence component'
        GoTo 1
     End If
     ypc(:,:,:) = y(1:mya,1:nuy,:)
     y(1:mya,1:nuy,0) = y(1:mya,1:nuy,0) + y(1:mya,1:nuy,1)
  End If
! - calculate climatologies -
  Call calc_climates (.false., lb, nsn, mya, y(1:mya,:,:), ave(:,:), sdev(:,:), ncu(:), clim(1:mya,:,:), ifail)
  If (lstop) GoTo 2
  If (ifail /= 0) GoTo 2
  If (igauss == 2) Call set_gamma (.true.)
!
! Isolate training period if climatological period extends outside it
  If (icex /= 0) Then
     Allocate (vcopy(mya,nuy,nsn), Stat=ifail)
     If (ifail /= 0) Then
        ifail = 4
        cproc = 'run_analysis'
        carg1 = 'climatological period outside of training period'
        GoTo 1
     End If
     vcopy(:,:,:) = y(1:mya,itrain1:itrain2,1:nsn)
     Deallocate (y)
     Allocate (y(mya,nuy,nsn), Stat=ifail)
     If (ifail /= 0) Then
        ifail = 4
        cproc = 'run_analysis'
        carg1 = 'climatological period outside of training period'
        GoTo 1
     End If
     y(:,:,:) = vcopy(:,:,:)
     Deallocate (vcopy)
  End If
  Call get_clim_limits (icex, isq_mn, &
       itrain1=it1)
  it2 = it1 + nt - 1
  Allocate (ktuse(nt))
  ktuse(:) = kcuse(it1:it2)
!
! Standardize
  If (.not.lpfv) Then ! - for PFV, istd = 0 so standardization is unnecessary-
! - training data -
     Do isn = 1, nsn ! - start at 1 to leave combined season unstandardized (nsn is always 1) -
        Call stdize (istd, igauss, isn, mya, nuy, ave(:,isn), sdev(:,isn), y(1:mya,1:nuy,isn))
     End Do
! - climatological data -
     Do isn = lb, nsn
        Call stdize (istd, igauss, isn, mya, ncu(isn), ave(:,isn), sdev(:,isn), clim(1:mya,1:ncu(isn),isn))
     End Do
  End If
  If (update_progress_meter(.false.) /= 0) GoTo 2
!
! Calculate thresholds
  Call calc_thresholds (igauss, mya, nsn, ncu(lb:), clim(:,:,lb:), tobs, tobst, ifail)
  If (ifail /= 0) GoTo 2
  If (update_progress_meter(.false.) /= 0) GoTo 2
! - temporarily remove persistence component -
  If (lprst) Then
     y(1:mya,1:nuy,:) = ypc(:,:,:)
     If (igauss /= 0) Call stdize (istd, igauss, 1, mya, nuy, ave(:,1), sdev(:,1), y(1:mya,1:nuy,1)) ! - restandardise -
  End If
!
! Check validity of EOF and CCA settings
  Call start_progress_meter ('Beginning analysis',.false.)
  If (.not.lcde) Then
     ifail = check_pcs(ianal)
     If (ifail == -1) GoTo 2
  Else ! - calculate mode numbers -
     Select Case (ianal)
      Case (ia_cca) ! - CCA -
        ifail = get_eofx_settings()
        ifail = get_eofy_settings()
        ifail = get_cca_settings()
      Case (ia_pcr) ! - PCR -
        ifail = get_eofx_settings()
     End Select
     ifail = 0
  End If
!
! Perform analysis
  istdo = istd
  lzero = (igauss == 0 .and. izero == 1)
  Select Case (ianal)
   Case (ia_cca, ia_pcr, ia_mlr, ia_gcm) ! - CCA, PCR, MLR, GCM -
     ifail = perform_analysis()
   Case (ia_pfv) ! - PFV -
     ifail = perform_pfv(ierr)
  End Select
  If (lstop) GoTo 2
  Select Case (ifail)
   Case (0)
     Continue
   Case (-1)
     GoTo 2
   Case (11:)
     ifail = ifail - 10
     cproc = 'match_fields'
     carg1 = 'PFV'
     carg2 = yfield(ierr)%var
     carg3 = yfield(ierr)%cssn
     iarg1 = ierr
     GoTo 1
   Case Default
     cproc = 'perform_'//ca_atypes_a(ianal)
     carg1 = ca_atypes_a(ianal)
     GoTo 1
  End Select
!
! Add persistence component to cross-validated hindcasts
  If (lprst) Then
     If (igauss == 2) Call set_gamma (.false.)
     Select Case (igauss)
      Case (0)
        yhat(:,:) = yhat(:,:) + y(1:mya,1:nuy,0)
        If (istd /= 0) Call stdize (istd, igauss, 0, mya, nuy, ave(:,0), sdev(:,0), yhat(:,:))
        If (iretro == 1) Then
           yret(:,:) = yret(:,:) + y(1:mya,1:nur,0)
           If (istd /= 0) Call stdize (istd, igauss, 0, mya, nur, ave(:,0), sdev(:,0), yret(:,:))
        End If
        y(1:mya,1:nuy,0) = y(1:mya,1:nuy,0) + y(1:mya,1:nuy,1) ! - replace persistence component -
        If (istd /= 0) Then
           Do isn = lb, nsn
              Call stdize (istd, igauss, isn, mya, nuy, ave(:,isn), sdev(:,isn), y(1:mya,1:nuy,isn))
           End Do
        End If
      Case (1, 2)
        y(1:mya,1:nuy,:) = ypc(:,:,:)
        If (istd /= 0) Call ustdize (istd, igauss, 1, mya, nuy, ave(:,1), sdev(:,1), yhat(:,:))
        yhat(:,:) = yhat(:,:) + y(1:mya,1:nu,0)
        If (istd /= 0) Call stdize (istd, igauss, 0, mya, nuy, ave(:,0), sdev(:,0), yhat(:,:))
        y(1:mya,1:nuy,0) = y(1:mya,1:nuy,0) + y(1:mya,1:nuy,1) ! - replace persistence component -
        If (istd /= 0) Then
           Do isn = lb, nsn
              Call stdize (istd, igauss, isn, mya, nuy, ave(:,isn), sdev(:,isn), y(1:mya,1:nuy,isn))
           End Do
        End If
     End Select
! - add predicted component to persistence season -
     yfield(0)%tprd%edate = yfield(1)%tprd%edate
     yfield(0)%cssn = get_cssn(yfield(0)%tprd,yfile%iseq,1)
  End If
  If (lstop) GoTo 2
!
! Allocate additional memory
  Call init_results (ianal, Max(MaxVal(xfield(:)%nlt), MaxVal(yfield(:)%nlt)), Max(MaxVal(xfield(:)%nlg), MaxVal(yfield(:)%nlg)), &
       ifail)
  If (ifail /= 0) Then
     ifail = 1
     cproc = 'perform_'//ca_atypes_a(ianal)
     carg1 = ca_atypes_a(ianal)
     GoTo 1
  End If
!
! Determine forecast categories
  Print *, 'Identifying categories ...'
  If (.not.lpfv) Then 
     If (iretro == 1) Then
        Call set_progress_increment ((one - prog)*Real(nu, Kind=rp)/Real(mya*(nu + nur), Kind=rp),.false.)
        is0=0
     Else
        Call set_progress_increment ((one - prog)/Real(mya, Kind=rp),.false.)
     End If
     If (set_cv_categories(lb) /= 0) GoTo 2
  Else
     If (nsn /= nmn) Then
        is0 = 0
     Else
        is0 = xfile%fdate%imn - 1
     End If
  End If
  If (iretro == 1) Then
     Call set_progress_increment ((one - prog)/Real(mya, Kind=rp),.false.)
     If (set_ra_categories(ianal, lb, nsn, is0, ncu(lb:), clim(:,:,lb:)) /= 0) GoTo 2
  End If
!
! Restructure seasonal data
! - compress seasonal data (not previously compressed) -
  Select Case (iseq)
   Case (isq_sn, isq_so)
     Call compress (nt, kcuse(:), nur, &
          ia=irobs(1:mya,:), aa=rfps(1:mya,:,:))
! - reduce seasons and year to one time dimension -
     Allocate (vcopy(mya,nuy,nsn), Stat=ifail)
     If (ifail /= 0) Then
        ifail = 4
        cproc = 'run_analysis'
        carg1 = 'PFV'
        GoTo 1
     End If
     vcopy(:,:,:) = y(1:mya,1:nuy,:)
     Deallocate (y)
     Allocate (y(mya,nss*nuy,1), Stat=ifail)
     kk = 0
     Do k = 1, nuy
        Do isn = 1, nss
           kk = kk + 1
           y(:,kk,1) = vcopy(:,k,isn)
        End Do
     End Do
     Deallocate (vcopy)
  End Select
!
! Permit forecasts if Z data file is set
  If (zfile%cfile%lset) Then
     ifcast = 1
  Else
     ifcast = 0
  End If
!
! Permit scree plots if number of modes is more than one
  ifail = set_scree()
!
#if GUI == 1
! Initialise variable limits
  Call set_ilims (xfield(:)%nv, xfield(:)%nva, iffx, ilimx0, ilimx1, ilimx2)
  Call set_ilims (yfield(1:)%nv, yfield(1:)%nva, iffy, ilimy0, ilimy1, ilimy2)
!
#endif
! End analysis
  irv = iretro
  Call window_update (irv)
! - set length of verification period -
  Select Case (iseq)
   Case (isq_sn, isq_so)
     ntr = nt
  End Select
! - switch off calculation flags -
  Call analysis_flags ('end')
  Call window_update (icalc)
  im_ov = iretro
! - reord change -
  ifail = record_results()
! - end progress meter -
  Call end_progress_meter (.true.,.true.)
  ifail = 0
  Return
!
! Terminate analysis
1 Call cpt_error (Trim(cproc), ifail, .false., &
       c_arg1=Trim(carg1), c_arg2=Trim(carg2), c_arg3=Trim(carg3), i_arg1=iarg1)
2 If (lstop) Then
     mxa = 0
     mya = 0
     Call window_update (mxa)
     Call window_update (mya)
  End If
  Call analysis_flags ('off')
  ifail = close_analysis() ! - clear memory -
#if GUI == 1
  Call clear_window$ (ih_con)
#endif
  Write (Unit=*, Fmt=*)
  Call stop_progress_meter (.true.)
!
  Return
 End Function run_analysis
!
!
!
 Function cpt_reset()
!
! Resets CPT ready for calculation
!
! Function type
  Integer :: cpt_reset
!
! Executable Statements
!
! Reset
  cpt_reset = reset('Resetting', .true.)
!
  Return
 End Function cpt_reset
!
!
!
 Function reset(msg,lflags)
!
! Clears memory for recalculation
!
! Modules
  Use data_cpt_constants, Only: ia_pfv
#if GUI == 1
  Use clrwin$,            Only: clear_window$
#endif
  Use data_numbers,       Only: zero
  Use screen,             Only: window_update
  Use gui,                Only: prompt_yesno
#if GUI == 1
  Use gui,                Only: ih_con, &
                                close_windows
#endif
  Use menus,              Only: im_ov
  Use progress_meter,     Only: set_progress_meter
  Use iofiles,            Only: mxa, mya
  Use settings,           Only: isave
  Use analysis,           Only: ianal, icalc, &
                                analysis_flags
!
! Function type
  Integer :: reset
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: msg ! - message -
!
  Logical, Intent(In) :: lflags ! - reset analysis flags? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim
!
! Executable Statements
!
! Check for reset
  If (icalc == 1) Then
     reset = 1
     If (Len_Trim(msg) > 0) Then
        If (prompt_yesno(.true., .false., msg//' will clear the current results.', &
                         msg2='Do you wish to Continue?') /= 1) Return
     End If
     isave = isave - 2
!
! Close existing windows and clear memory
#if GUI == 1
     reset = close_windows()
     Call clear_window$ (ih_con)
#endif
     reset = reset_mem()
     Select Case (ianal)
      Case Default
        im_ov = 0
      Case (ia_pfv)
        im_ov = 1
     End Select
!
! Reset numbers of used variables
     mxa = 0
     mya = 0
     Call window_update (mxa)
     Call window_update (mya)
!
! Reset progress meter
     Call set_progress_meter (zero)
     Print *, ' '
     If (lflags) Call analysis_flags ('off')
  End If
  reset = 2
!
  Return
 End Function reset
!
!
!
#if GUI == 1
 Function check_reset() &
          Result (check)
!
! Prompts for reset
!
! Modules
  Use screen,   Only: window_update
  Use gui,      Only: ic, ih_nt
  Use settings, Only: nt, nt_old, &
                      record_change
  Use analysis, Only: icalc
!
! Function type
  Integer :: check
!
! Executable Statements
!
! Check for reset
  check = 1
  If (icalc == 1) Then
     If (ic == ih_nt) Then ! - length of training period -
        check = reset('Modifying length of training period', .true.)
     End If
     If (check == 1) Then
        nt = nt_old
        Call window_update (nt)
        check = 2
        Return
     End If
  End If
!
! Update backups
  nt_old = nt
  check = record_change()
!
  Return
 End Function check_reset
#endif
!
!
!
 Function reset_mem()
!
! Clears memory
!
! Modules
  Use analysis,  Only: close_analysis
  Use bootstrap, Only: close_boot, close_pboot
  Use forecast,  Only: close_fcast
!
! Function type
  Integer :: reset_mem
!
! Executable Statements
!
! Clear memory
  reset_mem = close_pboot()
  reset_mem = close_boot()
  reset_mem = close_fcast()
  reset_mem = close_analysis()
!
  reset_mem = 2
!
  Return
 End Function reset_mem
End Module cpt_actions
