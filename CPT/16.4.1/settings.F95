! Author: Simon Mason
Module settings
!
! Modules
  Use data_numbers, Only: rp
  Use menus,        Only: switch
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: change_igauss, change_setting, check_n, check_nt, get_cv, get_maxn, get_retro, init_elevel, record_change,       &
            record_results, reset_settings, reset_limits, reset_limits_nf, reset_limits_nt, set_it1, set_nused, set_options, &
            set_zero, training_refresh
  Public :: proj_read_settings, proj_write_settings
#if GUI == 1
  Public :: switch_elev
#endif
!
! Scalars
!
! Integer scalars
  Integer, Public :: hcw       ! - half cross-validation window -
  Integer, Public :: ibi       ! - number of probability bins identifier -
  Integer, Public :: ibs       ! - binned probabilities for verification scores? -
  Integer, Public :: ieav      ! - ensemble-average flag -
  Integer, Public :: iensc     ! - ensemble counting flag -
  Integer, Public :: iensc_bk  ! - backup ensemble counting flag -
  Integer, Public :: iev       ! - error-variance flag -
  Integer, Public :: iexc      ! - exceedance (0) or non-exceedance (1) probabilities -
  Integer, Public :: ifpr      ! - forecast probability rounding -
  Integer, Public :: igauss    ! - transform to gaussian flag -
  Integer, Public :: igauss_bk ! - backup transform to gaussian flag -
  Integer, Public :: igcms     ! - GCM standardization option -
  Integer, Public :: igcms_bk  ! - backup GCM standardization option -
  Integer, Public :: igood     ! - goodness index flag -
  Integer, Public :: igto      ! - gaussian transformation option -
  Integer, Public :: ig0       ! - lower index limit -
  Integer, Public :: imc       ! - model combination option -
  Integer, Public :: intp      ! - interpolation option -
  Integer, Public :: iodds     ! - odds ratio flag -
  Integer, Public :: iprecf    ! - forecast precision -
  Integer, Public :: iretro    ! - retroactive forecast flag -
  Integer, Public :: isave     ! - change made to project settings flag -
!                                  (0 = no changes;
!                                   1 = changes made to settings;
!                               2,4,6 = results calculated, no changes to settings;
!                               3,5,7 = results calculated, changes to settings;
!                                  -1 = project file needs updating)
  Integer, Public :: isem      ! - sort ensemble members flag -
  Integer, Public :: istd      ! - standardization option -
!                                  0 = no standardization
!                                  1 = anomalies
!                                  2 = standardized anomalies
!                                  3 = SPI
!                                  4 = % average
  Integer, Public :: istdo     ! - original standardization option -
  Integer, Public :: istd_bk   ! - backup standardization option -
  Integer, Public :: isynch    ! - allow synchronous predictors -
  Integer, Public :: iv        ! - current series -
  Integer, Public :: iva       ! - current available series -
  Integer, Public :: ivf       ! - current series by field -
  Integer, Public :: ivf_old   ! - previous series by field -
  Integer, Public :: ivfa      ! - current available series by field -
  Integer, Public :: izero     ! - zero-bound flag -
  Integer, Public :: lcw       ! - length of cross-validation window -
  Integer, Public :: lcw_old   ! - backup length of cross-validation window -
  Integer, Public :: m         ! - number of variables -
  Integer, Public :: mntu      ! - minimum length of training period to use -
  Integer, Public :: mnu       ! - minimum length of training period -
  Integer, Public :: nb        ! - number of probability bins -
  Integer, Public :: ncv       ! - length of cross-validated training period -
  Integer, Public :: nenf      ! - number of ensemble forecasts -
  Integer, Public :: nf        ! - number of forecasts -
  Integer, Public :: nmf       ! - number of forecasts in monthly file -
  Integer, Public :: ntm       ! - number of cases in monthly file to read -
  Integer, Public :: nret      ! - number of retroactive iterations -
  Integer, Public :: nretro    ! - model update interval -
  Integer, Public :: nt        ! - number of training cases -
  Integer, Public :: nt_old    ! - previous number of training cases -
  Integer, Public :: ntr       ! - number of retroactive verification cases -
  Integer, Public :: nt1       ! - initial number of training cases -
  Integer, Public :: nu        ! - number of used training cases -
  Integer, Public :: nur       ! - number of used retroactive verification cases -
  Integer, Public :: nu1       ! - initial number of used training cases -
  Integer, Public :: nyr       ! - number of years -
  Integer, Public :: xfd_old   ! - backup X first year of interest -
  Integer, Public :: yfd_old   ! - backup Y first year of interest -
!
! Real scalars
  Real(Kind=rp), Public :: clf  ! - forecast confidence level -
  Real(Kind=rp), Public :: dofr ! - degrees of freedom for regression -
!
! Logical scalars
  Logical, Public :: lzero ! - apply zero-bound? -
!
! Derived-type scalars
  Type(switch), Public :: sw_elev ! - expertise level -
!
Contains
!
!
 Subroutine reset_settings (ianal)
!
! Resets program settings
!
! Modules
  Use data_numbers,       Only: zero
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv, mnt, nbss
  Use iofiles,            Only: mxa, mya, mza, nx, ny, nz
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis indicator
!
! Executable Statements
!
! Reset numbers of variables
  nx = 0
  mxa = 0
  ny = 0
  mya = 0
  nz = 0
  mza = 0
!
! Reset length of training period and number of forecasts
  nt = 0
  nt_old = 0
  nf = 1
  Select Case (ianal)
   Case (ia_cca, ia_pcr, ia_mlr, ia_gcm)
    mnu = mnt
   Case (ia_pfv)
    mnu = 1
   Case Default
    mnu = 1
  End Select
!
! Retroactive forecast settings
  iretro = 0
  nretro = 1
  nt1 = 0
!
! Set transformation flag
  igauss = 0
  igauss_bk = 0
!
! Set zero-bound flag
  izero = 0
!
! Allow synchronous predictors flag
  isynch = 0
!
! Set average ensemble members flag
  ieav = 0
!
! Set sort ensemble members flag
  isem = 1
!
! Initialise odds flag
  iodds = 0
!
! Initialise probability of exceedance
  iexc = 0
!
! Verification settings
  ibi = 1
  nb = nbss(ibi)
  ibs = 0
!
! Forecast probablity rounding
  ifpr = 1
!
! Regression
  dofr = zero
!
! Backups
  xfd_old = 0
  yfd_old = 0
!
  Return
 End Subroutine reset_settings
!
!
!
 Function record_change()
!
! Records a change to be saved in project file
!
! Function type
  Integer :: record_change
!
! Executable Statements
!
! Record change
  isave = 1
  record_change = 2
!
  Return
 End Function record_change
!
!
!
 Function record_results()
!
! Records calculation of results to be saved in project file
!
! Function type
  Integer :: record_results
!
! Executable Statements
!
! Record calculation of results
  isave = isave + 2
  record_results = 2
!
  Return
 End Function record_results
!
!
!
 Subroutine training_refresh (ncy)
!
! Refreshes training period settings
!
! Modules
  Use screen, Only: window_update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ncy ! - length of climatological period (in years) -
!
! Executable Statements
!
! Refresh file information
  Call window_update (nt)
  Call window_update (ncy)
  Call window_update (nf)
!
  Return
 End Subroutine training_refresh
!
!
!
 Subroutine init_elevel (lcde)
!
! Initialises expertise level
!
! Arguments
  Logical, Intent(Out) :: lcde ! - calculate EOF mode settings? -
!
! Set whether to calculate EOF mode settings
  lcde = (.not.sw_elev%on)
!
  Return
 End Subroutine init_elevel
!
!
!
#if GUI == 1
 Function switch_elev ()
!
! Toggles expertise level
!
! Modules
  Use menus, Only: Assignment(=)
!
! Function Type
  Integer :: switch_elev
!
! Executable Statements
!
! Toggle expertise level
  sw_elev = sw_elev%ion - 1
  switch_elev = 2
!
  Return
 End Function switch_elev
!
!
!
#endif
 Function change_setting(i, ip, mn, mx) &
          Result (change)
!
! Change a setting
!
! Modules
  Use maths,  Only: force_limits
  Use screen, Only: window_update
!
! Function type
  Integer :: change
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: mn ! - minimum value -
  Integer, Intent(In) :: mx ! - maximum value -
!
! Input/output scalars
  Integer, Intent(InOut) :: i  ! - value -
  Integer, Intent(InOut) :: ip ! - previous value -
!
! Locals
!
! Local scalars
  Integer :: iadj ! - adjustment indicator -
!
! Executable Statements
!
! Check for double call
  change = 1
  If (i == ip) Return
!
! Loop if necessary
  Call force_limits (i, mn, mx, iadj)
  If (iadj  /=  0) Call window_update (i)
  ip = i
  change = 2
!
  Return
 End Function change_setting
!
!
!
 Function change_igauss() &
          Result (change)
!
! Toggles Y transformation flag
!
! Function type
  Integer :: change
!
! Executable Statements
!
! Record change
  change = record_change()
!
! Toggle Y transformation flag
#if GUI == 1
  change = 2
#else
  change = 0
#endif
  Select Case (igauss)
   Case (0)
     igauss = igto
   Case Default
     igauss = 0
  End Select
  igauss_bk = igauss
!
  Return
 End Function change_igauss
!
!
!
 Function set_options(iw, i1, i2, i3, i4, i5, i6, i7, i8, r1, r2, r3, l1, &
          lrc)
!
! Sets options
! set_options = 0 if backing up
!               1 if no change is recorded
!               2 if change is identified (and recorded if lrc/ = .false.)
!               3 any changes are cancelled
!
! Function type
  Integer :: set_options
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Integer, Intent(In), Optional :: iw ! - windows response (iw == 1 if action is to be taken)  -
!
  Logical, Intent(In), Optional :: lrc ! - record change? (default = .true.) -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: i1 ! - first integer argument -
  Integer, Intent(InOut), Optional :: i2 ! - second integer argument -
  Integer, Intent(InOut), Optional :: i3 ! - third integer argument -
  Integer, Intent(InOut), Optional :: i4 ! - fourth integer argument -
  Integer, Intent(InOut), Optional :: i5 ! - fifth integer argument -
  Integer, Intent(InOut), Optional :: i6 ! - sixth integer argument -
  Integer, Intent(InOut), Optional :: i7 ! - seventh integer argument -
  Integer, Intent(InOut), Optional :: i8 ! - eighth integer argument -
!
  Real(Kind=rp), Intent(InOut), Optional :: r1 ! - first real argument -
  Real(Kind=rp), Intent(InOut), Optional :: r2 ! - second real argument -
  Real(Kind=rp), Intent(InOut), Optional :: r3 ! - third real argument -
!
  Logical, Intent(InOut), Optional :: l1 ! - first logical argument -
!
! Locals
!
! Local scalars
  Integer, Save :: i1_bk ! - backup first integer argument -
  Integer, Save :: i2_bk ! - backup second integer argument -
  Integer, Save :: i3_bk ! - backup third integer argument -
  Integer, Save :: i4_bk ! - backup fourth integer argument -
  Integer, Save :: i5_bk ! - backup fifth integer argument -
  Integer, Save :: i6_bk ! - backup sixth integer argument -
  Integer, Save :: i7_bk ! - backup seventh integer argument -
  Integer, Save :: i8_bk ! - backup eighth integer argument -
!
  Real(Kind=rp), Save :: r1_bk ! - first real argument -
  Real(Kind=rp), Save :: r2_bk ! - second real argument -
  Real(Kind=rp), Save :: r3_bk ! - third real argument -
!
  Logical, Save :: l1_bk ! - backup first logical argument -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Back up arguments
  If (.not.Present(iw)) Then
     If (Present(i1)) i1_bk = i1
     If (Present(i2)) i2_bk = i2
     If (Present(i3)) i3_bk = i3
     If (Present(i4)) i4_bk = i4
     If (Present(i5)) i5_bk = i5
     If (Present(i6)) i6_bk = i6
     If (Present(i7)) i7_bk = i7
     If (Present(i8)) i8_bk = i8
     If (Present(r1)) r1_bk = r1
     If (Present(r2)) r2_bk = r2
     If (Present(r3)) r3_bk = r3
     If (Present(l1)) l1_bk = l1
     set_options = 0
!
! Check for changes
  Else
     If (iw == 1) Then
        set_options = 1
        If (Present(lrc)) Then
           If (.not.lrc) Then
              set_options = 2
              Return
           End If
        End If
        If (Present(i1)) Then
           If (i1 /= i1_bk) set_options = record_change()
        End If
        If (Present(i2)) Then
           If (i2 /= i2_bk) set_options = record_change()
        End If
        If (Present(i3)) Then
           If (i3 /= i3_bk) set_options = record_change()
        End If
        If (Present(i4)) Then
           If (i4 /= i4_bk) set_options = record_change()
        End If
        If (Present(i5)) Then
           If (i5 /= i5_bk) set_options = record_change()
        End If
        If (Present(i6)) Then
           If (i6 /= i6_bk) set_options = record_change()
        End If
        If (Present(i7)) Then
           If (i7 /= i7_bk) set_options = record_change()
        End If
        If (Present(i8)) Then
           If (i8 /= i8_bk) set_options = record_change()
        End If
        If (Present(r1)) Then
           If (r1 /= r1_bk) set_options = record_change()
        End If
        If (Present(r2)) Then
           If (r2 /= r2_bk) set_options = record_change()
        End If
        If (Present(r3)) Then
           If (r3 /= r3_bk) set_options = record_change()
        End If
        If (Present(l1)) Then
           If (l1 .neqv. l1_bk) set_options = record_change()
        End If
     Else
        If (Present(i1)) i1 = i1_bk
        If (Present(i2)) i2 = i2_bk
        If (Present(i3)) i3 = i3_bk
        If (Present(i4)) i4 = i4_bk
        If (Present(i5)) i5 = i5_bk
        If (Present(i6)) i6 = i6_bk
        If (Present(i7)) i7 = i7_bk
        If (Present(i8)) i8 = i8_bk
        If (Present(r1)) r1 = r1_bk
        If (Present(r2)) r2 = r2_bk
        If (Present(r3)) r3 = r3_bk
        If (Present(l1)) l1 = l1_bk
        set_options = 3
     End If
  End If
!
  Return
 End Function set_options
!
!
!
 Function set_zero()
!
! Toggles zero-bound option, and resets standardization option if necessary
!
! Modules
  Use errors, Only: cpt_error
!
! Function type
  Integer :: set_zero
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Toggle zero-bound option
  izero = 1 - izero
!
! Check for invalid standardization setting
  set_zero = 2
  If (izero == 1) Return
  Select Case (istd)
   Case (3) ! - SPI -
     ifail = 1
   Case (4) ! - % average -
     ifail = 2
   Case Default
     set_zero = record_change()
     Return
  End Select
!
! Errors
  izero = 1
  Call cpt_error ('set_zero', ifail, .false.)
!
  Return
 End Function set_zero
!
!
!
 Function get_cv()
!
! Initialises settings for cross-validated analysis
!
! Function type
  Integer :: get_cv
!
! Executable Statements
!
! Set retroactive flag
  iretro = 0
  nretro = 1
  nt1 = nt
  nret = 1
  ntr = 0
  nur = 0
  get_cv = 0
!
  Return
 End Function get_cv
!
!
!
 Function get_retro()
!
! Initialises settings for retroactive analysis
!
! Modules
  Use gui, Only: init_win, prompt_integer, win_prompt
!
! Function type
  Integer :: get_retro
#if GUI == 0
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
  Character(Len=64) :: cprmpt ! - prompt -
#endif
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod
#if GUI == 0
  Intrinsic Ceiling, Real, Trim
#endif
!
! Executable Statements
!
! Calculate estimate of initial training period length
  If (nt1 > nt+1-mnu .or. nt1 < mnu) Then
     If (Mod(nt,2) == 0) Then
        nt1 = nt/2
     Else
        nt1 = nt/2 + 1
     End If
  End If
!
! Prompt for initial training period length
  iretro = 0
  Call init_win (title='Retroactive forecasts', &
       fclose=close_retro)
#if GUI == 1
  If (nt-mnu > mnu) Then
     Call prompt_integer ('Length of initial training period', nt1, &
          itab=35, ilow=1, ihgh=nt+1-mnu, lfocus = .true.)
  Else
     Call prompt_integer ('Length of initial training period', nt1, &
          itab=35)
  End If
  Call prompt_integer ('Training period update interval', nretro, &
       itab=35, ilow=1)
  get_retro = win_prompt(.true., &
              fok=check_nt1)
#else
  Write (Unit=cprmpt, Fmt='(A,I0,A)') 'Length of initial training period (suggested value ', nt1, ')'
  Call prompt_integer (Trim(cprmpt), nt1, &
       ifail=ifail, itab=35, ilow=1, ihgh=nt+1-mnu, lzero=.true., fcb1=check_nt1)
  If (ifail /= 0) GoTo 1
  Call prompt_integer ('Training period update interval', nretro, &
       ifail=ifail, itab=35, ilow=1, lzero=.true.)
  If (ifail /= 0) GoTo 1
  get_retro = win_prompt(.true.)
#endif
  nret = Ceiling((Real(ntr)-0.5)/Real(nretro))
  get_retro = 0
  Return
!
#if GUI == 0
! Errors
1 iretro = -1
  get_retro = 0
!
  Return
#endif
 End Function get_retro
!
!
!
 Function close_retro()
!
! Halts retroactive analysis
!
! Function type
  Integer :: close_retro
!
! Executable Statements
!
! Set retroactive flag
  If (iretro /= 1) iretro = -1
  close_retro = 0
!
  Return
 End Function close_retro
!
!
!
 Function check_nt1() &
          Result (check)
!
! Checks initial length of training period
!
! Modules
  Use screen, Only: window_update
  Use gui,    Only: set_cursor_waiting
  Use errors, Only: cpt_error
!
! Function type
  Integer :: check
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Ceiling
  Intrinsic Real
!
! Executable Statements
!
! Force nt1 to be at least lcw+mnu-1
1 If (nt1 < mntu) Then
     ifail = 1
     Call cpt_error ('check_nt1', ifail, .false., &
          i_arg1=mnu-1)
#if GUI == 1
     check = 1
#else
     check = ifail
     If (ifail == 1) Return
#endif
     nt1 = mntu
#if GUI == 0
     GoTo 1
#endif
     Call window_update (nt1)
     Call set_cursor_waiting (0)
     iretro = 0
  Else
     ntr = nt - nt1
     nret = Ceiling((Real(ntr)-0.5)/Real(nretro))
     iretro = 1
     check = 0
  End If
!
  Return
 End Function check_nt1
!
!
!
 Function check_nt(nt, it1, n) &
          Result (check)
!
! Checks for availability of data for training period
!
! Function type
  Integer :: check
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt  ! - length of training period -
  Integer, Intent(In) :: it1 ! - first case to use -
  Integer, Intent(In) :: n   ! - number of cases in file -
!
! Executable Statements
!
! Determine whether there are sufficient data
  If (it1+nt-1 <= n) Then
     check = 0
  Else
     check = n + 1 - it1
  End If
!
  Return
 End Function check_nt
!
!
!
 Function check_n(n, lcw, lacw, ianal, iretro, ifail) &
          Result (check)
!
! Checks for sufficient number of cases and that number of cases is sufficiently more than length of cross-validation window.
!
! On exit:
!    check = 0 All checks passed
!    check = 1 n < mnu, and iretro=0
!    check = 2 n < 2*mnu, and iretro=1
!    check = 3 lcw > n-3, and iretro=0; lcw reset to n-3
!    check = 4 lcw > n-mnu*2, and iretro=1; lcw reset to n-mnu*2
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv, mnt
!
! Function type
  Integer :: check
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n      ! - number of cases -
  Integer, Intent(In) :: ianal  ! - analysis indicator -
  Integer, Intent(In) :: iretro ! - retroactive flag -
!
  Logical, Intent(In) :: lacw ! - adjust cross-validation window? -
!
! Input/output scalars
  Integer, Intent(InOut) :: lcw    ! - length of cross-validation window -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Mod
!
! Executable Statements
!
! Set minimum number of cases
  Select Case (ianal)
   Case (ia_cca, ia_pcr, ia_mlr, ia_gcm)
    mnu = mnt
    mntu = mnt + lcw - 1
   Case (ia_pfv)
    mnu = 1
    mntu = 1
  End Select
!
! Determine whether there are sufficient data
  Select Case (iretro)
   Case (0) ! - cross-validation -
     If (n < mnu) Then
        ifail = mnu
        check = 1
        Return
     End If
   Case (1) ! - retroactive -
     If (n < 2*mnu-1) Then
        ifail = 2*mnu - 1
        check = 2
        Return
     End If
  End Select
!
! Force lcw to leave at least 3
  Select Case (iretro)
   Case (0) ! - cross-validation -
     If (lcw > n-(mnu-1)) Then
        If (lacw) Then
           lcw = n - (mnu - 1)
           If (Mod(lcw,2) == 0) lcw = lcw - 1
        End If
        ifail = 0
        check = record_change()
        check = 3
        Return
     End If
   Case (1) ! - retroactive -
     If (lcw > n-2*(mnu-1)) Then
        If (lacw) Then
           lcw = n - 2*(mnu - 1)
           If (Mod(lcw,2) == 0) lcw = lcw - 1
        End If
        ifail = 2*mnu - 1
        check = record_change()
        check = 4
        Return
     End If
  End Select
  check = 0
!
  Return
 End Function check_n
!
!
!
 Subroutine get_maxn (icsv, lsetx, lsety, lprst, iseqx, iseqy, iseqm, lensn, lsn0, lsnn, fdatex, fdatey, bdatem, &
            gdatex, gdatey, gdatem, nt, ntm, &
            ih, mn)
!
! Calculates number of remaining pairs of data given start dates
!
! Modules
#if GUI == 1
  Use clrwin$,             Only: set_integer_limits$
  Use data_numbers,        Only: ihuge
#endif
  Use data_gui_constants,  Only: i_handle
  Use data_time_constants, Only: isq_sn, isq_so, isq_yr, nmn
  Use time,                Only: pdate, &
                                 iseq, &
                                 date_diff
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icsv  ! - calculate seasonal values flag -
  Integer, Intent(In) :: iseqx ! - X-file time-sequencing -
  Integer, Intent(In) :: iseqy ! - Y-file time-sequencing -
  Integer, Intent(In) :: iseqm ! - M-file time-sequencing -
  Integer, Intent(In) :: lensn ! - length of season -
  Integer, Intent(In) :: lsn0  ! - length of persistence seasons -
  Integer, Intent(In) :: lsnn  ! - length of last season -
!
  Logical, Intent(In) :: lprst ! - persistence component? -
  Logical, Intent(In) :: lsetx ! - X-file set? -
  Logical, Intent(In) :: lsety ! - Y-file set? -
!
  Type(pdate), Intent(In) :: fdatex ! - start-at date in X-file -
  Type(pdate), Intent(In) :: fdatey ! - start-at date in Y-file -
  Type(pdate), Intent(In) :: bdatem ! - beginning date in M-file -
  Type(pdate), Intent(In) :: gdatex ! - last date in X-file -
  Type(pdate), Intent(In) :: gdatey ! - last date in Y-file -
  Type(pdate), Intent(In) :: gdatem ! - last date in M-file -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: mn ! - minimum number -
!
  Integer(Kind=i_handle), Optional :: ih ! - handle for prompt -
!
! Input/output scalars
  Integer, Intent(InOut) :: nt  ! - number of cases available -
  Integer, Intent(InOut) :: ntm ! - number of monthly cases available -
!
! Input scalars
!
! Locals
!
! Local scalars
  Integer :: lssn ! - length of season -
  Integer :: ntx  ! - number of remaining cases in X-file -
  Integer :: n    ! - number of cases available -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Floor, Min, Mod, Real
#if GUI == 1
  Intrinsic Present
#endif
!
! Executable Statements
!
! Calculate remaining number of cases given start dates 
  If (lsetx) Then
     ntx = 1 + date_diff(fdatex, gdatex, iseqx)
  Else
     ntx = 0
  End If
  Select Case (icsv)
   Case (0) ! - no seasonal calculations -
     If (lsety) Then
        n = 1 + date_diff(fdatey, gdatey, iseqy)
        If (lsetx) n = Min(ntx, n)
     Else
        n = ntx
     End If
     ntm = n
   Case Default
     If (.not.lprst) Then
        lssn = lensn
     Else
        If (iseq /= isq_so) Then
           lssn = lensn - lsn0
        Else
           lssn = lsnn
        End If
     End If
     ntm = 1 + date_diff(bdatem, gdatem, iseqm)
     Select Case (iseq)
      Case (isq_yr) ! - yearly sequencing -
        If (ntm > lssn) Then
           n = 1 + Floor(Real(ntm - lssn)/Real(nmn))
        Else
           n = 0
        End If
        n = Min(ntx, n)
      Case (isq_sn) ! - seasonal sequencing -
        n = 1 + date_diff(fdatey, gdatey, iseqy)
        n = Min(ntx, n)
      Case (isq_so) ! - seasonally overlapping sequencing -
        If (ntm > lssn) Then
           n = ntm - Mod(lssn, nmn)
        Else
           n = ntm
        End If
        n = Min(ntx, n)
      Case Default ! - monthly sequencing -
        If (ntx > ntm) Then
           n = ntm - Mod(lssn, nmn)
        Else
           n = Min(ntx + Mod(lssn, nmn), ntm)
        End If
     End Select
  End Select
!
#if GUI == 1
! Set limits
  If (Present(ih)) Then
     If (lsetx .or. lsety) Then
        If (n > 0) Then
           Call set_integer_limits$ (ih, mn, n)
        Else
           Call set_integer_limits$ (ih, 0, 0)
        End If
     Else
        Call set_integer_limits$ (ih, 0, ihuge)
     End If
  End If
!
#endif
! Reset if smaller than current nt
  If (nt /= 0) Then
     nt = Min(n, nt)
  Else
     nt = n
  End If
!
  Return
 End Subroutine get_maxn
!
!
!
 Subroutine reset_limits (icsv, lprst, lensn, lsn, lb, nsn)
!
! Resets training-period and forecast-period limits
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icsv  ! - calculate seasonal values flag -
  Integer, Intent(In) :: lb    ! - lower bound -
  Integer, Intent(In) :: lensn ! - length of season -
  Integer, Intent(In) :: nsn   ! - number of seasons -
!
  Logical, Intent(In) :: lprst ! - persistence component? -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsn ! - lengths of seasons -
!
! Executable Statements
!
! Reset training-period limits
  Call reset_limits_nt (icsv, lprst, lensn, lsn, lb, nsn)
!
! Reset forecast-period limits
  Call reset_limits_nf (icsv, lprst, lsn, lb, nsn)
!
  Return
 End Subroutine reset_limits
!
!
!
 Subroutine reset_limits_nf (icsv, lprst, lsn, lb, nsn)
!
! Resets forecast-period limits
!
! Modules
#if GUI == 1
  Use clrwin$,      Only: set_integer_limits$
  Use data_numbers, Only: ihuge
#endif
  Use time,         Only: Operator(+), &
                          date_diff
  Use screen,       Only: window_update
  Use gui,          Only: ih_nf
  Use iofiles,      Only: mfile, xfile, yfile, zfile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icsv ! - calculate seasonal values flag -
  Integer, Intent(In) :: lb   ! - lower bound -
  Integer, Intent(In) :: nsn  ! - number of seasons -
!
  Logical, Intent(In) :: lprst ! - persistence component? -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsn ! - lengths of seasons -
!
! Locals
!
! Local scalars
  Integer :: nfm ! - number of cases in monthly file to read -
  Integer :: d   ! - difference betwen X and Z start-at dates -
!
! Executable Statements
!
! Reset training-period limits
  If (zfile%cfile%lset) Then
     d = date_diff(xfile%fdate, zfile%fdate, xfile%iseq)
     If (.not.lprst) Then
        Call get_maxn (   0, zfile%cfile%lset,          .false.,   lprst, xfile%iseq, yfile%iseq, zfile%iseq,      0,       0,   &
                    0, zfile%fdate, yfile%fdate+d,   zfile%bdate, zfile%prdn%sdate, yfile%prdn%sdate, zfile%prdn%sdate, nf, nfm, &
             ih=ih_nf, mn=1)
     Else
        Call get_maxn (icsv, zfile%cfile%lset, yfile%cfile%lset, .false., xfile%iseq, yfile%iseq, mfile%iseq, lsn(0), lsn(lb),   &
             lsn(nsn), zfile%fdate, yfile%fdate+d, mfile%fdate+d, zfile%prdn%sdate, yfile%prdn%sdate, mfile%prdn%sdate, nf, nfm, &
             ih=ih_nf, mn=1)
     End If
  Else
#if GUI == 1
     Call set_integer_limits$ (ih_nf, 0, ihuge)
#endif
     nf = 0
  End If
  Call window_update (nf)
!
  Return
 End Subroutine reset_limits_nf
!
!
!
 Subroutine reset_limits_nt (icsv, lprst, lensn, lsn, lb, nsn)
!
! Resets training-period limits
!
! Modules
#if GUI == 1
  Use clrwin$,      Only: set_integer_limits$
  Use data_numbers, Only: ihuge
#endif
  Use screen,       Only: window_update
  Use gui,          Only: ih_nt
  Use iofiles,      Only: mfile, xfile, yfile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icsv  ! - calculate seasonal values flag -
  Integer, Intent(In) :: lb    ! - lower bound -
  Integer, Intent(In) :: lensn ! - length of season -
  Integer, Intent(In) :: nsn   ! - number of seasons -
!
  Logical, Intent(In) :: lprst ! - persistence component? -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsn ! - lengths of seasons -
!
! Executable Statements
!
! Reset training-period limits
  If (xfile%cfile%lset .or. yfile%cfile%lset) Then
     Call get_maxn (icsv, xfile%cfile%lset, yfile%cfile%lset, lprst, xfile%iseq, yfile%iseq, mfile%iseq, lensn, lsn(lb), lsn(nsn), &
          xfile%fdate, yfile%fdate, mfile%bdate, xfile%prdn%sdate, yfile%prdn%sdate, mfile%prdn%sdate, nt, ntm, &
          ih=ih_nt, mn=mnu)
  Else
#if GUI == 1
     Call set_integer_limits$ (ih_nt, 0, ihuge)
#endif
     nt = 0
  End If
  Call window_update (nt)
!
  Return
 End Subroutine reset_limits_nt
!
!
!
 Function set_nused(lset, kuse)
!
! Adjusts indices based on data availability
!
! Function type
  Integer :: set_nused
!
! Arguments
!
! Input arrays
  Logical, Intent(In) :: lset ! - set nur? -
!
! Input arrays
  Logical, Dimension(:), Intent(In) :: kuse ! - used case? -
!
! Locals
!
! Local scalars
  Integer :: k  ! - case index -
  Integer :: kk ! - available case index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count
!
! Executable Statements
!
! Check retroactive settings
  If (iretro == 1) Then
     If (mnu > 1) Then
        mntu = lcw + mnu - 1
     Else ! - PFV -
        nt1 = nt
        mntu = 1
     End If
     nu1 = nt1 - Count(.not.kuse(1:nt1))
     If (nu1 >= mntu) Then
        set_nused = 0
     Else
        nu1 = mntu
        kk = 0
        Do k = 1, nt
           If (kuse(k)) Then
              kk = kk + 1
              If (kk == nu1) Then
                 nt1 = k
                 Exit
              End If
           End If
        End Do
        Select Case (mntu)
         Case (2:)
           set_nused = 1
         Case (1)
           set_nused = 2
        End Select
     End If
     If (lset) nur = nu - nu1
  Else
     nu1 = 0
     set_nused = 0
  End If
!
  Return
 End Function set_nused
!
!
!
 Function set_it1(isq, sdate, fdate) &
          Result (it1)
!
! Checks that first date is available.
!
! Modules
  Use time, Only: pdate, &
                  date_diff
!
! Function type
  Integer :: it1
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - sequence -
!
  Type(pdate), Intent(In) :: sdate ! - first date available -
  Type(pdate), Intent(In) :: fdate ! - first date of interest -
!
! Executable Statements
!
! Force first date of interest not to be before first date available
  it1 = 1 + date_diff(sdate, fdate, isq)
!
  Return
 End Function set_it1
!
!
!
 Subroutine proj_read_settings (iin, icalc, pver, ifail)
!
! Reads program settings from project file
!
! Modules
  Use data_numbers,       Only: sp, tolv
  Use data_cpt_constants, Only: nbss
  Use time,               Only: itime, lead
  Use iofiles,            Only: mxa, mya, mza, nx, ny, nz
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin   ! - input file unit number -
  Integer, Intent(In) :: icalc ! - results available flag -
!
  Real(Kind=sp), Intent(In) :: pver ! - project version number -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: lxt  ! - minimum of nxa and nu (not needed) -
  Integer :: lxyt ! - Min(Max(nxa,nya),nu) (not needed) -
  Integer :: lyt  ! - minimum of nya and nu (not needed) -
  Integer :: ntc  ! - length of climatological period with length of training period (not needed) -
!
! Executable Statements
!
! Read program settings
  istd_bk = istd
  If (pver > 16.02_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi, ibs, ieav, iensc, iev, iexc, ifpr, igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo, isynch, itime, izero, lcw,      nenf, nf, nretro, nt, &
                                   nt_old, nt1, xfd_old, yfd_old
  Else If (pver > 16.00_sp) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi, ibs, ieav, iensc, iev,       ifpr, igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo, isynch, itime, izero, lcw,      nenf, nf, nretro, nt, &
                                   nt_old, nt1, xfd_old, yfd_old
  Else If (pver > 15.04_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi, ibs,       iensc, iev,       ifpr, igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo, isynch, itime, izero, lcw,      nenf, nf, nretro, nt, &
                                   nt_old, nt1, xfd_old, yfd_old
  Else If (pver > 15.03_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi,            iensc, iev,       ifpr, igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo, isynch, itime, izero, lcw,      nenf, nf, nretro, nt, &
                                   nt_old, nt1, xfd_old, yfd_old
  Else If (pver > 14.06_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi,            iensc, iev,       ifpr, igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo, isynch,        izero, lcw, ntc, nenf, nf, nretro, nt, &
                                   nt_old, nt1, nx, ny, nz, xfd_old, yfd_old
  Else If (pver > 14.05_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi,            iensc, iev,             igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo, isynch,        izero, lcw, ntc, nenf, nf, nretro, nt, &
                                   nt_old, nt1, nx, ny, nz, xfd_old, yfd_old
  Else If (pver > 14.04_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) clf, ibi,            iensc, iev,             igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                                   intp, iodds, iprecf, isem, istd, istdo,                izero, lcw, ntc, nenf, nf, nretro, nt, &
                                   nt_old, nt1, nx, ny, nz, xfd_old, yfd_old
  Else
     Read (Unit=iin, Err=1, End=2) clf, ibi,            iensc, iev,             igauss, igauss_bk, igcms, igcms_bk, igood,       &
                                   intp, iodds, iprecf, isem, istd, istdo,                izero, lcw, ntc, nenf, nf, nretro, nt, &
                                   nt_old, nt1, nx, ny, nz, xfd_old, yfd_old
  End If
  If (pver < 14.05_sp-tolv) isynch = 0
  If (pver < 14.06_sp-tolv) ifpr = 1
  If (pver < 15.04_sp-tolv) ibs = 0
  If (pver < 16.00_sp)      ieav = 0
  If (pver < 16.02_sp-tolv) iexc = 0
  nb = nbss(ibi)
!
! Read results settings if saved
  If (icalc == 1) Then
     lzero = (igauss == 0 .and. izero == 1)
     If (pver > 16.01_sp-tolv) Then
        Read (Unit=iin, Err=1, End=2) dofr, iretro, iv, iva, ivf, ivf_old, ivfa, lead, &
                                      nret, ntr, nu, nur, nu1
     Else If (pver > 16.00_sp) Then
        Read (Unit=iin, Err=1, End=2) dofr, iretro, iv, iva, ivf, ivf_old, ivfa, lead, lxt, lxyt, lyt, &
                                      nret, ntr, nu, nur, nu1
     Else
        Read (Unit=iin, Err=1, End=2) dofr, iretro, iv, iva, ivf, ivf_old, ivfa, lead, lxt, lxyt, lyt, &
                                      mxa, mya, mza, nret, ntr, nu, nur, nu1
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
 End Subroutine proj_read_settings
!
!
!
 Subroutine proj_write_settings (iout, icalc, ifail)
!
! Writes program settings to project file
!
! Modules
  Use time, Only: itime, lead
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout  ! - output file unit number -
  Integer, Intent(In) :: icalc ! - results available flag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write program settings
  Write (Unit=iout, Err=1) clf, ibi, ibs, ieav, iensc, iev, iexc, ifpr, igauss, igauss_bk, igcms, igcms_bk, igood, igto, &
                           intp, iodds, iprecf, isem, istd, istdo, isynch, itime, izero, lcw,      nenf, nf, nretro, nt, &
                           nt_old, nt1, xfd_old, yfd_old
!
! Write results settings if saved
  If (icalc == 1) Write (Unit=iout, Err=1) dofr, iretro, iv, iva, ivf, ivf_old, ivfa, lead, nret, ntr, nu, nur, nu1
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
!
 End Subroutine proj_write_settings
End Module settings
