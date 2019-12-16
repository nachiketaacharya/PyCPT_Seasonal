! Author: Simon Mason
Module cpt_options
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: change_language, get_boot_opts, get_cca_opts, get_climatology_opts, get_cv_opts, get_eofx_opts, get_eofy_opts, &
            get_fcast_opts, get_gcm_opts, get_goodness_opts, get_missing_opts, get_startat_opts, get_tailoring_opts,       &
            get_targetssn_opts, get_time_opts, get_transform_opts, get_verif_opts
#if GUI == 1
  Public :: get_console_opts, get_gfont_opts, get_gsize_opts, get_marker_opts, get_text_opts
#else
  Public :: get_error_opts, get_nt
#endif
!
Contains
!
!
 Function get_boot_opts() &
          Result (opts)
!
! Prompts for bootstrap settings
!
! Modules
  Use data_numbers, Only: zero, oneh
  Use labels,       Only: l_resamples, l_scoresp
  Use gui,          Only: iw, &
                          add_check_button, box_close, box_open, boxes_open, init_win, prompt_integer, prompt_real, &
                          win_prompt
#if GUI == 1
  Use gui,          Only: gui_creturn
#endif
  Use menus,        Only: im_ov
  Use settings,     Only: set_options
  Use analysis,     Only: ihind
  Use bootstrap,    Only: clb, ibcl, ipval, nboot, nperm
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup resampling settings
  opts = set_options(i1=nboot, i2=nperm, r1=clb)
!
! Prompt for resampling settings
  Call init_win (title=Trim(l_resamples%c), &
       lresize=.true.)
! - probabilistic score confidence intervals / skill map p-values -
  Call boxes_open (1, 1 + im_ov + 2*ihind, .false., .true.)
  If (im_ov == 1) Then
     Call box_open (.true., .false., &
          title=Trim(l_scoresp%c))
#if GUI == 1
     Call gui_creturn (.true.)
#endif
     Call add_check_button ('Calculate confidence intervals', ibcl, .true.)
     Call box_close (.true.)
     Call box_close (.false.)
  End If
! - skill map p-values -
  If (ihind == 1) Then
     Call box_open (.true., .false., &
          title='Skill maps')
#if GUI == 1
     Call gui_creturn (.true.)
#endif
     Call add_check_button ('Calculate p-values', ipval, .true.)
     Call box_close (.true.)
     Call box_close (.false.)
! - permutation settings -
     Call box_open (.true., .false., &
          title='Permutations')
     Call prompt_integer ('Number of permutations', nperm, &
          itab=32, ilow=100, lfocus=.true.)
     Call box_close (.true.)
     Call box_close (.false.)
  End If
! - bootstrap settings -
  Call box_open (.true., .false., &
       title='Bootstrapping')
  Call prompt_integer ('Number of bootstrap samples', nboot, &
       itab=32, ilow=100, lfocus=.true.)
  Call prompt_real ('Confidence level (%)', clb, &
       itab=32, elow=zero, ehgh=oneh, lfocus=.true.)
  Call box_close (.true.)
  Call box_close (.false.)
! - prompt to confirm -
  iw = win_prompt(.true.)
!
! Restore old settings If cancelled
  opts = set_options(iw=iw, i1=nboot, i2=nperm, r1=clb)
  opts = 2
!
  Return
 End Function get_boot_opts
!
!
!
 Function get_cca_opts() &
          Result (opts)
!
! Sets CCA options
!
! Modules
  Use pcs,         Only: lcde, &
                         get_cca_settings
  Use cpt_actions, Only: reset
!
! Function type
  Integer :: opts
!
! Executable Statements
!
! Check for completed calculations
  opts = 2
  If (reset('Changing the CCA settings', .true.) == 1) Return
!
! Prompt for CCA options
  lcde = .false.
  opts = get_cca_settings()
!
  Return
 End Function get_cca_opts
!
!
!
 Function get_climatology_opts() &
          Result (opts)
!
! Redefines climatological period for forecasts
!
! Modules
  Use data_cpt_constants,  Only: ia_pfv
  Use data_io_constants,   Only: ioutstd
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr, nmn
  Use labels,              Only: l_climper
  Use time,                Only: pdate, &
                                 Operator(==), Operator(+), Operator(-), Operator(<), &
                                 iseq, nss, &
                                 add_to_month, reset_iseq
  Use gui,                 Only: iw, &
                                 add_check_button, box_close, box_open, ganged_rbs, init_win, print_text, prompt_pdate, win_prompt
  Use arrays,              Only: get_flag, set_flag
  Use iofiles,             Only: xfile, yfile
  Use settings,            Only: nt, &
                                 record_change, set_options
  Use fields,              Only: yfield
  Use season,              Only: icsv, lprst, lsn, nsn
  Use climate,             Only: Operator(/=), &
                                 climate_bak, climate_per, icex, ieclim, &
                                 climatology, label_climate, print_climatology, set_climate, valid_climate
  Use analysis,            Only: ianal, icalc
  Use forecast,            Only: ifc
  Use cpt_actions,         Only: reset
!
! Function type
  Integer :: opts
!
! Locals
!
! Local scalars
  Integer :: iy2 ! - year in two-digit format -
!
  Character(Len=64) :: ctxt ! - text message -
!
! Local arrays
  Integer, Dimension(0:2) :: iscs ! - seasonal calculation options -
!
  Type(pdate), Dimension(2) :: dlimits ! - date limits -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim, Mod, Trim
!
! Executable Statements
!
! Set default climatology If necesssary
  opts = set_options(i1=ieclim, i2=icsv)
  climate_bak = climate_per
  If (lprst) Call reset_iseq (isq_new=isq_mn)
  If (ianal /= ia_pfv) Then
     If (climate_per%d1 == 0) climate_per%d1 = yfile%fdate - lsn(0)
     If (climate_per%d2 == 0) climate_per%d2 = yfile%fdate + ((nt - 1) - lsn(0))
  Else
     If (ieclim == 1 .and. yfile%prd1%sdate < xfile%prd1%sdate) Then
        If (climate_per%d1 == 0) climate_per%d1 = yfile%prd1%sdate
        If (climate_per%d2 == 0) climate_per%d2 = yfile%fdate - nss
     Else
        If (climate_per%d1 == 0) climate_per%d1 = yfile%prd1%sdate
        If (climate_per%d2 == 0) climate_per%d2 = yfile%prdn%sdate
     End If
  End If
  Call set_flag (icsv+1, iscs(0:))
!
! Prompt for climatological period
  Call init_win (title=Trim(l_climper%c))
  Call print_text ('Please specify climatological period:', &
       lbold=.true.)
  If (yfile%prd1%edate%iyr /= yfile%prd1%sdate%iyr) Then
     Call print_text ('The year should be for the first month of the season.')
     iy2 = Mod(yfile%fdate%iyr+1, 100)
     Write (Unit=ctxt, Fmt='(3A,I4,A,I2.2,A,I4)') &
          'For example, for ', Trim(yfile%cssn),' ', yfile%fdate%iyr,'/', iy2,' enter ', yfile%fdate%iyr
     Call print_text (Trim(ctxt))
  End If
  dlimits(1) = yfile%prd1%sdate - lsn(0)
  dlimits(2) = yfile%prdn%sdate
  Call print_text ('(Data limits are '//Trim(yfile%cprd1)//' to '//Trim(yfile%cprdn)//')')
  Call print_text (' ')
  Select Case (yfile%iseq)
   Case (isq_yr, isq_mn) ! - yearly, and monthly sequences -
     Call prompt_pdate ('First year', isq_yr, climate_per%d1, &
               itab=20, dlimits=dlimits, lfocus=.true.)
#if GUI == 0
     dlimits(1) = climate_per%d1
#endif
     Call prompt_pdate ('Last year', isq_yr, climate_per%d2, &
               itab=20, dlimits=dlimits, lfocus=.true.)
   Case (isq_sn, isq_so) ! - seasonal sequences -
     Call prompt_pdate ('First year ('//Trim(yfield(1)%cssn)//')', isq_yr, climate_per%d1, &
               itab=20, dlimits=dlimits, lfocus=.true.)
#if GUI == 0
     dlimits(1) = climate_per%d1
#endif
     Call prompt_pdate ('Last year ('//Trim(yfield(nss)%cssn)//')', isq_yr, climate_per%d2, &
               itab=20, dlimits=dlimits, lfocus=.true.)
   Case (1:) ! - daily, pentad, weekly, dekadal sequences -
     Call prompt_pdate ('First date', iseq, climate_per%d1, &
               itab=20, dlimits=dlimits, lfocus=.true.)
#if GUI == 0
     dlimits(1) = climate_per%d1
#endif
     Call prompt_pdate ('Last date', iseq, climate_per%d2, &
               itab=20, dlimits=dlimits, lfocus=.true.)
  End Select
! - prompt for extension beyond training period -
  Call print_text (' ')
  Call add_check_button ('Climatological period can extend beyond training period?', ieclim, .false.)
! - prompt for seasonal aggregation -
  If (icsv > 0 .and. icalc == 0) Then
     Call print_text (' ')
     Call print_text (' ')
     Call box_open (.true., .false., &
          title='Seasonal calculations')
     Call ganged_rbs (2, iscs(1:), &
          crbs=(/'Seasonal averages', 'Seasonal totals  '/))
     Call box_close (.false.)
  End If
  If (lprst) Call reset_iseq ()
  iw = win_prompt(.true., &
       fok=valid_climate)
!
! Reset climatological period
  If (iw == 1) icsv = get_flag(iscs(0:)) - 1
  opts = set_options(iw=iw, i1=ieclim, i2=icsv)
  If (iw == 1) Then
     If (climate_per /= climate_bak) Then
        opts = record_change() ! - record change for project file -
        Select Case (yfile%iseq)
         Case (isq_yr)
           climate_per%d1%imn = add_to_month(yfile%prd1%sdate%imn, -lsn(0))
           climate_per%d1%idy = 0
           climate_per%d2%imn = climate_per%d1%imn
           climate_per%d2%idy = 0
         Case (isq_sn)
           climate_per%d1%idy = 0
           climate_per%d2%imn = yfield(nsn)%tprd%sdate%imn
           climate_per%d2%idy = 0
         Case (isq_so)
           climate_per%d1%idy = 0
           climate_per%d2%imn = add_to_month(climate_per%d2%imn, nmn-1)
           climate_per%d2%idy = 0
        End Select
!
! Recalculate climatology
        If (icalc == 0) Then
           Call set_climate (icalc, lsn(0))
        Else
           If (icex == 0) Then
              opts = climatology(icalc, lsn(0), nsn, lprst)
              Select Case (opts)
               Case (0)
                 If (ifc == 2) ifc = 1
               Case (-1)
                 If (reset('Extending the climatological period outside the range of the analysed data', .true.) == 1) Then
                    climate_per = climate_bak
                    iw = 0
                 End If
              End Select
           Else
              If (reset('Adjusting the climatological period', .true.) == 1) Then
                 climate_per = climate_bak
                 iw = 0
              Else
                 Call set_climate (icalc, lsn(0))
                 Call label_climate (lsn(0))
              End If
           End If
        End If
     End If
!
! Restore old settings
  Else
     climate_per = climate_bak
  End If
  If (iw /= 1) Then
     If (Len_Trim(climate_per%clim1) > 0) Call print_climatology (.true., &
                                               iout=ioutstd)
  End If
  opts = 2
!
  Return
 End Function get_climatology_opts
!
!
!
#if GUI == 1
 Function get_console_opts() &
          Result (opts)
!
! Prompts for console size adjustments
!
! Modules
  Use labels,   Only: c_def=>c_lab, &
                      l_defaults, l_console, &
                      make_label
  Use gui,      Only: icon_h, icon_v, iw, &
                      init_win, print_text, prompt_integer, win_prompt
  Use settings, Only: set_options
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup current setting
  opts = set_options(i1=icon_h, i2=icon_v)
!
! Prompt for console scaling
  Call init_win (title=Trim(l_console%c))
  Call print_text ('Set console size adjustments:', &
       lbold=.true.)
  Call print_text (' ')
  Call prompt_integer ('Height adjustment', icon_v, &
       itab=25, lbold=.false.)
  Call prompt_integer ('Width adjustment', icon_h, &
       itab=25, lbold=.false.)
  c_def = make_label('T', l_defaults, &
          lak=.true.)
  iw = win_prompt(.true., &
       cextra=Trim(c_def), iset1=icon_h, ival1=0, iset2=icon_v, ival2=0)
!
! Confirm changes
  opts = set_options(iw=iw, i1=icon_h, i2=icon_v)
  opts = 1
!
  Return
 End Function get_console_opts
!
!
!
#endif
 Function get_cv_opts() &
          Result (opts)
!
! Prompts for length of cross-validation period
!
! Modules
  Use labels,      Only: cg_lcv, &
                         l_cvopts
  Use gui,         Only: iw, &
                         init_win, prompt_integer, win_prompt
  Use settings,    Only: lcw, lcw_old, &
                         set_options
  Use cpt_actions, Only: reset
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Check for reset
  opts = 2
  If (reset('Modifying length of cross-validation window', .true.) == 1) Return
!
! Get training data settings
  lcw_old = lcw
  opts = set_options(i1=lcw)
  Call init_win (title=Trim(l_cvopts%c))
  Call prompt_integer (Trim(cg_lcv)//' (must be odd)', lcw, &
       Ifly=2, ilow=1, fcb2=force_odd_lcw)
! - prompt to confirm -
  iw = win_prompt(.true.)
  opts = set_options(iw=iw, i1=lcw)
!
  Return
 End Function get_cv_opts
!
!
!
 Function force_odd_lcw()
!
! Forces length of cross-validation window to be odd
!
! Modules
  Use maths,    Only: force_odd
  Use screen,   Only: window_update
#if GUI == 0
  Use gui,      Only: print_warning
#endif
  Use settings, Only: lcw, lcw_old
!
! Function type
  Integer :: force_odd_lcw
!
#if GUI == 0
!
! Locals
!
! Local scalars
  Integer :: lcw_prv ! - previous value -
!
  Character(Len=40) :: cmsg ! - warning message -
!
#endif
! Executable Statements
!
! Force length of cross-validation window to be odd
#if GUI == 0
  lcw_prv = lcw
#endif
  Call force_odd (lcw, lcw_old)
  Call window_update (lcw)
  lcw_old = lcw
#if GUI == 0
  If (lcw /= lcw_prv) Then
     Write (Unit=cmsg, Fmt=*) 'Length has been changed to:', lcw
     Call print_warning ('Length of cross-validation window must be odd.', &
          msg2=cmsg)
  End If
#endif
  force_odd_lcw = 2
!
  Return
 End Function force_odd_lcw
!
!
!
 Function get_eofx_opts() &
          Result (opts)
!
! Calls function to set X EOF options
!
! Modules
  Use pcs, Only: ieofx, ieofx_bk, &
                 get_eofx_settings
!
! Function type
  Integer :: opts
!
! Executable Statements
!
! Prompt for EOF options
  opts = eof_opts(get_eofx_settings)
  ieofx_bk = ieofx
!
  Return
 End Function get_eofx_opts
!
!
!
 Function get_eofy_opts() &
          Result (opts)
!
! Calls function to set Y EOF options
!
! Modules
  Use pcs, Only: get_eofy_settings
!
! Function type
  Integer :: opts
!
! Executable Statements
!
! Prompt for EOF options
  opts = eof_opts(get_eofy_settings)
!
  Return
 End Function get_eofy_opts
!
!
!
 Function eof_opts(feof) &
          Result (opts)
!
! Sets EOF options
!
! Modules
  Use pcs,         Only: lcde
  Use cpt_actions, Only: reset
!
! Function type
  Integer :: opts
!
! Arguments
!
! Procedure arguments
  Integer, External :: feof ! - EOF settings function -
!
! Executable Statements
!
! Check for completed calculations
  opts = 2
  If (reset('Changing the modes settings', .true.) == 1) Return
  lcde = .false.
  opts = feof()
!
  Return
 End Function eof_opts
!
!
!
#if GUI == 0
 Function get_error_opts() &
          Result (opts)
!
! Prompts for user-error options
!
! Modules
  Use data_cpt_constants, Only: nueo
  Use gui,                Only: iprompt
  Use labels,             Only: cg_errhand_t, &
                                ls_uerropts
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for user-error handling
  iprompt = iprompt + 1
  Call get_option_list (Trim(cg_errhand_t), 'Select option for how to handle user input errors:', nueo, ls_uerropts(:)%c, iprompt)
  iprompt = iprompt - 1
  opts = 2
!
  Return
 End Function get_error_opts
!
!
!
#endif
 Function get_fcast_opts() &
          Result (opts)
!
! Prompts for forecast settings
!
! Modules
  Use data_numbers,       Only: zero, oneh
  Use data_cpt_constants, Only: nev
  Use labels,             Only: cg_pinterval_t, &
                                l_fcastsets
  Use gui,                Only: iw, &
                                add_check_button, box_close, box_open, ganged_rbs, init_win, print_subtitle, prompt_integer, &
                                prompt_real, win_prompt
  Use arrays,             Only: get_flag, set_flag
  Use settings,           Only: clf, iev, ifpr, iodds, iprec, iretro, nenf, &
                                set_options
  Use forecast,           Only: ifc
!
! Function type
  Integer :: opts
!
! Locals
!
! Local arrays
  Integer :: ievo(nev) ! - error-variance options flags -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup forecast settings
  opts = set_options(i1=iev, i2=iodds, i3=nenf, i4=iprec, i5=ifpr, r1=clf)
!
! Identify error variance option
! - reset invalid option -
  If (iretro == 0 .and. iev == 2) iev = 1
  Call set_flag (iev, ievo)
!
! Prompt for forecast settings
  Call init_win (title=Trim(l_fcastsets%c), &
       lresize=.true.)
! - prediction confidence interval -
  Call box_open (.false., .false.)
  Call box_open (.true., .false., &
       Trim(cg_pinterval_t))
  Call print_subtitle ('Prediction interval width')
  Call prompt_real ('Confidence level (%)', clf, &
       itab=32, elow=zero, ehgh=oneh, lfocus=.true., lbold=.false.)
! - error variance -
  Call print_subtitle ('Error variance')
  Select Case (iretro)
   Case (0)
     Call ganged_rbs (2, ievo, &
          crbs=(/'Cross-validated error variance', 'Fitted error variance         '/))
   Case (1)
     Call ganged_rbs (3, ievo, &
          crbs=(/'Cross-validated error variance', 'Retroactive error variance    ', 'Fitted error variance         '/))
  End Select
! - ensemble size -
  Call print_subtitle ('Ensemble size')
  Call prompt_integer ('Number of ensemble forecasts', nenf, &
       itab=32, ilow=1, lbold=.false.)
! - odds -
  Call print_subtitle ('Odds')
  Call add_check_button ('Show as odds relative to climatology', iodds, .false.)
! - precision -
  Call print_subtitle ('Precision')
  Call prompt_integer ('Number of decimal places (Max 8)', iprec, &
       itab=32, ilow=0, ihgh=8, lbold=.false.)
  Call prompt_integer ('Forecast probability rounding', ifpr, &
       itab=32, ilow=1, ihgh=10, lbold=.false.)
  Call box_close (.true.)
  Call box_close (.false.)
! - prompt to confirm -
  iw = win_prompt(.true.)
!
! Identify error variance type
  If (iw == 1) iev = get_flag(ievo)
  If (iretro == 0 .and. iev == 2) iev = 3
  opts = set_options(iw=iw, i1=iev, i2=iodds, i3=nenf, r1=clf)
  If (opts == 2) Then
     If (ifc == 2) ifc = 1
  End If
  opts = set_options(iw=iw, i4=iprec, i5=ifpr)
  opts = 2
!
  Return
 End Function get_fcast_opts
!
!
!
 Function get_gcm_opts() &
          Result (opts)
!
! Sets EOF options
!
! Modules
  Use data_cpt_constants, Only: ccorrects, cmcs, ncos, nensc, nmcs
  Use gui,                Only: iw, &
                                box_close, box_open, ganged_rbs, init_win, print_text, set_tabs, win_prompt
  Use arrays,             Only: get_flag, set_flag
  Use settings,           Only: iensc, igauss, igauss_bk, igcms, imc, intp, &
                                set_options
  Use pcs,                Only: ih_eof
  Use analysis,           Only: jm_gauss
  Use cpt_actions,        Only: reset
!
! Function type
  Integer :: opts
!
! Locals
!
! Local arrays
  Integer, Dimension(    2) :: intps  ! - interpolation options -
  Integer, Dimension(nensc) :: ienscs ! - ensemble forecasting flags -
  Integer, Dimension( nmcs) :: imcs   ! - model combination options flags -
  Integer, Dimension( ncos) :: icos   ! - correction options flags -
!
! Executable Statements
!
! Identify standardization option
  opts = set_options(i1=intp, i2=igcms, i3=imc, i4=iensc)
  Call set_flag (2-intp, intps)
  Call set_flag (igcms+1, icos)
  Call set_flag (imc, imcs)
  Call set_flag (iensc, ienscs)
!
! Check for completed calculations
  opts = 1
  If (reset('Changing the GCM options', .true.) == 1) Return
!
! Prompt for GCM options
  Call init_win (title='GCM options', &
       ihandle=ih_eof) ! - get window handle -
  Call set_tabs (1, (/37/))
! - interpolation option -
  Call box_open (.false., .false., &
       title='Model grid')
  Call print_text ('Select method:')
  Call ganged_rbs (2, intps, (/'Interpolate      ', 'Nearest gridpoint'/))
  Call box_close (.false.)
! - model standardization -
  Call box_open (.false., .false., &
       'Model climatology')
  Call print_text ('Select method:')
  Call ganged_rbs (ncos, icos, ccorrects)
  Call box_close (.false.)
! - model combination -
  Call box_open (.false., .false., &
       title='Model combination')
  Call print_text ('Select method:')
  Call ganged_rbs (nmcs, imcs, cmcs)
  Call box_close (.false.)
! - ensemble counting -
!$$$$$$   Call box_open ('Ensemble forecasting', .false.)
!$$$$$$   Call print_text ('Select method:')
!$$$$$$   Call ganged_rbs (nensc, ienscs, censc)
!$$$$$$   Call box_close (.false.)
  iw = win_prompt(.true.)
!
! Set option
  If (iw == 1) Then
     intp = 2 - get_flag(intps)
     igcms = get_flag(icos) - 1
     imc = get_flag(imcs)
     iensc = get_flag(ienscs)
  End If
  opts = set_options(iw=iw, i1=intp, i2=igcms, i3=imc, i4=iensc)
  If (iw == 1) Then
     If (igcms == 3) Then
        jm_gauss = 1
        igauss = igauss_bk
     Else
        jm_gauss = 0
        igauss_bk = igauss
        igauss = 0
     End If
  End If
!
  Return
 End Function get_gcm_opts
!
!
!
#if GUI == 1
 Function get_gfont_opts() &
          Result (opts)
!
! Prompts for graphics font scaling
!
! Modules
  Use data_numbers, Only: zero, tenth
  Use labels,       Only: c_def=>c_lab, &
                          l_defaults, l_fntscal, &
                          make_label
  Use screen,       Only: font1, font2, font3, &
                          reset_gfonts
  Use gui,          Only: iw, &
                          init_win, print_text, prompt_real, win_prompt
  Use settings,     Only: set_options
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup current setting
  opts = set_options(r1=font1%fscale, r2=font1%fscale, r3=font3%fscale)
!
! Prompt for font scaling
  Call init_win (title=Trim(l_fntscal%c))
  Call print_text ('Set font scaling:', &
       lbold=.true.)
  Call print_text (' ')
  Call prompt_real ('Default font', font1%fscale, &
       itab=25, rfly=tenth, elow=zero, lbold=.false.)
  Call prompt_real ('Title font',font2%fscale, &
       itab=25, rfly=tenth, elow=zero, lbold=.false.)
  Call prompt_real ('Numeric labels font', font3%fscale, &
       itab=25, rfly=tenth, elow=zero, lbold=.false.)
  c_def = make_label('T', l_defaults, &
          lak=.true.)
  iw = win_prompt(.true., &
       cextra=Trim(c_def), fextra=reset_gfonts)
!
! Confirm changes
  opts = set_options(iw=iw, r1=font1%fscale, r2=font1%fscale, r3=font3%fscale)
  opts = 1
!
  Return
 End Function get_gfont_opts
!
!
!
#endif
 Function get_goodness_opts() &
          Result (opts)
!
! Prompts for correlation coefficient to calculate goodness index
!
! Modules
  Use data_cpt_constants, Only: ncorr
  Use labels,             Only: ca_corrs, &
                                l_goodness
  Use settings,           Only: igood
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for goodness index correlation coefficient
  Call get_option_list (Trim(l_goodness%c), 'Select score for calculating goodness index:', ncorr, ca_corrs, igood)
  opts = 1
!
  Return
 End Function get_goodness_opts
!
!
!
#if GUI == 1
 Function get_gsize_opts() &
          Result (opts)
!
! Prompts for graphics scaling
!
! Modules
  Use data_numbers,  Only: zero,tenth
  Use labels,        Only: l_graphscal
  Use gui,           Only: iw, &
                           init_win, prompt_real, win_prompt
  Use settings,      Only: set_options
  Use graphics,      Only: gszs
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! External routines
  Interface
   Subroutine init_image_sizes ()
   End Subroutine init_image_sizes
  End Interface
!
! Executable Statements
!
! Backup current setting
  opts = set_options(r1=gszs)
  Call init_win (title=Trim(l_graphscal%c))
  Call prompt_real ('Graphics Scaling', gszs, &
       rfly=tenth, elow=zero, lbold=.true.)
  iw = win_prompt(.true.)
!
! Confirm changes
  opts = set_options(iw=iw, r1=gszs)
!
! Updaet graphics sizes
  If (iw == 1) Then
     Call init_image_sizes ()
  End If
  opts = 1
!
  Return
 End Function get_gsize_opts
!
!
!
#endif
#if GUI == 1
 Function get_marker_opts() &
          Result (opts)
!
! Prompts for station markers
!
! Modules
  Use data_numbers,       Only: zero, tenth
  Use data_gui_constants, Only: nmrk
  Use labels,             Only: ca_markers, &
                                l_stnmrks
  Use gui,                Only: iw, &
                                ganged_rbs, init_win, print_text, prompt_real, win_prompt
  Use arrays,             Only: get_flag
  Use settings,           Only: set_options
  Use markers,            Only: gmarker
!
! Function type
  Integer :: opts
!
! Locals
!
! Local arrays
  Integer, Dimension(nmrk) :: jmrk
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup current setting
  opts = set_options(i1=gmarker%imrk, r1=gmarker%rscale)
!
! Prompt for station marker
  jmrk(:) = 0
  jmrk(gmarker%imrk) = 1
  Call init_win (title=Trim(l_stnmrks%c))
  Call print_text ('Select station marker:', &
       lbold=.true.)
  Call print_text (' ')
  Call ganged_rbs (nmrk, jmrk, ca_markers)
!
! Prompt for scaling
  Call print_text (' ')
  Call prompt_real ('Scaling', gmarker%rscale, &
       rfly=tenth, elow=zero, lbold=.true.)
  iw = win_prompt(.true.)
!
! Confirm changes
  If (iw == 1) gmarker%imrk=get_flag(jmrk)
  opts = set_options(iw=iw, i1=gmarker%imrk, r1=gmarker%rscale)
  opts = 1
!
  Return
 End Function get_marker_opts
!
!
!
#endif
 Function get_missing_opts() &
          Result (opts)
!
! Prompts for missing value settings and options
!
! Modules
  Use data_cpt_constants, Only: ia_pfv, ipm, nmo
  Use labels,             Only: l_fcasts, l_missings, l_predictandy, l_predictorx                                
#if GUI == 1
  Use gui,                Only: gui_creturn, print_text
#endif
  Use gui,                Only: iw, &
                                boxes_open, ganged_rbs, init_win, win_prompt
  Use arrays,             Only: get_flag
  Use iofiles,            Only: xfile, yfile
  Use settings,           Only: set_options
  Use fields,             Only: xfield, yfield
  Use missing,            Only: immx, immy, ipmx, ipmy, ipvx, ipvy, nnsx, nnsy
  Use analysis,           Only: ianal
!
! Function type
  Integer :: opts
!
! Locals
!
! Local arrays
  Integer :: irmx(nmo) ! - replace missing X option flags -
  Integer :: irmy(nmo) ! - replace missing Y option flags -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Max, Trim
!
! Executable Statements
!
! Backup missing value settings
  opts = set_options(i1=immx, i2=immy, i3=ipmx, i4=ipmy, i5=ipvx, i6=ipvy, i7=nnsx, i8=nnsy, &
         r1=xfield(1)%rmiss, r2=yfield(1)%rmiss)
!
! Prompt for missing value options
  Call init_win (title=Trim(l_missings%c))
!
! X variables
! - identify missing value replacement option -
  irmx(:) = 0
  irmx(Max(1, immx)) = 1
! - prompt for missing value settings -
  Call boxes_open (2, 1, .true., .false.)
  If (ianal /= ia_pfv) Then
     Call missing_opts (Trim(l_predictorx%c), .true., xfile%idstr, xfield(1)%rmiss, ipmx, ipvx, nnsx, irmx)
  Else
     Call missing_opts (Trim(l_fcasts%c), .false., xfile%idstr, xfield(1)%rmiss, ipmx, ipvx, nnsx, irmx)
  End If
#if GUI == 1
  Call print_text (' ', &
       ladvance=.false.)
#endif
!
! Y variables
! - identify missing value replacement option -
  irmy(:) = 0
  irmy(Max(1, immy)) = 1
! - prompt for missing value settings -
  If (ianal /= ia_pfv) Then
     Call missing_opts (Trim(l_predictandy%c), .true., yfile%idstr, yfield(1)%rmiss, ipmy, ipvy, nnsy, irmy)
  Else
     Call missing_opts ('Observations', .true., yfile%idstr, yfield(1)%rmiss, ipmy, ipvy, nnsy, irmy)
  End If
  iw = win_prompt(.true.)
!
! Identify missing value replacement methods
! - X variables -
  If (iw == 1) Then
     immx = get_flag(irmx)
! - Y variables -
     immy = get_flag(irmy)
  End If
!
! Confirm settings
  opts = set_options(iw=iw, i1=immx, i2=immy, i3=ipmx, i4=ipmy, i5=ipvx, i6=ipvy, i7=nnsx, i8=nnsy, &
         r1=xfield(1)%rmiss, r2=yfield(1)%rmiss)
!
! Duplicate missing value flag for additional fields
  If (opts /= 3) Then
     If (xfile%nfl > 1) xfield(2:)%rmiss = xfield(1)%rmiss
     If (yfile%nfl > 1) yfield(2:)%rmiss = yfield(1)%rmiss
  End If
!
  Return
!
 Contains
!
!
  Subroutine missing_opts (title, lreplace, idstr, rmiss, ipmv, ipvv, nnsv, irmv)
!
! Pormpts for missing value options
!
! Modules
  Use data_numbers, Only: rp
  Use labels,       Only: ca_dstruct_l
  Use screen,       Only: iwdth
  Use gui,          Only: box_close, box_open, print_text, prompt_integer, prompt_real
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: idstr ! - data structure flag -
!
  Character(Len=*), Intent(In) :: title ! - title -
!
  Logical, Intent(In) :: lreplace ! - prompt for replacement options? -
!
! Input/output scalars
  Integer, Intent(InOut) :: ipmv ! - Maximum % of missing data -
  Integer, Intent(InOut) :: ipvv ! - Maximum % of missing variables -
  Integer, Intent(InOut) :: nnsv ! - number of near-neighbours -
!
  Real(KIND=rp), Intent(InOut) :: rmiss ! - missing value flag -
!
! Input/output scalars
  Integer, Dimension(:), Intent(InOut) :: irmv ! - missing value replacement option flags -
!
! Locals
!
! Local scalars
  Integer :: itab ! - tab position -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for missing value options
  Call box_open (.true., .false., &
       title=title)
  If (iwdth >= 1024) Then
     itab = 32
  Else
     itab = 24
  End If
  Call prompt_real ('Missing value flag', rmiss, &
       itab=itab, iwid=10, lfocus=.true.)
  Call prompt_integer ('Maximum % of missing values', ipmv, &
       itab=itab, ilow=0, ihgh=ipm, iwid=10, lfocus=.true.)
  Call prompt_integer ('Maximum % of missing '//Trim(ca_dstruct_l(idstr)), ipvv, &
       itab=itab, ilow=0, ihgh=ipm, iwid=10, lfocus=.true.)
  If (lreplace) Then
     Call prompt_integer ('Number of near-neighbours', nnsv, &
          itab=itab, ilow=1, iwid=10, lfocus=.true.)
#if GUI == 1
     Call gui_creturn (.true.)
#endif
     Call box_open (.false., .false., &
          'Missing value Replacement')
     Call print_text ('Select method:')
     Call ganged_rbs (4, irmv, &
          crbs=(/'Long-term means     ', 'Long-term medians   ', 'Random numbers      ', 'Best near-neighbours'/))
     Call box_close (.false.)
  End If
  Call box_close (.true.)
  Call box_close (.false.)
!
  Return
  End Subroutine missing_opts
 End Function get_missing_opts
!
!
!
#if GUI == 0
 Function get_nt()
!
! Prompts for length of training period
!
! Modules
  Use data_numbers,        Only: ihuge
  Use labels,              Only: cg_lentp
  Use gui,                 Only: prompt_integer
  Use iofiles,             Only: xfile, yfile
  Use settings,            Only: mnu, ntm, nt, &
                                 get_maxnt
  Use season,              Only: icsv, lensn, lb, lprst, lsn, nsn
!
! Function type
  Integer :: get_nt
!
! Locals
!
! Local scalars
  Integer :: mx ! - maximum length of training period -
  Integer :: mn ! - minimum length of training period -
!
  Character(Len=13) :: cfmt ! - format statement -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Determine maximum length of training period
  mx = xfile%nt
  If (xfile%nt > 0) Then
     If (yfile%nt > 0) Then
        mn = mnu
        Call get_maxnt (icsv, lprst, lensn, lsn(lb:), lb, nsn, mx, ntm)
     Else
        mn = 0
        mx = xfile%nt
     End If
  Else
     mn = 0
     If (yfile%nt > 0) Then
        mx = yfile%nt
     Else
        mx = 0
     End If
  End If
!
! Get training data settings
  Print *, ' '
  If (mx > 0) Then
     Write (Unit=cfmt, Fmt='(A,2(I1,A))') 
     Write (Unit=*, Fmt='(A,I0,A,I0,A)', Advance='no') Trim(cg_lentp)//' (minimum is ', mn, '; maximum is ', mx, '): '
  Else
     Write (Unit=cfmt, Fmt='(A,I1,A)') '(A,I0,A)'
     Write (Unit=*, Fmt=cfmt, Advance='no') Trim(cg_lentp)//' (minimum is ', mn, '): '
     mx = ihuge
  End If
  Call prompt_integer(' ', nt, &
       ilow=mn, ihgh=mx)
  get_nt = 0
!
  Return
 End Function get_nt
!
!
!
#endif
!
!
!
 Function get_startat_opts() &
          Result (opts)
!
! Prompts for start-at date
!
! Modules
#if GUI == 1
  Use clrwin$,            Only: set_integer_limits$
#endif
  Use data_cpt_constants, Only: c_file
  Use screen,             Only: window_update
#if GUI == 1
  Use gui,                Only: ih_nt
#endif
  Use iofiles,            Only: mfile, xfile, yfile, zfile
  Use settings,           Only: nt, ntm, &
                                get_maxnt
#if GUI == 1
  Use settings,           Only: mnu
#endif
  Use season,             Only: icsv, lb, lensn, lsn, lprst, nsn, &
                                len_ssn
  Use forecast,           Only: change_fcast
  Use cpt_file,           Only: ixyz
!
! Function type
  Integer :: opts
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Exeutable Statements
!
! Prompt for start-at date
  opts = 1
  Select Case (c_file(ixyz))
   Case ('X', 'x')
     Call get_startat ('X', 1, (/len_ssn(xfile%prd1%sdate, xfile%prd1%edate, xfile%iseq)/), xfile, ifail)
   Case ('Y', 'y')
     Call get_startat ('Y', lb, lsn(lb:), yfile, ifail, &
          cprd1=mfile%cprd1, cprdn=mfile%cprdn)
     If (ifail == 0) Then
        mfile%bdate = yfile%bdate
        mfile%fdate = yfile%fdate
     End If
   Case ('Z', 'z')
     Call get_startat ('Forecast', 1, (/len_ssn(zfile%prd1%sdate, zfile%prd1%edate, zfile%iseq)/), zfile, ifail, &
          itp=yfile%idstr, cbf=change_fcast)
   Case ('1')
     Call get_startat ('First', 1, (/len_ssn(xfile%prd1%sdate, xfile%prd1%edate, xfile%iseq)/), xfile, ifail)
   Case ('2')
     Call get_startat (' Second', 1, (/len_ssn(yfile%prd1%sdate, yfile%prd1%edate, yfile%iseq)/), yfile, ifail)
   Case ('3')
     Call get_startat ('Update', 1, (/len_ssn(zfile%prd1%sdate, zfile%prd1%edate, zfile%iseq)/), zfile, ifail)
  End Select
!
! Reset training period limits
  Select Case (c_file(ixyz))
   Case ('X', 'x', 'Y', 'y')
     If (xfile%cfile%lset .or. yfile%cfile%lset) Then
        Call get_maxnt (icsv, lprst, lensn, lsn(lb:), lb, nsn, nt, ntm)
        Call window_update (nt)
#if GUI == 1
        Call set_integer_limits$ (ih_nt, mnu, nt)
#endif
     End If
  End Select
!
  Return
!
 Contains
!
!
  Subroutine get_startat (cxyz, lb, lsn, afile, ifail, &
             cprd1, cprdn, itp, cbf)
!
! Modules
  Use data_io_constants, Only: lprd
  Use labels,            Only: cg_data1, cg_datan, &
                               l_startat
  Use time,              Only: iseq, &
                               set_bdate
  Use screen,            Only: window_update
  Use gui,               Only: iw, &
                               init_win, print_cvalue, prompt_pdate, win_prompt
  Use iofiles,           Only: ifile
  Use settings,          Only: record_change, set_options
  Use forecast,          Only: cdates, &
                               update_fprds
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lb ! - lower-bound -
!
  Character(Len=*), Intent(In) :: cxyz ! - file type -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itp ! - target period indicator -
!
  Character(Len=*), Intent(In), Optional :: cprd1 ! - first period of data -
  Character(Len=*), Intent(In), Optional :: cprdn ! - last period of data -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - input file -
!
! Input arrays
  Integer, Dimension(lb:), Intent(In) :: lsn ! - length of season -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Procedure arguments
  Integer, External, Optional :: cbf ! - call-back function -
!
  Character(Len=lprd) :: cp1 ! - first period of data -
  Character(Len=lprd) :: cpn ! - last period of data -
!
! Locals
!
! Local scalars
  Integer :: itab ! - tab setting -
  Integer :: iwid ! - width -
!
  Logical :: ltp ! - add target period? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim, Max, Present, Trim
!
! Exeutable Statements
!
! Backup current date
  ifail = 1
  opts = set_options(i1=afile%fdate%iyr, i2=afile%fdate%imn, i3=afile%fdate%idy)
!
! Prompt for start-at date
  Call init_win (title=Trim(l_startat%c)//' ('//cxyz//')')
  itab = 30
  If (Present(cprd1)) Then
     cp1 = cprd1
  Else
     cp1 = afile%cprd1
  End If
  If (Present(cprdn)) Then
     cpn = cprdn
  Else
     cpn = afile%cprdn
  End If
  iwid = Max(Len_Trim(cp1), Len_Trim(cpn))
  Call print_cvalue (Trim(cg_data1), cp1, &
       iwid=iwid, itab=itab)
  Call print_cvalue (Trim(cg_datan), cpn, &
       iwid=iwid, itab=itab)
  If (Present (itp)) Then
     If (itp > 0) Then
        ltp = .true.
     Else
        ltp = .false.
     End If    
  Else
     ltp = .false.
  End If    
  If (.not.ltp) Then
     Call prompt_pdate (Trim(l_startat%c), iseq, afile%fdate, &
          itab=itab, lbold=.true., lfocus=.true., dlimits=(/afile%prd1%sdate, afile%prdn%sdate/), &
          modify2=cbf)
  Else
     iw = update_fprds()
     Call prompt_pdate (Trim(l_startat%c), iseq, afile%fdate, &
          itab=itab, lbold=.true., lfocus=.false., dlimits=(/afile%prd1%sdate, afile%prdn%sdate/), &
          modify1=update_fprds, modify2=cbf)
     Call print_cvalue ('Target period', cdates(1))
  End If
! - prompt to confirm -
  iw = win_prompt(.true.)
!
! Restore old settings if cancelled
  opts = set_options(iw=iw, i1=afile%fdate%iyr, i2=afile%fdate%imn, i3=afile%fdate%idy)
!
! Reset beginning dates
  If (iw /= 1) Return
  Call set_bdate (lb, lsn(lb:), nsn, afile%prd1%sdate, afile%fdate, afile%bdate, &
       cprdb=afile%cprdb)
  Call window_update (afile%cprdb)
  ifail = record_change()
  ifail = 0
!
  Return
  End Subroutine get_startat
 End Function get_startat_opts
!
!
!
 Function get_tailoring_opts () &
          Result (opts)
!
! Prompts for threshold and other tailoring settings
!
! Modules
#if GUI == 1
  Use clrwin$,            Only: winio$
#endif
  Use data_numbers,       Only: rp, zero, one, onehth, sfMax
  Use data_cpt_constants, Only: cstds,  cthrs, ng, nstd, nthr, nts
  Use labels,             Only: c_lab, &
                                l_probs, l_tailor, l_threshs, &
                                ls_cats, ls_threshs, &
                                make_label
  Use time,               Only: pdate, &
                                Operator(+), Operator(==), &
                                iseq
  Use gui,                Only: iw, &
                                box_close, box_open, boxes_open, ganged_rbs, init_win, print_text, prompt_pdate, prompt_real, &
                                win_prompt
  Use menus,              Only: im_trns
  Use progress_meter,     Only: lcalc
  Use errors,             Only: cpt_error
  Use arrays,             Only: fcast, fpls, yrpls, &
                                get_flag, set_flag
  Use iofiles,            Only: yfile
  Use settings,           Only: igauss, igauss_bk, igto, istd, izero, nt, &
                                record_change
  Use season,             Only: lb, nsn
  Use categories,         Only: ithr, pthr, refs, thr, &
                                check_refdates
  Use analysis,           Only: ianal, icalc
  Use forecast,           Only: ifc
  Use verif,              Only: irv
  Use cpt_actions,        Only: reset
!
! Function type
  Integer :: opts
!
! Locals
!
! Local scalars
  Integer :: j        ! - threshold index -
  Integer :: istd_old ! - backup standardization option -
  Integer :: ithr_old ! - backup threshold type -
  Integer :: ifail    ! - error indicator -
!
  Character(Len=16) :: cref ! - reference year -
!
! Local arrays
  Integer, Dimension(nthr) :: jthr
  Integer, Dimension(nstd) :: istds ! - standardization options flags -
!
  Real(KIND=rp), Dimension(nts) :: pthr_old ! - backup percentile thresholds -
  Real(KIND=rp), Dimension(nts) :: thr_old  ! - backup absolute thresholds -
!
  Type(pdate), Dimension(  2) :: dlimits  ! - date limits -
  Type(pdate), Dimension(nts) :: refs_old ! - backup reference dates -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any, Sum, Trim
!
! Executable Statements
!
! Backup settings
  opts = 1
  istd_old = istd
  ithr_old = ithr
  thr_old = thr
  pthr_old(:) = pthr(:)
  If (yfile%cfile%lset) Then
     Do j = 1, nts
        If (refs(j)%iyr == 0) refs(j) = yfile%fdate + (j-1)
     End Do
  End If
  refs_old(:) = refs(:)
!
! Convert upper threshold to right-tail probability
  pthr(2) = one - pthr(2)
!
! Identify standardization option
  Call set_flag (istd+1, istds)
!
! Prompt for threshold settings
  Call set_flag (ithr, jthr)
1 Call init_win (title=Trim(l_tailor%c))
  Call boxes_open (1,2, .true., .false.)
! - standardization -
  Call box_open (.true., .false., &
       title='Standardization')
  Call print_text ('Select method:')
  Call ganged_rbs (nstd-(1-izero), istds, &
       crbs=cstds)
  Call box_close (.true.)
  Call box_close (.false.)
! - percentile thresholds -
  Call box_open (.true., .false., &
       title=Trim(l_threshs%c))
  Call boxes_open (1, 3, .true., .false.)
#if GUI == 1
  Call ganged_rbs (3, jthr)
#else
  Call print_text ('Threshold type')
  Call ganged_rbs (nthr, ithr, cthrs(:))
  Select Case (ithr)
   Case (1)
#endif
     Call box_open (.true., .false., &
          title=Trim(l_probs%c))
#if GUI == 1
     iw = winio$('%rb[]%bf %ws:%`bf%ff&', jthr(1), Trim(cthrs(1)))
#endif
     c_lab = make_label('L', ls_cats(ng))
     Call prompt_real ('     '//Trim(c_lab), pthr(2), &
          itab=35, iwid=10, rfly=onehth, rlow=zero, rhgh=one, lbold=.false.)
     c_lab = make_label('L', ls_cats(1))
     Call prompt_real ('     '//Trim(c_lab), pthr(1), &
          itab=35, iwid=10, rfly=onehth, rlow=zero, rhgh=one, lbold=.false.)
     Call box_close (.true.)
     Call box_close (.false.)
! - absolute thresholds -
#if GUI == 0
   Case (2)
#endif
     Call box_open (.true., .false., &
          title=Trim(cthrs(2)))
#if GUI == 1
     iw = winio$('%rb[]%bf %ws:%`bf%ff&', jthr(2), Trim(cthrs(2)))
#endif
     Do j = nts, 1, -1
        Call prompt_real ('     '//Trim(ls_threshs(j)%c), thr(j), &
             itab=35, iwid=10, rlow=-sfMax, rhgh=sfMax, lbold=.false.)
     End Do
     Call box_close (.true.)
     Call box_close (.false.)
! - references -
#if GUI == 0
   Case (3)
#endif
     Call box_open (.true., .false., &
          title=Trim(cthrs(3)))
#if GUI == 1
     iw = winio$('%rb[]%bf %ws:%`bf%ff&', jthr(3), Trim(cthrs(3)))
#else
6    Continue
#endif
     If (yfile%fdate%iyr < yfile%prdn%sdate%iyr) Then
        dlimits(1) = yfile%fdate
        dlimits(2) = yfile%prdn%sdate
        Do j = 1, nts
           Write (Unit=cref, Fmt='(A,I2)') '     Reference',j
           Call prompt_pdate (cref, iseq, refs(j), &
                itab=39, dlimits=dlimits, lbold=.false., lfocus=.true.)
        End Do
     Else
        Do j = 1, nts
           Write (Unit=cref, Fmt='(A,I2)') '     Reference',j
           Call prompt_pdate (cref, iseq, refs(j), &
                itab=39, lbold=.false., lfocus=.true.)
        End Do
     End If
     Call box_close (.true.)
     Call box_close (.false.)
#if GUI == 0
     If (check_refdates(yfile%fdate, nt) /= 0) GoTo 6
  End Select
#endif
  Call box_close (.true.)
  Call box_close (.false.)
  iw = win_prompt(.true.)
!
! Confirm changes
  If (iw == 1) Then
     If (Sum(pthr(:)) > one) Then ! - check for valid probability thresholds -
        Call cpt_error ('threshold_opts', ifail, .false.)
        GoTo 1
     End If
     istd = get_flag(istds) - 1
     If (istd /= istd_old) Then
        If (istd == 3) Then ! - switch on Gaussian flag using gamma distribution for SPI -
           If (icalc == 1 .and. igauss/=2) Then ! - reset If SPI selected -
              If (reset('Selecting SPI', .true.) == 1) Then
                 istd = istd_old
                 iw = 0
              End If
           End If
           igto = 2
           igauss_bk = igauss
           igauss = igto
           im_trns = 0
        Else
           im_trns = 1
        End If
     End If
#if GUI == 1
     ithr = get_flag(jthr)
#endif
     If (Any(pthr(:) /= pthr_old(:))) pthr(2) = one - pthr(2)
     If (ithr /= ithr_old .or. Any(thr /= thr_old) .or. Any(pthr(:) /= pthr_old(:)) .or.  &
         istd /= istd_old .or. .not.(refs(:) == refs_old(:))) Then
        If (ifc == 2) ifc=1
     Else
        Return
     End If
     opts = record_change()
     If (thr(1) > thr(2)) Then ! - ensure ascending order of thresholds -
        thr_old = thr
        thr(1) = thr_old(2)
        thr(2) = thr_old(1)
     End If
     If (icalc == 1) Then
        If (ithr == 3) Then
           ifail = check_refdates(yfile%fdate, nt)
           If (ifail /= 0) Call cpt_error ('check_refdates', ifail, .true.)
        End If
        lcalc = .true.
        ifail = tailoring()
        If (ifail /= 0) Then
           istd = istd_old
           ithr = ithr_old
           thr = thr_old
           pthr(:) = pthr_old(:)
           Return
        End If
        irv = 0
        lcalc = .false.
     End If
  Else
     istd = istd_old
     ithr = ithr_old
     thr = thr_old
     pthr(:) = pthr_old(:)
     Return
  End If
  opts = 2
!
  Return
!
 Contains
!
!
  Function tailoring()
!
! Restandardises data based on tailoring options
!
! Modules
  Use data_numbers,   Only: one
  Use progress_meter, Only: end_progress_meter, set_progress_increment, start_progress_meter, update_progress_meter
  Use arrays,         Only: y, yhat, yret
  Use iofiles,        Only: mya
  Use settings,       Only: igauss, iretro, istd, nenf, nf, nu, nur
  Use distribs,       Only: restdize
  Use categories,     Only: tobs, tobst, &
                            calc_thresholds, set_cv_categories, set_ra_categories
  Use climate,        Only: ave, clim, ncu, sdev
!
! Function type
  Integer :: tailoring
!
! Locals
!
! Local scalars
  Integer :: ie    ! - ensemble member index -
  Integer :: isn   ! - season index -
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
!
! Executable Statements
!
! Confirm settings
  tailoring=-1
  Call start_progress_meter ('Recalculating thresholds and categories', .true.)
  Select Case (igauss)
   Case (0)
     Call set_progress_increment (one/Real(3+mya*(iretro+1), Kind=rp), .false.)
   Case Default
     Call set_progress_increment (one/Real(4+mya*(iretro+1), Kind=rp), .false.)
  End Select
! - restandardise -
  If (istd /= istd_old) Then
     Do isn = lb, nsn
        Call restdize (istd_old, istd, igauss, isn, mya, nu, ave(:,isn), sdev(:,isn), y(1:mya,1:nu,isn))
        Call restdize (istd_old, istd, igauss, isn, mya, ncu(isn), ave(:,isn), sdev(:,isn), clim(1:mya,:,isn))
     End Do
     Call restdize (istd_old, istd, igauss, lb, mya, nu, ave(:,lb), sdev(:,lb), yhat(:,:))
     If (iretro == 1) Then
        Call restdize (istd_old, istd, igauss, lb, mya, nur, ave(:,lb), sdev(:,lb), yret(:,:))
        If (ifc == 2) Then
           Call restdize (istd_old, istd, igauss, lb, mya, nur, ave(:,lb), sdev(:,lb), yrpls(:,:)%lower)
           Call restdize (istd_old, istd, igauss, lb, mya, nur, ave(:,lb), sdev(:,lb), yrpls(:,:)%upper)
        End If
     End If
     If (ifc == 2) Then
        Do ie = 0, nenf
           Call restdize (istd_old, istd, igauss, lb, mya, nf, ave(:,lb), sdev(:,lb), fcast(:,:,ie))
        End Do
        Call restdize (istd_old, istd, igauss, lb, mya, nf, ave(:,lb), sdev(:,lb), fpls(:,:)%lower)
        Call restdize (istd_old, istd, igauss, lb, mya, nf, ave(:,lb), sdev(:,lb), fpls(:,:)%upper)
     End If
  End If
  If (update_progress_meter(.false.) /= 0) Return
! - redefine categories and thresholds -
  If (set_cv_categories(lb) /= 0) Return
  If (iretro == 1) Then
     If (set_ra_categories(ianal, lb, nsn, 0, ncu(lb:), clim(:,:,lb:)) /= 0) Return
  End If
  Call calc_thresholds (igauss, mya, nsn, ncu(lb:), clim(:,:,lb:), tobs, tobst, ifail)
  If (ifail /= 0) Return
  Call end_progress_meter (.true., .true.)
  tailoring = 0
!
  Return
  End Function tailoring
 End Function get_tailoring_opts
!
!
!
 Function get_targetssn_opts() &
          Result (opts)
!
! Prompts for target season
!
! Modules
  Use screen,      Only: window_update
  Use settings,    Only: nt, ntm, &
                         get_maxnt
  Use fields,      Only: yfield
  Use season,      Only: icsv, lb, lensn, lprst, lsn, nsn, &
                         prompt_season, set_seasons
  Use cpt_actions, Only: reset

!
! Function type
  Integer :: opts
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Exeutable Statements
!
! Check for completed calculations
  opts = 1
  If (reset('Changing the season settings', .true.) == 1) Return
!
! Identify season
  nsn = 1
1 Call prompt_season (nsn, lb, yfield, ifail)
  If (ifail/=0) Return
!
! Set seasons
  Call set_seasons (nsn, yfield(:)%cssn, yfield(:)%tprd)
  If (ifail /= 0) GoTo 1
!
! Set default length of training period
  Call get_maxnt (icsv, lprst, lensn, lsn(lb:), lb, nsn, nt, ntm)
  Call window_update (nt)
!
  Return
 End Function get_targetssn_opts
!
!
!
#if GUI == 1
 Function get_text_opts() &
          Result (opts)
!
! Prompts for text scaling
!
! Modules
  Use data_numbers, Only: zero, tenth, one
  Use labels,       Only: c_def=>c_lab, &
                          l_defaults, l_txtscal, &
                          make_label
  Use screen,       Only: dts_adj
  Use gui,          Only: iw, &
                          init_win, print_text, prompt_real, win_prompt
  Use settings,     Only: set_options
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup current setting
  opts = set_options(r1=dts_adj)
!
! Prompt for font scaling
  Call init_win (title=Trim(l_txtscal%c))
  Call print_text ('Set text scaling:', &
       lbold=.true.)
  Call print_text (' ')
  Call prompt_real ('Default text scaling', dts_adj, &
       itab=25, rfly=tenth, elow=zero, lbold=.false.)
  c_def = make_label('T', l_defaults, &
          lak=.true.)
  iw = win_prompt(.true., &
       cextra=Trim(c_def), rset1=dts_adj, rval1=one)
!
! Confirm changes
  opts = set_options(iw=iw, r1=dts_adj)
  opts = 1
!
  Return
 End Function get_text_opts
!
!
!
#endif
 Function get_time_opts() &
          Result (opts)
!
! Prompts for handling of lengths of months
!
! Modules
  Use data_cpt_constants, Only: ctls, ntl
  Use labels,             Only: l_timesets
  Use time,               Only: itime
  Use gui,                Only: iw, &
                                ganged_rbs, init_win, print_text, win_prompt
  Use arrays,             Only: get_flag, set_flag
!
! Function type
  Integer :: opts
!
! Locals
!
! Local arrays
  Integer, Dimension(ntl) :: itls ! - time options -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup time settings
  Call set_flag (itime, itls)
!
! Prompt for transformation option
  Call init_win (title=Trim(l_timesets%c))
  Call print_text ('Select option for lengths of months:')
  Call ganged_rbs (ntl, itls, ctls)
  iw = win_prompt(.true.)
!
! Set option
  If (iw == 1) itime = get_flag(itls)
!
! Restore old settings If cancelled
  opts = 2
!
  Return
 End Function get_time_opts
!
!
!
 Function get_transform_opts() &
          Result (opts)
!
! Prompts for Y-data transformation options
!
! Modules
  Use data_cpt_constants, Only: ctos, nto
  Use labels,             Only: l_transforms
  Use settings,           Only: igauss, igto
!
! Function type
  Integer :: opts
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for transformation option
  Call get_option_list (Trim(l_transforms%c), 'Select distribution for transformation:', nto, ctos, igto)
  If (igauss > 0) igauss = igto    
  opts = 2
!
  Return
 End Function get_transform_opts
!
!
!
 Function get_verif_opts() &
          Result (opts)
!
! Prompts for verification options
!
! Modules
  Use data_cpt_constants, Only: nbs, nbss
  Use labels,             Only: l_verifsets
#if GUI == 1
  Use gui,                Only: iw, &
                                win_prompt
#endif
  Use gui,                Only: init_win, print_text, prompt_boolean, prompt_list
  Use settings,           Only: ibi, ibs, nb, &
                                set_options
!
! Function type
  Integer :: opts
!
! Locals
!
! Local scalars
  Integer :: i ! - bin index -
  Integer :: n ! - size of bins -
!
! Local arrays
  Character(Len=16), Dimension(nbs) :: cnbs ! - number of probability bins -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Construct labels
  Do i = 1, nbs
     n = 100/(nbss(i) - 1)
     Write (Unit=cnbs(i), Fmt='(A,I0,A)') '(', n, '% bins)'
  End Do
!
! Backup verification settings
  opts = set_options(i1=ibi, i2=nb, i3=ibs)
!
! Prompt for number of probability bins
  Call init_win (title=Trim(l_verifsets%c))
! - number of probability bins -
  Call prompt_list ('Number of probability bins', .false., cnbs(:), nbs, ibi, &
       ns=nbss(:))
  Call print_text (' ')
! - binned probabilities -
  Call prompt_boolean ('Calculate scores using binned probabilities', ibs)
#if GUI == 1
! - prompt to confirm -
  iw = win_prompt(.true.)
  nb = nbss(ibi)
!
! Restore old settings If cancelled
  opts = set_options(iw=iw, i1=ibi, i2=nb, i3=ibs)
#endif
  opts = 2
!
  Return
 End Function get_verif_opts
!
!
!
 Function change_language()
!
! Changes language of CPT prompts
!
! Modules
#if GUI == 0
  Use data_text,         Only: nlang, &
                               clang
#endif
  Use data_text,         Only: ilang
  Use data_io_constants, Only: ddir
  Use labels,            Only: init_chars, init_labels
#if GUI == 0
  Use labels,            Only: l_exit
  Use gui,               Only: prompt_integer
#endif
  Use errors,            Only: cpt_error
!
! Function type
  Integer :: change_language
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for language preferences
1 change_language = 0
#if GUI == 0
  Print *, ' '
  Do ilang = 1, nlang
     Write (Unit=*, Fmt='(1X,I1,2A)') ilang, '.  ', clang(ilang)
  End Do
  Write (Unit=*, Fmt='(1X,A)') '0.  '//Trim(l_exit%c)
  Call prompt_integer (' ', ilang)
  Select Case (ilang)
! - quit -
   Case (0)
     Return
! - set language -
   Case (1:nlang)
#endif
     Call init_chars ()
     Call init_labels (ilang, Trim(ddir), ifail)
     If (ifail /= 0) Then
        Call cpt_error ('init_labels', ifail, .false.)
        change_language = 1
     End If
! - error -
#if GUI == 0
     If (ifail /= 0) GoTo 1
   Case Default
     GoTo 1
  End Select
!
#endif
  Return
 End Function change_language
!
!
!
 Subroutine get_option_list (ctitle, cdesc, nopt, copts, iopt)
!
! Prompts for option from a list
!
! Modules
  Use gui,                Only: iw, &
                                ganged_rbs, init_win, print_text, win_prompt
  Use arrays,             Only: get_flag, set_flag
  Use settings,           Only: set_options
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nopt ! - number of options -
!
  Character(Len=*), Intent(In) :: ctitle ! - title -
  Character(Len=*), Intent(In) :: cdesc  ! - description -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: copts ! - list of options -
!
! Input/output scalars
  Integer, Intent(InOut) :: iopt ! - selected options -
!
! Locals
!
! Local arrays
  Integer, Dimension(nopt) :: iopts ! - options -
!
! Executable Statements
!
! Backup settings
  iw = set_options(i1=iopt)
  Call set_flag (iopt, iopts)
!
! Prompt for option
  Call init_win (title=ctitle)
  Call print_text (cdesc)
  Call ganged_rbs (nopt, iopts, copts)
  iw = win_prompt(.true.)
!
! Set option
  If (iw == 1) iopt = get_flag(iopts)
!
! Restore old settings If cancelled
  iw = set_options(iw=iw, i1=iopt)
!
  Return
 End Subroutine get_option_list
End Module cpt_options
