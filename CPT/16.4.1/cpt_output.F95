! Author: Simon Mason
Module cpt_output
!
! Modules
  Use data_numbers,      Only: rp
  Use data_io_constants, Only: lvar
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: save_data, save_good, save_results, update_results, write_results
#if GUI == 0
  Public :: print_results
#endif
!
! Scalars
!
! Integer scalars
  Integer, Private :: ihp1 = 1 ! - property sheet 1 -
  Integer, Private :: ihp2 = 2 ! - property sheet 2 -
  Integer, Private :: ihp3 = 3 ! - property sheet 3 -
  Integer, Private :: ihp4 = 4 ! - property sheet 4 -
  Integer, Private :: ihp5 = 5 ! - property sheet 5 -
  Integer, Private :: ihp6 = 6 ! - property sheet 6 -
!
  Character(Len=lvar+1), Private :: cout ! - output field -
!
Contains
!
!
 Function save_data(ifc)
!
! Prompts for output files and saves requested results
!
! Modules
#if GUI == 1
  Use clrwin$,            Only: update_window$, winio$
  Use data_numbers,       Only: one, two
#endif
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use labels,             Only: cg_indata_t, &
                                l_data, l_fcasts, l_modes, l_outres
#if GUI == 1
  Use labels,             Only: l_close, l_edit, l_file, l_options
  Use gui,                Only: ih_con, &
                                set_text_size, win_prompt
#else
  Use labels,             Only: l_exit
  Use gui,                Only: print_text, print_warning, prompt_integer
#endif
  Use gui,                Only: iw, &
                                box_close, box_open, gui_creturn, init_psheet, init_win, window_end
#if GUI == 1
  Use menus,              Only: add_menu
  Use iofiles,            Only: itags,  &
                                files_reset, get_fmt
#else
  Use iofiles,            Only: new_ids, &
                                geto_gen
#endif
  Use iofiles,            Only: avfile, ccfile, cofile, fcfile, fefile, flfile, fofile, fpfile, fsfile, fvfile, fxfile, pbfile, &
                                rbfile, rcfile, rlfile, rpfile, thfile, xefile, xifile, xlfile, xmfile, xofile, xsfile, xtfile, &
                                ycfile, yefile, yhfile, ylfile, ymfile, yofile, yrfile, ysfile, ytfile, zifile,                 &
                                prompt_results_file
  Use settings,           Only: iretro
#if GUI == 1
  Use settings,           Only: igcms
#endif
  Use scores,             Only: lcorr
  Use analysis,           Only: ianal, ihind
!
! Function Type
  Integer :: save_data
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ifc ! - recalculate forecast flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
#if GUI == 0
  Intrinsic Ubound
!
! Executable Statements
!
! Get output files
1 Call init_win (title=Trim(l_outres%c), &
       lresize=.true.)
#endif
!
! Create input data property sheet
  Call init_psheet (ihp1, Trim(cg_indata_t))
  Call box_open (.false., .true.)
  Call box_open (.true., .false., & 
       title=Trim(l_data%c))
  Call prompt_results_file (xofile, .true.)
  If (ianal == ia_gcm .or. ianal == ia_pfv) Then
     Call prompt_results_file (xifile, .false.)
  End If
  Call prompt_results_file (yofile, .false.)
  Select Case (ianal)
   Case Default
     Call prompt_results_file (ycfile, .false.)
     If (iretro == 1) Then
        Call prompt_results_file (rcfile, .false.)
     End If
   Case (ia_pfv)
     Call prompt_results_file (rcfile, .false.)
  End Select
  Call box_close (.true.)
  Call box_close (.false.)
  Call window_end ()
!
! Create hindcasts property sheet
  If (ihind == 1) Then
     Call init_psheet (ihp2, 'Hindcasts')
     Call box_open (.false., .true.)
     Call box_open (.true., .false., &
          title='Predictions')
     Call prompt_results_file (yhfile, .true.)
     If (iretro == 1) Then
        Call prompt_results_file (yrfile, .false.)
        Call prompt_results_file (rpfile, .false.)
        Call prompt_results_file (rlfile, .false.)
     End If
     Call box_close (.true.)
     Call box_close (.false.)
     Call window_end ()
  End If
!
! Create EOF property sheet
  Select Case (ianal)
   Case (ia_cca, ia_pcr)
     Call init_psheet (ihp3, Trim(l_modes%c))
! - X EOF files -
     Call box_open (.false., .true.)
     Call box_open (.true., .false., &
          title='Prefiltering')
     Call box_open (.false., .false., &
          title='X Modes files:')
     Call prompt_results_file (xefile, .true.)
     Call prompt_results_file (xlfile, .false.)
     Call prompt_results_file (xsfile, .false.)
     Call box_close (.false.)
! - Y EOF files -
     If (ianal == ia_cca) Then
        Call gui_creturn (.false.)
        Call box_open (.false., .false., &
             title='Y Modes files:')
        Call prompt_results_file (yefile, .true.)
        Call prompt_results_file (ylfile, .false.)
        Call prompt_results_file (ysfile, .false.)
        Call box_close (.false.)
     End If
     Call box_close (.true.)
     Call box_close (.false.)
     Call window_end ()
  End Select
!
! Create CCA property sheet
  Select Case (ianal)
   Case (ia_cca)
     Call init_psheet (ihp4, 'CCA Results')
     Call box_open (.false., .true.)
     Call box_open (.true., .false., &
          title='CCA files')
     Call box_open (.false., .false., &
          title='Canonical correlations:')
     Call prompt_results_file (ccfile, .true.)
     Call box_close (.false.)
     Call gui_creturn (.false.)
     Call box_open (.false., .false., &
          title='X CCA maps:')
     Call prompt_results_file (xmfile, .true.)
     Call prompt_results_file (xtfile, .false.)
     Call box_close (.false.)
     Call gui_creturn (.false.)
     Call box_open (.false., .false., &
          title='Y CCA maps:')
     Call prompt_results_file (ymfile, .true.)
     Call prompt_results_file (ytfile, .false.)
     Call box_close (.false.)
     Call box_close (.true.)
     Call box_close (.false.)
     Call window_end ()
!
! Create PCR property sheet
   Case (ia_pcr)
     Call init_psheet (ihp4, 'PCR Results')
! - regression output files -
     Call box_open (.false., .true.)
     Call box_open (.true., .false., &
          title='Regression files:')
     Call prompt_results_file (rbfile, .true.)
     Call prompt_results_file (pbfile, .false.)
     Call box_close (.true.)
     Call box_close (.false.)
     Call window_end ()
!
! Create MLR/GCM property sheet
   Case (ia_mlr, ia_gcm)
     If (ianal == ia_mlr) Then
        Call init_psheet (ihp4, 'MLR Results')
     Else
        Call init_psheet (ihp4, 'GCM Results')
     End If
! - regression output files -
     Call box_open (.false., .true.)
     Call box_open (.true., .false., &
          title='Regression files:')
     Call prompt_results_file (rbfile, .true.)
     Call box_close (.true.)
     Call box_close (.false.)
     Call window_end ()
  End Select
!
! Create forecasts property sheet
  If (ifc == 2) Then
     Call init_psheet (ihp5, Trim(l_fcasts%c))
     Call box_open (.false., .true.)
     Call box_open (.true., .false., &
          title=Trim(l_fcasts%c))
! - forecast probabilities -
     Call box_open (.false., .false., &
          title='Probabilistic forecasts:')
     Call prompt_results_file (fpfile, .true.)
     Call prompt_results_file (fofile, .false.)
     Call box_close (.false.)
! - forecasts -
     If (ihind == 1) Then
        Call gui_creturn (.false.)
        Call box_open (.false., .false., &
             title='Deterministic forecasts:')
        Call prompt_results_file (fvfile, .true.)
        Call prompt_results_file (fsfile, .false.) ! - forecast ensembles -
        Call prompt_results_file (flfile, .false.) ! - prediction limits -
        Call prompt_results_file (fefile, .false.) ! - prediction error variances -
        Call box_close (.false.)
! - predictor PC scores -
        Call gui_creturn (.false.)
        Call box_open (.false., .false., &
             title='Predictors:')
        Select Case (ianal)
         Case Default
           Call prompt_results_file (fxfile, .true.)
         Case (ia_gcm)
           Call prompt_results_file (zifile, .true.)
        End Select
        If (ianal == ia_cca) Then
           Call prompt_results_file (fcfile, .false.)
        End If
        Call box_close (.false.)
     End If
     Call box_close (.true.)
     Call box_close (.false.)
     Call window_end ()
  End If
!
! Create climatologies property sheet
  Call init_psheet (ihp6, 'Climatologies')
  Call box_open (.false., .true.)
  Call box_open (.true., .false., &
       title='Climatologies')
  Call prompt_results_file (thfile, .true.)
  Call prompt_results_file (avfile, .false.)
  If (lcorr) Then
     Call prompt_results_file (cofile, .false.)
  End If
  Call box_close (.true.)
  Call box_close (.false.)
  Call window_end ()
!
#if GUI == 1
! Add menu bar
  Call init_win (title=Trim(l_outres%c), &
       lresize=.true.)
! - file menu -
  Call add_menu ('mn', 1,  -1, l_file)
  Call add_menu ('mn', 2,   0, l_close, &
       scbf='Exit')
! - edit menu -
  Call add_menu ('mn', 1,  -1, l_edit)
  Call add_menu ('mn', 2, 131, '&File formats', &
       cbf=get_fmt)
! - options menu -
  Call add_menu ('mn', 1,  -1, l_options)
  Call add_menu ('mn', 2, 133, 'Include CPT &tags', &
       icheck=itags, itoggle=itags)
!
! Get output files
  Call set_text_size (two)
  iw = winio$('%cn%bfOutput files%`bf&')
  Call set_text_size (one)
!
! Add property sheets
  Select Case (ianal)
   Case (ia_cca, ia_pcr) ! - CCA, PCR -
     If (ifc == 2) Then
        iw = winio$('%nl%nl%6ps&', ihp1, ihp6, ihp3, ihp4, ihp2, ihp5)
     Else
        iw = winio$('%nl%nl%5ps&', ihp1, ihp6, ihp3, ihp4, ihp2)
     End If
   Case (ia_mlr) ! - MLR -
     If (ifc == 2) Then
        iw = winio$('%nl%nl%5ps&', ihp1, ihp6, ihp4, ihp2, ihp5)
     Else
        iw = winio$('%nl%nl%4ps&', ihp1, ihp6, ihp4, ihp2)
     End If
   Case (ia_gcm) ! - GCM -
     If (igcms == 3) Then
        If (ifc == 2) Then
           iw = winio$('%nl%nl%5ps&', ihp1, ihp6, ihp4, ihp2, ihp5)
        Else
           iw = winio$('%nl%nl%4ps&', ihp1, ihp6, ihp4, ihp2)
        End If
     Else
        If (ifc == 2) Then
           iw = winio$('%nl%nl%4ps&', ihp1, ihp6, ihp2, ihp5)
        Else
           iw = winio$('%nl%nl%3ps&', ihp1, ihp6, ihp2)
        End If
     End If
   Case (ia_pfv) ! - PFV -
     If (igcms == 3) Then
        If (ifc == 2) Then
           iw = winio$('%nl%nl%4ps&', ihp1, ihp6, ihp4, ihp5)
        Else
           iw = winio$('%nl%nl%3ps&', ihp1, ihp6, ihp4)
        End If
     Else
        If (ifc == 2) Then
           iw = winio$('%nl%nl%3ps&', ihp1, ihp6, ihp5)
        Else
           iw = winio$('%nl%nl%2ps&', ihp1, ihp6)
        End If
     End If
  End Select
!
! Prompt
  iw = win_prompt(.true.)
#else
     Call print_text (' ')
     Write (Unit=*, Fmt='(A)') ' 0. '//Trim(l_exit%c)
     Call print_text (' ')
#endif
!
#if GUI == 1
! Save requested results
  Select Case (iw)
   Case (1)
     Call write_results ()
   Case Default
     Call files_reset ()
  End Select
!
! Restore console
  Call update_window$ (ih_con)
  save_data = 2
!
#else
! Identify file
  Call prompt_integer (' ', iw)
2 If (iw == 0) Then
     save_data = 0
     Return
  Else if (iw > 0 .and. iw <= Ubound(new_ids, Dim=1)) Then ! - convert old responses -
     iw = new_ids(iw)
     Call print_warning ('Old file identifier has been converted.', &
          nopause=.true.)
     GoTo 2
  Else If (iw == xofile%nfile) Then ! - X input data output file -
     Call geto_gen (xofile)
  Else If (iw == yofile%nfile) Then ! - Y input data output file -
     Call geto_gen (yofile)
  Else If (iw == xifile%nfile) Then ! - interpolated X output file -
     If (ianal == ia_gcm .or. ianal == ia_pfv) Call geto_gen (xifile)
  Else If (iw == ycfile%nfile) Then ! - Y categories output file -
     Call geto_gen (ycfile)
  Else If (iw == rcfile%nfile) Then ! - retroactive Y categories output file -
     If (iretro == 1) Call geto_gen (rcfile)
  Else If (iw == yhfile%nfile) Then ! - cross-validated predictions output file -
     If (ihind == 1) Call geto_gen (yhfile)
  Else If (iw == yrfile%nfile) Then ! - retroactive predictions output file -
     If (iretro == 1 .and. ihind == 1) Call geto_gen (yrfile)
  Else If (iw == rpfile%nfile) Then ! - retroactive forecast probabilities output file -
     If (iretro == 1 .and. ihind == 1) Call geto_gen (rpfile)
  Else If (iw == rlfile%nfile) Then ! - retroactive prediction limits output file -
     If (iretro == 1 .and. ihind == 1) Call geto_gen (rlfile)
  Else If (iw == xefile%nfile) Then ! - X eigenvalues output file -
     If (ianal == ia_cca .or. ianal == ia_pcr) Call geto_gen (xefile)
  Else If (iw == xlfile%nfile) Then ! - X spatial loadings output file -
     If (ianal == ia_cca .or. ianal == ia_pcr) Call geto_gen (xlfile)
  Else If (iw == xsfile%nfile) Then ! - X temporal scores output file -
     If (ianal == ia_cca .or. ianal == ia_pcr) Call geto_gen (xsfile)
  Else If (iw == yefile%nfile) Then ! - Y eigenvalues output file -
     If (ianal == ia_cca) Call geto_gen (yefile)
  Else If (iw == ylfile%nfile) Then ! - Y spatial loadings output file -
     If (ianal == ia_cca) Call geto_gen (ylfile)
  Else If (iw == ysfile%nfile) Then ! - Y temporal scores output file -
     If (ianal == ia_cca) Call geto_gen (ysfile)
  Else If (iw == ccfile%nfile) Then ! - canonical correlations output file -
     If (ianal == ia_cca) Call geto_gen (ccfile)
  Else If (iw == xmfile%nfile) Then ! - X homogeneous covariance maps output file -
     If (ianal == ia_cca) Call geto_gen (xmfile)
  Else If (iw == xtfile%nfile) Then ! - X homogeneous covariance maps time series output file -
     If (ianal == ia_cca) Call geto_gen (xtfile)
  Else If (iw == ymfile%nfile) Then ! - Y homogeneous covariance maps output file -
     If (ianal == ia_cca) Call geto_gen (ymfile)
  Else If (iw == ytfile%nfile) Then ! - Y homogeneous covariance maps time series output file -
     If (ianal == ia_cca) Call geto_gen (ytfile)
  Else If (iw == rbfile%nfile) Then ! - regression coefficients output file -
     If (ianal == ia_pcr .or. ianal == ia_mlr) Call geto_gen (rbfile)
  Else If (iw == pbfile%nfile) Then ! - PC regression coefficients output file -
     If (ianal == ia_pcr) Call geto_gen (pbfile)
  Else If (iw == fpfile%nfile) Then ! - forecast probabilities output file -
     If (ifc == 2) Call geto_gen (fpfile)
  Else If (iw == fofile%nfile) Then ! - forecast odds output file -
     If (ifc == 2) Call geto_gen (fofile)
  Else If (iw == fvfile%nfile) Then ! - forecasts output file -
     If (ifc == 2) Call geto_gen (fvfile)
  Else If (iw == fsfile%nfile) Then ! - forecast ensembles output file -
     If (ifc == 2) Call geto_gen (fsfile)
  Else If (iw == flfile%nfile) Then ! - prediction limits output file -
     If (ifc == 2) Call geto_gen (flfile)
  Else If (iw == fefile%nfile) Then ! - prediction error variances output file -
     If (ifc == 2) Call geto_gen (fefile)
  Else If (iw == fxfile%nfile) Then ! - predictor time scores output file -
     If (ifc == 2) Call geto_gen (fxfile)
  Else If (iw == zifile%nfile) Then ! - interpolated forecast predictors output file -
     If (ifc == 2) Call geto_gen (zifile)
  Else If (iw == fcfile%nfile) Then ! - predictor CCA time scores output file -
     If (ifc == 2) Call geto_gen (fcfile)
  Else If (iw == thfile%nfile) Then ! - thresholds output file -
     Call geto_gen (thfile)
  Else If (iw == avfile%nfile) Then ! - averages output file -
     Call geto_gen (avfile)
  Else If (iw == cofile%nfile) Then ! - correlations output file -
     If (lcorr) Call geto_gen (cofile)
  Else
     GoTo 1
  End If
!
! Save requested results
  Call write_results ()
  GoTo 1
!
#endif
  Return
 End Function save_data
!
!
!
 Function save_good()
!
! Prompts for goodness index output file
!
! Modules
  Use iofiles, Only: gifile, &
                     geto_gen
!
! Function Type
  Integer :: save_good
!
! Prompt for goodness index
  Call geto_gen (gifile)
  save_good = 2
!
  Return
 End Function save_good
!
!
!
#if GUI == 1
 Function save_results(outfile)
!
! Prompts for results files
!
! Modules
  Use clrwin$, Only: update_window$
  Use labels,  Only: l_close, l_edit, l_file, l_outres
  Use gui,     Only: ih_con, iw, &
                     init_win, prompt_file, win_prompt
  Use menus,   Only: add_menu
  Use iofiles, Only: ofile, &
                     nfile, &
                     files_reset, get_fmt, geto_file
!
! Function Type
  Integer :: save_results
!
! Arguments
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: outfile
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Open window
  Call init_win (title=Trim(l_outres%c), &
       lresize=.true.)
!
! Add menu bar
  Call add_menu ('mn', 1,  -1, l_file)
  Call add_menu ('mn', 2,   0, l_close, &
       scbf='EXIT')
  Call add_menu ('mn', 1,  -1, l_edit)
  Call add_menu ('mn', 2, 133, '&File format', &
       iset1=nfile, ival1=outfile%nfile, cbf=get_fmt)
!
! Get output file
  Call prompt_file (Trim(outfile%cfile%desc), outfile%cfile%fname, outfile%nfile, nfile, geto_file)
! - prompt -
  iw = win_prompt(.true.)
!
! Save requested results
  Select Case (iw)
   Case (2)
     Call write_results ()
   Case Default
     Call files_reset ()
  End Select
!
! Restore console
  Call update_window$ (ih_con)
  save_results = 2
!
  Return
 End Function save_results
!
!
!
#else
 Function save_results(outfile1, &
          outfile2)
!
! Prompts for results files
!
! Modules
  Use iofiles, Only: ofile, &
                     geto_gen
!
! Arguments
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: outfile1
!
! - Optional input/output scalars -
  Type(ofile), Intent(InOut), Optional :: outfile2
!
! Function Type
  Integer :: save_results
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Get output file
  Call geto_gen (outfile1)
  If (Present(outfile2)) Call geto_gen (outfile2)
!
! Save requested results
  Call write_results ()
  save_results = 0
!
  Return
 End Function save_results
!
!
!
 Function print_results (f_update, &
          f_calc)
!
! Prints results to the screen
!
! Modules
  Use data_io_constants, Only: ioutstd
  Use gui,               Only: print_text
  Use iofiles,           Only: yfile
  Use settings,          Only: ivfa
  Use fields,            Only: iffy, yfield, &
                               check_ivf, print_field
!
! Function Type
  Integer :: print_results
!
! Arguments
!
! Procedure arguments
  Integer, External :: f_update ! - update function -
!
! - Optional procedure arguments -
  Integer, External, Optional :: f_calc ! - calculation function -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Repeat for each field and variable
  Do iffy = 1, yfile%nfl
     Call print_field (ioutstd, yfield(iffy), yfile%nfs, yfile%nls, yfile%nms)
     Call print_text (' ')
     Do ivfa = 1, yfield(iffy)%nva
        print_results = check_ivf(iffy, ivfa)
        If (Present(f_calc)) print_results = f_calc()
        print_results = f_update()
     End Do
  End Do
  iffy = 1
  print_results = 0
!
  Return
 End Function print_results
!
!
!
#endif
 Function update_results(ioutr, ih_r, ir, chead, lclim, print_res) &
          Result (update)
!
! Prints updated results
!
! Modules
#if GUI == 1
  Use clrwin$,            Only: clear_window$
#endif
  Use data_gui_constants, Only: i_handle
  Use gui,                Only: print_text, set_cursor_waiting
  Use arrays,             Only: iuse=>iusey
  Use iofiles,            Only: yfile
  Use space,              Only: coor, rlaty, rlngdy, &
                                make_map_coor
  Use fields,             Only: iffy, yfield
  Use climate,            Only: print_climatology
!
! Function Type
  Integer :: update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ioutr ! - output unit number -
  Integer, Intent(In) :: ir    ! - array index -
!
  Integer(Kind=i_handle), Intent(In) :: ih_r ! - handle -
!
  Character(Len=*), Intent(In) :: chead ! - heading -
!
  Logical, Intent(In) :: lclim ! - print cliamtology -
!
! Procedure arguments
  Integer, External :: print_res
!
! Locals
!
! Local scalars
  Integer :: i       ! - latitude index -
  Integer :: j       ! - longitude index -
  Integer :: ij, ija ! - indices -
  Integer :: ilat    ! - latitude index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Sum
!
! Executable Statements
!
#if GUI == 1
! Clear window
  Call clear_window$ (ih_r)
!
#endif
! Print climatology
  If (lclim) Then
#if GUI == 1
     Call print_climatology (.false., &
          iout=ioutr)
#else
     Call print_climatology (.true., &
          iout=ioutr)
#endif
  End If
!
! Calculate offsets
  If (iffy == 1) Then
     ij = 0
     ija = 1
  Else
     ij = Sum(yfield(1:iffy-1)%nv)
     ija = Sum(yfield(1:iffy-1)%nva) + 1
  End If
!
! Print results
! - gridded data -
  Select Case (yfile%idstr)
   Case (1)
     Call print_text (chead, &
          iout=ioutr)
     Call print_text (' ', &
          iout=ioutr)
     Do i = 1, yfield(iffy)%region%nlts
        ilat = yfield(iffy)%region%nlt1 + i - 1
        Do j = 1, yfield(iffy)%region%nlgs
           Call set_cursor_waiting (1)
           ij = ij + 1
           If (iuse(ija) == ij) Then
              coor = make_map_coor(rlaty(ilat,iffy), rlngdy(j,iffy))
              update = print_res(ioutr, ija, ir, 0, iffy)
              ija = ija + 1
           End If
           Call set_cursor_waiting (0)
        End Do
     End Do
! - station and unreferenced data -
   Case (2, 3)
     Call print_text (chead, &
          iout=ioutr)
     Call print_text (' ', &
          iout=ioutr)
     Do j = 1, yfield(iffy)%nva
        Call set_cursor_waiting (1)
        update = print_res(ioutr, ija, ir, iuse(ija)-ij, iffy)
        ija = ija + 1
        Call set_cursor_waiting (0)
     End Do
  End Select
  update=0
!
  Return
 End Function update_results
!
!
!
 Subroutine write_results ()
!
! Writes output to files
!
! Modules
  Use data_numbers,        Only: one
  Use data_cpt_constants,  Only: ia_gcm, nep, ng, npall, nts, nwr
  Use data_io_constants,   Only: lstn
  Use data_time_constants, Only: isq_mn
  Use labels,              Only: ca_pvs_t, ca_scores, &
                                 l_correls, &
                                 ls_cats, ls_wrlts
  Use time,                Only: Assignment(=), Operator(+), &
                                 iseq, nss, &
                                 date_diff
  Use progress_meter,      Only: end_progress_meter, set_progress_increment, start_progress_meter, stop_progress_meter
  Use arrays,              Only: b, b0, c, e, f, fcast, fpls, iusex, iusey, lwk, pev, tex, x, xiny, y, yhat, yret, yrpls, ziny
  Use iofiles,             Only: ifile, &
                                 mxa, mya, nx, ny, nz, &
                                 atfile, avfile, ccfile, cofile, exfile, fcfile, fefile, flfile, fofile, fpfile, fsfile, fvfile, &
                                 fxfile, pbfile, psfile, pvfile, rbfile, rcfile, rhfile, rlfile, rofile, rpfile, rrfile, skfile, &
                                 thfile, wrfile, xefile, xifile, xlfile, xmfile, xofile, xsfile, xtfile, x2file, ycfile, yefile, &
                                 yhfile, ylfile, ymfile, yofile, yrfile, ysfile, ytfile, zifile, &
                                 mfile, xfile, yfile, zfile, &
                                 files_reset
  Use settings,            Only: clf, dofr, iexc, ig0, iodds, nb, nenf, nf, nt, ntr, nt1, nur
  Use space,               Only: cstndx, cstndy, rlatdx, rlatdy, rlngdx, rlngdy
  Use fields,              Only: sfield, xfield, yfield, zfield
  Use scores,              Only: afp, b0s, b1s, corrs, far, fars, hit, hits, ifq, ips, iskill, orf, pscores, pvalues, rhit, roca, &
                                 rocas, skills, tends, wrlts
  Use season,              Only: lb, lprst, nsn
  Use pcs,                 Only: bz, eofx, eofy, fsx, fsxc, hx_map, hy_map, hx_ser, hy_ser, lxt, lyt, mcc, mxe, mye, mu, &
                                 nco, svx, svy, tsx, tsy
  Use missing,             Only: kfuse, ktuse
  Use categories,          Only: fps, icatc, iobs, irobs, rfps, odds, tobs
  Use climate,             Only: ave, climate_per
  Use analysis,            Only: ianal
!
! Locals
!
! Local scalars
  Integer :: i     ! - index -
  Integer :: it0_f ! - time 0 for forecasts -
  Integer :: it0_r ! - time 0 for retractive predictands -
  Integer :: it0_x ! - time 0 for predictors -
  Integer :: it0_y ! - time 0 for predictands -
  Integer :: it0_z ! - time 0 for predictors in forecast period -
  Integer :: ifail ! - error indicator -
  Integer :: nsave ! - number of files to save -
!
  Character(Len= 1) :: ctype ! - output Type -
  Character(Len=64) :: cvar  ! - field name -
!
  Logical :: lout ! - output? -
!
  Type(ifile) :: dfile ! - descriptor file -
!
! Local arrays
  Integer, Dimension(nwr) :: iwk ! - workspace -
!
  Character(Len=lstn), Dimension(:,:), Allocatable :: cvars ! - variable names -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Nint, Real, Reshape, Transpose, Trim
!
! Executable Statements
!
! Update console
  Call start_progress_meter (' ', .true.)
! - count number of files to save -
  nsave = 0
  If (xofile%cfile%lset) nsave = nsave + 1
  If (yofile%cfile%lset) nsave = nsave + 1
  If (xifile%cfile%lset) nsave = nsave + 1
  If (ycfile%cfile%lset) nsave = nsave + 1
  If (rcfile%cfile%lset) nsave = nsave + 1
  If (x2file%cfile%lset) nsave = nsave + 1
  If (yhfile%cfile%lset) nsave = nsave + 1
  If (yrfile%cfile%lset) nsave = nsave + 1
  If (rpfile%cfile%lset) nsave = nsave + 1
  If (rlfile%cfile%lset) nsave = nsave + 1
  If (xefile%cfile%lset) nsave = nsave + 1
  If (yefile%cfile%lset) nsave = nsave + 1
  If (xlfile%cfile%lset) nsave = nsave + 1
  If (ylfile%cfile%lset) nsave = nsave + 1
  If (xsfile%cfile%lset) nsave = nsave + 1
  If (ysfile%cfile%lset) nsave = nsave + 1
  If (ccfile%cfile%lset) nsave = nsave + 1
  If (xmfile%cfile%lset) nsave = nsave + 1
  If (ymfile%cfile%lset) nsave = nsave + 1
  If (xtfile%cfile%lset) nsave = nsave + 1
  If (ytfile%cfile%lset) nsave = nsave + 1
  If (rbfile%cfile%lset) nsave = nsave + 1
  If (pbfile%cfile%lset) nsave = nsave + 1
  If (fpfile%cfile%lset) nsave = nsave + 1
  If (fofile%cfile%lset) nsave = nsave + 1
  If (fvfile%cfile%lset) nsave = nsave + 1
  If (fsfile%cfile%lset) nsave = nsave + 1
  If (flfile%cfile%lset) nsave = nsave + 1
  If (fefile%cfile%lset) nsave = nsave + 1
  If (exfile%cfile%lset) nsave = nsave + 1
  If (fxfile%cfile%lset) nsave = nsave + 1
  If (zifile%cfile%lset) nsave = nsave + 1
  If (fcfile%cfile%lset) nsave = nsave + 1
  If (thfile%cfile%lset) nsave = nsave + 1
  If (avfile%cfile%lset) nsave = nsave + 1
  If (cofile%cfile%lset) nsave = nsave + 1
  If (skfile%cfile%lset) nsave = nsave + 1
  If (pvfile%cfile%lset) nsave = nsave + 1
  If (rofile%cfile%lset) nsave = nsave + 1
  If (psfile%cfile%lset) nsave = nsave + 1
  If (rrfile%cfile%lset) nsave = nsave + 1
  If (atfile%cfile%lset) nsave = nsave + 1
  If (rhfile%cfile%lset) nsave = nsave + 1
  If (wrfile%cfile%lset) nsave = nsave + 1
  If (nsave == 0) Then
     Print *, 'No files to save.'
     Call files_reset ()
     Call stop_progress_meter (.true.)
     Return
  End If
  Call set_progress_increment (one/Real(nsave, Kind=rp), .false.)
!
! Calculate time-0
  it0_x = xfile%it1 - 1
  it0_y = yfile%it1 - 1
  it0_z = zfile%it1 - 1
  it0_r = yfile%it1 + nt1 - 1
  it0_f = yfile%it1 + date_diff(xfile%fdate, yfile%fdate, iseq) - date_diff(zfile%fdate, yfile%fdate, iseq) - 1
!
! Save input data with missing values estimated
! - X input data -
  If (xofile%cfile%lset) Then
     If (xfile%ngs == 1) Then
        Call write_cptdata (.true., .true., 'T', xofile, xfile, xfield(1:), xfield(1)%tprd%sdate, nt, nx, &
             x(:,:,1), iusex(:), ktuse(:), rlatdx(:,:), rlngdx(:,:), cstndx(:,:), ifail, &
             it0=it0_x)
     Else
        Call write_cptprobs (xofile, xfile, xfield(1:), xfield(1)%tprd%sdate, it0_x, nt, nx, &
             x(:,:,:), iusex(:), ktuse(:), rlatdx(:,:), rlngdx(:,:), cstndx(:,:), ifail)
     End If
     If (ifail < 0) GoTo 1
  End If
! - interpolated X data -
  If (xifile%cfile%lset) Then
     Allocate (sfield(xfile%nfl))
     sfield(:) = xfield(:)
     sfield(:)%nlt = yfield(:)%nlt
     sfield(:)%nlg = yfield(:)%nlg
     sfield(:)%region = yfield(:)%region
     If (xfile%ngs == 1) Then
        Call write_cptdata (.true., .true., 'T', xifile, xfile, sfield(1:), xfield(1)%tprd%sdate, nt, ny, &
             xiny(:,:,ig0), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
             it0=it0_x)
     Else
        dfile = xfile
        dfile%nfs = yfile%nfs
        dfile%nls = yfile%nls
        Call write_cptprobs (xifile, dfile, sfield(1:), xfield(1)%tprd%sdate, it0_x, nt, nx, &
             rfps(:,:,:), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     End If
     Deallocate (sfield)
     If (ifail < 0) GoTo 1
  End If
! - Y input data -
  If (yofile%cfile%lset) Then
     Call write_cptdata (.true., .true., 'T', yofile, yfile, yfield(1:), yfield(lb)%tprd%sdate, nt, ny, &
          y(:,:,1), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_y)
     If (ifail < 0) GoTo 1
  End If
! - Y categories -
  If (ycfile%cfile%lset) Then
     Call write_cptdata (.true., .true., 'T', ycfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, nt, ny, &
          Real(iobs(:,:), Kind=rp), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_y)
     If (ifail < 0) GoTo 1
  End If
! - retroactive categories -
  If (rcfile%cfile%lset) Then
     Call write_cptdata (.true., .true., 'T', rcfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, ntr, ny,    &
          Real(irobs(:,:), Kind=rp), iusey(:), ktuse(nt1+1:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_r)
     If (ifail < 0) GoTo 1
  End If
! - merged input data -
  If (x2file%cfile%lset) Then
     Call write_cptdata (.true., .false., 'T', x2file, xfile, xfield(1:), xfield(1)%tprd%sdate, nt, nx, &
          x(:,:,1), iusex(:), ktuse(:), rlatdx(:,:), rlngdx(:,:), cstndx(:,:), ifail, &
          it0=it0_x, nofs=xfile%nfs+yfile%nfs)
     If (ifail < 0) GoTo 1
     Call write_cptdata (.false., .true., 'T', x2file, yfile, yfield(1:), yfield(1)%tprd%sdate, nt, ny, &
          y(:,:,1), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_y)
     If (ifail < 0) GoTo 1
  End If
!
! Save hindcasts
! - cross-validated predictions -
  If (yhfile%cfile%lset) Then
     Call write_cptdata (.true., .true., 'T', yhfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, nt, ny, &
          yhat(:,:), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_y)
     If (ifail < 0) GoTo 1
  End If
! - retroactive predictions -
  If (yrfile%cfile%lset) Then
     Call write_cptdata (.true., .true., 'T', yrfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, ntr, ny, &
          yret(:,:), iusey(:), ktuse(nt1+1:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_r)
     If (ifail < 0) GoTo 1
  End If
! - retroactive forecast probabilities -
  If (rpfile%cfile%lset) Then
     Call write_cptprobs (rpfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, it0_r, ntr, ny, &
          rfps(:,:,:), iusey(:), ktuse(nt1+1:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - retroactive prediction limits -
  If (rlfile%cfile%lset) Then
     Call write_cptdata3 ('L', rlfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, it0_r, ntr, ny, 2,  &
          Reshape((/yrpls(:,:)%lower, yrpls(:,:)%upper/), (/mya, nur, 2/)), iusey(:), ktuse(nt1+1:), &
          rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          nofs=yfile%nfs*2, cl=clf)
     If (ifail < 0) GoTo 1
  End If
!
! Save PCA and CCA results
! - eigenvalues -
  If (xefile%cfile%lset) Then ! - X eigenvalues -
     Call write_cpteigs (xefile, lxt, nt, svx, ifail)
     If (ifail < 0) GoTo 1
  End If
  If (yefile%cfile%lset) Then ! - Y eigenvalues -
     Call write_cpteigs (yefile, lyt, nt, svy, ifail)
     If (ifail < 0) GoTo 1
  End If
! - spatial loadings -
  If (xlfile%cfile%lset) Then ! - X spatial loadings -
     Call write_cptload ('X_loadings', xlfile, xfile, xfield(:), mxe, nx, eofx(:,:), &
          iusex(:), rlatdx(:,:), rlngdx(:,:), cstndx(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
  If (ylfile%cfile%lset) Then ! - Y spatial loadings -
     Call write_cptload ('Y_loadings', ylfile, yfile, yfield(1:), mye, ny, eofy(:,:), &
          iusey(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - temporal scores -
  If (xsfile%cfile%lset) Then ! - X scores -
     Call write_cptscor ('X scores', xsfile, xfile%prd1, it0_x, nt, mxe, xfield(1)%rmiss, tsx(:,:), ktuse(:), ifail)
     If (ifail < 0) GoTo 1
  End If
  If (ysfile%cfile%lset) Then ! - Y scores -
     Call write_cptscor ('Y scores', ysfile, yfile%prd1, it0_y, nt, mye, yfield(1)%rmiss, tsy(:,:), ktuse(:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - canonical correlations -
  If (ccfile%cfile%lset) Then
     dfile%idstr = 3
     dfile%nfl = 1
     dfile%nfs = 1
     dfile%nls = 1
     Allocate (sfield(1))
     sfield(1)%nv = 1
     sfield(1)%nva = 1
     sfield(1)%var = 'correlation'
     sfield(1)%unitc = 'N/A'
     sfield(1)%rmiss = -999.0_rp
     sfield(1)%tprd%sdate%iyr = 1
     sfield(1)%tprd%sdate%imn = 0
     sfield(1)%tprd%sdate%idy = 0
     sfield(1)%tprd%edate = sfield(1)%tprd%sdate
     Allocate (lwk(mcc))
     lwk(:) = .true.
     Allocate (cvars(1,1))
     cvars(:,:) = 'correlation'
     Call write_cptdata (.true., .true., 'Mode', ccfile, dfile, sfield(:), sfield(1)%tprd%sdate, mcc, 1, &
          Reshape(mu(1:mcc), (/1,mcc/)), (/1/), lwk, rlatdy(:,:), rlatdy(:,:), cvars(:,:), ifail)
     Deallocate (cvars)
     Deallocate (lwk)
     Deallocate (sfield)
     If (ifail < 0) GoTo 1
  End If
! - X homogeneous maps -
  If (xmfile%cfile%lset) Then
     Call write_cptload ('X_CCA_loadings', xmfile, xfile, xfield(:), nco, nx, &
          hx_map(:,:), iusex(:), rlatdx(:,:), rlngdx(:,:), cstndx(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - Y homogeneous maps -
  If (ymfile%cfile%lset) Then
     Call write_cptload ('Y_CCA_loadings', ymfile, yfile, yfield(1:), nco, ny, &
          hy_map(:,:), iusey(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - X homogeneous map series -
  If (xtfile%cfile%lset) Then
     Call write_cptscor ('X_CCA_scores', xtfile, xfile%prd1, it0_x, nt, nco, xfield(1)%rmiss, hx_ser(:,:), ktuse(:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - Y homogeneous map series -
  If (ytfile%cfile%lset) Then
     Call write_cptscor ('Y_CCA_scores', ytfile, yfile%prd1, it0_y, nt, nco, yfield(1)%rmiss, hy_ser(:,:), ktuse(:), ifail)
     If (ifail < 0) GoTo 1
  End If
!
! Save regression coefficients
! - regression coefficients -
  If (rbfile%cfile%lset) Then
     Select Case (ianal)
      Case Default
        Call write_cptregr (rbfile, b(:,:), b0(:,:), ifail)
      Case (ia_gcm)
        Call write_cptload ('constants', rbfile, yfile, yfield(1:), 1, ny, b0(:,:), iusey, rlatdy, rlngdy, cstndy, ifail, &
             cload2='coefficients', v2=b(:,:))
     End Select
     If (ifail < 0) GoTo 1
  End If
! - PC regression coefficients -
  If (pbfile%cfile%lset) Then
     Call write_cptload ('coefficients', pbfile, yfile, yfield(1:), mxe, ny, Transpose(bz), iusey, rlatdy, rlngdy, cstndy, ifail)
     If (ifail < 0) GoTo 1
  End If
!
! Save forecasts
! - forecast probabilities -
  If (fpfile%cfile%lset) Then
     Call write_cptprobs (fpfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, it0_f, nf, ny, &
          fps(:,:,:), iusey(:), kfuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - forecast odds -
  If (fofile%cfile%lset) Then
     Select Case (iodds)
      Case (0)
        ctype='O'
      Case (1)
        ctype='R'
     End Select
     Call write_cptdata3 (ctype, fofile, yfile, yfield(1:), yfield(lb)%tprd%sdate, it0_f, nf, ny, ng, &
          odds(:,:,:,iodds+1), iusey(:), kfuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          nofs=yfile%nfs*ng)
     If (ifail < 0) GoTo 1
  End If
! - deterministic forecasts -
  If (fvfile%cfile%lset) Then
     Call write_cptdata (.true., .true., 'T', fvfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, nf, ny, &
          fcast(:,:,0), iusey(:), kfuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_f)
     If (ifail < 0) GoTo 1
  End If
! - forecast ensembles -
  If (fsfile%cfile%lset) Then
     Call write_cptdata3 ('E', fsfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, it0_f, nf, ny, nenf, &
          fcast(:,:,1:), iusey(:), kfuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          nofs=yfile%nfs*nenf)
     If (ifail < 0) GoTo 1
  End If
! - prediction limits -
  If (flfile%cfile%lset) Then
     Call write_cptdata3 ('L', flfile, yfile, yfield(1:), yfield(lb)%tprd%sdate, it0_f, nf, ny, 2, &
          Reshape((/fpls(:,:)%lower, fpls(:,:)%upper/), (/mya, nf, 2/)), iusey(:), kfuse(:),       &
          rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          nofs=yfile%nfs*2, cl=clf)
     If (ifail < 0) GoTo 1
  End If
! - prediction error variances -
  If (fefile%cfile%lset) Then
     Allocate (sfield(yfile%nfl))
     sfield(:) = yfield(1:)
     Do i = 1, yfile%nfl
        sfield(i)%var = Trim(sfield(i)%var)//' pred_err_var'
        If (Trim(sfield(i)%unitc) /= 'N/A') sfield(i)%unitc = Trim(sfield(i)%unitc)//' (squared)'
     End Do
     Call write_cptdata (.true., .true.,'T', fefile, yfile, sfield(1:), yfield(lb)%tprd%sdate, nf, ny, &
          pev(:,:)**2, iusey(:), kfuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
          it0=it0_f, ctagi='dof', itag=Nint(dofr))
     Deallocate (sfield)
     If (ifail < 0) GoTo 1
  End If
! - exceedence probabilities -
  If (exfile%cfile%lset) Then
     Call write_cptexcs (exfile, iexc, nep, c, e, f, tex(1), tex(2), ifail)
     If (ifail < 0) GoTo 1
  End If
! - predictor time scores -
  If (fxfile%cfile%lset) Then
     Call write_cptscor ('X_scores', fxfile, zfile%prd1, it0_z, nf, mxe, xfield(1)%rmiss, fsx(:,:), kfuse(:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - interpolated Z data -
  If (zifile%cfile%lset) Then
     Allocate (sfield(zfile%nfl))
     sfield(:) = zfield(:)
     sfield(:)%nlt = yfield(:)%nlt
     sfield(:)%nlg = yfield(:)%nlg
     sfield(:)%region = yfield(:)%region
     If (zfile%ngs == 1) Then
        Call write_cptdata (.true., .true., 'T', zifile, zfile, sfield(1:), zfield(1)%tprd%sdate, nf, ny, &
             ziny(:,:,ig0), iusey(:), kfuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail, &
             it0=it0_z)
     Else
        dfile = zfile
        dfile%nfs = yfile%nfs
        dfile%nls = yfile%nls
        Call write_cptprobs (zifile, dfile, sfield(1:), zfield(1)%tprd%sdate, it0_z, nf, nz, &
             rfps(:,:,:), iusey(:), ktuse(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     End If
     Deallocate (sfield)
     If (ifail < 0) GoTo 1
  End If
! - predictor CCA time scores -
  If (fcfile%cfile%lset) Then
     Call write_cptscor ('X_CCA_scores', fcfile, zfile%prd1, it0_z, nf, nco, xfield(1)%rmiss, fsxc(:,:), kfuse(:), ifail)
     If (ifail < 0) GoTo 1
  End If
!
! Save climatologies
  lout = (thfile%cfile%lset .or. avfile%cfile%lset)
  If (lout) Then
     dfile = yfile
     If (mfile%iseq == isq_mn .and. nss > 1) Then
        dfile%nfl = dfile%nfl*nss
     Else
        dfile%nls = dfile%nls/nss
     End If
     Allocate (sfield(dfile%nfl))
     sfield(:) = yfield(:)
! - set dates to climatological dates -
     sfield(:)%tprd%sdate%iyr = climate_per%d1%iyr
     sfield(:)%tprd%edate%iyr = climate_per%d2%iyr
     If (dfile%nfl > 1) Then
        Do i = 2, dfile%nfl
           If (yfield(i)%tprd%sdate%imn < yfield(i-1)%tprd%sdate%imn) sfield(i:)%tprd%sdate%iyr = sfield(i:)%tprd%sdate%iyr + 1
        End Do
        Do i = dfile%nfl-1, 1, -1
           If (yfield(i)%tprd%sdate%imn > yfield(i+1)%tprd%sdate%imn) sfield(:i)%tprd%sdate%iyr = sfield(:i)%tprd%sdate%iyr - 1
        End Do
     End If
     Allocate (lwk(nsn))
     lwk(:) = .true.
  End If
! - thresholds -
  If (thfile%cfile%lset) Then
     Call write_cptdata3 ('T', thfile, dfile, sfield(1:), sfield(lb)%tprd%sdate, 0, nsn, ny, nts, &
          Reshape(tobs, (/mya, nsn, nts/)), iusey(:), lwk(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - averages -
  If (avfile%cfile%lset) Then
     Do i = 1, nss
        sfield(i)%var = Trim(sfield(i)%var)//' '//Trim(sfield(i)%cssn)//' average'
        Call write_cptdata ((i==1), (i==nss), 'T', avfile, dfile, sfield(1+i-1:), sfield(lb+i-1)%tprd%sdate, 1, ny, &
             ave(:,i:), iusey(:), lwk(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
        If (ifail < 0) GoTo 1
     End Do
  End If
  If (lout) Then
     Deallocate (lwk)
     Deallocate (sfield)
  End If
! - correlations -
  If (cofile%cfile%lset) Then
     Allocate (sfield(xfile%nfs))
     sfield(:) = xfield(:)
     sfield(:)%tprd%sdate%iyr = 1
     sfield(:)%tprd%sdate%imn = 0
     sfield(:)%tprd%sdate%idy = 0
     sfield(:)%tprd%edate = sfield(:)%tprd%sdate
     Call write_cptdata (.true., .true., Trim(l_correls%c), cofile, xfile, sfield(:), yfield(lb)%tprd%sdate, 1, nx, &
          Reshape(corrs, (/mxa,1/)), iusex(:), (/.true./), rlatdx(:,:), rlngdx(:,:), cstndx(:,:), ifail)
     Deallocate (sfield)
     If (ifail < 0) GoTo 1
  End If
!
! Save scores
  lout=(skfile%cfile%lset .or. pvfile%cfile%lset .or. psfile%cfile%lset)
  If (lout) Then
     Allocate (sfield(yfile%nfs))
     sfield(:) = yfield(1:yfile%nfs)
     If (lprst) sfield(1:)%tprd%sdate = yfield(0)%tprd%sdate
     sfield(:)%unit  = 'unitless'
     sfield(:)%unitc = 'unitless'
  End If
! - deterministic skill scores -
  If (skfile%cfile%lset) Then
     Call write_cptload (Trim(ca_scores(iskill)), skfile, yfile, sfield(:), 1, ny, &
          Reshape(skills, (/mya, 1/)), iusey(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - p-values -
  If (pvfile%cfile%lset) Then
     Call write_cptload ('p-values', pvfile, yfile, sfield(:), 1, ny, &
          Reshape(pvalues,(/mya, 1/)), iusey(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - ROC results -
  If (rofile%cfile%lset) Then ! - individual point -
     Call write_cptroc (rofile, nt, 2, hit, far, roca, ifail)
     If (ifail < 0) GoTo 1
  End If
! - probabilistic scores -
  If (psfile%cfile%lset) Then
     cvar = ca_pvs_t(ips)
     If (ips > npall) cvar = Trim(cvar)//' ('//Trim(ls_cats(icatc)%c)//')'
     Call write_cptload (Trim(cvar), psfile, yfile, sfield(:), 1, ny, &
          Reshape(pscores,(/mya, 1/)), iusey(:), rlatdy(:,:), rlngdy(:,:), cstndy(:,:), ifail)
     If (ifail < 0) GoTo 1
  End If
  If (lout) Deallocate (sfield)
! - ROC results -
  If (rrfile%cfile%lset) Then ! - all points -
     Call write_cptroc (rrfile, nb, ng, hits, fars, rocas, ifail)
     If (ifail < 0) GoTo 1
  End If
! - attributes diagram -
  If (atfile%cfile%lset) Then
     Call write_cptrel (atfile, ng, nb, afp, orf, ifq, tends, b0s, b1s, ifail)
     If (ifail < 0) GoTo 1
  End If
! - ranked-hits -
  If (rhfile%cfile%lset) Then
     Call write_cptscor ('Ranked_hits', rhfile, yfile%prd1+nt1, it0_r, ntr, ng, yfield(1)%rmiss, Transpose(rhit), &
          ktuse(nt1+1:), ifail)
     If (ifail < 0) GoTo 1
  End If
! - weather roulette results -
  If (wrfile%cfile%lset) Then
     dfile%idstr = 3
     dfile%nfl = 1
     dfile%nfs = 1
     dfile%nls = 1
     Allocate (sfield(1))
     sfield(1) = yfield(1)
     sfield(1)%nv = nwr
     sfield(1)%nva = nwr
     sfield(1)%var = 'weather_roulette'
     sfield(1)%unitc = 'N/A'
     sfield(1)%tprd = yfile%prd1 + (yfile%it1 - 2)
     Do i = 1, nwr
        iwk(i) = i
     End Do
     Allocate (lwk(ntr+1))
     lwk(1) = .true.
     lwk(2:ntr+1) = ktuse(nt1+1:)
     Allocate (cvars(nwr,1))
     cvars(:,1) = ls_wrlts(:)%c
     Call write_cptdata (.true., .true., 'T', wrfile, dfile, sfield(:), yfield(lb)%tprd%sdate+(yfile%it1-2), ntr+1, nwr, &
          wrlts(:,0:), iwk(:), lwk(:), rlatdy(:,:), rlngdy(:,:), cvars(:,:), ifail)
     Deallocate (lwk)
     Deallocate (sfield)
     If (ifail < 0) GoTo 1
  End If
!
! Update console
1 Call end_progress_meter (.true., (.true. .and. ianal /= 0))
!
  Return
!
 End Subroutine write_results
!
!
!
 Subroutine write_cptdata (lopen, lclose, ccol, afile, dfile, afield, sdate, nt, nv, v, iuse, kuse, rlat, rlng, cstn, ifail, &
            it0, nofs, itag, ctagi)
!
! Outputs data in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_io_constants, Only: iout, lstn
  Use time,              Only: pdate, &
                               Operator(/=)
  Use iofiles,           Only: ifile, ofile
  Use fields,            Only: field, &
                               tfield
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: nv  ! - number of variables -
!
  Logical, Intent(In) :: lopen  ! - open file? -
  Logical, Intent(In) :: lclose ! - close file? -
!
  Character(Len=*), Intent(In) :: ccol ! - column indicator -
!
  Type(ifile), Intent(In) :: dfile ! - descriptor file -
!
  Type(pdate), Intent(In) :: sdate ! - start date -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: it0  ! - time 0 -
  Integer, Intent(In), Optional :: nofs ! - number of output fields -
  Integer, Intent(In), Optional :: itag ! - extra integer tag value -
!
  Character(Len=*), Intent(In), Optional :: ctagi ! - extra integer tag -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v    ! - data -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Size
!
! Executable Statements
!
! Open output file
  ifail = 0
  If (lopen) Then
     If (Present(nofs)) Then
        Call open_output (iout, afile, nofs, dfile%idstr, nv, ifail)
     Else
        Call open_output (iout, afile, dfile%nfs, dfile%idstr, nv, ifail)
     End If
  End If
  If (ifail == 0) Then
!
! Set field information
     Allocate (tfield(Size(afield, Dim=1)))
     tfield(:) = afield(:)
     If (sdate /= tfield(1)%tprd%sdate) tfield(:)%tprd%sdate = sdate ! - set start date for possible persistence component -
!
! Print data
     Select Case (dfile%idstr)
      Case (1) ! - gridded -
        Call write_grid (iout, afile, tfield(:), dfile%nfs, dfile%nls, nt, nv, &
             v(:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), ifail, &
             it0=it0, itag=itag, ctagi=ctagi)
      Case (2) ! - station -
        If (Present(it0)) Then
           Call write_stns (iout, afile, tfield(:), dfile%nfs, dfile%nls, it0, nt, nv, &
                v(:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), cstn(:,:), ifail, &
                itag=itag, ctagi=ctagi)
        Else
           Call write_stns (iout, afile, tfield(:), dfile%nfs, dfile%nls, 0, nt, nv, &
                v(:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), cstn(:,:), ifail, &
                itag=itag, ctagi=ctagi)
        End If
      Case (3) ! - unreferenced -
        If (Present(it0)) Then
           Call write_unrf (iout, afile, tfield(:), dfile%nfs, dfile%nls, it0, nt, nv, &
                v(:,:), iuse(:), kuse(:), cstn(:,:), ccol, ifail, &
                itag=itag, ctagi=ctagi)
        Else
           Call write_unrf (iout, afile, tfield(:), dfile%nfs, dfile%nls, 0, nt, nv, &
                v(:,:), iuse(:), kuse(:), cstn(:,:), ccol, ifail, &
                itag=itag, ctagi=ctagi)
        End If
     End Select
     Deallocate (tfield)
  End If
!
! Close file and reset
  If (lclose) Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptdata
!
!
!
 Subroutine write_cptdata3 (ctype, afile, dfile, afield, sdate, it0, nt, nv, n3, v, iuse, kuse, rlat, rlng, cstn, ifail, &
            nofs, cl)
!
! Outputs 3-Dimensional data in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_io_constants, Only: iout, lstn
  Use time,              Only: pdate
  Use iofiles,           Only: ifile, ofile
  Use fields,            Only: field, &
                               tfield
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: it0 ! - time 0 -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: nv  ! - number of stations -
  Integer, Intent(In) :: n3  ! - third Dimension -
!
  Character(Len=1), Intent(In) :: ctype ! - output Type -
!
  Type(ifile), Intent(In) :: dfile ! - descriptor file -
!
  Type(pdate), Intent(In) :: sdate ! - start date -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: nofs ! - number of output fields -
!
  Real(Kind=rp), Intent(In), Optional :: cl ! - confidence level -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: v    ! - data -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Open output file
  ifail = 0
  If (Present(nofs)) Then
     Call open_output (iout, afile, nofs, dfile%idstr, nv, ifail)
  Else
     Call open_output (iout, afile, dfile%nfs, dfile%idstr, nv, ifail)
  End If
  If (ifail == 0) Then
!
! Set field information
     Allocate (tfield(dfile%nfl))
     tfield(:) = afield(:)
     tfield(:)%tprd%sdate = sdate ! - set start date for possible persistence component -
!
! Print data
     Select Case (dfile%idstr)
      Case (1)
        Call write_grids (ctype, iout, afile, tfield(1:), dfile%nfs, dfile%nls, it0, nt, nv, n3, &
             v(:,:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), ifail, &
             cl=cl)
      Case (2)
        Call write_stnss (ctype, iout, afile, tfield(1:), dfile%nfs, dfile%nls, it0, nt, nv, n3, &
             v(:,:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), cstn(:,:), ifail, &
             cl=cl)
      Case (3)
        Call write_unrfs (ctype, iout, afile, tfield(1:), dfile%nfs, dfile%nls, it0, nt, nv, n3, &
             v(:,:,:), iuse(:), kuse(:), cstn(:,:), ifail, &
             cl=cl)
     End Select
     Deallocate (tfield)
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptdata3
!
!
!
 Subroutine write_cptscor (cscore, afile, prd1, it0, nt, nv, rmiss, v, kuse, ifail)
!
! Outputs temporal scores in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_io_constants, Only: iout
  Use time,              Only: pprd
  Use iofiles,           Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: it0 ! - time 0 -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: nv  ! - number of grids -
!
  Real(Kind=rp), Intent(In) :: rmiss ! - missing values -
!
  Character(Len=*), Intent(In) :: cscore ! - score description -
!
  Type(pprd), Intent(In) :: prd1 ! - first period -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v ! - data -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, 1, 3, nt, ifail)
!
! Print data
  If (ifail == 0) Call write_scor (iout, cscore, afile, nt, it0, prd1, nv, kuse(:), rmiss, v(:,:), ifail)
!
! CLose file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptscor
!
!
!
 Subroutine write_cptload (cload, afile, dfile, afield, ne, nv, v, iuse, rlat, rlng, cstn, ifail, &
                           cload2, v2)
!
! Outputs loadings in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_io_constants, Only: iout, lstn
  Use iofiles,           Only: ifile, ofile
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ne ! - number of loadings -
  Integer, Intent(In) :: nv ! - number of variables -
!
  Character(Len=*), Intent(In) :: cload ! - loadings description -
!
  Type(ifile), Intent(In) :: dfile ! - descriptor file -
!
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: cload2 ! - second loadings description -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v    ! - data -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! - optional input arrays -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: v2 ! - additional data -
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, dfile%nfs, dfile%idstr, nv, ifail)
  If (ifail == 0) Then
!
! Print data
     Select Case (dfile%idstr)
      Case (1)
        Call write_load (iout, afile, afield(:), dfile%idstr, dfile%nfs, dfile%nls, ne, nv, iuse(:), v(:,:), cload, ifail, &
             rlat=rlat(:,:), rlng=rlng(:,:))
        If (Present(cload2) .and. Present (v2)) Then
           Call write_load (iout, afile, afield(:), dfile%idstr, dfile%nfs, dfile%nls, ne, nv, iuse(:), v2(:,:), cload2, ifail, &
                rlat=rlat(:,:), rlng=rlng(:,:))
        End If
      Case (2)
        Call write_load (iout, afile, afield(:), dfile%idstr, dfile%nfs, dfile%nls, ne, nv, iuse(:), v(:,:), cload, ifail, &
             rlat=rlat(:,:), rlng=rlng(:,:), cstn=cstn(:,:))
        If (Present(cload2) .and. Present (v2)) Then
           Call write_load (iout, afile, afield(:), dfile%idstr, dfile%nfs, dfile%nls, ne, nv, iuse(:), v2(:,:), cload2, ifail, &
                rlat=rlat(:,:), rlng=rlng(:,:), cstn=cstn(:,:))
        End If
      Case (3)
        Call write_load (iout, afile, afield(:), dfile%idstr, dfile%nfs, dfile%nls, ne, nv, iuse(:), v(:,:), cload, ifail, &
             cstn=cstn(:,:))
        If (Present(cload2) .and. Present (v2)) Then
           Call write_load (iout, afile, afield(:), dfile%idstr, dfile%nfs, dfile%nls, ne, nv, iuse(:), v2(:,:), cload2, ifail, &
                cstn=cstn(:,:))
        End If
     End Select
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptload
!
!
!
 Subroutine write_cpteigs (afile, ne, nt, v, ifail)
!
! Outputs eigenvalues in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_io_constants, Only: iout
  Use iofiles,           Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ne ! - total number of eigenvalues -
  Integer, Intent(In) :: nt ! - number of cases -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: v ! - data -
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, 1, 3, 3, ifail)
!
! Print data
  If (ifail == 0) Then
     Call write_eigs (iout, afile, ne, v(:), nt, ifail)
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cpteigs
!
!
!
 Subroutine write_cptprobs (afile, dfile, afield, sdate, it0, nt, nv, v, iuse, kuse, rlat, rlng, cstn, ifail)
!
! Outputs probabilities in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_numbers,       Only: one
  Use data_cpt_constants, Only: ng, nts
  Use data_io_constants,  Only: iout, lstn
  Use time,               Only: pdate, &
                                nss
  Use iofiles,            Only: ifile, ofile
  Use fields,             Only: field, &
                                tfield
  Use categories,         Only: pthr
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: it0 ! - time 0 -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: nv  ! - number of variables -
!
  Type(ifile), Intent(In) :: dfile ! - descriptor file -
!
  Type(pdate), Intent(In) :: sdate ! - start date -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: v ! - data -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! Locals
!
! Local arrays
  Real(Kind=rp), Dimension(ng) :: probs ! - climatological probabilities -
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, dfile%nfs, dfile%idstr, nv, ifail, &
       ng=ng)
  If (ifail == 0) Then
!
! Calculate category probabilities
     probs(1) = pthr(1)
     probs(2) = pthr(2) - pthr(1)
     probs(ng) = one - pthr(nts)
!
! Set field information
     Allocate (tfield(dfile%nfl))
     tfield(:) = afield(:)
     tfield(:)%rmiss = -one
     tfield(:)%tprd%sdate = sdate ! - set start date for possible persistence component -
!
! Print data
     Select Case (dfile%idstr)
      Case (1)
        Call write_grids ('P', iout, afile, tfield(:), dfile%nfs, dfile%nls/nss, it0, nt, nv, ng, &
             v(:,:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), ifail, &
             prob=probs)
       Case (2)
        Call write_stnss ('P', iout, afile, tfield(:), dfile%nfs, dfile%nls/nss, it0, nt, nv, ng, &
             v(:,:,:), iuse(:), kuse(:), rlat(:,:), rlng(:,:), cstn(:,:), ifail, &
             prob=probs)
       Case (3)
        Call write_unrfs ('P', iout, afile, tfield(:), dfile%nfs, dfile%nls/nss, it0, nt, nv, ng, &
             v(:,:,:), iuse(:), kuse(:), cstn(:,:), ifail, &
             prob=probs)
     End Select
     Deallocate (tfield)
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptprobs
!
!
!
 Subroutine write_cptroc (afile, nb, ng, hit, far, roca, ifail)
!
! Prints ROC results
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: iout
  Use iofiles,           Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nb ! - number of bins -
  Integer, Intent(In) :: ng ! - number of categories -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: roca ! - ROC areas -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: hit ! - hit rates -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: far ! - false-alarm rates -
!
! Locals
!
! Local scalars
  Integer :: j ! - category index -
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, 1, 3, nb*ng, ifail)
!
! Print data
  If (ifail == 0) Then
     Do j = 1, ng
        Call write_cols (iout, afile, nb, 'Hit rates', hit(:,j), ifail, &
             cpt_field='ROC', cpt_nrow=nb, cpt_ncol=2, cpt_row='threshold', cpt_col='index', &
             c2='False-alarm rates', v2=far(:,j),                                            &
             cpt_tagr='area', cpt_rtag=roca(j))
        If (ifail /= 0) Exit
     End Do
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptroc
!
!
!
 Subroutine write_cptexcs (afile, iexc, nep, c, e, f, tmin, tmax, ifail)
!
! Prints exceedence probabilities
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: iout
  Use labels,            Only: cg_thresh
  Use iofiles,           Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iexc ! - exceedance (0) or non-exceedance (1) probabilities -
  Integer, Intent(In) :: nep  ! - number of exceedance probabilities -
!
  Real(Kind=rp), Intent(In) :: tmin ! - minimum threshold -
  Real(Kind=rp), Intent(In) :: tMax ! - Maximum threshold -
!
! INput/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(0:), Intent(In) :: c ! - climatological probabilities -
  Real(Kind=rp), Dimension(0:), Intent(In) :: e ! - empirical probabilities -
  Real(Kind=rp), Dimension(0:), Intent(In) :: f ! - forecast probabilities -
!
! Locals
!
! Local scalars
  Integer :: i ! - threshold index -
!
  Real(Kind=rp) :: tinc ! - threshold increment -
!
  Character(Len=32) :: cfield ! - CPT field -
!
! Local arrays
  Real(Kind=rp), Dimension(0:nep) :: t ! - thresholds -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, 1, 3, nep+1, ifail)
!
! Set appropriate X-increment
  If (ifail == 0) Then
     tinc = (tmax - tmin)/Real(nep, Kind=rp)
     t(0) = tmin
     Do i = 1, nep
        t(i) = t(i-1) + tinc
     End Do
!
! Select title
  Select Case (iexc)
   Case (0)
     cfield = 'exceedance-probs'
   Case (1)
     cfield = 'nonexceedance-probs'
  End Select
!
! Print data
     Call write_cols (iout, afile, nep+1, Trim(cg_thresh), t(0:), ifail, &
          cpt_field=Trim(cfield), cpt_nrow=nep, cpt_ncol=4, cpt_row='threshold', cpt_col='prob', &
          c2='Climatological', v2=c(0:), c3='Empirical', v3=e(0:), c4='Forecast', v4=f(0:))
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptexcs
!
!
!
 Subroutine write_cptrel (afile, ng, nb, afp, orf, ifq, tends, b0s, b1s, ifail)
!
! Prints reliability diagram results
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,      Only: half, one, fifty, oneh
  Use data_io_constants, Only: iout
  Use labels,            Only: c_lab, ls_cats, &
                               make_label
  Use iofiles,           Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ng ! - number of categories -
  Integer, Intent(In) :: nb ! - number of probability bins -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:,0:), Intent(In) :: ifq ! - frequency -
!
  Real(Kind=rp), Dimension(0:), Intent(In) :: b0s ! - intercepts -
  Real(Kind=rp), Dimension(0:), Intent(In) :: b1s ! - slopes -
!
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: afp   ! - bin-averaged forecast probability -
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: orf   ! - observed relative frequency -
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: tends ! - tendencies -
!
! Locals
!
! Local scalars
  Integer :: j   ! - category index -
  Integer :: k   ! - probability bin index -
  Integer :: kmn ! - lowest populated probability bin -
  Integer :: kmx ! - largest populated probability bin -
!
  Real(Kind=rp) :: b0  ! - intercept for no-skill line -
  Real(Kind=rp) :: tfq ! - total frequency -
!
! Local arrays
  Real(Kind=rp), Dimension(nb) :: ffq ! - fitted frequencies -
  Real(Kind=rp), Dimension(nb) :: rfq ! - relative frequencies -
  Real(Kind=rp), Dimension(nb) :: sk0 ! - zero-skill -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real,Sum, Trim
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, 1, 3, nb*(ng+1), ifail)
!
! Calculate outputs
  If (ifail == 0) Then
     Do j = 0, ng
! - total frequency -
        tfq = Real(Sum(ifq(:,j)), Kind=rp)
        b0 = half*tends(1,j)
! - minimum frequency -
        Do k = 1, nb
           If (ifq(k,j) > 0) Then
              kmn = k
              Exit
           End If
        End Do
! - Maximum frequency -
        Do k = nb, kmn, -1
           If (ifq(k,j) > 0) Then
              kmx = k
              Exit
           End If
        End Do
! - relative frequencies -
        Do k = 1, nb
           rfq(k) = Real(100*ifq(k,j), Kind=rp)/tfq
! - fitted frequencies -
           If (k >= kmn .and. k <= kmx) Then
              ffq(k) = b0s(j) + b1s(j)*afp(k,j)
           Else
              ffq(k) = -one
           End If
! - zero-skill -
           If (ifq(k,j) > 0) Then
              sk0(k) = b0 + fifty*afp(k,j)/oneh
           Else
              sk0(k) = b0 + fifty*Real(k-1)/Real(nb-1, Kind=rp)
           End If
        End Do
!
! Print data
! - print data -
        c_lab = make_label('L', ls_cats(j))
        Call write_cols (iout, afile, nb, 'Forecast probability (x)', afp(:,j), ifail, &
             cpt_field='attributes', cpt_nrow=nb, cpt_ncol=4, cpt_row='probability', cpt_col='attribute', &
             c2='Observed relative frequency (y)', v2=orf(:,j),                                           &
             c3='Frequency', i3=ifq(:,j),                                                                 &
             c4='Relative frequency (%)', v4=rfq(:),                                                      &
             c5='Fitted frequency (%)', v5=ffq(:),                                                        &
             c6='No skill (%)', v6=sk0(:),                                                                &
             cpt_c=j, cpt_tagc='category', cpt_ctag=Trim(c_lab))
! - print scores -
        Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) 'Observed frequency'
        If (output_value(iout, r=tends(1,j), ldelim=.true.) /= 0) GoTo 1
        Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) 'Average forecast  '
        If (output_value(iout, r=tends(2,j), ldelim=.true.) /= 0) GoTo 1
        Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) 'Intercept (%)     '
        If (output_value(iout, r=b0s(j), ldelim=.true.) /= 0) GoTo 1
        Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) 'Slope (%)         '
        If (output_value(iout, r=b1s(j), ldelim=.true.) /= 0) GoTo 1
     End Do
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_cptrel
!
!
!
 Subroutine write_cptregr (afile, b, b0, ifail)
!
! Outputs regression coefficients in CPT format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!   -1 If output should be halted for all subsequent files
!
! Modules
  Use data_io_constants, Only: iout
  Use arrays,            Only: iusex, iusey
  Use iofiles,           Only: ofile, &
                               nx, ny, xfile, yfile
  Use space,             Only: cstndx, cstndy, rlatdx, rlatdy, rlngdx, rlngdy
  Use fields,            Only: xfield, yfield
!
! Arguments
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: b0 ! - intercepts -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: b  ! - coefficients -
!
! Executable Statements
!
! Open output file
  ifail = 0
  Call open_output (iout, afile, xfile%nfs, yfile%idstr, ny, ifail)
!
! Print data
  If (ifail == 0) Then
     Select Case (xfile%idstr)
      Case (1)
        Select Case (yfile%idstr)
         Case (1)
           Call write_regr (iout, afile, xfile, xfield(:), nx, iusex(:), yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                rlatx=rlatdx(:,:), rlngx=rlngdx(:,:), rlaty=rlatdy(:,:), rlngy=rlngdy(:,:))
         Case (2)
           Call write_regr (iout, afile, xfile, xfield(:), nx, iusex(:), yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                rlatx=rlatdx(:,:), rlngx=rlngdx(:,:), rlaty=rlatdy(:,:), rlngy=rlngdy(:,:), cstny=cstndy(:,:))
         Case (3)
           Call write_regr (iout, afile, xfile, xfield(:), nx, iusex(:), yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                rlatx=rlatdx(:,:), rlngx=rlngdx(:,:), cstny=cstndy(:,:))
        End Select
      Case (2)
        Select Case (yfile%idstr)
         Case (1)
           Call write_regr (iout, afile, xfile, xfield, nx, iusex, yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                rlatx=rlatdx(:,:), rlngx=rlngdx(:,:), cstnx=cstndx(:,:), rlaty=rlatdy(:,:), rlngy=rlngdy(:,:))
         Case (2)
           Call write_regr (iout, afile, xfile, xfield, nx, iusex, yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                rlatx=rlatdx(:,:), rlngx=rlngdx(:,:), cstnx=cstndx(:,:), rlaty=rlatdy(:,:), rlngy=rlngdy(:,:), cstny=cstndy(:,:))
         Case (3)
           Call write_regr (iout,afile, xfile, xfield, nx, iusex, yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                rlatx=rlatdx(:,:), rlngx=rlngdx(:,:), cstnx=cstndx(:,:), cstny=cstndy(:,:))
        End Select
      Case (3)
        Select Case (yfile%idstr)
         Case (1)
           Call write_regr (iout, afile, xfile, xfield, nx, iusex(:), yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                cstnx=cstndx(:,:), rlaty=rlatdy(:,:), rlngy=rlngdy(:,:), b0=b0(:,1))
         Case (2)
           Call write_regr (iout, afile, xfile, xfield(:), nx, iusex(:), yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                cstnx=cstndx(:,:), rlaty=rlatdy(:,:), rlngy=rlngdy(:,:), cstny=cstndy(:,:), b0=b0(:,1))
         Case (3)
           Call write_regr (iout, afile, xfile, xfield(:), nx, iusex(:), yfile, yfield(1:), ny, iusey(:), b(:,:), ifail, &
                cstnx=cstndx(:,:), cstny=cstndy(:,:), b0=b0(:,1))
        End Select
     End Select
  End If
!
! Close file and reset
  Call close_output (iout, afile, ifail)
!
  Return
 End Subroutine write_cptregr
!
!
!
 Subroutine write_grid (iout, afile, afield, nfs, nls, nt, nv, v, iuse, kuse, rlat, rlng, ifail, &
            ctagi, itag, it0)
!
! Outputs data in gridded format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: ffmts, lprd
  Use time,              Only: pdate, &
                               Assignment(=), Operator(==), Operator(+), &
                               get_cdate
  Use iofiles,           Only: ofile, &
                               itags
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of grids -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: it0  ! - time 0 -
  Integer, Intent(In), Optional :: itag ! - extra integer tag value -
!
  Character(Len=*), Intent(In), Optional :: ctagi ! - extra integer tag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v    ! - data -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! Locals
!
! Local scalars
  Integer :: k   ! - time index -
  Integer :: l   ! - field / lagged field index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged field index -
  Integer :: ii  ! - available series index -
  Integer :: ij  ! - current variable -
  Integer :: kk  ! - available time index -
!
  Character(Len=lprd) :: csp ! - output label for latitude-longitude space -
!
  Type(pdate) :: mdate ! - start date -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustL, Present, Sum, Trim
!
! Executable Statements
!
! Print gridded data
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Do ifd = 1, nfs
        Do ilf = 1, nls
           l = (ifd - 1)*nls + ilf
           If (l > 1) Then
              ij = Sum(afield(1:l-1)%nv)
              ii = Sum(afield(1:l-1)%nva) + 1
           Else
              ij = 0
              ii = 1
           End If
           If (afield(l)%mdate == 0) Then
              mdate = 0
           Else
              If (.not.Present(it0)) mdate=0
           End If
           kk = 0
           Do k = 1, nt
              If (kuse(k)) kk = kk + 1
              If ((.not.kuse(k) .and. k > 2) .or. kk == 0) Cycle
              If (itags == 1) Then
                 If (Present(it0)) Then
                    If (.not.(afield(l)%mdate == 0)) mdate = afield(l)%mdate + (it0 + k - 1)
                    If (k == 1) Then
                       If (ilf == 1) Then
                          Call write_tag (iout, ifail, &
                                          cpt_field=afield(l)%var, cpt_z=afield(l)%z, cpt_m=afield(l)%member,     &
                                          cpt_model=afield(l)%Model, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1), &
                                          cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs,         &
                                          cpt_row='Y', cpt_col='X', cpt_units=afield(l)%unitc,                    &
                                          cpt_missing=afield(l)%rmiss,                                            &
                                          cpt_tagi=ctagi, cpt_itag=itag)
                       Else
                          Call write_tag (iout, ifail, &
                                          cpt_model=afield(l)%Model, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1), &
                                          cpt_tagi=ctagi, cpt_itag=itag)
                       End If
                    Else
                       Call write_tag (iout, ifail, &
                                       cpt_model=afield(l)%Model, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1), &
                                       cpt_tagi=ctagi, cpt_itag=itag)
                    End If
                 Else
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_model=afield(l)%Model,             &
                                    cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                                    cpt_row='Y', cpt_col='X', cpt_units=afield(l)%unitc,            &
                                    cpt_missing=afield(l)%rmiss,                                    &
                                    cpt_tagi=ctagi, cpt_itag=itag)
                 End If
                 csp = ' '
                 If (ifail /= 0) GoTo 1
              Else
                 If (Present(it0)) Then
                    csp = get_cdate(afield(l)%tprd+(it0+k-1), 1, 1)
                 Else
                    Write (Unit=csp, Fmt=*) afield(l)%rmiss
                    csp = AdjustL(csp)
                 End If
              End If
              If (kuse(kk)) Then
                 Call write_gridded_block_v10 (iout, csp, afield(l), rlat(:,l), rlng(:,l), ij, ifail, &
                      v=v(ii:,kk), iuse=iuse(ii:))
              Else
                 Call write_gridded_block_v10 (iout, csp, afield(l), rlat(:,l), rlng(:,l), ij, ifail)
              End If
              If (ifail /= 0) GoTo 1
           End Do
        End Do
     End Do
   Case ('unformatted')
     Call write_unform (iout, afile, nt, nv, v, iuse, kuse, afield(1)%rmiss, ifail)
     If (ifail /= 0) Return
   Case ('GrADS')
     Call write_unform (iout, afile, nt, nv, v, iuse, kuse, afield(1)%rmiss, ifail)
     If (ifail /= 0) Return
     Call write_grads_grid_ctrl (iout+1, afile, afield(:), nfs, nls, nt, 1, rlat(:,:), rlng(:,:), ifail)
     If (ifail /= 0) Return
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_grid
!
!
!
 Subroutine write_stns (iout, afile, afield, nfs, nls, it0, nt, nv, v, iuse, kuse, rlat, rlng, cstn, ifail, &
            ctagi, itag)
!
! Outputs data in station format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: ffmts, lstn
  Use time,              Only: Operator(+), &
                               get_cdate
  Use iofiles,           Only: ofile, &
                               itags
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: it0  ! - time 0 -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of stations -
!
  Type(ofile), Intent(In) :: afile  ! - output file -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: itag ! - extra integer tag value -
!
  Character(Len=*), Intent(In), Optional :: ctagi ! - extra integer tag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v    ! - data -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! Locals
!
! Local scalars
  Integer :: i   ! - station index -
  Integer :: k   ! - time index -
  Integer :: l   ! - field / lagged field index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged field index -
  Integer :: ii  ! - available series index -
  Integer :: ij  ! - current variable -
  Integer :: kk  ! - available time index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Sum, Trim
!
! Executable Statements
!
! Print tag line
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Do ifd = 1, nfs
        l = (ifd - 1)*nls + 1
        If (itags == 1) Then
           Call write_tag (iout, ifail, &
                           cpt_field=afield(l)%var, cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, &
                           cpt_row='T', cpt_col='station', cpt_units=Trim(afield(l)%unitc), &
                           cpt_missing=afield(l)%rmiss,                                     &
                           cpt_tagi=ctagi, cpt_itag=itag)
           If (ifail /= 0) GoTo 1
        End If
!
! Print station information
        Do i = 1, afield(l)%nv
           If (output_value(iout, c=cstn(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
        End Do
        If (itags == 1) Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) 'cpt:Y'
        Do i = 1, afield(l)%nv
           If (output_value(iout, r=rlat(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
        End Do
        If (itags == 1) Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) 'cpt:X'
        Do i = 1, afield(l)%nv
           If (output_value(iout, r=rlng(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
        End Do
!
! Print station data
        kk = 0
        Do k = 1, nt
           If (kuse(k)) kk = kk + 1
           If (ifd == 1) Then
              ii = 1
              ij = 0
           Else
              ii = Sum(afield(1:l-1)%nva) + 1
              ij = Sum(afield(1:l-1)%nv)
           End If
           Do ilf = 1, nls
              cout = get_cdate(afield(l)%tprd+(it0+k-1), 1, 1)
              If (output_value(iout, c=cout, ladv=.false.) /= 0) GoTo 1
              Do i = 1, afield(l)%nv
                 If (kuse(k)) Then
                    ij = ij + 1
                    If (iuse(ii) == ij) Then
                       If (output_value(iout, r=v(ii,kk), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                       ii = ii + 1
                    Else
                       If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                    End If
                 Else
                    If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                 End If
                 If (ii > nv) Exit
              End Do
              Write (Unit=iout, Fmt='()', Advance='yes', Err=1)
           End Do
        End Do
     End Do
   Case ('unformatted')
     Call write_unform (iout, afile, nt, nv, v, iuse, kuse, afield(1)%rmiss, ifail)
     If (ifail /= 0) Return
   Case ('GrADS')
     Call write_grads_stns (iout, afile, nt, nls, afield(1)%nva, v, afield(:)%tprd%sdate, &
          cstn(:,1), rlat(:,1), rlng(:,1), iuse, kuse, afield(1)%rmiss,ifail)
     If (ifail /= 0) Return
     Call write_grads_stns_ctrl (iout+1, afile, afield(:), nfs,nls, nt, 1, ifail)
     If (ifail /= 0) Return
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_stns
!
!
!
 Subroutine write_unrf (iout, afile, afield, nfs, nls, it0, nt, nv, v, iuse, kuse, cstn, cpt_row, ifail, &
            b0, ctagi, itag)
!
! Outputs data in unreferenced format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: ffmts
  Use time,              Only: Operator(+), &
                               get_cdate
  Use iofiles,           Only: ofile, &
                               itags
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: it0  ! - time 0 -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of series -
!
  Character(Len=*), Intent(In) :: cpt_row ! - CPT row -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: itag ! - extra integer tag value -
!
  Real(Kind=rp), Intent(In), Optional :: b0 ! - regression constant -
!
  Character(Len=*), Intent(In), Optional :: ctagi ! - extra integer tag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v ! - data -
!
  Character(Len=*), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! Locals
!
! Local scalars
  Integer :: i   ! - series index -
  Integer :: k   ! - time index -
  Integer :: l   ! - field / lagged field index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged field index -
  Integer :: ii  ! - available series index -
  Integer :: ij  ! - current variable -
  Integer :: kk  ! - available time index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Sum, Trim
!
! Executable Statements
!
! Print tag line
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Do ifd = 1, nfs
        l = (ifd - 1)*nls + 1
        If (itags == 1) Then
           Call write_tag (iout, ifail, &
                           cpt_field=afield(l)%var, cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, &
                           cpt_row=cpt_row, cpt_col='index', cpt_units=afield(l)%unitc,     &
                           cpt_missing=afield(l)%rmiss,                                     &
                           cpt_tagi=ctagi, cpt_itag=itag)
           If (ifail /= 0) GoTo 1
        End If
!
! Print station information
        If (Present(b0)) Then
           If (output_value(iout, c='Constant', ldelim=.true., ladv=.false.) /= 0) GoTo 1
        End If
        Do i = 1, afield(l)%nv
           If (output_value(iout, c=cstn(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
        End Do
!
! Print unreferenced data
        kk = 0
        Do k = 1, nt
           If (kuse(k)) kk = kk + 1
           If (ifd == 1) Then
              ii = 1
              ij = 0
           Else
              ii = Sum(afield(1:l-1)%nva) + 1
              ij = Sum(afield(1:l-1)%nv)
           End If
           Do ilf = 1, nls
              cout = get_cdate(afield(l)%tprd+(it0+k-1), 1, 1)
              If (output_value(iout, c=cout, ladv=.false.) /= 0) GoTo 1
              Do i = 1, afield(l)%nv
                 If (kuse(k)) Then
                    ij = ij + 1
                    If (iuse(ii) == ij) Then
                       If (output_value(iout, r=v(ii,kk), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                       ii = ii + 1
                    Else
                       If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                    End If
                 Else
                    If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                 End If
                 If (ii > nv) Exit
              End Do
              Write (Unit=iout, Fmt='()', Advance='yes', Err=1)
           End Do
        End Do
     End Do
   Case ('unformatted')
     Call write_unform (iout, afile, nt, nv, v, iuse, kuse, afield(1)%rmiss, ifail)
     If (ifail /= 0) Return
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_unrf
!
!
!
 Subroutine write_grids (ctype, iout, afile, afield, nfs, nls, it0, nt, nv, n3, v, iuse, kuse, rlat, rlng, ifail, &
            prob, cl)
!
! Outputs data in gridded format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: ffmts, lprd
  Use labels,            Only: c_thresh=>c_lab, &
                               ls_cats, ls_threshs, &
                               make_label
  Use time,              Only: pdate, &
                               Assignment(=), Operator(==), Operator(+), &
                               get_cdate
  Use iofiles,           Only: ofile, &
                               itags
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: it0  ! - time 0 -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of grids -
  Integer, Intent(In) :: n3   ! - third Dimension -
!
  Character(Len=1), Intent(In) :: ctype ! - output Type -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! - Optional input scalars -
  Real(Kind=rp), Intent(In), Optional :: cl ! - confidence level -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: v ! - data -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! - Optional input arrays -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: prob ! - climatological probabilities -
!
! Locals
!
! Local scalars
  Integer :: k   ! - time index -
  Integer :: l   ! - field / lagged field index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged field index -
  Integer :: i3  ! - third Dimension index -
  Integer :: ii  ! - available series index -
  Integer :: ij  ! - current variable -
  Integer :: kk  ! - available time index -
!
  Character(Len=lprd) :: csp ! - output label for latitude-longitude space -
!
  Type(pdate) :: mdate ! - start date -
!
! Local arrays
  Character(Len=16), Dimension(n3) :: c_thrs ! - thresholds -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Sum, Trim
!
! Executable Statements
!
! Print gridded data
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Do ifd = 1, nfs
        Do ilf = 1, nls
           l = (ifd - 1)*nls + ilf
           If (l > 1) Then
              ij = Sum(afield(1:l-1)%nv)
              ii = Sum(afield(1:l-1)%nva) + 1
           Else
              ij = 0
              ii = 1
           End If
           If (afield(l)%mdate == 0) mdate = 0
           kk = 0
           Do k = 1, nt
              If (kuse(k)) kk = kk + 1
              If ((.not.kuse(k) .and. k > 2) .or. kk == 0) Cycle
              If (.not.afield(l)%mdate == 0) mdate = afield(l)%mdate + (it0 + k - 1)
              Do i3 = 1, n3
                 If (itags == 1) Then
                    If (i3 == 1) Then
                       If (k == 1) Then
                          If (ilf == 1) Then
                             Select Case (ctype)
                              Case ('E') ! - ensemble -
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(l)%var, cpt_z=afield(l)%z, cpt_m=i3,           &
                                                cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1),                    &
                                                cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                                                cpt_row='Y', cpt_col='X', cpt_units=afield(l)%unitc,            &
                                                cpt_missing=afield(l)%rmiss)
                              Case ('L') ! - prediction limits -
                                c_thresh = make_label('L', ls_threshs(i3))
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(l)%var, cpt_z=afield(l)%z,                     &
                                                cpt_m=afield(l)%member, cpt_clev=cl, cpt_limit=c_thresh,        &
                                                cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1),                    &
                                                cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                                                cpt_row='Y', cpt_col='X', cpt_units='%',                        &
                                                cpt_missing=afield(l)%rmiss)
                              Case ('M') ! - Models -
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(i3)%var, cpt_model=afield(i3)%Model,             &
                                                cpt_z=afield(i3)%z, cpt_m=afield(i3)%member,                      &
                                                cpt_s=mdate, cpt_t=afield(i3)%tprd+(it0+k-1),                     &
                                                cpt_nrow=afield(i3)%region%nlts, cpt_ncol=afield(i3)%region%nlgs, &
                                                cpt_row='Y', cpt_col='X', cpt_units=afield(i3)%unitc,             &
                                                cpt_missing=afield(i3)%rmiss)
                              Case ('O') ! - odds -
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(l)%var, cpt_c=i3, cpt_z=afield(l)%z,                &
                                                cpt_m=afield(l)%member, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1), &
                                                cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs,      &
                                                cpt_row='Y', cpt_col='X', cpt_units='odds',                          &
                                                cpt_missing=afield(l)%rmiss)
                              Case ('P') ! - probabilities -
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(l)%var, cpt_c=i3, cpt_prob=prob(i3),           &
                                                cpt_z=afield(l)%z, cpt_m=afield(l)%member, cpt_s=mdate,         &
                                                cpt_t=afield(l)%tprd+(it0+k-1),                                 &
                                                cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                                                cpt_row='Y', cpt_col='X', cpt_units='probability (%)',          &
                                                cpt_missing=afield(l)%rmiss)
                              Case ('R') ! - relative odds -
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(l)%var, cpt_c=i3, cpt_z=afield(l)%z,                &
                                                cpt_m=afield(l)%member, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1), &
                                                cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs,      &
                                                cpt_row='Y', cpt_col='X', cpt_units='relative odds',                 &
                                                cpt_missing=afield(l)%rmiss)
                              Case ('T') ! - thresholds -
                                c_thresh = make_label('L', ls_threshs(i3))
                                Call write_tag (iout, ifail, &
                                                cpt_field=afield(l)%var, cpt_z=afield(l)%z,                     &
                                                cpt_m=afield(l)%member, cpt_limit=c_thresh,                     &
                                                cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1),                    &
                                                cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                                                cpt_row='Y', cpt_col='X', cpt_units=afield(l)%unitc,            &
                                                cpt_missing=afield(l)%rmiss)
                             End Select
                          Else
                             Select Case (ctype)
                              Case ('E') ! - ensemble -
                                Call write_tag (iout, ifail, &
                                                cpt_m=i3, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                              Case ('L', 'T') ! - prediction limits, and thresholds -
                                c_thresh = make_label('L', ls_threshs(i3))
                                Call write_tag (iout, ifail, &
                                                cpt_limit=c_thresh, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                              Case ('M') ! - Model -
                                Call write_tag (iout, ifail, &
                                                cpt_model=afield(i3)%Model, &
                                                cpt_s=mdate, cpt_t=afield(i3)%tprd+(it0+k-1))
                              Case ('O', 'R') ! - odds, and relative odds -
                                Call write_tag (iout, ifail, &
                                                cpt_c=i3, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                              Case ('P') ! - probabilities -
                                Call write_tag (iout, ifail, &
                                                cpt_c=i3, cpt_prob=prob(i3), cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                             End Select
                          End If
                       Else
                          Select Case (ctype)
                           Case ('E') ! - ensemble -
                             Call write_tag (iout, ifail, &
                                             cpt_m=i3, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                           Case ('L', 'T') ! - prediction limits, and thresholds -
                             c_thresh = make_label('L', ls_threshs(i3))
                             Call write_tag (iout, ifail, &
                                             cpt_limit=c_thresh, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                           Case ('M') ! - Model -
                             Call write_tag (iout, ifail, &
                                             cpt_field=Trim(afield(i3)%var), cpt_model=Trim(afield(i3)%Model), &
                                             cpt_s=mdate, cpt_t=afield(i3)%tprd+(it0+k-1))
                           Case ('O', 'R') ! - odds, and relative odds -
                             Call write_tag (iout, ifail, &
                                             cpt_c=i3, cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                           Case ('P') ! - probabilities -
                             Call write_tag (iout, ifail, &
                                             cpt_c=i3, cpt_prob=prob(i3), cpt_s=mdate, cpt_t=afield(l)%tprd+(it0+k-1))
                          End Select
                       End If
                    Else
                       Select Case (ctype)
                        Case ('E') ! - ensemble -
                          Call write_tag (iout, ifail, &
                                          cpt_m=i3)
                        Case ('L', 'T') ! - prediction limits, and thresholds -
                          c_thresh = make_label('L', ls_threshs(i3))
                          Call write_tag (iout, ifail, &
                                          cpt_limit=c_thresh)
                        Case ('M') ! - Model -
                          Call write_tag (iout, ifail, &
                                          cpt_field=Trim(afield(i3)%var), cpt_model=Trim(afield(i3)%Model), &
                                          cpt_s=mdate, cpt_t=afield(i3)%tprd+(it0+k-1))
                        Case ('O', 'R') ! - odds, and relative odds -
                          Call write_tag (iout, ifail, &
                                          cpt_c=i3)
                        Case ('P') ! - probabilities -
                          Call write_tag (iout, ifail, &
                                          cpt_c=i3, cpt_prob=prob(i3))
                       End Select
                    End If
                    csp = ' '
                    If (ifail /= 0) GoTo 1
                 Else
                    csp = get_cdate(afield(l)%tprd+(it0+k-1), 1, 1)
                 End If
                 If (kuse(kk)) Then
                    Call write_gridded_block_v10 (iout, csp, afield(l), rlat(:,l), rlng(:,l), ij, ifail, &
                         v=v(ii:,kk,i3), iuse=iuse(ii:))
                 Else
                    Call write_gridded_block_v10 (iout, csp, afield(l), rlat(:,l), rlng(:,l), ij, ifail)
                 End If
                 If (ifail /= 0) GoTo 1
              End Do
           End Do
        End Do
     End Do
   Case ('unformatted')
     Do i3 = 1, n3
        Call write_unform (iout, afile, nt, nv, v(:,:,i3), iuse, kuse, afield(1)%rmiss, ifail)
        If (ifail /= 0) Return
     End Do
   Case ('GrADS')
     Do i3 = 1, n3
        Call write_unform (iout, afile, nt, nv, v(:,:,i3), iuse, kuse, afield(1)%rmiss, ifail)
        If (ifail /= 0) Return
     End Do
     Select Case (ctype)
      Case ('E') ! - ensemble -
      Case ('L', 'T') ! - prediction limits, and thresholds -
        Do i3 = 1, n3
           c_thrs(i3) = Trim(make_label('L', ls_threshs(i3)))
        End Do
        Call write_grads_grid_ctrl (iout+1, afile, afield(:), nfs, nls, nt, n3, rlat(:,:), rlng(:,:), ifail, &
             cvar=c_thrs(:))
      Case ('M') ! - Models -
      Case ('O', 'P', 'R') ! - odds, probabilities, relative odds -
        Call write_grads_grid_ctrl (iout+1, afile, afield(:), nfs, nls, nt, n3, rlat(:,:), rlng(:,:), ifail, &
             cvar=ls_cats(1:)%c)
     End Select
     If (ifail /= 0) Return
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_grids
!
!
!
 Subroutine write_stnss (ctype, iout, afile, afield, nfs, nls, it0, nt, nv, n3, v, iuse, kuse, rlat, rlng, cstn, ifail, &
            prob, cl)
!
! Outputs data in station format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: ffmts, lstn
  Use labels,            Only: c_thresh=>c_lab, &
                               ls_cats, ls_threshs, &
                               make_label
  Use time,              Only: Operator(+), &
                               get_cdate
  Use iofiles,           Only: ofile, &
                               itags
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: it0  ! - time 0 -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of stations -
  Integer, Intent(In) :: n3   ! - third Dimension -
!
  Character(Len=1), Intent(In) :: ctype ! - output Type -
!
  Type(ofile), Intent(In) :: afile  ! - output file -
!
! - Optional input scalars -
  Real(Kind=rp), Intent(In), Optional :: cl ! - confidence level -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rlng ! - longitudes -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: v ! - data -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! - Optional input arrays -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: prob ! - climatological probabilities -
!
! Locals
!
! Local scalars
  Integer :: i   ! - station index -
  Integer :: k   ! - time index -
  Integer :: l   ! - field / lagged field index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged field index -
  Integer :: i3  ! - third Dimension index -
  Integer :: ii  ! - available series index -
  Integer :: ij  ! - current variable -
  Integer :: kk  ! - available time index -
!
! Local arrays
  Character(Len=16), Dimension(n3) :: c_thrs ! - thresholds -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Sum, Trim
!
! Executable Statements
!
! Print tag line
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Do ifd = 1, nfs
        l = (ifd - 1)*nls + 1
        Do i3 = 1, n3
           If (itags == 1) Then
              If (i3 == 1) Then
                 Select Case (ctype)
                  Case ('E') ! - ensemble -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_m=i3,                                      &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units=afield(l)%unitc, cpt_missing=afield(l)%rmiss)
                  Case ('L') ! - prediction limits -
                    c_thresh = make_label('L', ls_threshs(i3))
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var,                                                &
                                    cpt_clev=cl, cpt_limit=c_thresh,                                        &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units=afield(l)%unitc, cpt_missing=afield(l)%rmiss)
                  Case ('M') ! - Models -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(i3)%var,                                                &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(i3)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units=afield(i3)%unitc, cpt_missing=afield(i3)%rmiss)
                  Case ('O') ! - odds -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_c=i3,                                      &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units='odds', cpt_missing=afield(l)%rmiss)
                  Case ('P') ! - probabilities -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_c=i3, cpt_prob=prob(i3),                   &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units='probability (%)', cpt_missing=afield(l)%rmiss)
                  Case ('R') ! - relative odds -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_c=i3,                                      &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units='relative odds', cpt_missing=afield(l)%rmiss)
                  Case ('T') ! - thresholds -
                    c_thresh = make_label('L', ls_threshs(i3))
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var,                                                &
                                    cpt_limit=c_thresh,                                                     &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units=afield(i3)%unitc, cpt_missing=afield(i3)%rmiss)
                 End Select
              Else
                 Select Case (ctype)
                  Case ('E') ! - ensemble -
                    Call write_tag (iout, ifail, &
                                    cpt_m=i3)
                  Case ('L', 'T') ! - prediction limits, and thresholds -
                    c_thresh = make_label('L', ls_threshs(i3))
                    Call write_tag (iout, ifail, &
                                    cpt_limit=c_thresh)
                  Case ('M') ! - Models -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(i3)%var,                                                &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(i3)%nv, cpt_row='T', cpt_col='station', &
                                    cpt_units=afield(i3)%unitc, cpt_missing=afield(i3)%rmiss)
                  Case ('O', 'R') ! - odds, and relative odds -
                    Call write_tag (iout, ifail, &
                                    cpt_c=i3)
                  Case ('P') ! - probabilities -
                    Call write_tag (iout, ifail, &
                                    cpt_c=i3, cpt_prob=prob(i3))
                 End Select
              End If
              If (ifail /= 0) GoTo 1
           End If
!
! Print station information
           Do i = 1, afield(l)%nv
              If (output_value(iout, c=cstn(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
           End Do
           If (itags == 1) Then
              If (output_value(iout, c='cpt:Y', ladv=.false.) /= 0) GoTo 1
           End If
           Do i = 1, afield(l)%nv
              If (output_value(iout, r=rlat(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
           End Do
           If (itags == 1) Then
              If (output_value(iout, c='cpt:X', ladv=.false.) /= 0) GoTo 1
           End If
           Do i = 1, afield(l)%nv
              If (output_value(iout, r=rlng(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
           End Do
!
! Print station data
           kk = 0
           Do k = 1, nt
              If (kuse(k)) kk = kk + 1
              If (ifd == 1) Then
                 ii = 1
                 ij = 0
              Else
                 ii = Sum(afield(1:l-1)%nva) + 1
                 ij = Sum(afield(1:l-1)%nv)
              End If
              Do ilf = 1, nls
                 cout = get_cdate(afield(l)%tprd+(it0+k-1), 1, 1)
                 If (output_value(iout, c=cout, ladv=.false.) /= 0) GoTo 1
                 Do i = 1, afield(l)%nv
                    If (kuse(k)) Then
                       ij = ij + 1
                       If (iuse(ii) == ij) Then
                          If (output_value(iout, r=v(ii,kk,i3), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                          ii = ii + 1
                       Else
                          If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                       End If
                    Else
                       If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                    End If
                 End Do
                 Write (Unit=iout, Fmt='()', Advance='yes', Err=1)
              End Do
           End Do
        End Do
     End Do
   Case ('unformatted')
     Do i3 = 1, n3
        Call write_unform (iout, afile, nt, nv, v(:,:,i3), iuse, kuse, afield(1)%rmiss, ifail)
        If (ifail /= 0) Return
     End Do
   Case ('GrADS')
     Call write_grads_stns3 (iout, afile, nt, nls, afield(1)%nva, n3, v(:,:,:), afield(:)%tprd%sdate, &
          cstn(:,1), rlat(:,1), rlng(:,1), iuse, kuse, afield(1)%rmiss, ifail)
     If (ifail /= 0) Return
     Select Case (ctype)
      Case ('E') ! - ensemble -
      Case ('L', 'T') ! - prediction limits, and thresholds -
        Do i3 = 1, n3
           c_thrs(i3) = Trim(make_label('L', ls_threshs(i3)))
        End Do
        Call write_grads_stns_ctrl (iout+1, afile, afield(:), nfs, nls, nt, n3, ifail, &
             cvar=c_thrs(:))
      Case ('M') ! - Models -
      Case ('O', 'P', 'R') ! - odds, probabilities, relative odds -
        Call write_grads_stns_ctrl (iout+1, afile, afield(:), nfs, nls, nt, n3, ifail, &
             cvar=ls_cats(1:)%c)
     End Select
     If (ifail /= 0) Return
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_stnss
!
!
!
 Subroutine write_unrfs (ctype, iout, afile, afield, nfs, nls, it0, nt, nv, n3, v, iuse, kuse, cstn, ifail, &
             prob, cl)
!
! Outputs data in unreferenced format
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants, Only: ffmts, lstn
  Use labels,            Only: c_thresh=>c_lab, &
                               ls_threshs, &
                               make_label
  Use time,              Only: Operator(+), &
                               get_cdate
  Use iofiles,           Only: ofile, &
                               itags
  Use fields,            Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: it0  ! - time 0 -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of series -
  Integer, Intent(In) :: n3   ! - third Dimension -
!
  Character(Len=1), Intent(In) :: ctype ! - output Type -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! - Optional input scalars -
  Real(Kind=rp), Intent(In), Optional :: cl ! - confidence level -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: v ! - data -
!
  Character(Len=lstn), Dimension(:,:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! - Optional input arrays -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: prob ! - climatological probabilities -
!
! Locals
!
! Local scalars
  Integer :: i   ! - series index -
  Integer :: k   ! - time index -
  Integer :: l   ! - lagged field index -
  Integer :: i3  ! - third Dimension index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged-field index -
  Integer :: ii  ! - available series index -
  Integer :: kk  ! - available time index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Print tag line
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Do ifd = 1, nfs
        l = (ifd - 1)*nfs + 1
        Do i3 = 1, n3
           If (itags == 1) Then
              If (i3 == 1) Then
                 Select Case (ctype)
                  Case ('E') ! - ensemble -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_m=i3,                                    &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units=afield(l)%unitc, cpt_missing=afield(l)%rmiss)
                  Case ('L') ! - prediction limits -
                    c_thresh = make_label('L', ls_threshs(i3))
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var,                                              &
                                    cpt_clev=cl, cpt_limit=c_thresh,                                      &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units=Trim(afield(l)%unitc), cpt_missing=afield(l)%rmiss)
                  Case ('M') ! - Models -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(i3)%var,                                              &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(i3)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units=afield(i3)%unitc, cpt_missing=afield(i3)%rmiss)
                  Case ('O') ! - odds -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_c=i3,                                    &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units='odds', cpt_missing=afield(l)%rmiss)
                  Case ('P') ! - probabilities -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_c=i3, cpt_prob=prob(i3),                 &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units='probability (%)', cpt_missing=afield(l)%rmiss)
                  Case ('R') ! - relative odds -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var, cpt_c=i3,                                    &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units='relative odds', cpt_missing=afield(l)%rmiss)
                  Case ('T') ! - thresholds -
                    c_thresh = make_label('L', ls_threshs(i3))
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(l)%var,                                              &
                                    cpt_limit=c_thresh,                                                   &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(l)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units=afield(i3)%unitc, cpt_missing=afield(i3)%rmiss)
                 End Select
              Else
                 Select Case (ctype)
                  Case ('E') ! - ensemble -
                    Call write_tag (iout, ifail, &
                                    cpt_m=i3)
                  Case ('L', 'T') ! - prediction limits, and thresholds -
                    c_thresh = make_label('L', ls_threshs(i3))
                    Call write_tag (iout, ifail, &
                                    cpt_limit=c_thresh)
                  Case ('M') ! - Models -
                    Call write_tag (iout, ifail, &
                                    cpt_field=afield(i3)%var,                                              &
                                    cpt_nrow=nt*nls, cpt_ncol=afield(i3)%nv, cpt_row='T', cpt_col='index', &
                                    cpt_units=afield(i3)%unitc, cpt_missing=afield(i3)%rmiss)
                  Case ('O', 'R') ! - odds, and relative odds -
                    Call write_tag (iout, ifail, &
                                    cpt_c=i3)
                  Case ('P') ! - probabilities -
                    Call write_tag (iout, ifail, &
                                    cpt_c=i3, cpt_prob=prob(i3))
                 End Select
              End If
              If (ifail /= 0) GoTo 1
           End If
!
! Print index names
           Do i = 1, afield(l)%nv
              If (output_value(iout, c=cstn(i,l), ldelim=.true., ladv=(i==afield(l)%nv)) /= 0) GoTo 1
           End Do
!
! Print unreferenced data
           kk = 0
           Do k = 1, nt
              If (kuse(k)) kk = kk + 1
              ii = 1
              Do ilf = 1, nls
                 l = (ifd - 1)*nfs + ilf
                 cout = get_cdate(afield(l)%tprd+(it0+k-1), 1, 1)
                 If (output_value(iout, c=cout, ladv=.false.) /= 0) GoTo 1
                 Do i = 1, afield(l)%nv
                    If (kuse(k)) Then
                       If (iuse(ii) == i) Then
                          If (output_value(iout, r=v(ii,kk,i3), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                          ii = ii + 1
                       Else
                          If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                       End If
                    Else
                       If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                    End If
                 End Do
                 Write (Unit=iout, Fmt='()', Advance='yes', Err=1)
              End Do
           End Do
        End Do
     End Do
   Case ('unformatted')
     Do i3 = 1, n3
        Call write_unform (iout, afile, nt, nv, v(:,:,i3), iuse, kuse, afield(1)%rmiss, ifail)
        If (ifail /= 0) Return
     End Do
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_unrfs
!
!
!
 Subroutine write_eigs (iout, afile, ne, sv, nt, ifail)
!
! Prints eigenvalues
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,      Only: ap, dp, sp, zero, oneh
  Use data_io_constants, Only: cprcs, faccs,ffmts
  Use iofiles,           Only: ofile, &
                               itags
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: ne   ! - total number of eigenvalues -
  Integer, Intent(In) :: nt   ! - number of cases -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: sv ! - singular values -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: i  ! - EOF index -
  Integer :: up ! - precision -
!
  Real(Kind=rp) :: df     ! - number of cases -
  Real(Kind=rp) :: evalue ! - eigenvalue -
  Real(Kind=rp) :: tvar   ! - total variance -
  Real(Kind=rp) :: pvar   ! - percentage variance -
  Real(Kind=rp) :: svar   ! - cumulative variance -
  Real(Kind=rp) :: spvar  ! - cumulative percentage variance -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Sum, Trim
!
! Executable Statements
!
! Calculate total variance
  df = Real(nt-1, Kind=rp)
  tvar = Sum(sv(1:ne)**2/df)
!
! Calculate and print eigenvalues
  svar = zero
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     If (itags == 1) Then
        Call write_tag (iout, ifail, &
                        cpt_field='eigenvalues', cpt_nrow=ne, cpt_ncol=3, cpt_row='Mode', cpt_col='index')
        If (ifail /= 0) GoTo 1
     End If
     If (output_value(iout, c='Mode', ladv=.false.) /= 0) GoTo 1
     If (output_value(iout, c='Eigenvalue', ldelim=.true., ladv=.false.) /= 0) GoTo 1
     If (output_value(iout, c='% variance', ldelim=.true., ladv=.false.) /= 0) GoTo 1
     If (output_value(iout, c='Cum. % variance', ldelim=.true., ladv=.false.) /= 0) GoTo 1
     Do i = 1, ne
        evalue = sv(i)**2/df
        svar = svar + evalue
        pvar = oneh*evalue/tvar
        spvar = oneh*svar/tvar
        Write (Unit=iout, Fmt='(I0)', Advance='no', Err=1) i
        If (output_value(iout, r=evalue, ladv=.false.) /= 0) GoTo 1
        If (output_value(iout, r=pvar,   ladv=.false.) /= 0) GoTo 1
        If (output_value(iout, r=spvar,  ladv=.true. ) /= 0) GoTo 1
     End Do
   Case ('unformatted')
     Select Case (Trim(faccs(afile%cfile%ffmt%iacc)))
      Case ('sequential', 'transparent', 'stream')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           Do i = 1, ne
              evalue = sv(i)**2/df
              svar = svar + evalue
              pvar = oneh*evalue/tvar
              spvar = oneh*svar/tvar
              Write (Unit=iout, Err=1) evalue, pvar, spvar
           End Do
        Else
           Do i = 1, ne
              evalue = sv(i)**2/df
              svar = svar + evalue
              pvar = oneh*evalue/tvar
              spvar = oneh*svar/tvar
              Write (Unit=iout, Err=1) Real(evalue, Kind=ap), Real(pvar, Kind=ap), Real(spvar, Kind=ap)
           End Do
        End If
      Case ('direct')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           Do i = 1, ne
              evalue = sv(i)**2/df
              svar = svar + evalue
              pvar = oneh*evalue/tvar
              spvar = oneh*svar/tvar
              Write (Unit=iout, Err=1, Rec=i) evalue, pvar, spvar
           End Do
        Else
           Do i = 1, ne
              evalue = sv(i)**2/df
              svar = svar + evalue
              pvar = oneh*evalue/tvar
              spvar = oneh*svar/tvar
              Write (Unit=iout, Err=1, Rec=i) Real(evalue, Kind=ap), Real(pvar, Kind=ap), Real(spvar, Kind=ap)
           End Do
        End If
     End Select
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_eigs
!
!
!
 Subroutine write_load (iout, afile, afield, idstr, nfs, nls, ne, nv, iuse, eof, cfld, ifail, rlat, rlng, cstn)
!
! Prints spatial loadings
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants,   Only: ffmts, lstn
  Use data_time_constants, Only: nmn
  Use iofiles,             Only: ofile, &
                                 itags
  Use fields,              Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout  ! - output unit number -
  Integer, Intent(In) :: idstr ! - data structure -
  Integer, Intent(In) :: nfs   ! - number of fields -
  Integer, Intent(In) :: nls   ! - number of lagged fields -
  Integer, Intent(In) :: ne    ! - number of EOF Modes -
  Integer, Intent(In) :: nv    ! - number of variables -
!
  Character(Len=*), Intent(In) :: cfld ! - field description -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: eof ! - spatial loadings -
!
  Type(field), Dimension(:), Intent(In) :: afield ! - field -
!
! - Optional input arrays -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: rlng ! - longitudes -
!
  Character(Len=lstn), Dimension(:,:), Intent(In), Optional :: cstn ! - station names -
!
! Locals
!
! Local arrays
  Logical, Dimension(ne) :: luse ! - dummy array -
!
! Local scalars
  Integer :: j   ! - longitude/station/series index -
  Integer :: ie  ! - EOF index -
  Integer :: l   ! - field index -
  Integer :: ii  ! - grid index -
  Integer :: ij  ! - available grid index -
  Integer :: ilf ! - lagged field index -
  Integer :: ifd ! - field index -
!
  Character(Len=nmn+10) :: cMode ! - Mode string -
  Character(Len=     8) :: csp   ! - output label for latitude-longitude space -
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Print formatted data
  luse(:) = .true.
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Select Case (idstr)
! - gridded data -
      Case (1)
        Do l = 1, nfs*nls
           Do ie = 1, ne
              If (itags == 1) Then
                 If (Trim(afield(l)%var)/='undefined') Then
                    If (nls > 1) Then
                       Write (Unit=cMode, Fmt='(I0,A)') ie,'_('//Trim(afield(l)%cssn)//')'
                       Call write_tag (iout, ifail, &
                            cpt_field=Trim(afield(l)%var)//'_'//Trim(cfld), cpt_cMode=cMode, cpt_nrow=afield(l)%region%nlts, &
                            cpt_ncol=afield(l)%region%nlgs, cpt_row='Y', cpt_col='X')
                    Else
                       Call write_tag (iout, ifail, &
                            cpt_field=Trim(afield(l)%var)//'_'//Trim(cfld), cpt_mode=ie, cpt_nrow=afield(l)%region%nlts, &
                            cpt_ncol=afield(l)%region%nlgs, cpt_row='Y', cpt_col='X')
                    End If
                 Else
                    If (nls > 1) Then
                       Write (Unit=cMode, Fmt='(I0,A)') ie,'_('//Trim(afield(l)%cssn)//')'
                       Call write_tag (iout, ifail, &
                            cpt_field=Trim(cfld), cpt_cMode=cMode, cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                            cpt_row='Y', cpt_col='X')
                    Else
                       Call write_tag (iout, ifail, &
                            cpt_field=Trim(cfld), cpt_mode=ie, cpt_nrow=afield(l)%region%nlts, cpt_ncol=afield(l)%region%nlgs, &
                            cpt_row='Y', cpt_col='X')
                    End If
                 End If
                 If (ifail /= 0) GoTo 1
                 csp = ' '
              Else
                 Write (Unit=csp, Fmt='(I0)') ie
              End If
              If (l>1) Then
                 ij = Sum(afield(1:l-1)%nv)
                 ii = Sum(afield(1:l-1)%nva) + 1
              Else
                 ij = 0
                 ii = 1
              End If
              Call write_gridded_block_v10 (iout, csp, afield(l), rlat(:,l), rlng(:,l), ij, ifail, &
                   v=eof(ii:,ie),iuse=iuse(ii:))
           End Do
        End Do
! - station data -
      Case (2)
        ii = 0
        ij = 1
        Do ifd = 1, nfs
           l = (ifd - 1)*nls + 1
           If (itags == 1) Then
              If (Trim(afield(l)%var)/='undefined') Then
                 Call write_tag (iout, ifail, &
                      cpt_field=Trim(afield(l)%var)//'_'//Trim(cfld), cpt_nrow=afield(l)%nv*nls, cpt_ncol=ne, cpt_row='station', &
                      cpt_col='Mode')
              Else
                 Call write_tag (iout, ifail, &
                      cpt_field=Trim(cfld), cpt_nrow=afield(l)%nv*nls, cpt_ncol=ne, cpt_row='station', cpt_col='Mode')
              End If
           End If
           Do ilf = 1, nls
              l = (ifd - 1)*nls + ilf
              If (output_value(iout, c='Station', ladv=.false.) /= 0) GoTo 1
              If (output_value(iout, c='Latitude', ldelim=.true., ladv=.false.) /= 0) GoTo 1
              If (output_value(iout, c='Longitude', ldelim=.true., ladv=.false.) /= 0) GoTo 1
              If (nls > 1) Then
                 If (output_value(iout, c='Season', ldelim=.true., ladv=.false.) /= 0) GoTo 1
              End If
              Do ie = 1, ne
                 If (output_value(iout, i=ie, ldelim=.true., ladv=(ie==ne)) /= 0) GoTo 1
              End Do
              Do j = 1, afield(l)%nv
                 ii = ii + 1
                 If (output_value(iout, c=cstn(j,l), ladv=.false.) /= 0) GoTo 1
                 If (output_value(iout, r=rlat(j,l), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                 If (output_value(iout, r=rlng(j,l), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                 If (nls > 1) Then
                    If (output_value(iout, c=afield(l)%cssn, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                 End If
                 If (iuse(ij) == ii) Then
                    Do ie = 1, ne
                       If (output_value(iout, r=eof(ij,ie), ldelim=.true., ladv=(ie==ne)) /= 0) GoTo 1
                    End Do
                    ij = ij + 1
                 Else
                    Do ie = 1, ne
                       If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=(ie==ne)) /= 0) GoTo 1
                    End Do
                 End If
              End Do
           End Do
        End Do
! - unreferenced data -
      Case (3)
        If (itags == 1) Then
           If (Trim(afield(1)%var) /= 'undefined') Then
              Call write_tag (iout, ifail, &
                   cpt_field=Trim(afield(1)%var)//'_'//Trim(cfld), cpt_nrow=afield(1)%nv, cpt_ncol=ne, cpt_row='index', &
                   cpt_col='Mode')
           Else
              Call write_tag (iout, ifail, &
                   cpt_field=Trim(cfld), cpt_nrow=nv, cpt_ncol=ne, cpt_row='index', cpt_col='Mode')
           End If
        End If
        ii = 0
        ij = 1
        Do l = 1, nls
           If (output_value(iout, c='Index', ladv=.false.) /= 0) GoTo 1
           If (nls > 1) Then
              If (output_value(iout, c='Season', ldelim=.true., ladv=.false.) /= 0) GoTo 1
           End If
           Do ie = 1, ne
              If (output_value(iout, i=ie, ldelim=.true., ladv=(ie==ne)) /= 0) GoTo 1
           End Do
           Do j = 1, afield(l)%nv
              ii = ii + 1
              If (output_value(iout, c=cstn(j,l), ladv=.false.) /= 0) GoTo 1
              If (nls > 1) Then
                 If (output_value(iout, c=afield(l)%cssn, ldelim=.true., ladv=.false.) /= 0) GoTo 1
              End If
              If (iuse(ij) == ii) Then
                 Do ie = 1, ne
                    If (output_value(iout, r=eof(ij,ie), ldelim=.true., ladv=(ie==ne)) /= 0) GoTo 1
                 End Do
                 ij = ij + 1
              Else
                 Do ie = 1, ne
                    If (output_value(iout, r=afield(l)%rmiss, ldelim=.true., ladv=(ie==ne)) /= 0) GoTo 1
                 End Do
              End If
           End Do
        End Do
     End Select
   Case ('unformatted')
     Call write_unform (iout, afile, ne, nv, eof, iuse, luse, afield(1)%rmiss, ifail)
     If (ifail /= 0) Return
   Case ('GrADS')
     Select Case (idstr)
      Case (1) ! - gridded data -
        Call write_unform (iout, afile, ne, nv, eof, iuse, luse, afield(1)%rmiss, ifail)
        If (ifail /= 0) Return
        Call write_grads_grid_ctrl (iout+1, afile, afield(:), nfs, nls, ne, 1, rlat(:,:), rlng(:,:), ifail)
        If (ifail /= 0) Return
      Case (2) ! - station data -
        Call write_grads_stns (iout, afile, ne, nls, afield(1)%nva, eof, afield(:)%tprd%sdate, &
             cstn(:,1), rlat(:,1), rlng(:,1), iuse(:), luse(:), afield(1)%rmiss, ifail)
        If (ifail /= 0) Return
        Call write_grads_stns_ctrl (iout+1, afile, afield(:), nfs, nls, ne, 1, ifail)
        If (ifail /= 0) Return
     End Select
     If (ifail /= 0) Return
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_load
!
!
!
 Subroutine write_scor (iout, cscore, afile, nt, it0, prd1, ne, kuse, miss, ts, ifail)
!
! Prints temporal scores
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,      Only: ap, dp, sp
  Use data_io_constants, Only: cprcs, faccs, ffmts
  Use time,              Only: pprd, &
                               Operator(+), &
                               get_cdate
  Use arrays,            Only: awk, rwk
  Use iofiles,           Only: ofile, &
                               itags
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: ne   ! - number of EOF Modes -
  Integer, Intent(In) :: it0  ! - index of 0th case -
!
  Real(Kind=rp), Intent(In) :: miss ! - missing value flag -
!
  Character(Len=*), Intent(In) :: cscore ! - score description -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
  Type(pprd), Intent(In) :: prd1 ! - first period -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: ts ! - time scores -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
! Locals
!
! Local scalars
  Integer :: i  ! - EOF index -
  Integer :: k  ! - case index -
  Integer :: kk ! - available case index -
  Integer :: up ! - precision -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Trim
!
! Executable Statements
!
! Print headers
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     If (itags == 1) Then
        Call write_tag (iout, ifail, &
                        cpt_field=cscore, cpt_nrow=nt, cpt_ncol=ne, cpt_row='T', cpt_col='index')
        If (ifail /= 0) GoTo 1
     End If
!
! Print temporal scores
     Do i = 1, ne
        If (output_value(iout, i=i, ldelim=.true., ladv=(i==ne)) /= 0) GoTo 1
     End Do
     kk = 0
     Do k = 1, nt
        If (kuse(k)) kk = kk + 1
        cout = get_cdate(prd1+(it0+k-1), 1, 1)
        If (output_value(iout, c=cout, ladv=.false.) /= 0) GoTo 1
        If (kuse(k)) Then
           Do i = 1, ne
              If (output_value(iout, r=ts(i,kk), ldelim=.true., ladv=(i==ne)) /= 0) GoTo 1
           End Do
        Else
           Do i = 1, ne
              If (output_value(iout, r=miss, ldelim=.true., ladv=(i==ne)) /= 0) GoTo 1
           End Do
        End If
     End Do
   Case ('unformatted')
     Select Case (Trim(faccs(afile%cfile%ffmt%iacc)))
      Case ('sequential', 'transparent', 'stream')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           Do i = 1, ne
              kk = 0
              Do k = 1, nt
                 If (kuse(k)) Then
                    kk = kk + 1
                    rwk(k) = ts(i,kk)
                 Else
                    rwk(k) = miss
                 End If
              End Do
              Write (Unit=iout, Err=1) (rwk(k), k=1,nt)
           End Do
        Else
           Do i = 1, ne
              kk = 0
              Do k = 1, nt
                 If (kuse(k)) Then
                    kk = kk + 1
                    awk(k) = Real(ts(i,kk), Kind=ap)
                 Else
                    awk(k) = Real(miss, Kind=ap)
                 End If
              End Do
              Write (Unit=iout, Err=1) (awk(k), k=1,nt)
           End Do
        End If
      Case ('direct')
        Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           Do i = 1, ne
              kk = 0
              Do k = 1, nt
                 If (kuse(k)) Then
                    kk = kk + 1
                    rwk(k) = ts(i,kk)
                 Else
                    rwk(k) = miss
                 End If
              End Do
              Write (Unit=iout, Rec=i, Err=1) (rwk(k), k=1,nt)
           End Do
        Else
           Do i = 1, ne
              kk = 0
              Do k = 1, nt
                 If (kuse(k)) Then
                    kk = kk + 1
                    awk(k) = Real(ts(i,kk), Kind=ap)
                 Else
                    awk(k) = Real(miss, Kind=ap)
                 End If
              End Do
              Write (Unit=iout, Rec=i, Err=1) (awk(k),k=1,nt)
           End Do
        End If
     End Select
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_scor
!
!
!
 Subroutine write_cols (iout, afile, nrow, c1, v1, ifail, &
                        cpt_field, cpt_nrow, cpt_ncol, cpt_row, cpt_col,            &
                        c2, c3, c4, c5, c6, i2, i3, i4, i5, i6, v2, v3, v4, v5, v6, &
                        cpt_c, cpt_tagi, cpt_itag, cpt_tagr, cpt_rtag, cpt_tagc, cpt_ctag)
!
! Prints columns
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,      Only: ap, dp, sp
  Use data_io_constants, Only: cprcs, faccs, ffmts
  Use iofiles,           Only: ofile, &
                               itags
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout     ! - output unit number -
  Integer, Intent(In) :: nrow     ! - number of rows -
  Integer, Intent(In) :: cpt_nrow ! - number of rows -
  Integer, Intent(In) :: cpt_ncol ! - number of columns -
!
  Character(Len=*), Intent(In) :: c1        ! - first column description -
  Character(Len=*), Intent(In) :: cpt_field ! - field -
  Character(Len=*), Intent(In) :: cpt_row   ! - rows -
  Character(Len=*), Intent(In) :: cpt_col   ! - columns -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: cpt_c    ! - current category -
  Integer, Intent(In), Optional :: cpt_itag ! - extra integer tag value -
!
  Real(Kind=rp), Intent(In), Optional :: cpt_rtag ! - extra Real tag value -
!
  Character(Len=*), Intent(In), Optional :: c2       ! - second column description -
  Character(Len=*), Intent(In), Optional :: c3       ! - third column description -
  Character(Len=*), Intent(In), Optional :: c4       ! - fourth column description -
  Character(Len=*), Intent(In), Optional :: c5       ! - fifth column description -
  Character(Len=*), Intent(In), Optional :: c6       ! - sixth column description -
  Character(Len=*), Intent(In), Optional :: cpt_tagi ! - extra integer tag -
  Character(Len=*), Intent(In), Optional :: cpt_tagr ! - extra Real tag -
  Character(Len=*), Intent(In), Optional :: cpt_tagc ! - extra character tag -
  Character(Len=*), Intent(In), Optional :: cpt_ctag ! - extra character tag value -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: v1 ! - first column -
!
! - Optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: i2 ! - second column -
  Integer, Dimension(:), Intent(In), Optional :: i3 ! - third column -
  Integer, Dimension(:), Intent(In), Optional :: i4 ! - fourth column -
  Integer, Dimension(:), Intent(In), Optional :: i5 ! - fifth column -
  Integer, Dimension(:), Intent(In), Optional :: i6 ! - sixth column -
!
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: v2 ! - second column -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: v3 ! - third column -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: v4 ! - fourth column -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: v5 ! - fifth column -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: v6 ! - sixth column -
!
! Locals
!
! Local scalars
  Integer :: i  ! - row index -
  Integer :: up ! - precision -
!
  Logical :: ladv1 ! - advance after first columm? -
  Logical :: ladv2 ! - advance after second columm? -
  Logical :: ladv3 ! - advance after third columm? -
  Logical :: ladv4 ! - advance after fourth columm? -
  Logical :: ladv5 ! - advance after fifth columm? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Determine when to advance
  If (Present(c2) .or. Present(i2)) Then
     ladv1 = .false.
  Else
     ladv1 = .true.
  End If
  If (Present(c3) .or. Present(i3)) Then
     ladv1 = .false.
     ladv2 = .false.
  Else
     ladv2 = .true.
  End If
  If (Present(c4) .or. Present(i4)) Then
     ladv1 = .false.
     ladv2 = .false.
     ladv3 = .false.
  Else
     ladv3 = .true.
  End If
  If (Present(c5) .or. Present(i5)) Then
     ladv1 = .false.
     ladv2 = .false.
     ladv3 = .false.
     ladv4 = .false.
  Else
     ladv4 = .true.
  End If
  If (Present(c6) .or. Present(i6)) Then
     ladv1 = .false.
     ladv2 = .false.
     ladv3 = .false.
     ladv4 = .false.
     ladv5 = .false.
  Else
     ladv5 = .true.
  End If  
!
! Print tag line
  If (itags == 1) Then
     Call write_tag (iout, ifail, &
                     cpt_field=cpt_field, cpt_nrow=cpt_nrow, cpt_ncol=cpt_ncol, cpt_row=cpt_row, cpt_col=cpt_col, &
                     cpt_c=cpt_c, cpt_tagi=cpt_tagi, cpt_itag=cpt_itag, cpt_tagr=cpt_tagr, cpt_rtag=cpt_rtag,     &
                     cpt_tagc=cpt_tagc, cpt_ctag=cpt_ctag)
     If (ifail /= 0) GoTo 1
  End If
!
! Print columns
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     If (output_value(iout, c=c1, ladv=ladv1) /= 0) GoTo 1
     If (Present(c2)) Then
        If (output_value(iout, c=c2, ldelim=.true., ladv=ladv2) /= 0) GoTo 1
     End If
     If (Present(c3)) Then
        If (output_value(iout, c=c3, ldelim=.true., ladv=ladv3) /= 0) GoTo 1
     End If
     If (Present(c4)) Then
        If (output_value(iout, c=c4, ldelim=.true., ladv=ladv4) /= 0) GoTo 1
     End If
     If (Present(c5)) Then
        If (output_value(iout, c=c5, ldelim=.true., ladv=ladv5) /= 0) GoTo 1
     End If
     If (Present(c6)) Then
        If (output_value(iout, c=c6, ldelim=.true., ladv=.true.) /= 0) GoTo 1
     End If
     Do i = 1, nrow
        If (output_value(iout, r=v1(i), ldelim=.true., ladv=ladv1) /= 0) GoTo 1
        If (Present(i2)) Then
           If (output_value(iout, i=i2(i), ldelim=.true., ladv=ladv2) /= 0) GoTo 1
        End If
        If (Present(v2)) Then
           If (output_value(iout, r=v2(i), ldelim=.true., ladv=ladv2) /= 0) GoTo 1
        End If
        If (Present(i3)) Then
           If (output_value(iout, i=i3(i), ldelim=.true., ladv=ladv3) /= 0) GoTo 1
        End If
        If (Present(v3)) Then
           If (output_value(iout, r=v3(i), ldelim=.true., ladv=ladv3) /= 0) GoTo 1
        End If
        If (Present(i4)) Then
           If (output_value(iout, i=i4(i), ldelim=.true., ladv=ladv4) /= 0) GoTo 1
        End If
        If (Present(v4)) Then
           If (output_value(iout, r=v4(i), ldelim=.true., ladv=ladv4) /= 0) GoTo 1
        End If
        If (Present(i5)) Then
           If (output_value(iout, i=i5(i), ldelim=.true., ladv=ladv4) /= 0) GoTo 1
        End If
        If (Present(v5)) Then
           If (output_value(iout, r=v5(i), ldelim=.true., ladv=ladv5) /= 0) GoTo 1
        End If
        If (Present(i6)) Then
           If (output_value(iout, i=i6(i), ldelim=.true., ladv=.true.) /= 0) GoTo 1
        End If
        If (Present(v6)) Then
           If (output_value(iout, r=v6(i), ldelim=.true., ladv=.true.) /= 0) GoTo 1
        End If
     End Do
   Case ('unformatted')
     Select Case (Trim(faccs(afile%cfile%ffmt%iacc)))
      Case ('sequential', 'transparent', 'stream')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           Write (Unit=iout, Err=1) v1(:)
           If (Present(i2)) Write (Unit=iout, Err=1) i2(:)
           If (Present(v2)) Write (Unit=iout, Err=1) v2(:)
           If (Present(i3)) Write (Unit=iout, Err=1) i3(:)
           If (Present(v3)) Write (Unit=iout, Err=1) v3(:)
           If (Present(i4)) Write (Unit=iout, Err=1) i4(:)
           If (Present(v4)) Write (Unit=iout, Err=1) v4(:)
           If (Present(i5)) Write (Unit=iout, Err=1) i5(:)
           If (Present(v5)) Write (Unit=iout, Err=1) v5(:)
           If (Present(i6)) Write (Unit=iout, Err=1) i6(:)
           If (Present(v6)) Write (Unit=iout, Err=1) v6(:)
        Else
           Write (Unit=iout, Err=1) Real(v1(:), Kind=ap)
           If (Present(i2)) Write (Unit=iout, Err=1) Real(i2(:), Kind=ap)
           If (Present(v2)) Write (Unit=iout, Err=1) Real(v2(:), Kind=ap)
           If (Present(i3)) Write (Unit=iout, Err=1) Real(i3(:), Kind=ap)
           If (Present(v3)) Write (Unit=iout, Err=1) Real(v3(:), Kind=ap)
           If (Present(i4)) Write (Unit=iout, Err=1) Real(i4(:), Kind=ap)
           If (Present(v4)) Write (Unit=iout, Err=1) Real(v4(:), Kind=ap)
           If (Present(i5)) Write (Unit=iout, Err=1) Real(i5(:), Kind=ap)
           If (Present(v5)) Write (Unit=iout, Err=1) Real(v5(:), Kind=ap)
           If (Present(i6)) Write (Unit=iout, Err=1) Real(i6(:), Kind=ap)
           If (Present(v6)) Write (Unit=iout, Err=1) Real(v6(:), Kind=ap)
        End If
      Case ('direct')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           Write (Unit=iout, Err=1, Rec=2) v1(:)
           If (Present(i2)) Write (Unit=iout, Err=1, Rec=3) i2(:)
           If (Present(v2)) Write (Unit=iout, Err=1, Rec=3) v2(:)
           If (Present(i3)) Write (Unit=iout, Err=1, Rec=3) i3(:)
           If (Present(v3)) Write (Unit=iout, Err=1, Rec=3) v3(:)
           If (Present(i4)) Write (Unit=iout, Err=1, Rec=3) i4(:)
           If (Present(v4)) Write (Unit=iout, Err=1, Rec=3) v4(:)
           If (Present(i5)) Write (Unit=iout, Err=1, Rec=3) i5(:)
           If (Present(v5)) Write (Unit=iout, Err=1, Rec=3) v5(:)
           If (Present(i6)) Write (Unit=iout, Err=1, Rec=3) i6(:)
           If (Present(v6)) Write (Unit=iout, Err=1, Rec=3) v6(:)
        Else
           Write (Unit=iout, Err=1, Rec=2) Real(v1(:), Kind=ap)
           If (Present(i2)) Write (Unit=iout, Err=1, Rec=3) Real(i2(:), Kind=ap)
           If (Present(v2)) Write (Unit=iout, Err=1, Rec=3) Real(v2(:), Kind=ap)
           If (Present(i3)) Write (Unit=iout, Err=1, Rec=3) Real(i3(:), Kind=ap)
           If (Present(v3)) Write (Unit=iout, Err=1, Rec=3) Real(v3(:), Kind=ap)
           If (Present(i4)) Write (Unit=iout, Err=1, Rec=3) Real(i4(:), Kind=ap)
           If (Present(v4)) Write (Unit=iout, Err=1, Rec=3) Real(v4(:), Kind=ap)
           If (Present(i5)) Write (Unit=iout, Err=1, Rec=3) Real(i5(:), Kind=ap)
           If (Present(v5)) Write (Unit=iout, Err=1, Rec=3) Real(v5(:), Kind=ap)
           If (Present(i6)) Write (Unit=iout, Err=1, Rec=3) Real(i6(:), Kind=ap)
           If (Present(v6)) Write (Unit=iout, Err=1, Rec=3) Real(v6(:), Kind=ap)
        End If
     End Select
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_cols
!
!
!
 Subroutine write_regr (iout, afile, xfile, xfield, nx, iusex, yfile, yfield, ny, iusey, b, ifail, &
            rlatx, rlngx, cstnx, rlaty, rlngy, cstny, b0)
!
! Prints regression coefficients
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,      Only: ap, dp, sp
  Use data_io_constants, Only: cprcs, faccs, ffmts, lstn
  Use arrays,            Only: awk, rwk
  Use iofiles,           Only: ifile, ofile, &
                               itags
  Use space,             Only: make_coors
  Use fields,            Only: field, tfield
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout   ! - output unit number -
  Integer, Intent(In) :: nx     ! - number of X variables -
  Integer, Intent(In) :: ny     ! - number of Y variables -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
  Type(ifile), Intent(In) :: xfile ! - X input file -
  Type(ifile), Intent(In) :: yfile ! - Y input file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iusex ! - X variable flags -
  Integer, Dimension(:), Intent(In) :: iusey ! - Y variable flags -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: b ! - regression coefficients -
!
  Type(field), Dimension(:), Intent(In) :: xfield ! - X field -
  Type(field), Dimension(:), Intent(In) :: yfield ! - Y field -
!
! - Optional input arrays -
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: b0 ! - regression constants -
!
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: rlatx ! - X latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: rlngx ! - X longitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: rlaty ! - Y latitudes -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: rlngy ! - Y longitudes -
!
  Character(Len=lstn), Dimension(:,:), Intent(In), Optional :: cstnx ! - X station names -
  Character(Len=lstn), Dimension(:,:), Intent(In), Optional :: cstny ! - Y station names -
!
! Locals
!
! Local scalars
  Integer :: i, j   ! - indices -
  Integer :: ib0    ! - b0 Present indicator -
  Integer :: ifdx   ! - X field index -
  Integer :: ifdy   ! - Y field index -
  Integer :: ilfx   ! - X lagged field index -
  Integer :: ilfy   ! - Y lagged field index -
  Integer :: lx     ! - X field / lagged field index -
  Integer :: ly     ! - Y field / lagged field index -
  Integer :: i1, j1 ! - indices -
  Integer :: i2, j2 ! - indices -
  Integer :: iy, jy ! - indices -
  Integer :: up     ! - precision -
!
  Real(Kind=ap) :: miss ! - missing value flag -
!
  Character(Len=16) :: cxcol ! - X column -
  Character(Len=16) :: cyrow ! - Y row -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Real, Trim
!
! Executable Statements
!
! Print regression coefficients
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Allocate (tfield(xfile%nfl))
     tfield(:) = xfield(:)
     tfield(:)%Unit = 'none'
     tfield(:)%unitc = 'none'
     Select Case (xfile%idstr)
      Case (1) ! - gridded X data -
        i2 = 1
        i1 = 0
        Select Case (yfile%idstr)
         Case (1) ! - gridded Y data -
           Do ifdy = 1, yfile%nfs
              Do ilfy = 1, yfile%nls
                 ly = (ifdy - 1)*yfile%nls + ilfy
                 Do iy = 1, yfield(ly)%region%nlts
                    Do jy = 1, yfield(ly)%region%nlgs
                       i1 = i1 + 1
                       If (yfile%nfs > 1) Then
                          If (yfile%nls > 1) Then
                             cout = Trim(make_coors(rlaty(iy,ly), rlngy(jy,ly)))//'; '//Trim(yfield(ly)%var)//&
                                                                                 &'; '//Trim(yfield(ly)%cssn)
                          Else
                             cout = Trim(make_coors(rlaty(iy,ly), rlngy(jy,ly)))//'; '//Trim(yfield(ly)%var)
                          End If
                       Else
                          If (yfile%nls > 1) Then
                             cout = Trim(make_coors(rlaty(iy,ly), rlngy(jy,ly)))//'; '//Trim(yfield(ly)%cssn)
                          Else
                             cout = Trim(make_coors(rlaty(iy,ly), rlngy(jy,ly)))
                          End If
                       End If
                       Call write_regr_grid (ifail)
                       If (ifail /= 0) GoTo 1
                    End Do
                 End Do
              End Do
           End Do
         Case (2, 3) ! - nongridded Y data -
           Do ifdy = 1, yfile%nfs
              Do ilfy = 1, yfile%nls
                 ly = (ifdy - 1)*yfile%nls + ilfy
                 Do i = 1, yfield(ly)%nv
                    i1 = i1 + 1
                    If (yfile%nfs > 1) Then
                       If (yfile%nls > 1) Then
                          cout = Trim(cstny(i,ly))//'; '//Trim(yfield(ly)%var)//'; '//Trim(yfield(ly)%cssn)
                       Else
                          cout = Trim(cstny(i,ly))//'; '//Trim(yfield(ly)%var)
                       End If
                    Else
                       If (yfile%nls > 1) Then
                          cout = Trim(cstny(i,ly))//'; '//Trim(yfield(ly)%cssn)
                       Else
                          cout = Trim(cstny(i,ly))
                       End If
                    End If
                    Call write_regr_grid (ifail)
                    If (ifail /= 0) GoTo 1
                 End Do
              End Do
           End Do
        End Select
      Case (2, 3) ! - nongridded X data -
        Select Case (xfile%idstr)
         Case (2) ! - station X data -
           cxcol = 'X station'
         Case (3) ! - unreferenced X data -
           cxcol = 'X index'
        End Select
        If (Present(b0)) Then
           ib0 = 1
        Else
           ib0 = 0
        End If
        Select Case (yfile%idstr)
         Case (1) ! - gridded Y data -
           cyrow = 'Y grid'
           Do ifdx = 1, xfile%nfs
              Do ilfx = 1, xfile%nls
                 lx = (ifdx - 1)*xfile%nls + ilfx
                 Call write_regr_headers (ifail)
                 If (ifail /= 0) GoTo 1
                 Select Case (xfile%idstr)
                  Case (2) ! - station X data -
                    If (output_value(iout, c='cpt:Y', ldelim=.true., ladv=.false.) /= 0) GoTo 1
                    Do i = 1, xfield(lx)%nv
                       If (output_value(iout, r=rlatx(i,lx), ldelim=.true., ladv=(i==xfield(lx)%nv)) /= 0) GoTo 1
                    End Do
                    If (output_value(iout, c='cpt:X', ldelim=.true., ladv=.false.) /= 0) GoTo 1
                    Do i = 1, xfield(lx)%nv
                       If (output_value(iout, r=rlngx(i,lx), ldelim=.true., ladv=(i==xfield(lx)%nv)) /= 0) GoTo 1
                    End Do
                 End Select
                 iy = 1
                 jy = 0
                 Do ifdy = 1, yfile%nfs
                    Do ilfy = 1, yfile%nls
                       ly = (ifdy - 1)*yfile%nls + ilfy
                       Do i = 1, yfield(ly)%region%nlts
                          Do j = 1, yfield(ly)%region%nlgs
                             jy = jy + 1
                             If (output_value(iout, r=rlaty(i,ly), ladv=.false.) /= 0) GoTo 1
                             If (output_value(iout, r=rlngy(i,ly), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                             If (yfile%nfs > 1) Then
                                If (output_value(iout, c=yfield(ly)%var, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                             End If
                             If (yfile%nls > 1) Then
                                If (output_value(iout, c=yfield(ly)%cssn, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                             End If
                             Call write_regr_nongrid (ifail)
                             If (ifail /= 0) GoTo 1
                          End Do
                       End Do
                    End Do
                 End Do
              End Do
           End Do
         Case (2, 3) ! - nongridded Y data -
           Select Case (yfile%idstr)
            Case (2) ! - station X data -
              cyrow = 'Y station'
            Case (3) ! - unreferenced X data -
              cyrow = 'Y index'
           End Select
           Do ifdx = 1, xfile%nfs
              Do ilfx = 1, xfile%nls
                 lx = (ifdx - 1)*xfile%nls + ilfx
                 Call write_regr_headers (ifail)
                 If (ifail /= 0) GoTo 1
                 iy = 1
                 jy = 0
                 Do ifdy = 1, yfile%nfs
                    Do ilfy = 1, yfile%nls
                       ly = (ifdy - 1)*yfile%nls + ilfy
                       Do j = 1, yfield(ly)%nv
                          jy = jy + 1
                          If (output_value(iout, c=cstny(j,ly), ladv=.false.) /= 0) GoTo 1
                          Select Case (yfile%idstr)
                           Case (2) ! - station Y data -
                             If (output_value(iout, r=rlaty(j,ly), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                             If (output_value(iout, r=rlngy(j,ly), ldelim=.true., ladv=.false.) /= 0) GoTo 1
                          End Select
                          If (yfile%nfs > 1) Then
                             If (output_value(iout, c=yfield(ly)%var, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                          End If
                          If (yfile%nls > 1) Then
                             If (output_value(iout, c=yfield(ly)%cssn, ldelim=.true., ladv=.false.) /= 0) GoTo 1
                          End If
                          Call write_regr_nongrid (ifail)
                          If (ifail /= 0) GoTo 1
                       End Do
                    End Do
                 End Do
              End Do
           End Do
        End Select
     End Select
     Deallocate (tfield)
   Case ('unformatted')
     Select Case (Trim(faccs(afile%cfile%ffmt%iacc)))
      Case ('sequential', 'transparent', 'stream')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           i2 = 1
           Do i1 = 1, ny
              If (iusey(i2) == i1) Then
                 j2 = 1
                 Do j1 = 1, nx
                    If (iusex(j2) == j1) Then
                       rwk(j1) = b(j2,i2)
                       j2 = j2 + 1
                    Else
                       rwk(j1) = xfield(1)%rmiss
                    End If
                 End Do
                 Write (Unit=iout, Err=1) (rwk(j1), j1=1,nx)
                 i2 = i2 + 1
              Else
                 Write (Unit=iout, Err=1) (xfield(1)%rmiss, j1=1,nx)
              End If
           End Do
        Else
           miss = Real(xfield(1)%rmiss, Kind=ap)
           i2 = 1
           Do i1 = 1, ny
              If (iusey(i2) == i1) Then
                 j2 = 1
                 Do j1 = 1, nx
                    If (iusex(j2) == j1) Then
                       awk(j1) = Real(b(j2,i2), Kind=ap)
                       j2 = j2 + 1
                    Else
                       awk(j1) = miss
                    End If
                 End Do
                 Write (Unit=iout, Err=1) (awk(j1), j1=1,nx)
                 i2 = i2 + 1
              Else
                 Write (Unit=iout, Err=1) (miss, j1=1,nx)
              End If
           End Do
        End If
      Case ('direct')
        Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
         Case ('double')
           up = dp
         Case ('single')
           up = sp
        End Select
        If (up == rp) Then
           i2 = 1
           Do i1 = 1, ny
              If (iusey(i2) == i1) Then
                 j2 = 1
                 Do j1 = 1, nx
                    If (iusex(j2) == j1) Then
                       rwk(j1) = b(j2,i2)
                       j2 = j2 + 1
                    Else
                       rwk(j1) = xfield(1)%rmiss
                    End If
                 End Do
                 Write (Unit=iout, Rec=i1, Err=1) (rwk(j1), j1=1,nx)
                 i2 = i2 + 1
              Else
                 Write (Unit=iout, Rec=i1, Err=1) (xfield(1)%rmiss, j1=1,nx)
              End If
           End Do
        Else
           miss = Real(xfield(1)%rmiss, Kind=ap)
           i2 = 1
           Do i1 = 1, ny
              If (iusey(i2) == i1) Then
                 j2 = 1
                 Do j1 = 1, nx
                    If (iusex(j2) == j1) Then
                       awk(j1) = Real(b(j2,i2), Kind=ap)
                       j2 = j2 + 1
                    Else
                       awk(j1) = miss
                    End If
                 End Do
                 Write (Unit=iout, Rec=i1, Err=1) (awk(j1), j1=1,nx)
                 i2 = i2 + 1
              Else
                 Write (Unit=iout, Rec=i1, Err=1) (miss, j1=1,nx)
              End If
           End Do
        End If
     End Select
  End Select
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
!
 Contains
!
!
  Subroutine write_regr_headers (ifail)
!
! Prints headers for regression data
!
! Arguments
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Print headers for regression data
  If (itags == 1) Then
     Call write_tag (iout, ifail, &
          cpt_field=xfield(lx)%var, cpt_nrow=ny, cpt_ncol=xfield(lx)%nv+ib0, &
          cpt_row=Trim(cyrow), cpt_col=Trim(cxcol), cpt_units='none', cpt_missing=xfield(lx)%rmiss)
     If (ifail /= 0) GoTo 1
  End If
  Select Case (yfile%idstr)
   Case (1, 2) ! - gridded and station Y data -
     If (yfile%idstr == 2) Then
        If (output_value(iout, c='Station',   ldelim=.true., ladv=.false.) /= 0) GoTo 1
     End If
     If (output_value(iout, c='Latitude',  ldelim=.true., ladv=.false.) /= 0) GoTo 1
     If (output_value(iout, c='Longitude', ldelim=.true., ladv=.false.) /= 0) GoTo 1
   Case (3) ! - unreferenced Y data -
     If (output_value(iout, c='Index', ldelim=.true., ladv=.false.) /= 0) GoTo 1
  End Select
  If (yfile%nfs > 1) Then
     If (output_value(iout, c='Field', ldelim=.true., ladv=.false.) /= 0) GoTo 1
  End If
  If (yfile%nls > 1) Then
     If (output_value(iout, c='Season', ldelim=.true., ladv=.false.) /= 0) GoTo 1
  End If
  If (Present(b0)) Then
     If (output_value(iout, c='Constant', ldelim=.true., ladv=.false.) /= 0) GoTo 1
  End If
  Do i = 1, xfield(lx)%nv
     If (output_value(iout, c=cstnx(i,lx), ldelim=.true., ladv=(i==xfield(lx)%nv)) /= 0) GoTo 1
  End Do
  Select Case (xfile%idstr)
   Case (2) ! - station X data -
     If (output_value(iout, c='cpt:Y', ldelim=.true., ladv=.false.) /= 0) GoTo 1
     Do i = 1, xfield(lx)%nv
        If (output_value(iout, r=rlatx(i,lx), ldelim=.true., ladv=(i==xfield(lx)%nv)) /= 0) GoTo 1
     End Do
     If (output_value(iout, c='cpt:X', ldelim=.true., ladv=.false.) /= 0) GoTo 1
     Do i = 1, xfield(lx)%nv
        If (output_value(iout, r=rlngx(i,lx), ldelim=.true., ladv=(i==xfield(lx)%nv)) /= 0) GoTo 1
     End Do
  End Select
!
  ifail = 0
  Return
!
! Errors
1 ifail = 1
!
  Return
  End Subroutine write_regr_headers
!
!
!
  Subroutine write_regr_grid (ifail)
!
! Prints gridded regression data
!
! Arguments
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Print gridded regression data
  If (itags == 1) Call write_tag (iout, ifail, &
                       cpt_name=cout) 
  If (i1 == iusey(i2)) Then
     Call write_grid (iout, afile, tfield, xfile%nfs, xfile%nls, 1, nx, b(:,i2:i2), iusex, (/.true./), &
          rlatx, rlngx, ifail)
     i2 = i2 + 1
  Else
     Call write_grid (iout, afile, tfield, xfile%nfs, xfile%nls, 1, nx, b(:,1:1), iusex, (/.false./), &
          rlatx, rlngx, ifail)
  End If
!
  Return
  End Subroutine write_regr_grid
!
!
!
  Subroutine write_regr_nongrid (ifail)
!
! Prints nongridded regression data
!
! Arguments
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: ij ! - indices -
  Integer :: ix ! - indices -
  Integer :: jx ! - indices -
!
! Executable Statements
!
! Print nongridded regression data
  If (iusey(iy) == jy) Then
     If (Present(b0)) Then
        If (output_value(iout, r=b0(iy), ldelim=.true., ladv=.false.) /= 0) GoTo 1
     End If
     ix = 1
     jx = 0
     If (lx > 1) Then
        ix = ix + Sum(xfield(1:lx-1)%nva)
        jx = jx + Sum(xfield(1:lx-1)%nv)
     End If
     Do ij = 1, xfield(lx)%nv
        jx = jx + 1
        If (iusex(ix) == jx) Then
           If (output_value(iout, r=b(ix,iy), ldelim=.true., ladv=(ij==xfield(lx)%nv)) /= 0) GoTo 1
           ix = ix + 1
        Else
           If (output_value(iout, r=xfield(lx)%rmiss, ldelim=.true., ladv=(ij==xfield(lx)%nv)) /= 0) GoTo 1
        End If
     End Do
     iy = iy + 1
  Else
     Do ix = 1, xfield(lx)%nv + ib0
        If (output_value(iout, r=xfield(lx)%rmiss, ldelim=.true., ladv=(ix==xfield(lx)%nv+ib0)) /= 0) GoTo 1
     End Do
  End If
!
  ifail = 0
  Return
!
! Errors
1 ifail = 1
!
  Return
  End Subroutine write_regr_nongrid
 End Subroutine write_regr
!
!
!
 Subroutine write_gridded_block_v10 (iout, csp, afield, rlat, rlng, ioff, ifail, &
             v, iuse)
!
! Outputs a block of gridded data in v10 format without the tag line
!
! On exit Returns:
!    0 If no errors
!    1 If problem writing to file
!
! Modules
  Use data_numbers, Only: zero, r180, r360
  Use fields,       Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: ioff ! - offset -
!
  Character(Len=*), Intent(In) :: csp ! - output label for latitude-longitude space -
!
  Type(field), Intent(In) :: afield ! - field -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:), Intent(In) :: rlng ! - longitudes -
!
! - Optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:), Intent(In), Optional :: v ! - data -
!
! Locals
!
! Local scalars
  Integer :: i  ! - latitude index -
  Integer :: j  ! - longitude index -
  Integer :: ii ! - available series index -
  Integer :: ij ! - current variable -
!
  Real(Kind=rp) :: roff ! - longitude offset -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic All, Present
!
! Executable Statements
!
! Print gridded data
  ij = ioff
  ii = 1
  If (All(rlng(1:afield%region%nlgs) < -r180)) Then
     roff = -r360
  Else If (All(rlng(1:afield%region%nlgs) < r360)) Then
     roff = zero
  Else
     roff = r360
  End If
  If (output_value(iout, c=csp, ladv=.false.) /= 0) GoTo 1
  Do j = 1, afield%region%nlgs
     If (output_value(iout, r=rlng(j)-roff, ldelim=.true., ladv=(j==afield%region%nlgs)) /= 0) GoTo 1
  End Do
  Do i = 1, afield%region%nlts
     If (output_value(iout, r=rlat(i), ldelim=.false., ladv=.false.) /= 0) GoTo 1
     Do j = 1, afield%region%nlgs
        If (Present(v)) Then
           ij = ij + 1
           If (iuse(ii) == ij) Then
              If (output_value(iout, r=v(ii), ldelim=.true., ladv=(j==afield%region%nlgs)) /= 0) GoTo 1
              ii = ii + 1
           Else
              If (output_value(iout, r=afield%rmiss, ldelim=.true., ladv=(j==afield%region%nlgs)) /= 0) GoTo 1
           End If
        Else
           If (output_value(iout, r=afield%rmiss, ldelim=.true., ladv=(j==afield%region%nlgs)) /= 0) GoTo 1
        End If
     End Do
  End Do
  ifail = 0
!
  Return
!
! Error
1 ifail = 1
!
  Return
 End Subroutine write_gridded_block_v10
!
!
!
 Subroutine write_unform (iout, afile, nt, nv, v, iuse, kuse, miss, ifail)
!
! Outputs unformatted data
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,      Only: ap, dp, sp
  Use data_io_constants, Only: cprcs, faccs
  Use arrays,            Only: rwk
  Use iofiles,           Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nv   ! - number of series -
!
  Real(Kind=rp), Intent(In) :: miss ! - missing value flag -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v ! - data -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
! Locals
!
! Local scalars
  Integer :: i  ! - series index -
  Integer :: j  ! - available series index -
  Integer :: k  ! - time index -
  Integer :: kk ! - available time index -
  Integer :: up ! - precision -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Identify precision
  Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
   Case ('double')
     up = dp
   Case ('single')
     up = sp
  End Select
!
! Print unformatted data
  kk = 0
  Do k = 1, nt
     If (kuse(k)) kk = kk + 1
     j = 1
     Do i = 1, nv
        If (kuse(k)) Then
           If (iuse(j) == i) Then
              rwk(i) = v(j,kk)
              j = j + 1
           Else
              rwk(i) = miss
           End If
        Else
           rwk(i) = miss
        End If
     End Do
     Select Case (Trim(faccs(afile%cfile%ffmt%iacc)))
      Case ('sequential', 'transparent', 'stream')
        If (up == rp) Then
           Write (Unit=iout, Err=1) (rwk(i),i=1,nv)
        Else
           Write (Unit=iout, Err=1) (Real(rwk(i), Kind=ap),i=1,nv)
        End If
      Case ('direct')
        If (up == rp) Then
           Write (Unit=iout, Rec=k, Err=1) (rwk(i),i=1,nv)
        Else
           Write (Unit=iout, Rec=k, Err=1) (Real(rwk(i), Kind=ap),i=1,nv)
        End If
     End Select
  End Do
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_unform
!
!
!
 Subroutine write_grads_stns (iout, afile, nt, nls, nst, v, sdate, cstn, rlat, rlng, iuse, kuse, miss, ifail)
!
! Outputs station data in unformatted GrADS format (assumes that the stations are identical for each field and lagged field)
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,        Only: ap, dp, sp, zero
  Use data_io_constants,   Only: cprcs
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr, nmn
  Use time,                Only: pdate, &
                                 iseq, &
                                 date_diff
  Use iofiles,             Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nls  ! - number of lagged-fields -
  Integer, Intent(In) :: nst  ! - number of stations -
!
  Real(Kind=rp), Intent(In) :: miss ! - missing value flag -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:), Intent(In) :: rlng ! - longitudes -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v ! - data -
!
  Character(Len=*), Dimension(:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(pdate), Dimension(:), Intent(In) :: sdate ! - start dates -
!
! Locals
!
! Local scalars
  Integer :: i   ! - series index -
  Integer :: ii  ! - series index -
  Integer :: ilf ! - lagged-field index -
  Integer :: j   ! - available series index -
  Integer :: k   ! - time index -
  Integer :: kk  ! - available time index -
  Integer :: up  ! - precision -
!
  Real(Kind=rp) :: d    ! - data -
  Real(Kind=rp) :: rlag ! - lag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Trim
!
! Executable Statements
!
! Identify precision
!
! Identify precision
  Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
   Case ('double')
     up = dp
   Case ('single')
     up = sp
  End Select
!
! Print unformatted data
  kk = 0
  Do k = 1, nt
     If (kuse(k)) kk = kk + 1
     rlag = zero
     Do ilf = 1, nls
        If (ilf > 1) Then
           Select Case (iseq)
            Case (isq_yr) ! - yearly -
              rlag = date_diff(sdate(1), sdate(ilf), isq_mn)/Real(nmn, Kind=rp)
            Case (isq_sn) ! - seasonal -
              rlag = date_diff(sdate(1), sdate(ilf), isq_mn)/Real(3, Kind=rp)
            Case (isq_mn, isq_so) ! - monthly -
              rlag = date_diff(sdate(1), sdate(ilf), isq_mn)
            Case (1:) ! - daily and sub-seasonal -
              rlag = date_diff(sdate(1), sdate(ilf), 1)
           End Select
        End If
        Do i = 1, nst
           j = iuse(i)
           If (up == rp) Then
              Write (Unit=iout) cstn(j)(1:8), rlat(j), rlng(j), zero, 1, 1
           Else
              Write (Unit=iout) cstn(j)(1:8), Real(rlat(j), Kind=ap), Real(rlng(j), Kind=ap), Real(zero, Kind=ap), 1, 1
           End If
           If (kuse(k)) Then
              ii = nst*(ilf - 1) + i
              d = v(ii,kk)
           Else
              d = miss
           End If
           If (up == rp) Then
              Write (Unit=iout, Err=1) d
           Else
              Write (Unit=iout, Err=1) Real(d, Kind=ap)
           End If
        End Do
     End Do
     If (up == rp) Then
        Write (Unit=iout) cstn(1)(1:8), rlat(1), rlng(1), rlag, 0, 1
     Else
        Write (Unit=iout) cstn(1)(1:8), Real(rlat(1), Kind=ap), Real(rlng(1), Kind=ap), Real(rlag, Kind=ap), 0, 1
     End If
  End Do
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_grads_stns
!
!
!
 Subroutine write_grads_stns3 (iout, afile, nt, nls, nst, n3, v, sdate, cstn, rlat, rlng, iuse, kuse, miss, ifail)
!
! Outputs station data in unformatted GrADS format (assumes that the stations are identical for each field and lagged field)
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_numbers,        Only: ap, dp, sp, zero
  Use data_io_constants,   Only: cprcs
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr, nmn
  Use time,                Only: pdate, &
                                 iseq, &
                                 date_diff
  Use iofiles,             Only: ofile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nls  ! - number of lagged-fields -
  Integer, Intent(In) :: nst  ! - number of stations -
  Integer, Intent(In) :: n3   ! - third dimension -
!
  Real(Kind=rp), Intent(In) :: miss ! - missing value flag -
!
  Type(ofile), Intent(In) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: rlat ! - latitudes -
  Real(Kind=rp), Dimension(:), Intent(In) :: rlng ! - longitudes -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: v ! - data -
!
  Character(Len=*), Dimension(:), Intent(In) :: cstn ! - station names -
!
  Logical, Dimension(:), Intent(In) :: kuse ! - cases used? -
!
  Type(pdate), Dimension(:), Intent(In) :: sdate ! - start dates -
!
! Locals
!
! Local scalars
  Integer :: i   ! - series index -
  Integer :: ii  ! - series index -
  Integer :: i3  ! - third dimension index -
  Integer :: ilf ! - lagged-field index -
  Integer :: j   ! - available series index -
  Integer :: k   ! - time index -
  Integer :: kk  ! - available time index -
  Integer :: up  ! - precision -
!
  Real(Kind=rp) :: d    ! - data -
  Real(Kind=rp) :: rlag ! - lag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Trim
!
! Executable Statements
!
! Identify precision
!
! Identify precision
  Select Case (Trim(cprcs(afile%cfile%ffmt%iprc)))
   Case ('double')
     up = dp
   Case ('single')
     up = sp
  End Select
!
! Print unformatted data
  kk = 0
  Do k = 1, nt
     If (kuse(k)) kk = kk + 1
     rlag = zero
     Do ilf = 1, nls
        If (ilf > 1) Then
           Select Case (iseq)
            Case (isq_yr) ! - yearly -
              rlag = date_diff(sdate(1), sdate(ilf), isq_mn)/Real(nmn, Kind=rp)
            Case (isq_sn) ! - seasonal -
              rlag = date_diff(sdate(1), sdate(ilf), isq_mn)/Real(3, Kind=rp)
            Case (isq_mn, isq_so) ! - monthly -
              rlag = date_diff(sdate(1), sdate(ilf), isq_mn)
            Case (1:) ! - daily and sub-seasonal -
              rlag = date_diff(sdate(1), sdate(ilf), 1)
           End Select
        End If
        Do i = 1, nst
           j = iuse(i)
           If (up == rp) Then
              Write (Unit=iout) cstn(j)(1:8), rlat(j), rlng(j), zero, 1, 1
           Else
              Write (Unit=iout) cstn(j)(1:8), Real(rlat(j), Kind=ap), Real(rlng(j), Kind=ap), Real(zero, Kind=ap), 1, 1
           End If
           Do i3 = 1, n3
              If (kuse(k)) Then
                 ii = nst*(ilf - 1) + i
                 d = v(ii,kk,i3)
              Else
                 d = miss
              End If
              If (up == rp) Then
                 Write (Unit=iout, Err=1) d
              Else
                 Write (Unit=iout, Err=1) Real(d, Kind=ap)
              End If
          End Do
        End Do
     End Do
     If (up == rp) Then
        Write (Unit=iout) cstn(1)(1:8), rlat(1), rlng(1), rlag, 0, 1
     Else
        Write (Unit=iout) cstn(1)(1:8), Real(rlat(1), Kind=ap), Real(rlng(1), Kind=ap), Real(rlag, Kind=ap), 0, 1
     End If
  End Do
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_grads_stns3
!
!
!
 Subroutine write_grads_grid_ctrl (iout, cfile, cfield, nfs, nls, nt, nvs, rlat, rlng, ifail, &
            cvar)
!
! Outputs a GrADS control file for gridded data
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants,   Only: cendian, grads_ext
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr
  Use time,                Only: iseq, &
                                 get_cdate
  Use iofiles,             Only: ofile
  Use fields,              Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nvs  ! - number of variables -
!
  Type(ofile), Intent(In) :: cfile ! - control file -
!
! Input arrays
  Real(Kind=rp), Intent(In) :: rlat(:,:) ! - latitudes -
  Real(Kind=rp), Intent(In) :: rlng(:,:) ! - longitudes -
!
  Type(field), Intent(In) :: cfield(:) ! - field specifications -
!
! - Optional input arrays -
  Character(Len=*), Dimension(:), Intent(In), Optional :: cvar ! - variable names -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Character(Len=18) :: cdate ! - format statement -
  Character(Len=32) :: cout1 ! - output field -
  Character(Len=32) :: cout2 ! - output field -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustL, Trim
!
! Executable Statements
!
! Open control file
  Open (Unit=iout, File=Trim(cfile%cfile%fdir)//Trim(cfile%cfile%fname)//grads_ext, &
        Action='write', IOstat=ifail, Form='formatted', Status='unknown')
  If (ifail /= 0) Then
     ifail = 1
     Return
  End If
!
! Write control file
  Write (Unit=iout, Fmt='(A)', Err=1) 'DSET    '//Trim(cfile%cfile%ffile)
  Write (Unit=iout, Fmt='(A)', Err=1) 'TITLE   '//Trim(cfile%cfile%desc)
  Write (Unit=cout1, Fmt=*) cfield(1)%rmiss
  Write (Unit=iout, Fmt='(2A)', Err=1) 'UNDEF   ', Trim(AdjustL(cout1))
  Write (Unit=iout, Fmt='(A)', Err=1) 'OPTIONS yrev sequential '//cendian//'_endian'
  Write (Unit=cout1, Fmt=*) rlng(1,1)
  Write (Unit=cout2, Fmt=*) rlng(2,1) - rlng(1,1)
  Write (Unit=iout, Fmt='(A,I0,4A)', Err=1) &
        'XDEF    ', cfield(1)%region%nlgs, ' LINEAR ', Trim(AdjustL(cout1)), ' ', Trim(AdjustL(cout2))
  Write (Unit=cout1, Fmt=*) rlat(cfield(1)%region%nlts,1)
  Write (Unit=cout2, Fmt=*) rlat(1,1) - rlat(2,1)
  Write (Unit=iout, Fmt='(A,I0,4A)', Err=1) &
        'YDEF    ', cfield(1)%region%nlts, ' LINEAR ', Trim(AdjustL(cout1)), ' ', Trim(AdjustL(cout2))
  cdate = get_cdate(cfield(1)%tprd%sdate, 3)
  Select Case (iseq)
   Case (isq_yr) ! - yearly -
     Write (Unit=iout, Fmt=     '(A,I0,3A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), ' 1yr'
   Case (isq_sn) ! - seasonal -
     Write (Unit=iout, Fmt=     '(A,I0,3A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), ' 3mn'
   Case (isq_mn, isq_so) ! - monthly -
     Write (Unit=iout, Fmt=     '(A,I0,3A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), ' 1mn'
   Case (1:) ! - daily and sub-seasonal -
     Write (Unit=iout, Fmt='(A,I0,2A,I0,A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate)//' ', iseq, 'dy'
  End Select
  Write (Unit=iout, Fmt='(A)', Err=1) 'ZDEF    1 LINEAR 1 1'
!
! Output variables list
  Call write_grads_vars_ctrl (iout, cfield, nfs, nls, nvs, ifail, &
       cvar=cvar)
  GoTo 2
!
! Error
1 ifail = 3
!
2 Close (Unit=iout)
!
  Return
 End Subroutine write_grads_grid_ctrl
!
!
!
 Subroutine write_grads_stns_ctrl (iout, cfile, cfield, nfs, nls, nt, nvs, ifail, &
            cvar)
!
! Outputs a GrADS control file for station data
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_io_constants,   Only: cendian, grads_ext
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr
  Use time,                Only: iseq, &
                                 get_cdate
  Use iofiles,             Only: ofile
  Use fields,              Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: nt   ! - number of cases -
  Integer, Intent(In) :: nvs  ! - number of variables -
!
  Type(ofile), Intent(In) :: cfile ! - data file -
!
! Input arrays
  Type(field), Intent(In) :: cfield(:) ! - field specifications -
!
! - Optional input arrays -
  Character(Len=*), Dimension(:), Intent(In), Optional :: cvar ! - variable names -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Character(Len=18) :: cdate ! - format statement -
  Character(Len=32) :: cout1 ! - output field -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustL, Trim
!
! Executable Statements
!
! Open control file
  Open (Unit=iout, File=Trim(cfile%cfile%fdir)//Trim(cfile%cfile%fname)//grads_ext, &
        Action='write', Form='formatted', Status='unknown')
  If (ifail /= 0) Then
     ifail = 1
     Return
  End If
!
! Write control file
  Write (Unit=iout, Fmt='(A)', Err=1) 'DSET    '//Trim(cfile%cfile%ffile)
  Write (Unit=iout, Fmt='(A)', Err=1) 'DType   station'
  Write (Unit=iout, Fmt='(A)', Err=1) 'STNMAP  '//Trim(cfile%cfile%fdir)//Trim(cfile%cfile%fname)//'.map'
  Write (Unit=iout, Fmt='(A)', Err=1) 'TITLE   '//Trim(cfile%cfile%desc)
  Write (Unit=cout1, Fmt=*) cfield(1)%rmiss
  Write (Unit=iout, Fmt='(2A)', Err=1) 'UNDEF   ', Trim(AdjustL(cout1))
  Write (Unit=iout, Fmt='(A)', Err=1) 'OPTIONS '//cendian//'_endian'
  cdate = get_cdate(cfield(1)%tprd%sdate, 3)
  Select Case (iseq)
   Case (isq_yr) ! - yearly -
     Write (Unit=iout, Fmt=     '(A,I0,3A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), ' 1yr'
   Case (isq_sn) ! - seasonal -
     Write (Unit=iout, Fmt=     '(A,I0,3A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), ' 3mn'
   Case (isq_mn, isq_so) ! - monthly -
     Write (Unit=iout, Fmt=     '(A,I0,3A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), ' 1mn'
   Case (1:) ! - daily and sub-seasonal -
     Write (Unit=iout, Fmt='(A,I0,2A,I0,A)', Err=1) 'TDEF    ', nt, ' LINEAR ', Trim(cdate), iseq, 'dy'
  End Select
!
! Output variables list
  Call write_grads_vars_ctrl (iout, cfield, nfs, nls, nvs, ifail, &
       cvar=cvar)
  GoTo 2
!
! Error
1 ifail = 3
!
2 Close (Unit=iout)
!
  Return
 End Subroutine write_grads_stns_ctrl
!
!
!
 Subroutine write_grads_vars_ctrl (iout, cfield, nfs, nls, nvs, ifail, &
            cvar)
!
! Outputs variables list in a GrADS control file
!
! On exit Returns:
!    0 If no errors
!    3 If problem writing to file
!
! Modules
  Use data_text, Only: iachara_l, nlet
  Use fields,    Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: nfs  ! - number of fields -
  Integer, Intent(In) :: nls  ! - number of lagged fields -
  Integer, Intent(In) :: nvs  ! - number of variables -
!
! Input arrays
  Type(field), Intent(In) :: cfield(:) ! - field specifications -
!
! - Optional input arrays -
  Character(Len=*), Dimension(:), Intent(In), Optional :: cvar ! - variable names -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: ia, ib ! - ASCII collating sequence indicators -
  Integer :: ifd    ! - field index -
  Integer :: ilf    ! - lagged-field index -
  Integer :: iv     ! - variable index -
  Integer :: l      ! - field index -
  Integer :: nfl    ! - total number of fields and lagged-fields -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Achar, Mod, Present, Trim
!
! Executable Statements
!
! Output variables list
  nfl = nfs*nls*nvs
  Write (Unit=iout, Fmt='(A,I0)', Err=1) 'VARS ', nfl
  If (nfl <= nlet) Then
     ia = iachara_l
     Do iv = 1, nvs
        Do ifd = 1, nfs
           Do ilf = 1, nls
              l = (ifd - 1)*nls + ilf
              If (Present(cvar)) Then
                 If (Trim(cfield(l)%unitc) == 'N/A') Then
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia), 0, 99, Trim(cfield(l)%var)//'_'//Trim(cfield(l)%cssn)//'_'//Trim(cvar(iv))
                 Else
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia), 0, 99, &
                          Trim(cfield(l)%var)//'_'//Trim(cfield(l)%cssn)//'_'//Trim(cvar(iv))//' '//Trim(cfield(l)%unitc)
                 End If
              Else
                 If (Trim(cfield(l)%unitc) == 'N/A') Then
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia), 0, 99, Trim(cfield(l)%var)//'_'//cfield(l)%cssn
                 Else
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia), 0, 99, Trim(cfield(l)%var)//'_'//cfield(l)%cssn//' '//Trim(cfield(l)%unitc)
                 End If
              End If
              ia = ia + 1
           End Do
        End Do
     End Do
  Else
     ia = iachara_l
     ib = iachara_l
     Do iv = 1, nvs
        Do ifd = 1, nfs
           Do ilf = 1, nls
              l = (ifd - 1)*nls + ilf
              If (Present(cvar)) Then
                 If (Trim(cfield(l)%unitc) == 'N/A') Then
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia)//Achar(ib), 0, 99, Trim(cfield(l)%var)//'_'//cfield(l)%cssn//'_'//Trim(cvar(iv))
                 Else
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia)//Achar(ib), 0, 99, &
                          Trim(cfield(l)%var)//'_'//cfield(l)%cssn//'_'//Trim(cvar(iv))//' '//Trim(cfield(l)%unitc)
                 End If
              Else
                 If (Trim(cfield(l)%unitc) == 'N/A') Then
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia)//Achar(ib), 0, 99, Trim(cfield(l)%var)//'_'//cfield(l)%cssn
                 Else
                    Write (Unit=iout, Fmt='(2X,A,I2,I3,1X,A)', Err=1) &
                          Achar(ia)//Achar(ib), 0, 99, Trim(cfield(l)%var)//'_'//cfield(l)%cssn//' '//Trim(cfield(l)%unitc)
                 End If
              End If
              ib = ib + 1
              If (Mod(ib-iachara_l, nlet) == 0) Then
                 ib = iachara_l
                 ia = ia + 1
              End If
           End Do
        End Do
     End Do
  End If
  Write (Unit=iout, Fmt='(A)', Err=1) 'ENDVARS'
!
  ifail = 0
!
  Return
!
! Error
1 ifail = 3
!
  Return
 End Subroutine write_grads_vars_ctrl
!
!
!
 Subroutine open_output (iout, afile, nfs, idstr, lrec, ifail, &
            ng)
!
! Opens CPT output file and prints XMLNS header
!
! On exit Returns:
!    0 If no errors
!    3 If problem creating output file
!
! Modules
  Use data_io_constants, Only: cxmlns_cpt, faccs, ffmts
#ifdef NAGFOR
  Use data_io_constants, Only: lfli
#endif
  Use iofiles,           Only: ofile, &
                               itags, &
                               set_fmt
  Use analysis,          Only: ianal
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout  ! - output unit number -
  Integer, Intent(In) :: nfs   ! - number of fields -
  Integer, Intent(In) :: idstr ! - data structure identifier -
  Integer, Intent(In) :: lrec  ! - record length -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: ng ! - number of categories -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Character(Len=11) :: cformat ! - format -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Set format
  Call set_fmt (lrec, afile%cfile%ffmt)
  If (Trim(ffmts(afile%cfile%ffmt%ifmt))/='GrADS') Then
     cformat = ffmts(afile%cfile%ffmt%ifmt)
  Else
     cformat = 'unformatted'
     If (idstr == 2) afile%cfile%ffmt%iacc=3
  End If
!
! Open output file
  Select Case (Trim(faccs(afile%cfile%ffmt%iacc)))
   Case ('sequential')
     If (ianal /= 0) Print *, 'Saving: '//Trim(afile%cfile%ffile)
     Open (Unit=iout, File=Trim(afile%cfile%ffile), Access='sequential', Action='Write', Form=Trim(cformat), &
#ifdef NAGFOR
           IOstat=ifail, Recl=lfli, Status='unknown')
#else
           IOstat=ifail, Status='unknown')
#endif
   Case ('direct')
     If (afile%cfile%ffmt%lrec == 0) afile%cfile%ffmt%lrec = lrec*afile%cfile%ffmt%iprc*4
     Write (Unit=*, Fmt='(3A,I0,A)') 'Saving: ', Trim(afile%cfile%ffile),' (record length ',afile%cfile%ffmt%lrec,')'
     Open (Unit=iout, File=Trim(afile%cfile%ffile), Access='direct', Action='Write', Form=Trim(cformat), &
           IOstat=ifail, Recl=afile%cfile%ffmt%lrec, Status='unknown')
   Case ('transparent', 'stream')
     If (ianal /= 0) Print *, 'Saving: '//Trim(afile%cfile%ffile)
     Open (Unit=iout, File=Trim(afile%cfile%ffile), Access=Trim(faccs(afile%cfile%ffmt%iacc)), Action='Write', Form=Trim(cformat), &
           IOstat=ifail, Status='unknown')
  End Select
!
! Error
  If (ifail /= 0) Then
     ifail = 1
     Return
  End If
!
! Omit CPT tags
  If (itags == 0) Return
!
! Print XMLNS header
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     Write (Unit=iout, Fmt='(A)',IOstat=ifail) 'xmlns:cpt='//cxmlns_cpt
     If (ifail /= 0) GoTo 1
!
! Print number of fields
     If (nfs > 1) Then
        Call write_tag (iout, ifail, &
                        cpt_nfields=nfs)
     End If
     If (ifail /= 0) GoTo 1
     If (Present(ng)) Call write_tag (iout, ifail, &
                           cpt_ncats=ng)
  End Select
!
! Error
1 If (ifail /= 0) ifail = 3
!
  Return
 End Subroutine open_output
!
!
!
 Subroutine close_output (iout, afile, ifail)
!
! Closes CPT output file and prints XMLNS header
!
! On exit Returns:
!    0 If no errors
!   -1 If problem updating progress meter
!
! Modules
  Use progress_meter, Only: update_progress_meter
  Use errors,         Only: cpt_error
  Use iofiles,        Only: ofile, &
                            file_reset
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout  ! - output unit number -
!
! Input/output scalars
  Integer, Intent(InOut) :: ifail ! - error indicator -
!
  Type(ofile), Intent(InOut) :: afile ! - output file -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Close file
  Close (Unit=iout)
!
! Report errors 
  If (ifail /= 0) Call cpt_error ('write_results', ifail, .true., &
                       c_arg1=Trim(afile%cfile%desc), c_arg2=Trim(afile%cfile%ffile))
  Call file_reset (afile%cfile%lset, afile%cfile%fname)
!
! Update progress meter
  If (update_progress_meter(.false.) /= 0) ifail = -1
!
  Return
 End Subroutine close_output
!
!
!
 Subroutine write_tag (iout, ifail, &
                       cpt_nfields, cpt_ncats, cpt_name, cpt_field, cpt_model, cpt_c, cpt_prob, cpt_cMode, cpt_mode,       &
                       cpt_missing, cpt_t, cpt_s, cpt_z, cpt_m, cpt_clev, cpt_limit, cpt_nrow, cpt_ncol, cpt_row, cpt_col, &
                       cpt_units, cpt_tagi, cpt_itag, cpt_tagr, cpt_rtag, cpt_tagc, cpt_ctag)
!
! Prints a CPT tag line
!
! On exit Returns:
!    0 If no errors
!    1 If problem writing to file
!
! Modules
  Use data_numbers,      Only: one
  Use data_io_constants, Only: lvar
  Use maths,             Only: iprec
  Use time,              Only: pdate, pprd, &
                               Operator(==), &
                               get_cdate
  Use space,             Only: level
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
!
! - Optional input scalars -
  Integer, Intent(In), Optional :: cpt_nfields ! - number of fields -
  Integer, Intent(In), Optional :: cpt_ncats   ! - number of categories -
  Integer, Intent(In), Optional :: cpt_c       ! - current category -
  Integer, Intent(In), Optional :: cpt_mode    ! - current Mode -
  Integer, Intent(In), Optional :: cpt_m       ! - ensemble member -
  Integer, Intent(In), Optional :: cpt_nrow    ! - number of rows -
  Integer, Intent(In), Optional :: cpt_ncol    ! - number of columns -
  Integer, Intent(In), Optional :: cpt_itag    ! - extra integer tag value -
!
  Real(Kind=rp), Intent(In), Optional :: cpt_prob    ! - climatological probability -
  Real(Kind=rp), Intent(In), Optional :: cpt_clev    ! - confidence level -
  Real(Kind=rp), Intent(In), Optional :: cpt_missing ! - missing values flag -
  Real(Kind=rp), Intent(In), Optional :: cpt_rtag    ! - extra Real tag value -
!
  Character(Len=*), Intent(In), Optional :: cpt_field ! - field -
  Character(Len=*), Intent(In), Optional :: cpt_model ! - Model -
  Character(Len=*), Intent(In), Optional :: cpt_name  ! - name -
  Character(Len=*), Intent(In), Optional :: cpt_cMode ! - current Mode -
  Character(Len=*), Intent(In), Optional :: cpt_row   ! - rows -
  Character(Len=*), Intent(In), Optional :: cpt_limit ! - confidence limit -
  Character(Len=*), Intent(In), Optional :: cpt_col   ! - columns -
  Character(Len=*), Intent(In), Optional :: cpt_units ! - units -
  Character(Len=*), Intent(In), Optional :: cpt_tagi  ! - extra integer tag -
  Character(Len=*), Intent(In), Optional :: cpt_tagr  ! - extra Real tag -
  Character(Len=*), Intent(In), Optional :: cpt_tagc  ! - extra character tag -
  Character(Len=*), Intent(In), Optional :: cpt_ctag  ! - extra character tag value -
!
!
  Type(level), Intent(In), Optional :: cpt_z ! - level -
!
  Type(pprd), Intent(In), Optional :: cpt_t ! - date -
!
  Type(pdate), Intent(In), Optional :: cpt_s ! - start date -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Character(Len=    15) :: cfmt ! - format statement -
  Character(Len=lvar+1) :: cout ! - output field -
!
  Logical :: lfirst ! - first output field? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustL, Len_Trim, Max, Present, Trim
!
! Executable Statements
!
! Print nfields
  ifail = 0
  If (Present(cpt_nfields)) Then
     Write (Unit=iout, Fmt='(A,I0)', Err=1) 'cpt:nfields=', cpt_nfields
     Return
  End If
!
! Print ncats
  If (Present(cpt_ncats)) Then
     Write (Unit=iout, Fmt='(A,I0)', Err=1) 'cpt:ncats=', cpt_ncats
     Return
  End If
!
! Print name
  If (Present(cpt_name)) Then
     Write (Unit=iout, Fmt='(A)', Err=1) 'cpt:Name='//Trim(cpt_name)
     Return
  End If
!
! Print tags line
  lfirst = .true.
  If (Present(cpt_field)) Then
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:field='//Trim(cpt_field)
     lfirst = .false.
  End If
  If (Present(cpt_tagi) .and. Present(cpt_itag)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A,I0)', Err=1, Advance='no') 'cpt:'//Trim(cpt_tagi)//'=', cpt_itag
     lfirst = .false.
  End If
  If (Present(cpt_tagr) .and. Present(cpt_rtag)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=cout, Fmt=*) cpt_rtag
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:'//Trim(cpt_tagr)//'='//Trim(AdjustL(cout))
     lfirst = .false.
  End If
  If (Present(cpt_tagc) .and. Present(cpt_ctag)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:'//Trim(cpt_tagc)//'='//cpt_ctag
     lfirst = .false.
  End If
  If (Present(cpt_model)) Then
     If (Len_Trim(cpt_model) > 0 .and. Trim(cpt_model) /= 'N/A') Then
        If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
        Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:Model='//Trim(cpt_model)
        lfirst = .false.
     End If
  End If
  If (Present(cpt_c)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A,I0)', Err=1, Advance='no') 'cpt:C=', cpt_c
     lfirst = .false.
  End If
  If (Present(cpt_prob)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=cout, Fmt=*) cpt_prob
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:clim_prob='//Trim(AdjustL(cout))
     lfirst = .false.
  End If
  If (Present(cpt_cMode)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:Mode='//Trim(cpt_cMode)
     lfirst = .false.
  End If
  If (Present(cpt_mode)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A,I0)', Err=1, Advance='no') 'cpt:Mode=', cpt_mode
     lfirst = .false.
  End If
  If (Present(cpt_m)) Then
     If (cpt_m > 0) Then
        If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
        Write (Unit=iout, Fmt='(A,I0)', Err=1, Advance='no') 'cpt:M=', cpt_m
        lfirst = .false.
     End If
  End If
  If (Present(cpt_limit)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(2A)', Err=1, Advance='no') 'cpt:climit=', Trim(cpt_limit)
     lfirst = .false.
  End If
  If (Present(cpt_clev)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:clev='
     Write (Unit=cfmt, Fmt='(A,2(I1,A))') '(F', iprec(cpt_clev,3)+3, '.', iprec(cpt_clev,3), ')'
     Write (Unit=iout, Fmt=cfmt, Err=1, Advance='no') cpt_clev
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') '%'
     lfirst = .false.
  End If
  If (Present(cpt_z)) Then
     If (Trim(cpt_z%unit)/='none') Then
        If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
        If (cpt_z%hght > one) Then
           Write (Unit=cfmt, Fmt='(A,I1,A)') '(F24.', Max(1, iprec(cpt_z%hght, 3)), ')'
        Else
           Write (Unit=cfmt, Fmt='(A,I1,A)') '(F24.', iprec(cpt_z%hght, 9), ')'
        End If
        Write (Unit=cout, Fmt=cfmt) cpt_z%hght
        Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:Z='//Trim(AdjustL(cout))//' '//Trim(cpt_z%unit)
        lfirst = .false.
     End If
  End If
  If (Present(cpt_t)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     cout = get_cdate(cpt_t, 1, 1)
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:T='//Trim(cout)
     lfirst = .false.
  End If
  If (Present(cpt_s)) Then
     If (.not.(cpt_s == 0)) Then
        If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
        cout = get_cdate(cpt_s, 1)
        Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:S='//Trim(cout)
        lfirst = .false.
     End If
  End If
  If (Present(cpt_nrow)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A,I0)', Err=1, Advance='no') 'cpt:nrow=', cpt_nrow
     lfirst = .false.
  End If
  If (Present(cpt_ncol)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A,I0)', Err=1, Advance='no') 'cpt:ncol=', cpt_ncol
     lfirst = .false.
  End If
  If (Present(cpt_row)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:row='//cpt_row
     lfirst = .false.
  End If
  If (Present(cpt_col)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:col='//cpt_col
     lfirst = .false.
  End If
  If (Present(cpt_units)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:units='//Trim(cpt_units)
     lfirst = .false.
  End If
  If (Present(cpt_missing)) Then
     If (.not.lfirst) Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') ', '
     Write (Unit=cout, Fmt=*) cpt_missing
     Write (Unit=iout, Fmt='(A)', Err=1, Advance='no') 'cpt:missing='//AdjustL(cout)
     lfirst = .false.
  End If
  Write (Unit=iout, Fmt='()', Err=1, Advance='yes')
!
  Return
!
1 ifail = 1
  Return
 End Subroutine write_tag
!
!
!
 Function output_value (iout, &
          c, i, r, ladv, ldelim)
!
! Modules
  Use iofiles,  Only: cdelim, ipreco
!
! Prints a delimited output
!
! Function type
  Integer :: output_value
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output unit number -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: i ! - output value -
!
  Real(Kind=rp), Intent(In), Optional :: r ! - output value -
!
  Character(Len=*), Intent(In), Optional :: c ! - output value -
!
  Logical, Intent(In), Optional :: ladv   ! - advance? -
  Logical, Intent(In), Optional :: ldelim ! - include delimitor? -
!
! Locals
!
! Local scalars
  Character(Len=4) :: cadv ! - advance? -
  Character(Len=8) :: cfmt ! - output format -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustL, Present, Trim
!
! Executable Statements
!
! Identify whether to advance
  If (Present(ladv)) Then
     If (ladv) Then
        cadv = 'yes'
     Else
        cadv = 'no'
     End If
  Else
     cadv = 'yes'
  End If
!
! Identify whether to include delimitor
  If (Present(ldelim)) Then
     If (ldelim) Write (Unit=iout, Fmt='(A)', Advance='no', Err=1) cdelim
  End If
!
! Print delimited output
! - integers -
  If (Present(i)) Then
     Write (Unit=iout, Fmt='(I0)', Advance=Trim(cadv), Err=1) i
! - reals -
  Else If (Present(r)) Then
     If (ipreco >= 0) Then
        Write (Unit=cfmt, Fmt='(A,I0,A)') '(F0.', ipreco, ')'
        Write (Unit=iout, Fmt=cfmt, Advance=Trim(cadv), Err=1) r
     Else
        Write (cout, Fmt=*) r
        Write (Unit=iout, Fmt='(A)', Advance=Trim(cadv), Err=1) Trim(AdjustL(cout))
     End If
! - characters -
  Else If (Present(c)) Then
     Write (Unit=iout, Fmt='(A)', Advance=Trim(cadv), Err=1) Trim(AdjustL(c))
  End If
!
  output_value = 0
  Return
!
! Errors
1 output_value = 1
  Return
 End Function output_value
End Module cpt_output
