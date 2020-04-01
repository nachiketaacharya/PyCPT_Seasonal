! Author: Simon Mason
Module gui
!
! Modules
  Use data_numbers,       Only: rp
#if GUI == 1
  Use data_numbers,       Only: one
  Use data_gui_constants, Only: mcol
#endif
  Use data_gui_constants, Only: i_handle
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: add_check_button, box_close, box_open, boxes_open, ganged_rbs, get_handle, get_output_unit, gui_creturn, init_psheet, &
            init_win, open_clearwin, print_advisory, print_cvalue, print_error, print_item, print_ivalue, print_rvalue,           &
            print_subtitle, print_text, print_warning, prompt_boolean, prompt_bullet, prompt_cvalue, prompt_file, prompt_integer, &
            prompt_list, prompt_pdate, prompt_real, prompt_select, prompt_yesno, set_cursor_waiting, set_tabs, set_win_title,     &
            use_url, win_prompt, window_end
#if GUI == 1
  Public :: add_bitmap, add_button, add_graphic_area, add_mdi_frame, add_prog_bar, add_slider_bar, adjust_y, close_windows, &
            get_cols, gui_centre, ic_current, init_cols, init_cw, leave_window, reset_jpg, set_colour_direction, set_font,  &
            set_text_size
#endif
!
! Scalars
!
! Integer scalars
  Integer(i_handle), Private :: ihc = 1001 ! - window handle counter -
!
  Integer, Public :: ic_par = 0 ! - parent window control variable -
  Integer, Public :: iw         ! - window argument -
  Integer, Public :: mcw        ! - maximum number of results windows -
#if GUI == 1
  Integer, Public :: ic          ! - current control handle -
  Integer, Public :: ic_con = 0  ! - console window control variable -
  Integer, Public :: icbg        ! - background colour -
  Integer, Public :: icbw = 0    ! - colour / black and white flag (0 = colour; 1 = black and white) -
  Integer, Public :: icoldom     ! - domain colour -
  Integer, Public :: icolh       ! - highlighted text colour -
  Integer, Public :: icoli = 0   ! - colour scale identifier -
  Integer, Public :: icon_h      ! - console window width adjustment -
  Integer, Public :: icon_v      ! - console window height adjustment -
  Integer, Public :: icrev = 0   ! - colour reversal flag -
  Integer, Public :: nclr = mcol ! - number of shading colours -
#else
  Integer, Public :: iprompt = 0 ! - action on invalid response to prompt: 0 = reprompt; 1 = ignore; 2 = stop -
#endif
!
  Integer(Kind=i_handle), Public :: ih_nf  ! - handle for nf prompt -
  Integer(Kind=i_handle), Public :: ih_nt  ! - handle for nt prompt -
#if GUI == 1
  Integer(Kind=i_handle), Public :: ih_con ! - console window handle -
!
! Real scalars
  Real(Kind=rp), Public :: jq ! - JPEG quality -
!
#endif
! Character scalars
  Character(Len=128), Public :: cwtitle ! - CPT window title -
!
! Arrays
!
#if GUI == 1
! Integer arrays
  Integer, Dimension(0:mcol+20), Public :: icol ! - colours -
!
  Integer, Dimension(:), Allocatable, Private :: ich ! - child window handles -
#endif
! Interfaces
!
! Generic interfaces
#if GUI == 1
  Interface add_button
   Module Procedure add_button_lab
   Module Procedure add_button_txt
  End Interface add_button
!
#else
  Interface ganged_rbs
   Module Procedure ganged_rbs_a
   Module Procedure ganged_rbs_s
  End Interface ganged_rbs
!
#endif
Contains
!
!
#if GUI == 1
 Subroutine reset_jpg ()
!
! Resets JPEG graphics quality to default setting
!
! Modules
  Use clrwin$,            Only: set_jpeg_quality$
  Use data_gui_constants, Only: djqual
!
! Executable Statements
!
! Initialise JPEG quality
  jq = set_jpeg_quality$(djqual)
!
  Return
 End Subroutine reset_jpg
!
!
!
 Subroutine init_cw (ifail)
!
! Initialises child windows
!
! Arguments
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Set maximum number of child windows
  Allocate (ich(mcw), Stat=ifail)
  If (ifail /= 0) Then
     ifail = 1
     Return
  End If
  ich(:) = 0
  ifail = 0
!
  Return
 End Subroutine init_cw
!
!
!
#endif
 Function get_output_unit() &
          Result (unit)
!
! Assigns an output unit number
!
#if GUI == 0
! Modules
  Use data_io_constants, Only: ioutstd
!
#endif
! Function type
  Integer :: unit
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Int
!
! Executable Statements
!
! Assign an output unit number
  unit = Int(get_handle())
#else
  unit = ioutstd
#endif
!
  Return
 End Function get_output_unit
!
!
!
 Function get_handle()
!
! Assigns a handle
!
! Function type
  Integer(Kind=i_handle) :: get_handle
!
! Executable Statements
!
! Assign a handle
  get_handle = ihc
  ihc = ihc + 1
!
  Return
 End Function get_handle
!
!
!
#if GUI == 1
 Subroutine init_cols ()
!
! Initialises colour indices
!
! Modules
  Use clrwin$,            Only: rgb$
  Use data_gui_constants, Only: col_bg, col_black, col_domain, col_htext, col_white
!
! Executable Statements
!
! Background colour
  icbg = rgb$(col_bg(1), col_bg(2), col_bg(3))
!
! Define basic colours
! - white -
  icol(0)=rgb$(col_white(1), col_white(2), col_white(3))
! - black -
  icol(1) = rgb$(col_black(1), col_black(2), col_black(3))
!
! Highlighted text colour
  icolh = rgb$(col_htext(1), col_htext(2), col_htext(3))
!
! Domain colour
  icoldom = rgb$(col_domain(1), col_domain(2), col_domain(3))
!
  Return
 End Subroutine init_cols
!
!
!
 Function get_cols(n)
!
! Initialises colour indices
!
! Modules
  Use clrwin$,            Only: rgb$
  Use data_gui_constants, Only: icol_gr, icol_mc, icol_p2, icol_p3, icol_sk
  Use data_io_constants,  Only: iin, ddir
!
! Function type
  Integer :: get_cols
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of colours -
!
! Locals
!
! Local scalars
  Integer :: i     ! - colour index -
  Integer :: irgb  ! - RGB component -
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Open colours file
  Open (Unit=iin, File=Trim(ddir)//'colours.txt', IOstat=ifail, Action='read', Form='formatted')
  If (ifail /= 0) GoTo 1
!
! Adjust number of colours for availability
  Select Case (icoli)
   Case (0) ! - simple graph colours -
     nclr = 3
   Case (icol_gr) ! - graph colours -
     nclr = 9
   Case (icol_p2) ! - forecast probabilities shading (blue-red) -
     nclr = 14
   Case (icol_p3) ! - forecast probabilities shading (blue-red with green) -
     nclr = 22
   Case (icol_sk, icol_mc) ! - skill map shading (blue-red) and scores shading (white-red) -
     nclr = n
     Select Case (n)
      Case (18:)
        nclr = 20
      Case (14:17)
        nclr = 18
      Case (12:13)
        nclr = 14
      Case Default
        nclr = 12
     End Select
  End Select
!
! Black and white colour scheme
  Select Case (icbw)
   Case (1)
     Call read_cols (icbw, icrev, 0, 0, icol(1:), ifail)
     If (ifail /= 0) GoTo 1
! - graph background colours -
     icol(11:16) = icol(0)
! - shading colours -
     Select Case (icrev)
      Case (0)
        Do i = 21, 20+nclr
           irgb = ((i - 20)*255)/(nclr + 1)
           icol(i) = rgb$(irgb, irgb, irgb)
        End Do
      Case (1)
        Do i = 21, 20+nclr
           irgb = ((nclr - (i - 21))*255)/(nclr + 1)
           icol(i) = rgb$(irgb, irgb, irgb)
        End Do
     End Select
!
! Colour schemes
   Case (0)
     Call read_cols (icbw, 0, 0, 0, icol(1:), ifail) ! - get primary colours -
     If (ifail /= 0) GoTo 1
! - skill shading colours -
     Select Case (icoli)
      Case Default
        Call read_cols (icbw, icrev, icoli, nclr, icol(1:), ifail) ! - get shading colours -
        If (ifail /= 0) GoTo 1
      Case (icol_mc)
! - monochrome colours -
        Select Case (icrev)
         Case (0)
           Do i = 21, 20+nclr
              irgb = ((nclr - (i - 21))*255)/(nclr + 1)
              icol(i) = rgb$(255, irgb, irgb)
           End Do
         Case (1)
           Do i = 21, 20+nclr
              irgb = ((nclr - (i - 21))*255)/(nclr + 1)
              icol(i) = rgb$(irgb, irgb, 255)
           End Do
        End Select
     End Select
     If (nclr > mcol) icol(21+nclr:) = 0
  End Select
!
! Close file
  Close (Unit=iin)
  get_cols = 0
  Return
!
! Errors
1 Close (Unit = iin)
  get_cols = 1
  Return
!
 Contains
!
!
  Subroutine read_cols (icbw, icrev, icoli, nclr, icol, ifail)
!
! Gets colour definitions
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icbw  ! - colour / black and white flag -
  Integer, Intent(In) :: icrev ! - colour reversal flag -
  Integer, Intent(In) :: icoli ! - colour scale identifier -
  Integer, Intent(In) :: nclr  ! - number of colour definitions -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input/output arrays
  Integer, Dimension(:), Intent(InOut) :: icol ! - colour definitions -
!
! Locals
!
! Local scalars
  Integer :: i      ! - colour index -
  Integer :: ibw    ! - black and white index -
  Integer :: icr    ! - colour reversal index -
  Integer :: ics    ! - colour scale identifier -
  Integer :: ic     ! - colour index -
  Integer :: ired   ! - red component -
  Integer :: igreen ! - green component -
  Integer :: iblue  ! - blue component -
  Integer :: n      ! - number of colours -
!
  CHARACTER :: c ! - tag -
!
! Executable Statements
!
! Get colour definitions
  Do
     Read (Unit=iin, Fmt=*, IOstat=ifail, End=1) c, ibw, icr, ics, n
     If (ifail /= 0) Cycle
     If (c == 'C') Then
        If (ibw == icbw .and. icr == icrev) Then
           If (icoli > 0) Then
              If (icoli == ics) Then
                 If (nclr > 0) Then
                    If (n == nclr) Exit
                 Else
                    Exit
                 End If
              End If
           Else
              Exit
           End If
        End If
        Do i = 1, n
           Read (Unit=iin, Fmt=*, Err=1)
        End Do
     End If
  End Do
  Do i = 1, n
     Read (Unit=iin, Fmt=*, Err=1) ic, ired, igreen, iblue
     icol(ic) = rgb$(ired, igreen, iblue)
  End Do
  ifail = 0
  Return
!
1 ifail = 1
  Return
  End Subroutine read_cols
 End Function get_cols
!
!
!
 Subroutine set_colour_direction (cvar)
!
! Sets initial colour direction based on the Y variable
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cvar ! - variable -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Set initial colour direction
  Select Case (Trim(cvar))
   Case ('enso', 'nino', 'sst', 'ssta', 'temp', 'tmax', 'tmin', 't2m', &
         'ENSO', 'NINO', 'SST', 'SSTA', 'TEMP', 'TMax', 'TMin') ! - temperatures -
     icrev = 1
   Case ('prec', 'precipitation', 'prate', 'prcp', 'rain', 'rainfall', &
         'PREC', 'PRECIPITATION', 'PRATE', 'PRCP', 'RAIN', 'RAINFAll') ! - rainfall -
     icrev = 0
   Case Default
     Return
  End Select
!
  Return
 End Subroutine set_colour_direction
!
!
!
#endif
 Subroutine init_psheet (ihp, title)
!
! Initialises a Clearwin+ property sheet
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ihp ! - property sheet handle -
!
  Character(Len=*), Intent(In) :: title ! - window title -
!
! Executable Statements
!
#if GUI == 1
! Create property sheet
! - set property sheet handle -
  iw = winio$('%sh&', ihp)
! - set window title -
  iw = winio$('%ca@&', title)
! - set default font -
  Call set_font ()
#else
! Indicate property
  Call print_text (' ')
  Call print_text (title)
#endif
!
  Return
 End Subroutine init_psheet
!
!
!
 Subroutine init_win &
            (title, idisable, ihandle, ipar, lch, lresize, c_opts, fclose, fclose2, fclose3)
!
#if GUI == 1
! Initialises ClearWin+ window
!
! Modules
  Use clrwin$, Only: winio$
  Use screen,  Only: window_update
#else
! Prints a centred title
!
#endif
! Arguments
!
! Input scalars
! - optional input scalars -
  Integer, Intent(In), Optional :: ipar ! - parent window handle -
!
  Integer(i_handle), Intent(In), Optional :: idisable ! - disable window handle -
!
  Character(Len=*), Intent(In), Optional :: title  ! - window title -
  Character(Len=*), Intent(In), Optional :: c_opts ! - additional options -
!
  Logical, Intent(In), Optional :: lch     ! - child window? -
  Logical, Intent(In), Optional :: lresize ! - resizeable window? -
!
! Output scalars
! - optional output scalars -
  Integer(i_handle), Intent(Out), Optional :: ihandle ! - window handle -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: fclose  ! - callback function on windows closure -
  Integer, External, Optional :: fclose2 ! - second callback function on windows closure -
  Integer, External, Optional :: fclose3 ! - third callback function on windows closure -
!
#if GUI == 1
! Locals
!
! Local scalars
  Integer :: i  ! - window index -
  Integer :: ic ! - current window index -
!
#endif
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
#if GUI == 1
  Intrinsic Trim
!
! Executable Statements
!
! Check number of open child windows
  If (Present(lch)) Then
     If (lch) Then
        ic = 0
        Do i = 1, mcw
           If (ich(i) >= 0) Then
              ic = i
              Exit
           End If
        End Do
        If (ic == 0) Then
           ich(mcw) = 1
           Call window_update (ich(mcw))
           ic = mcw
        End If
     End If
  End If
!
! Create window
! - set background colour -
  iw = winio$('%bg&', icbg)
! - set window title -
  If (Present(title)) Then
     iw = winio$('%ca@&', title)
     If (Present(lresize)) Then ! - make window resizeable -
        If (lresize) Then
           If (Present(c_opts)) Then
              iw = winio$('%ww['//Trim(c_opts)//']&')
           Else
              iw = winio$('%ww&')
           End If
        End If
     End If
  Else
     If (Present(c_opts)) Then
        iw = winio$('%ww['//Trim(c_opts)//']&')
     Else
        iw = winio$('%ww&')
     End If
  End If
! - set default font -
  Call set_font ()
! - set text size -
  Call set_text_size (one)
! - window handle -
  If (Present(ihandle)) iw = winio$('%hw&', ihandle)
! - child window -
  If (Present(ipar)) Then
     iw = winio$('%aw&', ipar) ! - create as child window -
     If (Present(lch)) Then
        If (lch) iw = winio$('%lw&', ich(ic)) ! - attach control variable -
     End If
  End If
! - disable host window -
  If (Present(idisable)) iw = winio$('%de&', idisable)
! - add callback function for window closure -
  If (Present(fclose)) iw = winio$('%cc&', fclose)
  If (Present(fclose2)) iw = winio$('%cc&', fclose2)
  If (Present(fclose3)) iw = winio$('%cc&', fclose3)
!
#else
! Executable Statements
!
! Print centred title
  Call print_text (' ')
  Call print_text (' ')
  Call print_text (title, &
       lcentre=.true.)
  Call print_text (' ')
  If (Present(ihandle)) ihandle = 0
#endif
!
  Return
 End Subroutine init_win
!
!
!
#if GUI == 1
 Subroutine leave_window (ic, &
            lowned)
!
! Leaves a window open
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Logical, Intent(In), Optional :: lowned ! - owned? -
!
! Output scalars
  Integer, Intent(Out) :: ic ! - window control -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! End window
  If (Present(lowned)) Then
     If (lowned) Then
        iw = winio$('%lw[owned]', ic)
     Else
        iw = winio$('%lw', ic)
     End If
  Else
     iw = winio$('%lw', ic)
  End If
!
  Return
 End Subroutine leave_window
!
!
!
#endif
 Subroutine window_end ()
!
! Complete a window
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
!
! Executable Statements
!
! End window
  iw = winio$('')
#endif
!
  Return
 End Subroutine window_end
!
!
!
#if GUI == 1
 Subroutine set_font ( &
            cfont, ts)
!
! Sets default font
!
! Modules
  Use clrwin$, Only: winio$
  Use screen,  Only: cdfont
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: cfont ! - font -
!
  Real(Kind=rp), Intent(In), Optional :: ts ! - text size -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Initialise font
  If (Present(cfont)) Then
     iw = winio$('%fn@&', Trim(cfont))
  Else
     iw = winio$('%fn@&', Trim(cdfont))
  End If
  If (Present(ts)) Then
     Call set_text_size (ts)
  Else
     Call set_text_size (one)
  End If
!
  Return
 End Subroutine set_font
!
!
!
 Subroutine set_text_size (ts)
!
! Sets text size
!
! Modules
  Use clrwin$, Only: winio$
  Use screen,  Only: dts
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In) :: ts ! - text-size -
!
! Executable Statements
!
! Set text size
  iw = winio$('%ts&', ts*dts)
!
  Return
 End Subroutine set_text_size
!
!
!
#endif
 Subroutine gui_centre ()
!
! Centres the next line of input
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
! Executable Statements
  iw = winio$('%cn&')
!
#endif
  Return
 End Subroutine gui_centre
!
!
!
 Subroutine gui_creturn (lhard)
!
! Advances GUI to a new line
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lhard ! - hard return? (.true. forces a newline; .false. only advances to the beyond the last output) -
!
! Executable Statements
!
! Advance to a new line
  If (lhard) Then
#if GUI == 0
     Call print_text (' ')
#else
     iw = winio$('%nl&')
  Else
     iw = winio$('%ff&')
#endif
  End If
!
  Return
 End Subroutine gui_creturn
!
!
!
 Subroutine set_win_title (ianal, &
            cproj)
!
! Sets CPT window title
!
! Modules
  Use data_version, Only: cpt_title
  Use labels,       Only: ca_atypes_a, &
                          ls_analyses
  Use screen,       Only: window_update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis type -
!
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: cproj ! - project name -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Set CPT window title
  If (Present(cproj) .and. ianal > 0) Then
     cwtitle = cpt_title//' - '//Trim(cproj)//' ('//ca_atypes_a(ianal)//')'
  Else
     Select Case (ianal)
      Case (0)
        cwtitle = cpt_title
      Case Default
        cwtitle = cpt_title//' - '//Trim(ls_analyses(ianal)%c)
     End Select
  End If
  Call window_update (cwtitle)
!
  Return
 End Subroutine set_win_title
!
!
!
 Function win_prompt(lcancel, &
          fok, cextra, iset1, ival1, iset2, ival2, rset1, rval1, fextra)
!
! Adds prompts to close window
!
#if GUI == 1
! Modules
  Use labels, Only: l_cancel
#endif
!
! Function type
  Integer :: win_prompt
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lcancel ! - add cancel button? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: ival1 ! - value to set to iset1 -
  Integer, Intent(In), Optional :: ival2 ! - additional value to set to iset2 -
!
  Real(Kind=rp), Intent(In), Optional :: rval1 ! - value to set to rset -
!
  Character(Len=*), Intent(In), Optional :: cextra ! - additional button text -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: iset1   ! - variable to set -
  Integer, Intent(InOut), Optional :: iset2   ! - additional variable to set -
!
  Real(Kind=rp), Intent(InOut), Optional :: rset1 ! - variable to set -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: fok    ! - OK button Callback function -
  Integer, External, Optional :: fextra ! - additional button Callback function -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Check number of open child windows
#if GUI == 1
  Call gui_creturn (.false.)
  Call gui_creturn (.true.)
  If (Present(cextra)) Then
     Call add_button ('tt', 1, '&OK', &
          ldef=.true., cbf=fok)
     Call add_button ('tt', 2, cextra, &
          iset1=iset1, ival1=ival1, iset2=iset2, ival2=ival2, rset1=rset1, rval1=rval1, cbf=fextra)
     If (lcancel) Call add_button ('tt', 3, l_cancel)
  Else
     Call add_button ('bt', 1, '&OK', &
          ldef=.true., iwdth=5, cbf=fok)
     If (lcancel) Call add_button ('bt', 3, l_cancel, &
                       iwdth=5)
  End If
  Call gui_creturn (.false.)
  Call window_end ()
#else
  iw = 1
  If (Present(fok)) Then
     win_prompt = fok()
  End If
#endif
  win_prompt = iw
!
  Return
 End Function win_prompt
!
!
!
 Subroutine open_clearwin (iout, ncol, nrow, lscale, &
            cbox, ih, i_out, maxl, lhs, lvs)
!
! Opens a Clearwin+ window
!
#if GUI == 1
! Modules
  Use clrwin$, Only: set_max_lines$, winio$
  Use screen,  Only: cffont, dts, dts_adj, f_scale
!
#endif
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ncol ! - number of columns (width) -
  Integer, Intent(In) :: nrow ! - number of rows -
!
  Logical, Intent(In), Optional :: lscale ! - scale window dimensions? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: i_out ! - specified output unit number -
  Integer, Intent(In), Optional :: maxl  ! - maximum number of lines -
!
  Integer(Kind=i_handle), Intent(In), Optional :: ih ! - handle -
!
  Character(Len=*), Intent(In), Optional :: cbox ! - box title -
!
  Logical, Intent(In), Optional :: lhs ! - add horizontal scroll bar? -
  Logical, Intent(In), Optional :: lvs ! - add vertical scroll bar? -
!
! Output scalars
  Integer, Intent(Out) :: iout ! - output unit number -
!
#if GUI == 1
! Locals
!
! Local scalars
  Real(Kind=rp) :: rscale ! - scaling -
!
  Character :: copt ! - option separator -
!
  Character(Len=64) :: ccfmt ! - ClearWin+ format statement -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Ceiling, Present, Real, Trim
!
! Executable Statements
!
! Set local font
  Call set_font (cfont=cffont)
!
#endif
! Assign output unit number
  If (Present(i_out)) Then
     iout = i_out
  Else
     iout = get_output_unit()
  End If
!
! Open Clearwin+ window
#if GUI == 1
  If (Present(cbox)) Call box_open (.true., .false., &
                          title=cbox)
! - add any handle -
  ccfmt = '%'
  If (Present(ih)) ccfmt = Trim(ccfmt)//'`'
! - add generic size -
  ccfmt = Trim(ccfmt)//'*.*cw'
! - add any options -
  copt = '['
  If (Present(lhs)) Then ! - horizontal scrolling -
     If (lhs) Then
        ccfmt = Trim(ccfmt)//'[hscroll'
        copt = ','
     Endif
  Endif
  If (Present(lvs)) Then ! - vertical scrolling -
     If (lvs) Then
        ccfmt = Trim(ccfmt)//copt//'vscroll'
        copt = ','
     Endif
  Endif
  ccfmt = Trim(ccfmt)//copt//'local_font'
  ccfmt = Trim(ccfmt)//']'
! - set scaling -
  If (lscale) Then
     rscale = dts*f_scale
  Else
     rscale = dts_adj
  End If
! - open window -
  If (Present(ih)) Then
     iw = winio$(Trim(ccfmt)//'&', Ceiling(Real(ncol,Kind=rp)*rscale), Ceiling(Real(nrow,Kind=rp)*rscale), iout, ih)
  Else
     iw = winio$(Trim(ccfmt)//'&', Ceiling(Real(ncol,Kind=rp)*rscale), Ceiling(Real(nrow,Kind=rp)*rscale), iout)
  End If
! - close box -
  If (Present(cbox)) Call box_close (.true.)
!
! Set maximum number of lines
  If (Present(maxl)) Call set_max_lines$ (ih, maxl)
#endif
!
  Return
 End Subroutine open_clearwin
!
!
!
#if GUI == 1
 Subroutine add_button_lab (cbt, ibt, lab, &
            ldef, igrey, iwdth, iset1, ival1, iset2, ival2, rset1, rval1, scbf, cbf)
!
! Adds a button
!
! Modules
  Use labels, Only: label, &
                    make_label
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ibt ! - button number -
!
  Character(Len=*), Intent(In) :: cbt   ! - button type -
!
  Type(label), Intent(In) :: lab ! - button label -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: igrey ! - greyed button flag -
  Integer, Intent(In), Optional :: iwdth ! - width of button -
  Integer, Intent(In), Optional :: ival1 ! - value to set to iset1 -
  Integer, Intent(In), Optional :: ival2 ! - value to set to iset2 -
!
  Real(Kind=rp), Intent(In), Optional :: rval1 ! - value to set to rset1 -
!
  Character(Len=*), Intent(In), Optional :: scbf ! - Silverfrost call-back function -
!
  Logical, Intent(In), Optional :: ldef ! - set as default? -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: iset1 ! - variable to set -
  Integer, Intent(InOut), Optional :: iset2 ! - variable to set -
!
  Real(Kind=rp), Intent(InOut), Optional :: rset1 ! - variable to set -
!
! Procedure arguments
  Integer, External, Optional :: cbf ! - call-back function -
!
! Locals
!
! Local scalars
  Character(Len=32) :: c_lab ! - formatted label text -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Construct label
  c_lab = Trim(make_label('T', lab, &
          lak=.true.))
!
! Add menu item
  Call add_button_txt (cbt, ibt, Trim(c_lab), &
       ldef=ldef, igrey=igrey, iwdth=iwdth, iset1=iset1, ival1=ival1, iset2=iset2, ival2=ival2, rset1=rset1, rval1=rval1, &
       scbf=scbf, cbf=cbf)
!
  Return
 End Subroutine add_button_lab
!
!
!
 Subroutine add_button_txt (cbt, ibt, c_lab, &
            ldef, igrey, iwdth, iset1, ival1, iset2, ival2, rset1, rval1, scbf, cbf)
!
! Adds a button
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ibt ! - button number -
!
  Character(Len=*), Intent(In) :: cbt   ! - button type -
  Character(Len=*), Intent(In) :: c_lab ! - button label -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: igrey ! - greyed button flag -
  Integer, Intent(In), Optional :: iwdth ! - width of button -
  Integer, Intent(In), Optional :: ival1 ! - value to set to iset1 -
  Integer, Intent(In), Optional :: ival2 ! - value to set to iset2 -
!
  Real(Kind=rp), Intent(In), Optional :: rval1 ! - value to set to rset1 -
!
  Character(Len=*), Intent(In), Optional :: scbf    ! - Silverfrost call-back function -
!
  Logical, Intent(In), Optional :: ldef ! - set as default? -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: iset1 ! - variable to set -
  Integer, Intent(InOut), Optional :: iset2 ! - variable to set -
!
  Real(Kind=rp), Intent(InOut), Optional :: rset1 ! - variable to set -
!
! Procedure arguments
  Integer, External, Optional :: cbf ! - call-back function -
!
! Locals
!
! Local scalars
  Character(Len=64) :: cbtfmt ! - button format -
  Character(Len= 8) :: cwdth  ! - button width -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Construct button
  If (ibt == 1) THen
     cbtfmt = '%cn%'                             ! - centre at the first button -
  Else
     cbtfmt = '  %'
  End If
  If (Present(igrey)) cbtfmt = Trim(cbtfmt)//'~' ! - add grey flag -
  If (Present(cbf) .or. Present(scbf) .or. &
     (Present(iset1) .and. Present(ival1)) .or. &
     (Present(iset1) .and. Present(ival1)) .or. &
     (Present(rset1) .and. Present(rval1))) Then
     cbtfmt = Trim(cbtfmt)//'^'                  ! - add callback -
  End If
  If (Present(iwdth)) Then
     Write (Unit=cwdth, Fmt='(I0)') iwdth
     cbtfmt = Trim(cbtfmt)//cwdth                ! - add button width -
  End If
  If (Present(ldef)) Then
     If (ldef) cbtfmt = Trim(cbtfmt)//'`'        ! - add default -
  End If
  cbtfmt = Trim(cbtfmt)//Trim(cbt)//'@'
  If (ibt /= 0) cbtfmt = Trim(cbtfmt)//'&'
!
! Construct button
  If (Present(igrey)) Then
     If (Present(rset1) .and. Present(rval1)) Then
        If (Present(iset1) .and. Present(ival1)) Then
           If (Present(iset2) .and. Present(ival2)) Then
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2)
                 End If
              End If
           Else ! - no second integer value to set -
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', 'SETF', rset1, rval1, 'SET', iset1, ival1)
                 End If
              End If
           End If
        Else ! - no first integer value to set -
           If (Present(scbf)) Then
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, '+', '+', 'SETF', rset1, rval1, scbf, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, '+', 'SETF', rset1, rval1, scbf)
              End If
           Else ! - no Silverfrost callback -
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, '+', 'SETF', rset1, rval1, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, 'SETF', rset1, rval1)
              End If
           End If
        End If
     Else ! - no real value to set -
        If (Present(iset1) .and. Present(ival1)) Then
           If (Present(iset2) .and. Present(ival2)) Then
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', '+', 'SET', iset1, ival1, 'SET', iset2, ival2, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', 'SET', iset1, ival1, 'SET', iset2, ival2, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', 'SET', iset1, ival1, 'SET', iset2, ival2, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', 'SET', iset1, ival1, 'SET', iset2, ival2)
                 End If
              End If
           Else ! - no second integer value to set -
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', '+', 'SET', iset1, ival1, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', 'SET', iset1, ival1, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, '+', 'SET', iset1, ival1, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         igrey, 'SET', iset1, ival1)
                 End If
              End If
           End If
        Else ! - no first integer value to set -
           If (Present(scbf)) Then
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, '+', scbf, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, scbf)
              End If
           Else ! - no Silverfrost callback -
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      igrey)
              End If
           End If
        End If
     End If
  Else ! - no grey flag -
     If (Present(rset1) .and. Present(rval1)) Then
        If (Present(iset1) .and. Present(ival1)) Then
           If (Present(iset2) .and. Present(ival2)) Then
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, 'SET', iset2, ival2)
                 End If
              End If
           Else ! - no second integer value to set -
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', 'SETF', rset1, rval1, 'SET', iset1, ival1, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', 'SETF', rset1, rval1, 'SET', iset1, ival1)
                 End If
              End If
           End If
        Else ! - no first integer value to set -
           If (Present(scbf)) Then
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      '+', '+', 'SETF', rset1, rval1, scbf, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      '+', 'SETF', rset1, rval1, scbf)
              End If
           Else ! - no Silverfrost callback -
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      '+', 'SETF', rset1, rval1, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      'SETF', rset1, rval1)
              End If
           End If
        End If
     Else ! - no real value to set -
        If (Present(iset1) .and. Present(ival1)) Then
           If (Present(iset2) .and. Present(ival2)) Then
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', '+', 'SET', iset1, ival1, 'SET', iset2, ival2, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', 'SET', iset1, ival1, 'SET', iset2, ival2, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', 'SET', iset1, ival1, 'SET', iset2, ival2, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', 'SET', iset1, ival1, 'SET', iset2, ival2)
                 End If
              End If
           Else ! - no second integer value to set -
              If (Present(scbf)) Then
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', '+', 'SET', iset1, ival1, scbf, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', 'SET', iset1, ival1, scbf)
                 End If
              Else ! - no Silverfrost callback -
                 If (Present(cbf)) Then
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         '+', 'SET', iset1, ival1, cbf)
                 Else ! - no callback -
                    iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                         'SET', iset1, ival1)
                 End If
              End If
           End If
        Else ! - no first integer value to set -
           If (Present(scbf)) Then
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      '+', scbf, cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      scbf)
              End If
           Else ! - no Silverfrost callback -
              If (Present(cbf)) Then
                 iw = winio$(Trim(cbtfmt), Trim(c_lab), &
                      cbf)
              Else ! - no callback -
                 iw = winio$(Trim(cbtfmt), Trim(c_lab))
              End If
           End If
        End If
     End If
  End If
!
  Return
 End Subroutine add_button_txt
!
!
!
 Subroutine adjust_y (dy)
!
! Adjusts the vertical position of the next control
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In) :: dy ! - ertical adjustment -
!
! Executable Statements
!
! Adjust the vertical position of the next control
  iw = winio$('%dy&', dy)
!
  Return
 End Subroutine adjust_y
!
!
!
#endif
 Subroutine box_open (lbold, lshade, &
            title, c_opts)
!
! Opens a window box
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lbold  ! - bold box text? -
  Logical, Intent(In) :: lshade ! - shaded background? -
!
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: c_opts ! - options -
  Character(Len=*), Intent(In), Optional :: title  ! - box text -
!
#if GUI == 1
! Locals
!
! Local scalars
  Character(Len=64) :: cbox ! - box command -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Create box command
  If (lshade) Then
     cbox = '%`ob['
  Else
     cbox = '%ob['
  End If
  If (Present(c_opts)) cbox = Trim(cbox)//c_opts//','
  If (Present(title)) Then
     cbox = Trim(cbox)//'named_c]['//title//']&'
     iw = winio$('%tc&', icolh)
  Else
     cbox = Trim(cbox)//'invisible]&'
  End If
  If (lbold) iw = winio$('%bf&')
!
! Open box
  iw = winio$(Trim(cbox))
!
! Restore default font
  If (lbold) iw = winio$('%`bf&')
  If (Present(title)) iw = winio$('%tc&', -1)
#else
  Call print_text (' ')
  If (Present(title)) Call print_text (title//':')
#endif
!
  Return
 End Subroutine box_open
!
!
!
 Subroutine boxes_open (ncol, nrow, linvisible, lshade)
!
! Opens a structured box
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ncol ! - number of columns -
  Integer, Intent(In) :: nrow ! - number of rows -
!
  Logical, Intent(In) :: linvisible ! - invisible? (.false. implies no_border) -
  Logical, Intent(In) :: lshade     ! - shaded background? -
!
#if GUI == 1
!
! Locals
!
! Local scalars
  Character(Len=32) :: cfmt ! - format -
!
! Executable Statements
!
! Open box
  cfmt = '%'
  IF (lshade) cfmt = Trim(cfmt)//'`'
  cfmt = Trim(cfmt)//'*.*ob['
  If (linvisible) Then
     cfmt = Trim(cfmt)//'invisible'
  Else
     cfmt = Trim(cfmt)//'no_border'
  End If
  cfmt = Trim(cfmt)//']&'
  iw = winio$(Trim(cfmt), ncol, nrow)
#else
  If (.not.linvisible) Call print_text (' ')
#endif
!
  Return
 End Subroutine boxes_open
!
!
!
 Subroutine box_close (lbold)
!
! Closes a window box
!
#if GUI == 1
! Modules
  Use clrwin$,      Only: winio$
  Use data_numbers, Only: zero
!
#endif
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lbold ! - bold box text? -
!
#if GUI == 1
! Executable Statements
!
! Close box
  If (lbold) Then
     iw = winio$('%bf%tc%cb%tc%`bf&', icolh, -1)
  Else
     iw = winio$('%tc%cb%tc&', icolh, -1)
  End If
  Call set_text_size (zero)
  Call gui_creturn (.true.)
  Call set_text_size (one)
#else
  If (lbold) Call print_text (' ')
#endif
!
  Return
 End Subroutine box_close
!
!
!
#if GUI == 1
 Subroutine ganged_rbs (nrb, irbs, &
            crbs, frb)
!
! Constructs a set of ganged radio buttons
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nrb ! - number of radio buttons to gang -
!
! Input arrays
! - optional input arrays -
  Character(Len=*), Dimension(:), Intent(In), Optional :: crbs ! - radio button names -
!
! Input/output arrays
  Integer, Dimension(:), Intent(InOut) :: irbs ! - radio button identifiers -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: frb ! - radio button Callback function -
!
! Locals
!
! Local scalars
  Integer :: i ! - radio button index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Gang items
  iw = winio$('%*`~ga&', nrb, irbs(:))
!
! Construct radio-buttons
  If (Present(crbs)) Then
     If (Present(frb)) Then
        iw = winio$('%^rb@&', crbs(1), irbs(1), frb)
        Do i = 2, nrb
           Call gui_creturn (.true.)
           iw = winio$('%^rb@&', crbs(i), irbs(i), frb)
        End Do
     Else
        iw = winio$('%rb@&', crbs(1), irbs(1))
        Do i = 2, nrb
           Call gui_creturn (.true.)
           iw = winio$('%rb@&', crbs(i), irbs(i))
        End Do
     End If
     Call gui_creturn (.false.)
  End If
!
  Return
 End Subroutine ganged_rbs
#else
!
!
!
 Subroutine ganged_rbs_a (nrb, irbs, crbs, &
            lprompt, ns)
!
! Prompts for one of a set of ganged options (array version)
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nrb ! - number of radio buttons to gang -
!
! - optional input scalars -
  Logical, Intent(In), Optional :: lprompt ! - prompt required? -
!
! Input/output arrays
  Integer, Dimension(:), Intent(InOut) :: irbs ! - radio button identifiers -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: crbs ! - radio button names -
!
! - optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: ns ! - non-sequential numbers -
!
! Locals
!
! Local scalars
  Integer :: irb ! - radio button index -
!
! Executable Statements
!
! Gang items
  Call prompt_list (' ', .false., crbs, nrb, irb, &
       lprompt=lprompt, ns=ns)
  If (irb > 0) Then
     irbs(:) = 0
     irbs(irb) = 1
  End If
!
  Return
 End Subroutine ganged_rbs_a
!
!
!
 Subroutine ganged_rbs_s (nrb, irb, crbs, &
            lprompt, ns)
!
! Prompts for one of a set of ganged options (scalar version)
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nrb ! - number of radio buttons to gang -
!
! - optional input scalars -
  Logical, Intent(In), Optional :: lprompt ! - prompt required? -
!
! Input/output scalars
  Integer, Intent(InOut) :: irb ! - radio button identifier -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: crbs ! - radio button names -
!
! - optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: ns ! - non-sequential numbers -
!
! Executable Statements
!
! Gang items
  Call prompt_list (' ', .false., crbs, nrb, irb, &
       lprompt=lprompt, ns=ns)
!
  Return
 End Subroutine ganged_rbs_s
#endif
!
!
!
#if GUI == 1
 Subroutine prompt_integer (cprmpt, i, &
            ifly, itab, iwid, ilow, ihgh, fcb1, fcb2, lbold, lfocus, ihandle)
#else
 Subroutine prompt_integer (cprmpt, i, &
            ifail, ifly, itab, iwid, ilow, ihgh, n, ns, altns, fcb1, fcb2, lbold, lfocus, lzero)
#endif
!
! Prompts for integer
!
! Modules
#if GUI == 1
  Use clrwin$,      Only: winio$
  Use data_numbers, Only: ihuge
#endif
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itab  ! - tab point -
  Integer, Intent(In), Optional :: iwid  ! - width of prompt -
  Integer, Intent(In), Optional :: ilow  ! - lower limit -
  Integer, Intent(In), Optional :: ihgh  ! - upper limit -
  Integer, Intent(In), Optional :: ifly  ! - flywheel increment -
!
  Logical, Intent(In), Optional :: lbold  ! - bold? -
  Logical, Intent(In), Optional :: lfocus ! - check on focus loss? -
#if GUI == 0
  Logical, Intent(In), Optional :: lzero  ! - return if response is 0? -
#endif
!
! Input-output scalars
  Integer, Intent(InOut) :: i ! - integer -
!
! Output scalars
#if GUI == 1
! - optional output scalars -
  Integer(Kind=i_handle), Intent(Out), Optional :: ihandle ! - control handle -
!
#else
  Integer, Intent(Out), Optional :: n     ! - value from a list (ns must also be present) -
  Integer, Intent(Out), Optional :: ifail ! - error indicator -
!
! Input arrays
! - optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: ns    ! - list of valid responses -
  Integer, Dimension(:), Intent(In), Optional :: altns ! - alternative non-sequential numbers -
!
#endif
! Procedure arguments
! - optional procedures -
  Integer, External, Optional :: fcb1 ! - Callback function -
  Integer, External, Optional :: fcb2 ! - Callback function -
!
! Locals
!
! Local scalars
#if GUI == 1
  Integer :: iwd ! - width -
#else
  Integer :: j    ! - index -
  Integer :: ierr ! - error indicator -
!
  Logical :: lns ! - perform search for non-sequential numbers? -
#endif
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
#if GUI == 0
  Intrinsic All, Max, Min
#endif
!
! Executable Statements
!
! Prompt for integer
#if GUI == 1
! - prompt -
  Call print_prompt (cprmpt, &
       itab=itab, lbold=lbold, lfocus=lfocus)
! - set limits -
  If (Present(ilow)) Then
     If (Present(ihgh)) Then
        iw = winio$('%il&', ilow, ihgh)
     Else
        iw = winio$('%il&', ilow, ihuge)
     End If
  Else
     If (Present(ihgh)) Then
        iw = winio$('%il&', -ihuge, ihgh)
     Else
        iw = winio$('%il&', -ihuge, ihuge)
     End If
  End If
! - flywheel increment -
  If (Present(ifly)) Then
     iw = winio$('%dd&', ifly)
  Else
     iw = winio$('%dd&', 1)
  End If
! - width -
  If (Present(iwid)) Then
     iwd=iwid
  Else
     iwd=6
  End If
! - prompt for integer -
  If (Present(fcb1) .and. Present(fcb2)) Then
     iw = winio$('%^*rd&', iwd, i, '+', fcb1, fcb2)
  Else If (Present(fcb1)) Then
     iw = winio$('%^*rd&', iwd, i, fcb1)
  Else If (Present(fcb2)) Then
     iw = winio$('%^*rd&', iwd, i, fcb2)
  Else
     iw = winio$('%*rd&', iwd, i)
  End If
  Call gui_creturn (.false.)
! - get control handle -
  If (Present(ihandle)) iw = winio$('%lc&', ihandle)
#else
1 Call print_prompt (cprmpt)
  Read (Unit=*, Fmt=*, Err=2) i
! - return if response is 0 -
  If (Present(ifail)) ifail = 0
  If (Present(lzero)) Then
     If (lzero .and. i == 0) Then
        If (Present(ifail)) ifail = 1
        Return
     End If
  End If
! - check whether value is within range -
  If (Present(ilow)) Then
     If (i < ilow) GoTo 2
  End If
  If (Present(ihgh)) Then
     If (i > ihgh) GoTo 2
  End If
! - check against alternative non-sequential numbers -
  If (Present(ns)) Then
     If (Present(altns)) Then
        If (Any(altns(:) == i)) Then
           If (Present(n)) n = i
           j = 1
           Do
              If (i == altns(j)) Exit
              j = j + 1
           End Do
           i = j
           lns = .false.
        Else
           lns = .true.
        End If
     Else
        lns = .true.
     End If
! - check against non-sequential numbers -
     If (lns) Then
        If (All(ns(:) /= i)) GoTo 2
        If (Present(n)) n = i
        j = 1
        Do
           If (i == ns(j)) Exit
           j = j + 1
        End Do
        i = j
     End If
  End If
! - call-back functions -
  If (Present(fcb1)) Then
     ierr = fcb1()
     If (Present(ifail)) ifail = ierr
     If (ierr /= 2 .and. ierr /= 0) GoTo 2
  End If
  If (Present(fcb2)) Then
     ierr = fcb2()
     If (Present(ifail)) ifail = ierr
     If (ierr /= 2 .and. ierr /= 0) GoTo 2
  End If
  Return
!
2 Select Case (iprompt)
   Case (0)
     GoTo 1
   Case (1)
     If (Present(ilow)) i = Max(i,ilow)
     If (Present(ihgh)) i = Min(i,ihgh)
   Case (2)
     Stop
  End Select
#endif
!
  Return
 End Subroutine prompt_integer
!
!
!
 Subroutine prompt_boolean (cprmpt, iyn)
!
! Prompts for logical
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
! Input-output scalars
  Integer, Intent(InOut) :: iyn ! - boolean -
!
! Executable Statements
!
! Prompt for boolean
#if GUI == 1
  iw = winio$('%`rb@&', cprmpt, iyn)
  Call gui_creturn (.false.)
#else
  iyn = prompt_yesno(.true., .false., cprmpt)
#endif
!
  Return
 End Subroutine prompt_boolean
!
!
!
 Subroutine prompt_cvalue (ctxt, cvalue, &
            itab, iwid)
!
! Prompts for a character value
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
#endif
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctxt   ! - text -
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: cvalue ! - character value -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itab ! - tab position -
  Integer, Intent(In), Optional :: iwid ! - width of value -
!
! Executable Statements
!
! Print text
#if GUI == 1
  Call print_prompt (ctxt, &
       itab=itab, lbold=.true.)
  iw = winio$('%*rs&', iwid, cvalue)
  Call gui_creturn (.false.)
#else
  Call print_prompt (ctxt)
  Read (Unit=*, Fmt=*) cvalue
#endif
!
  Return
 End Subroutine prompt_cvalue
!
!
!
 Subroutine prompt_real (cprmpt, r, &
            itab, iwid, rfly, rlow, rhgh, elow, ehgh, lfocus, lbold)
!
! Prompts for real
!
#if GUI == 1
! Modules
  Use clrwin$,      Only: winio$
  Use data_numbers, Only: bignum, eps
!
#endif
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itab ! - tab point -
  Integer, Intent(In), Optional :: iwid ! - width of prompt -
!
  Real(Kind=rp), Intent(In), Optional :: rfly ! - flywheel increment -
  Real(Kind=rp), Intent(In), Optional :: rlow ! - inclusive lower limit -
  Real(Kind=rp), Intent(In), Optional :: rhgh ! - inclusive upper limit -
  Real(Kind=rp), Intent(In), Optional :: elow ! - exclusive lower limit -
  Real(Kind=rp), Intent(In), Optional :: ehgh ! - exclusive upper limit -
!
  Logical, Intent(In), Optional :: lfocus ! - check on focus loss? -
  Logical, Intent(In), Optional :: lbold  ! - bold? -
!
! Input-output scalars
  Real(Kind=rp), Intent(InOut) :: r ! - real -
!
#if GUI == 1
! Locals
!
! Local scalars
  Integer :: iwd ! - width -
!
  Real(Kind=rp) :: rmn ! - lower limit -
  Real(Kind=rp) :: rmx ! - upper limit -
!
#endif
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
#if GUI == 0
  Intrinsic Max, Min
#endif
!
! Executable Statements
!
! Prompt for real
#if GUI == 1
! - prompt -
  Call print_prompt (cprmpt, &
       itab=itab, lbold=lbold, lfocus=lfocus)
! - set limits -
  If (Present(elow)) Then
     rmn = elow + eps
  Else
     rmn = -bignum
  End If
  If (Present(ehgh)) Then
     rmx = ehgh - eps
  Else
     rmx = bignum
  End If
  If (Present(rlow)) Then
     rmn = rlow
  Else
     rmn = -bignum
  End If
  If (Present(rhgh)) Then
     rmx = rhgh
  Else
     rmx = bignum
  End If
  iw = winio$('%fl&', rmn, rmx)
! - flywheel increment -
  If (Present(rfly)) Then
     iw = winio$('%df&', rfly)
  Else
     iw = winio$('%df&', one)
  End If
! - width -
  If (Present(iwid)) Then
     iwd = iwid
  Else
     iwd = 6
  End If
! - prompt for real -
  iw = winio$('%*rf&', iwd, r)
  Call gui_creturn (.false.)
#else
1 Call print_prompt (cprmpt)
  Read (Unit=*, Fmt=*, Err=2) r
  If (Present(rlow)) Then
     If (r < rlow) GoTo 2
  End If
  If (Present(elow)) Then
     If (.not.(r > elow)) GoTo 2
  End If
  If (Present(rhgh)) Then
     If (r > rhgh) GoTo 2
  End If
  If (Present(ehgh)) Then
     If (.not.(r < ehgh)) GoTo 2
  End If
  Return
!
2 Select Case (iprompt)
   Case (0)
     GoTo 1
   Case (1)
     If (Present(rlow)) r = Max(r, rlow)
     If (Present(rhgh)) r = Min(r, rhgh)
   Case (2)
     Stop
  End Select
#endif
!
  Return
 End Subroutine prompt_real
!
!
!
 Subroutine prompt_pdate (cprmpt, isq, adate, &
            itab, dlimits, lfocus, lbold, lyr, modify1, modify2)
!
! Prompts for a parsed date
!
! Modules
#if GUI == 1
  Use clrwin$,             Only: winio$
  Use data_numbers,        Only: ihuge
  Use data_time_constants, Only: cmon,isq_mn,isq_sn,isq_so,mdm,nmn
  Use time,                Only: pdate
#else
  Use data_time_constants, Only: isq_yr, nmn
  Use labels,              Only: ca_seq_l
  Use time,                Only: pdate, &
                                 Operator(<), Operator(>), &
                                 ndays
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isq ! - sequence -
!
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itab ! - tab point -
!
  Logical, Intent(In), Optional :: lbold  ! - bold? -
  Logical, Intent(In), Optional :: lfocus ! - focus? -
  Logical, Intent(In), Optional :: lyr ! - prompt for year? -
!
! Input/output scalars
  Type(pdate), Intent(InOut) :: adate ! - parsed date -
!
! Input arrays
! - optional input arrays -
  Type(pdate), Dimension(2), Intent(In), Optional :: dlimits ! - date limits -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: modify1 ! - first modify function -
  Integer, External, Optional :: modify2 ! - second modify function -
!
! Locals
!
! Local scalars
  Logical :: lpyr ! - prompt for year? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
#if GUI == 0
  Intrinsic Trim
#endif
!
! Executable Statements
!
! Identify whether to prompt for year
  If (Present(lyr)) Then
     lpyr = lyr
  Else
     lpyr = .true.
  End If
#if GUI == 1
!
! Prompt for parsed date
! - prompt -
  Call print_prompt (cprmpt, &
       itab=itab, lbold=lbold, lfocus=lfocus)
  Select Case (isq)
   Case (1:)
     iw = winio$('%il&', 1, mdm)
     If (Present(modify1)) Then
        iw = winio$('%dd%^2rd - &', 1, adate%idy, modify1)
        iw = winio$('%^`5.12ls - &', cmon, nmn, adate%imn, modify1)
     Else
        iw = winio$('%dd%2rd - &', 1, adate%idy)
        iw = winio$('%`5.12ls - &', cmon, nmn, adate%imn)
     End If
   Case (isq_mn, isq_sn, isq_so)
     iw = winio$('%il&', 1, mdm)
     If (Present(modify1)) Then
        iw = winio$('%^`5.12ls - &', cmon, nmn, adate%imn, modify1)
     Else
        iw = winio$('%`5.12ls - &', cmon, nmn, adate%imn)
     End If
  End Select
  If (lpyr) Then
     If (Present(dlimits)) Then
        iw = winio$('%il&', dlimits(1)%iyr, dlimits(2)%iyr)
     Else
        iw = winio$('%il&', 0, ihuge)
     End If
     If (Present(modify2)) Then
        If (Present(modify1)) Then
           iw = winio$('%dd%^5rd&', 1, adate%iyr, '+', modify1, modify2)
        Else
           iw = winio$('%dd%^5rd&', 1, adate%iyr, modify2)
        End If
     Else
        If (Present(modify1)) Then
           iw = winio$('%dd%^5rd&', 1, adate%iyr, modify1)
        Else
           iw = winio$('%dd%5rd&', 1, adate%iyr)
        End If
     End If
  End If
  Call gui_creturn (.false.)
#else
!
! Prompt for parsed date
  Call print_text (cprmpt)
! - prompt for year -
  If (lpyr) Then
     If (Present(dlimits)) Then
        Call prompt_integer (Trim(ca_seq_l(1)), adate%iyr, &
             ilow=dlimits(1)%iyr, ihgh=dlimits(2)%iyr)
     Else
        Call prompt_integer (Trim(ca_seq_l(1)), adate%iyr)
     End If
  End If
! - prompt for month -
  If (isq /= isq_yr) Then
     Call prompt_integer (Trim(ca_seq_l(2)), adate%imn, &
          ilow=1, ihgh=nmn)
     If (isq > 0) Then
        Call prompt_integer (Trim(ca_seq_l(3)), adate%imn, &
             ilow=1, ihgh=ndays(adate%iyr, adate%imn))
     End If
  End If
  If (Present(modify2)) iw = modify2()
#endif
!
  Return
 End Subroutine prompt_pdate
!
!
!
 Function prompt_bullet(ctitle, nbull, ibulls, cbulls, f_update, &
          lprompt, ns) &
          Result (prompt)
!
! Prompts for a selection from a bulleted list
!
! Function type
  Integer :: prompt
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nbull ! - number of bullets -
!
  Character(Len=*), Intent(In) :: ctitle ! - bulleted list title -
!
! - optional input scalars -
  Logical, Intent(In), Optional :: lprompt ! - prompt required? -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: cbulls ! - bullet labels -
!
! - optional input arrays -
  Integer, Dimension(:), Intent(In), Optional :: ns ! - non-sequential numbers -
!
#if GUI == 1
! Input/output arrays
  Integer, Dimension(:), Intent(InOut) :: ibulls ! - bullet identifiers -
#else
! Output arrays
  Integer, Dimension(:), Intent(Out) :: ibulls ! - bullet identifiers -
#endif
!
! Procedure arguments
  Integer, External :: f_update
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for bulleted item
  Call box_open (.true., .false., &
       title=Trim(ctitle))
#if GUI == 1
  Call ganged_rbs (nbull, ibulls, &
       crbs=cbulls, frb=f_update)
#else
  ibulls(:) = 0
  Call ganged_rbs (nbull, ibulls, cbulls, &
       lprompt=lprompt, ns=ns)
#endif
! - close box -
  Call box_close (.true.)
  prompt = 1
!
  Return
 End Function prompt_bullet
!
!
!
 Subroutine prompt_list (cprmpt, ldrop, clist, nlist, ilist, &
            igrey, itab, iwdth, ns, altns, ilistn, lprompt, cbf1, cbf2)
!
! Prompts for an item from a list
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
#else
  Use maths,  Only: magnitude
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nlist ! - number of items in list -
!
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
  Logical, Intent(In) :: ldrop ! - drop-down list? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: igrey ! - grey flag -
  Integer, Intent(In), Optional :: itab  ! - tab point -
  Integer, Intent(In), Optional :: iwdth ! - width of list -
!
  Logical, Intent(In), Optional :: lprompt ! - prompt required? -
!
! Input/output scalars
  Integer, Intent(InOut) :: ilist ! - selected item -
!
! Output scalars
! - optional output scalars -
  Integer, Intent(Out), Optional :: ilistn ! - associated number for given selection -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: clist ! - list items -
!
! - optional input scalars -
  Integer, Dimension(:), Intent(In), Optional :: ns    ! - non-sequential numbers -
  Integer, Dimension(:), Intent(In), Optional :: altns ! - alternative non-sequential numbers -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: cbf1 ! - first call-back function -
  Integer, External, Optional :: cbf2 ! - second call-back function -
!
! Locals
!
! Local scalars
#if GUI == 0
  Integer :: i ! - index -
!
  Logical :: lp ! - prompt required? -
!
#endif
  Character(Len=32) :: cfmt  ! - format statement -
  Character(Len= 8) :: cwdth ! - width -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 0
  Intrinsic MaxVal
#endif
!
! Executable Statements
!
! Prompt for listed item
  If (Present(iwdth)) Then
     Write (Unit=cwdth, Fmt='(I0)') iwdth
  Else
     cwdth = ' '
  End If
#if GUI == 1
  cfmt = '%'
  If (ldrop) cfmt = Trim(cfmt)//'`'                            ! - drop-down menu -
  If (Present(igrey)) cfmt = Trim(cfmt)//'~'                   ! - add grey flag -
  If (Present(cbf1) .or. Present(cbf2)) cfmt = Trim(cfmt)//'^' ! - add callback function
  cfmt = Trim(cfmt)//Trim(cwdth)//'ls'
  Call gui_creturn (.true.)
  Call print_prompt (cprmpt, &
       lbold=.true., itab=itab)
  If (Present(cbf1)) Then
     If (Present(cbf2)) Then
        If (Present(igrey)) Then
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, igrey, '+', cbf1, cbf2)
        Else ! - no greying -
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, '+', cbf1, cbf2)
        End If
     Else ! - no callback function -
        If (Present(igrey)) Then
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, igrey, cbf1)
        Else ! - no greying -
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, cbf1)
        End If
     End If
  Else
     If (Present(cbf2)) Then
        If (Present(igrey)) Then
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, igrey, cbf2)
        Else ! - no greying -
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, cbf2)
        End If
     Else ! - no callback function -
        If (Present(igrey)) Then
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist, igrey)
        Else ! - no greying -
           iw = winio$(Trim(cfmt)//'&', clist(:), nlist, ilist)
        End If
     End If
  End If
  Call gui_creturn (.false.)
  If (Present(ilistn)) ilistn = ilist
#else
  If (Present(lprompt)) Then
     lp = lprompt
  Else
     lp = .true.
  End If
  Call print_text (Trim(cprmpt))
  If (Present(ns)) Then
     Write (Unit=cfmt, Fmt='(A,I1,A)') '(I', magnitude(MaxVal(ns(:))), ',A)'
     Do i = 1, nlist
        Write (Unit=*, Fmt=cfmt) ns(i), '. '//clist(i)
     End Do
     If (lp) Then
        Call prompt_integer (' ', ilist, &
             ns=ns, altns=altns, n=ilistn, lzero=.true.)
     Else
        ilist = 0
     End If
  Else
     If (Present(iwdth)) Then
        Write (Unit=cfmt, Fmt='(A,I1,A)') '(I', magnitude(nlist)+1, ',A'//Trim(cwdth)//')'
     Else
        Write (Unit=cfmt, Fmt='(A,I1,A)') '(I', magnitude(nlist)+1, ',A)'
     End If
     Do i = 1, nlist
        Write (Unit=*, Fmt=cfmt) i, '. '//clist(i)
     End Do
     If (lp) Then
        Call prompt_integer (' ',ilist, &
             ilow=1, ihgh=nlist)
     Else
        ilist = 0
     End If
  End If
#endif
!
  Return
 End Subroutine prompt_list
!
!
!
 Subroutine prompt_select (cprmpt, nlist, isel, clist)
!
! Prompts for selection of items from a list
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
#else
  Use maths,  Only: magnitude
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nlist ! - number of items in list -
!
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
! Input/output arrays
  Integer, Dimension(:), Intent(InOut) :: isel ! - selected items -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: clist ! - list items -
!
! Locals
!
! Local scalars
  Integer :: i ! - list item index -
!
#if GUI == 0
  Character(Len=8) :: cfmt ! - format statement -
!
#endif
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for listed items
  Call gui_creturn (.true.)
  iw = winio$('%bf'//Trim(cprmpt)//':&')
  Do i = 1, nlist
     Call gui_creturn (.false.)
     iw = winio$('%`rb@&', Trim(clist(i)), isel(i))
  End Do
#else
  Write (Unit=*, Fmt='(A)') ' Enter selected items (enter 0 to quit)'
  Write (Unit=cfmt, Fmt='(A,I1,A)') '(I', magnitude(nlist)+1, ',A)'
  Do i = 1, nlist
     Write (Unit=*, Fmt=cfmt) i, '. '//clist(i)
  End Do
  Do
     Call prompt_integer (cprmpt,i, &
          ilow=0, ihgh=nlist)
     If (i == 0) Exit
     isel(i) = 1
  End Do
#endif
!
  Return
 End Subroutine prompt_select
!
!
!
 Subroutine prompt_file (cdesc, cname, ival, iset, cbf, &
            itab)
!
! Prompts for a file
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
  Use labels, Only: cg_browse_l
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ival ! - file identifier -
!
  Character(Len=*), Intent(In) :: cdesc ! - file description -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itab ! - tab point -
!
! Input/output scalars
  Integer, Intent(InOut) :: iset ! - file identifier -
!
  Character(Len=*), Intent(InOut) :: cname ! - file name -
!
! Procedure arguments
  Integer, External :: cbf ! - call-back function -
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
#endif
!
! Executable Statements
!
! Print file description
#if GUI == 1
  Call set_text_size (1.2_rp)
  iw = winio$('%bf&')         ! - use bold enlarged font -
  iw = winio$('%ws:&', cdesc) ! - print file description -
  iw = winio$('%`bf&')        ! - revert to default font -
  Call set_text_size (one)
!
! Indicate current file
  If (Present(itab)) Then
     Call set_tabs (1, (/itab/))
     iw = winio$('%ta&')
  Else
     Call print_text ('   ', &
          ladvance=.false.)
  End If
  iw = winio$('%ob%`25rs%cb&', cname)
!
! Add browse button
  Call add_button ('bt', 2, cg_browse_l, &
       iwdth=5, iset1=iset, ival1=ival, cbf=cbf)
#else
  Write (Unit=*, Fmt='(I3,2A)') ival, '. ', cdesc
#endif
!
  Return
 End Subroutine prompt_file
!
!
!
 Function prompt_yesno(lny, lcancel, msg1, &
          msg2, msg3, msg4, msg5) &
          Result (prompt)
!
! Creates a warning message window with a query of whether to continue
! Returns 1 for yes and 0 for no.
!
! Modules
#if GUI == 1
  Use labels, Only: l_cancel, l_no, l_yes
#endif
!
! Function type
  Integer :: prompt
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lny     ! - 'no' first? -
  Logical, Intent(In) :: lcancel ! - add cancel button? -
!
  Character(Len=*), Intent(In) :: msg1 ! - warning message -
!
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: msg2 ! - warning message -
  Character(Len=*), Intent(In), Optional :: msg3 ! - warning message -
  Character(Len=*), Intent(In), Optional :: msg4 ! - warning message -
  Character(Len=*), Intent(In), Optional :: msg5 ! - warning message -
!
! Locals
!
! Local scalars
#if GUI == 1
  Integer :: iyes ! - yes button number -
#else
  Character(Len=1) :: yn ! - yes/no response -
#endif
!
! Functions and Subroutines
!
! Intrinsic functions
#if GUI == 0
  Intrinsic Present
#endif
!
! Executable Statements
!
! Print warning message
#if GUI == 1
  Call print_warning (msg1, &
       msg2=msg2, msg3=msg3, msg4=msg4, msg5=msg5, nopause=.true.)
#else
  If (Present(msg5)) Then
     Call print_warning (msg1, &
          msg2=msg2, msg3=msg3, msg4=msg4, msg5=msg5//' (Y/N)', nopause=.true.)
  Else If (Present(msg4)) Then
     Call print_warning (msg1, &
          msg2=msg2, msg3=msg3, msg4=msg4//' (Y/N)', nopause=.true.)
  Else If (Present(msg3)) Then
     Call print_warning (msg1, &
          msg2=msg2, msg3=msg3//' (Y/N)', nopause=.true.)
  Else If (Present(msg2)) Then
     Call print_warning (msg1, &
          msg2=msg2//' (Y/N)', nopause=.true.)
  Else
     Call print_warning (msg1//' (Y/N)', &
          nopause=.true.)
  End If
#endif
!
! Prompt for yes/no response
#if GUI == 1
  Call gui_creturn (.false.)
  Call gui_creturn (.true.)
  If (lny) Then ! - no first -
     Call add_button ('bt', 1, l_no, &
          iwdth=5, ldef=.true.)
     Call add_button ('bt', 2, l_yes, &
          iwdth=5)
     If (lcancel) Call add_button ('bt', 3, l_cancel, &
                       iwdth=5)
     iyes = 2
     Call window_end ()
     If (iw == iyes) Then
        prompt = 1
     Else
        prompt = 0
     End If
  Else ! - yes first -
     Call add_button ('bt', 1, l_yes, &
          iwdth=5, ldef=.true.)
     Call add_button ('bt', 2, l_no, &
          iwdth=5)
     If (lcancel) Call add_button ('bt', 3, l_cancel, &
                       iwdth=5)
     Call window_end ()
     prompt = iw
  End If
#else
1 Read (Unit=*, Fmt='(A)') yn
  Select Case (yn)
   Case ('Y', 'y')
     prompt = 1
   Case ('N', 'n')
     prompt = 0
   Case Default
     Select Case (iprompt)
      Case (0)
        GoTo 1
      Case (1)
        prompt = 1
      Case (2)
        Stop
     End Select
  End Select
#endif
!
  Return
 End Function prompt_yesno
!
!
!
 Subroutine print_prompt (cprmpt, &
            cqual, itab, lbold, lfocus)
!
! Prints a prompt
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cprmpt ! - prompt -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: itab ! - tab point -
!
  Character(Len=*), Intent(In), Optional :: cqual ! - qualifying text -
!
  Logical, Intent(In), Optional :: lfocus ! - check on focus loss? -
  Logical, Intent(In), Optional :: lbold  ! - bold? -
!
#if GUI == 0
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim
!
#endif
! Executable Statements
!
! Prompt
#if GUI == 1
! - bold -
  If (Present(lbold)) Then
     If (lbold) iw = winio$('%bf&')
  End If
! - print description -
  iw = winio$('%ws&', cprmpt)
! - add qualifier -
  If (Present(cqual)) Then
     iw = winio$(' %`rs&', cqual)
  End If
! - restore standard font -
  If (Present(lbold)) Then
     If (lbold) iw = winio$('%`bf&')
  End If
! - add prompt separator -
  If (.not.Present(cqual)) iw = winio$(':&')
! - tab -
  If (Present(itab)) Then
     If (itab > 0) Call set_tabs (1, (/itab/))
     iw = winio$('%ta&')
  Else
     Call print_text ('  ', &
          ladvance=.false.)
  End If
! - focus -
  If (Present(lfocus)) Then
     If (lfocus) Then
        iw = winio$('%co[check_on_focus_loss]&')
     Else
        iw = winio$('%co[full_check]&')
     End If
  End If
#else
  If (Len_Trim(cprmpt) > 0) Then
     Write (Unit=*, Fmt='(1X,A)', Advance='no') cprmpt//':  '
  Else
     Write (Unit=*, Fmt='(1X,A)', Advance='no') '>  '
  End If
#endif
!
  Return
 End Subroutine print_prompt
!
!
!
 Subroutine print_cvalue (ctxt, cvalue, &
            iout, itab, iwid, ldepth)
!
! Prints a character value
!
! Modules
#if GUI == 1
  Use clrwin$,           Only: winio$
#else
  Use data_io_constants, Only: ioutstd
#endif
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctxt   ! - text -
  Character(Len=*), Intent(In) :: cvalue ! - character value -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
  Integer, Intent(In), Optional :: itab ! - tab position -
  Integer, Intent(In), Optional :: iwid ! - width of value -
!
  Logical, Intent(In), Optional :: ldepth ! - additional depth? -
!
! Locals
!
! Local scalars
  Integer :: ioutu ! - used output unit number -
!
  Character(Len=16) :: cfmt ! - format statement -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Print text
  If (Present(iout)) Then
     ioutu = iout
  Else
#if GUI == 0
     ioutu = ioutstd
#else
     ioutu = -1
     Call print_prompt (ctxt, &
          itab=itab, lbold=.true.)
     If (Present(ldepth)) Then
        If (ldepth) Then
           If (Present(iwid)) Then
              iw = winio$('%`*rs&', iwid, cvalue)
           Else
              iw = winio$('%`rs&', cvalue)
           End If
        Else
           If (Present(iwid)) Then
              iw = winio$('%`*rs[no_additional_depth]&', iwid, cvalue)
           Else
              iw = winio$('%`rs[no_additional_depth]&', cvalue)
           End If
        End If
     Else
        If (Present(iwid)) Then
           iw = winio$('%`*rs&', iwid, cvalue)
        Else
           iw = winio$('%`rs&', cvalue)
        End If
     End If
     Call gui_creturn (.false.)
#endif
  End If
  If (ioutu /= -1) Then
     If (Present(itab)) Then
        Write (Unit=cfmt, Fmt='(A,I0,A)') '(1X,2A,T', itab, ',A)'
        Write (Unit=ioutu, Fmt=cfmt) ctxt, ':', cvalue
     Else
        Call print_text (ctxt//':  '//cvalue, &
             iout=ioutu)
     End If
  End If
!
  Return
 End Subroutine print_cvalue
!
!
!
 Subroutine print_ivalue (ctxt, ivalue, &
            iout, itab, iwid, cqual)
!
! Prints an integer value
!
! Modules
#if GUI == 1
  Use clrwin$,           Only: winio$
#else
  Use data_io_constants, Only: ioutstd
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ivalue ! - integer value -
!
  Character(Len=*), Intent(In) :: ctxt ! - text -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
  Integer, Intent(In), Optional :: itab ! - tab position -
  Integer, Intent(In), Optional :: iwid ! - width of value -
!
  Character(Len=*), Intent(In), Optional :: cqual ! - qualifying text -
!
! Locals
!
! Local scalars
  Integer :: ioutu ! - used output unit number -
!
  Character(Len=16) :: cfmt   ! - format statement -
  Character(Len=16) :: cvalue ! - integer as a character string -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Print text
  If (Present(iout)) Then
     ioutu = iout
  Else
#if GUI == 0
     ioutu = ioutstd
#else
     ioutu = -1
     Call print_prompt (ctxt, &
          cqual=cqual, itab=itab, lbold=.true.)
     If (Present(iwid)) Then
        iw = winio$('%`*rd&', iwid, ivalue)
     Else
        iw = winio$('%`rd&', ivalue)
     End If
     Call gui_creturn (.false.)
#endif
  End If
  If (ioutu /= -1) Then
     If (Present(itab)) Then
        Write (Unit=cfmt, Fmt='(A,I0,A)') '(1X,2A,T', itab, ',I0)'
        Write (Unit=ioutu, Fmt=cfmt) ctxt, ':', ivalue
     Else
        Write (Unit=cvalue, Fmt='(I0)') ivalue
        Call print_text (ctxt//':  '//Trim(cvalue), &
             iout=ioutu)
     End If
  End If
!
  Return
 End Subroutine print_ivalue
!
!
!
 Subroutine print_rvalue (ctxt, rvalue, &
            iout, itab, iwid, idec, cqual)
!
! Prints a real value
!
! Modules
#if GUI == 1
  Use clrwin$,           Only: winio$
#else
  Use data_io_constants, Only: ioutstd
#endif
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In) :: rvalue ! - integer value -
!
  Character(Len=*), Intent(In) :: ctxt ! - text -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
  Integer, Intent(In), Optional :: itab ! - tab position -
  Integer, Intent(In), Optional :: iwid ! - width of value -
  Integer, Intent(In), Optional :: idec ! - number of decimal places -
!
  Character(Len=*), Intent(In), Optional :: cqual ! - qualifying text -
!
! Locals
!
! Local scalars
  Integer :: ioutu ! - used output unit number -
!
  Character(Len=16) :: cfmt   ! - format statement -
  Character(Len=16) :: cvalue ! - real as a character string -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustL, Present, Trim
!
! Executable Statements
!
! Print text
  If (Present(iout)) Then
     ioutu = iout
  Else
#if GUI == 0
     ioutu = ioutstd
#else
     ioutu = -1
     Call print_prompt (ctxt, &
          cqual=cqual, itab=itab, lbold=.true.)
     If (Present(iwid)) Then
        If (Present(idec)) Then
           iw = winio$('%`*.*rf&', iwid, idec, rvalue)
        Else
           iw = winio$('%`*rf&', iwid, rvalue)
        End If
     Else
        iw = winio$('%`rf&', rvalue)
     End If
     Call gui_creturn (.false.)
#endif
  End If
  If (ioutu /= -1) Then
     If (Present(iwid) .and. Present(idec)) Then
        Write (Unit=cfmt, Fmt='(A,I0,A,I0,A)') '(F', iwid, '.', idec, ')'
        Write (Unit=cvalue, Fmt=Trim(cfmt)) rvalue
     Else
        Write (Unit=cvalue, Fmt=*) rvalue
     End If
     If (Present(itab)) Then
        Write (Unit=cfmt, Fmt='(A,I0,A)') '(1X,2A,T', itab, ',A)'
        Write (Unit=ioutu, Fmt=cfmt) ctxt, ':', Trim(AdjustL(cvalue))
     Else
        Call print_text (ctxt//':  '//Trim(AdjustL(cvalue)), &
             iout=ioutu)
     End If
  End If
!
  Return
 End Subroutine print_rvalue
!
!
!
 Subroutine print_item (inum, ctxt, &
            iout, itab, iwid)
!
! Prints a numbered item
!
! Modules
#if GUI == 1
  Use clrwin$,           Only: winio$
#else
  Use data_io_constants, Only: ioutstd
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: inum ! - item number -
!
  Character(Len=*), Intent(In) :: ctxt ! - text -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
  Integer, Intent(In), Optional :: itab ! - tab position -
  Integer, Intent(In), Optional :: iwid ! - width of number -
!
! Locals
!
! Local scalars
  Integer :: ioutu ! - used output unit number -
!
  Character(Len=16) :: cfmt ! - format statement -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Print item
  If (Present(iout)) Then
     ioutu = iout
  Else
#if GUI == 0
     ioutu = ioutstd
#else
     ioutu = -1
     iw = winio$('%wd.&', inum)
     If (Present(itab)) Then
        Call set_tabs (1, (/itab/))
        iw = winio$('%ta&')
     Else
        Call print_text ('  ', &
             ladvance=.false.)
     End If
     iw = winio$('%`rs&', ctxt)
     Call gui_creturn (.false.)
#endif
  End If
  If (ioutu /= -1) Then
     If (Present(itab)) Then
        If (Present(iwid)) Then
           Write (Unit=cfmt, Fmt='(A,I0,A,I0,A)') '(1X,I', iwid+1, ',A,T', itab, ',A)'
        Else
           Write (Unit=cfmt, Fmt='(A,I0,A,I0,A)') '(1X,I0,A,T', itab, ',A)'
        End If
     Else
        If (Present(iwid)) Then
           Write (Unit=cfmt, Fmt='(A,I0,A)') '(1X,I', iwid+1, ',A,1X,A)'
        Else
           Write (Unit=cfmt, Fmt='(A)') '(1X,I0,A,1X,A)'
        End If
     End If
     Write (Unit=ioutu, Fmt=cfmt) inum, '.', ctxt
  End If
!
  Return
 End Subroutine print_item
!
!
!
 Subroutine print_subtitle (stitle, &
            iout)
!
! Prints a sub-title
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: stitle ! - subtitle -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Print subtitle
  If (Present(iout)) Then
     Call print_text (' ', &
          iout=iout)
     Call print_text (stitle, &
          iout=iout)
  Else
#if GUI == 1
     Call gui_creturn (.true.)
     Call print_text (stitle//':', &
          lbold=.true.)
     Call gui_creturn (.true.)
#else
     Call print_text (' ')
     Call print_text (stitle)
#endif
  End If
!
  Return
 End Subroutine print_subtitle
!
!
!
 Subroutine print_text (ctxt, &
            iout, itab, ladvance, lbold, lcentre, lcolour, lpad, llower, lupper)
!
! Prints text
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
#endif
  Use labels, Only: lowcase, upcase
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctxt ! - text -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iout ! - output unit number -
  Integer, Intent(In), Optional :: itab ! - indent -
!
  Logical, Intent(In), Optional :: ladvance ! - advance to new line? -
  Logical, Intent(In), Optional :: lbold    ! - bold? -
  Logical, Intent(In), Optional :: lcolour  ! - coloured? -
  Logical, Intent(In), Optional :: lcentre  ! - centre? -
  Logical, Intent(In), Optional :: lpad     ! - pad? -
  Logical, Intent(In), Optional :: llower   ! - lower case? -
  Logical, Intent(In), Optional :: lupper   ! - upper case? -
!
! Locals
!
! Local scalars
  Character(Len=8) :: cfmt ! - centring format -
  Character(Len=4) :: cadv ! - advance indicator -
!
  Character(Len=Len(ctxt)) :: c ! - copy of original text message -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len, Present
#if GUI == 1
  Intrinsic Len_Trim
#endif
!
! Executable Statements
!
! Create a copy of the text
  c = ctxt
!
! Convert case
  If (Present(lupper)) Then
     If (lupper) Call upcase (c)
  End If
  If (Present(llower)) Then
     If (llower) Call lowcase (c)
  End If
!
! Print text
  If (Present(ladvance)) Then
     If (ladvance) Then
        cadv = 'yes'
     Else
        cadv = 'no'
     End If
  Else
     cadv = 'yes'
  End If
! - calculate centring format -
  If (Present(iout)) Then
     If (Present(lcentre)) Then
        cfmt = get_format(lcentre, c)
     Else If (Present(itab)) Then
        Write (Unit=cfmt, Fmt='(A,I3,A)') '(',itab,'X,A)'
     Else
        cfmt = '(1X,A)'
     End If
     Write (Unit=iout, Fmt=cfmt, Advance=cadv) c ! - print centred text -
  Else
#if GUI == 1
     If (Present(itab)) Then
        Call set_tabs (1, (/itab/))
        iw = winio$('%ta&')
     End If
     If (Present(lbold)) Then
        If (lbold) iw = winio$('%bf&')
     End If
     If (Present(lcentre)) Then
        If (lcentre) iw = winio$('%cn&')
     End If
     If (Present(lcolour)) Then
        If (lcolour) iw = winio$('%tc&', icolh)
     End If
     If (Present(lpad)) Then
        If (lpad) Then
           iw = winio$('    %ws    &', c)
        Else
           iw = winio$('%ws&', c)
        End If
     Else
        If (Len_Trim(c) > 0 .or. cadv == 'yes') Then
           iw = winio$('%ws&', c)
        Else
           iw = winio$('%*ws&', Int(Len(c),Kind(0)), c)
        End If
     End If
     If (Present(lcolour)) Then
        If (lcolour) iw = winio$('%tc&', 1)
     End If
     If (Present(lbold)) Then
        If (lbold) iw = winio$('%`bf&')
     End If
     If (Trim(cadv) == 'yes') Call gui_creturn (.true.)
#else
     If (Present(lcentre)) Then
        cfmt=get_format(lcentre, c)
     Else If (Present(itab)) Then
        Write (Unit=cfmt, Fmt='(A,I3,A)') '(',itab,'X,A)'
     Else
        cfmt='(1X,A)'
     End If
     Write (Unit=*, Fmt=cfmt, Advance=cadv) c
#endif
  End If
!
  Return
!
 Contains
!
!
  Function get_format(lcentre, ctxt)
!
! Generates format constructor
!
! Modules
  Use screen, Only: iwdth
!
! Function type
  Character(Len=8) :: get_format
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctxt ! - text -
!
  Logical, Intent(In) :: lcentre ! - centre? -
!
! Locals
!
! Local scalars
  Integer :: i ! - centring offset -
!
! Executable Statements
!
! Get format
  If (lcentre) Then
     i = (iwdth - Len(ctxt))/2
     If (i > 0) Then
        Write (get_format, Fmt='(A,I3,A)') '(',i, 'X,A)'
     Else
        get_format='(1X,A)'
     End If
  End If
!
  Return
  End Function get_format
 End Subroutine print_text
!
!
!
 Subroutine print_error (msg1, msg2, msg3, email)
!
! Creates an error message window
!
! Modules
#if GUI == 1
  Use clrwin$,      Only: winio$
#endif
  Use data_version, Only: cpt_email
  Use labels,       Only: cg_error
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: msg1 ! - error message -
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: msg2 ! - error message -
  Character(Len=*), Intent(In), Optional :: msg3 ! - error message -
!
  Logical, Intent(In), Optional :: email ! - print email? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 0
  Intrinsic Len_Trim
#endif
!
! Executable Statements
!
! Create window
  Call init_win (title=Trim(cg_error))
!
! Print error message
#if GUI == 1
  iw = winio$('%si#&')
  Call print_text (Trim(cg_error), &
       lcentre=.true., lbold=.true., lupper=.true.)
  Call gui_creturn (.false.)
  Call print_text (msg1)
  If (Present(msg2)) Call print_text (msg2)
  If (Present(msg3)) Call print_text (msg3)
  If (Present(email)) Then
     Call gui_creturn (.true.)
     Call print_text ('If problems persist, contact '//cpt_email)
  End If
  iw = win_prompt(.false.)
#else
  Call print_text (Trim(cg_error)//': ', &
       ladvance=.false., lbold=.true., lupper=.true.)
  Call print_text (msg1)
  If (Present(msg2) .or. Present(msg3)) Then
     If (Present(msg2)) Call print_text (msg2, &
                             itab=Len_Trim(cg_error)+4)
     If (Present(msg3)) Call print_text (msg3, &
                             itab=Len_Trim(cg_error)+4)
  End If
  If (Present(email)) Call print_text ('If problems persist, contact '//cpt_email)
#endif
!
  Return
 End Subroutine print_error
!
!
!
 Subroutine print_warning (msg1, &
            msg2, msg3, msg4, msg5, nopause, lcursor)
!
! Creates a warning message window
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
#endif
  Use labels, Only: cg_warning
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: msg1 ! - warning message -
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: msg2 ! - warning message -
  Character(Len=*), Intent(In), Optional :: msg3 ! - warning message -
  Character(Len=*), Intent(In), Optional :: msg4 ! - warning message -
  Character(Len=*), Intent(In), Optional :: msg5 ! - warning message -
!
  Logical, Intent(In), Optional :: lcursor ! - cursor waiting? -
  Logical, Intent(In), Optional :: nopause ! - no pause? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 0
  Intrinsic Len_Trim
#endif
!
! Executable Statements
!
! Switch off cursor waiting
  If (Present(lcursor)) Call set_cursor_waiting (0)
!
! Create window
  Call init_win (title=Trim(cg_warning))
!
#if GUI == 1
! Print warning message
  iw = winio$('%si!&')
  Call print_text (Trim(cg_warning), &
       lcentre=.true., lbold=.true., lupper=.true.)
  Call gui_creturn (.false.)
  Call print_text (msg1)
  If (Present(msg2)) Call print_text (msg2)
  If (Present(msg3)) Call print_text (msg3)
  If (Present(msg4)) Call print_text (msg4)
  If (Present(msg5)) Call print_text (msg5)
  If (Present(nopause)) Then
     If (.not.nopause) iw = win_prompt(.false.)
  Else
     iw = win_prompt(.false.)
  End If
#else
!
! Print warning message
  Call print_text (Trim(cg_warning)//': ', &
       ladvance=.false., lbold=.true., lupper=.true.)
  Call print_text (msg1)
  If (Present(msg2) .or. Present(msg3) .or. Present(msg4) .or. Present(msg5)) Then
     If (Present(msg2)) Call print_text (msg2, &
                             itab=Len_Trim(cg_warning)+4)
     If (Present(msg3)) Call print_text (msg3, &
                             itab=Len_Trim(cg_warning)+4)
     If (Present(msg4)) Call print_text (msg4, &
                             itab=Len_Trim(cg_warning)+4)
     If (Present(msg5)) Call print_text (msg5, &
                             itab=Len_Trim(cg_warning)+4)
  End If
  If (Present(nopause)) Then
     If (nopause) Return
  End If
  Call gui_creturn (.true.)
#endif
!
! Set cursor waiting
  If (Present(lcursor)) Call set_cursor_waiting (1)
!
  Return
 End Subroutine print_warning
!
!
!
 Subroutine print_advisory (msg1, &
            msg2)
!
! Creates an advsory message window
!
! Modules
#if GUI == 1
  Use clrwin$, Only: winio$
#endif
  Use labels, Only: cg_advisory
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: msg1 ! - advisory message -
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: msg2 ! - error message -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 0
  Intrinsic Len_Trim
#endif
!
! Executable Statements
!
#if GUI == 1
! Create window
  Call init_win (title=Trim(cg_advisory))
!
! Print advisory message
  iw = winio$('%si*&')
  Call gui_creturn (.false.)
  Call print_text (msg1)
  If (Present(msg2)) Call print_text (msg2)
  iw=win_prompt(.false.)
#else
  Call print_text (Trim(cg_advisory)//': ', &
       ladvance=.false., lbold=.true., lupper=.true.)
  Call print_text (msg1)
  If (Present(msg2)) Then
     Call print_text (msg2, &
          itab=Len_Trim(cg_advisory)+4)
  End If
  Call gui_creturn (.true.)
#endif
!
  Return
 End Subroutine print_advisory
!
!
!
#if GUI == 1
 Subroutine add_graphic_area (img_wd, img_hg, &
            lbox, lmouse, img_handle, iset, ival, cbf, cbf2)
!
! Creates a drawing surface
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: img_wd ! - image width -
  Integer, Intent(In) :: img_hg ! - image height -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: ival ! - value to assign to iset -
!
  Integer(i_handle), Intent(In), Optional :: img_handle ! - image handle -
!
  Logical, Intent(In), Optional :: lbox   ! - box selection? -
  Logical, Intent(In), Optional :: lmouse ! - full mouse input? -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: iset ! - item to be set -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: cbf  ! - first callback function -
  Integer, External, Optional :: cbf2 ! - second callback function -
!
! Locals
!
! Local scalars
  Character(Len=64) :: cfmt ! - format -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Construct format statement
  cfmt = '%'
  If ((Present(iset) .and. Present(ival)) .or. Present(cbf) .or. Present(cbf2)) Then
     cfmt = Trim(cfmt)//'^'
  End If
  If (Present(img_handle)) cfmt = Trim(cfmt)//'`'
  cfmt = Trim(cfmt)//'gr[white,rgb_colours'
  If (Present(lbox)) Then
     If (lbox) cfmt = Trim(cfmt)//',box_selection'
  End If
  If (Present(lmouse)) Then
     If (lmouse) cfmt = Trim(cfmt)//',full_mouse_input'
  End If
  cfmt = Trim(cfmt)//']&'
!
! Create drawing surface
  If (Present(cbf2)) Then
     If (Present(cbf)) Then
        If (Present(iset) .and. Present(ival)) Then
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   '+', '+', 'SET', iset, ival, cbf, cbf2)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   '+', '+', 'SET', iset, ival, cbf, cbf2)
           Endif
        Else ! - no value to set -
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   '+', cbf, cbf2)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   '+', cbf, cbf2)
           Endif
        Endif
     Else ! - no callback function -
        If (Present(iset) .and. Present(ival)) Then
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   '+', 'SET', iset, ival, cbf2)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   '+', 'SET', iset, ival, cbf2)
           Endif
        Else ! - no value to set -
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   cbf2)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   cbf2)
           Endif
        Endif
     Endif
  Else ! - no second callback function -
     If (Present(cbf)) Then
        If (Present(iset) .and. Present(ival)) Then
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   '+', 'SET', iset, ival, cbf)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   '+', 'SET', iset, ival, cbf)
           Endif
        Else ! - no value to set -
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   cbf)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   cbf)
           Endif
        Endif
     Else ! - no callback function -
        If (Present(iset) .and. Present(ival)) Then
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle, &
                   'SET', iset, ival)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg, &
                   'SET', iset, ival)
           Endif
        Else ! - no value to set -
           If (Present(img_handle)) Then
              iw = winio$(Trim(cfmt), img_wd, img_hg, img_handle)
           Else ! - no image handle -
              iw = winio$(Trim(cfmt), img_wd, img_hg)
           Endif
        Endif
     Endif
  Endif
!
  Return
 End Subroutine add_graphic_area
!
!
!
 Subroutine add_bitmap (c_bm)
!
! Adds a bitmap
!
#if FTN95 == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: c_bm ! - bitmap -
!
#if FTN95 == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Add bitmap
  iw = winio$('%bm['//Trim(c_bm)//']&')
!
#endif
  Return
 End Subroutine add_bitmap
!
!
!
 Subroutine add_mdi_frame (iwid, ihgt)
!
! Adds a MDI frame
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iwid ! - width -
  Integer, Intent(In) :: ihgt ! - height -
!
! Executable Statements
!
! Add bitmap
  iw = winio$('%pv%fr&', iwid, ihgt)
!
  Return
 End Subroutine add_mdi_frame
!
!
!
 Subroutine add_prog_bar (iwid, prog)
!
! Adds a progress bar
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iwid ! - width of bar -
!
  Real(Kind=rp), Intent(In) :: prog ! - progress -
!
! Executable Statements
!
! Add progress bar
  iw = winio$('%*br[percentage]&', iwid, prog, icoldom)
!
  Return
 End Subroutine add_prog_bar
!
!
!
 Subroutine add_slider_bar (iwid, x, xmin, xmax, &
            c_min, c_max)
!
! Adds a slider bar
!
! Modules
  Use clrwin$, Only: winio$
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iwid ! - width of bar -
!
  Real(Kind=rp), Intent(In) :: x    ! - value -
  Real(Kind=rp), Intent(In) :: xmin ! - minimum value -
  Real(Kind=rp), Intent(In) :: xmax ! - maximum value -
!
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: c_min ! - text at minimum -
  Character(Len=*), Intent(In), Optional :: c_max ! - text at maximum -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Add slider bar
  If (Present(c_min)) iw = winio$(c_min//' &')
  iw = winio$('%*sl&', iwid, x, xmin, xmax)
  If (Present(c_max)) iw = winio$(' '//c_max//'&')
  iw = winio$('%ff&')
!
  Return
 End Subroutine add_slider_bar
!
!
!
#endif
 Subroutine add_check_button (ctxt, icheck, lbold)
!
! Adds a check button
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
!
#endif
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctxt ! - button text -
!
  Logical, Intent(In) :: lbold ! - bold text? -
!
! Input/output scalars
  Integer, Intent(InOut) :: icheck ! - on / off switch -
!
! Executable Statements
!
! Add check button
#if GUI == 1
  If (lbold) iw = winio$('%bf&')
  iw = winio$('%`rb@&', ctxt, icheck)
  If (lbold) iw = winio$('%`bf&')
#else
  icheck = prompt_yesno(.true., .false., ctxt)
#endif
!
  Return
 End Subroutine add_check_button
!
!
!
 Subroutine set_tabs (ntab, itabs)
!
! Sets windows tabs
!
#if GUI == 1
! Modules
  Use clrwin$, Only: winio$
  Use screen,  Only: dts
!
#endif
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ntab ! - number of tabs -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: itabs ! - tab positions -
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Nint
!
! Executable Statements
!
! Sets windows tabs
  Select Case (ntab)
   Case (1)
     iw = winio$('%1tl&', Nint(itabs(1)*dts))
!     iw = winio$('%*tl&', ntab, Nint(itabs(1)*dts))
   Case (2)
     iw = winio$('%2tl&', Nint(itabs(1)*dts), Nint(itabs(2)*dts))
!     iw = winio$('%*tl&', ntab, Nint(itabs(1)*dts), Nint(itabs(2)*dts))

  End Select
!  iw = winio$('%~*tl&', ntab, Nint(itabs(:)*dts))
!
#endif
  Return
 End Subroutine set_tabs
!
!
!
 Subroutine use_url (url)
!
! Opens a browser at given page (accesses the URL in batch version)
!
! Modules
#if GUI == 1
  Use clrwin$, Only: use_url$
#endif
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: url ! - URL -
!
! Executable Statements
!
#if GUI == 1
! Open a browser at given page
  Call use_url$ (Trim(url))
#else
! Print URL
  Print *, 'URL: ', Trim(url)
#endif
!
  Return
 End Subroutine use_url
!
!
!
 Subroutine set_cursor_waiting (icursor)
!
! Sets cursor
!
#if GUI == 1
! Modules
  Use clrwin$, Only: set_cursor_waiting$
!
#endif
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icursor ! - cursor waiting flag -
!
! Executable Statements
!
! Set cursor state
#if GUI == 1
  Call set_cursor_waiting$ (icursor)
#else
  Continue
#endif
!
  Return
 End Subroutine set_cursor_waiting
!
!
!
#if GUI == 1
#if FTN95 == 1
 Function ic_current()
#else
 Function ic_current() Bind(C)
#endif
!
! Identifies current control
!
! Modules
  Use clrwin$, Only: clearwin_info$
!
! Function type
  Integer :: ic_current
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Int
!
! Executable Statements
!
! Identify current control
  ic = Int(clearwin_info$('CURSOR_WINDOW'))
  ic_current = 2
!
  Return
 End Function ic_current
!
!
!
 Function close_windows()
!
! Closes all child windows
!
! Modules
  Use screen, Only: window_update
!
! Function type
  Integer :: close_windows
!
! Locals
!
! Local scalars
  Integer :: i ! - window index -
!
! Executable Statements
!
! Close existing windows
  Do i = 1, mcw
     If (ich(i) < 0) Then
        ich(i) = 0
        Call window_update (ich(i))
     End If
  End Do
  close_windows = 2
!
  Return
 End Function close_windows
#endif
End Module gui
