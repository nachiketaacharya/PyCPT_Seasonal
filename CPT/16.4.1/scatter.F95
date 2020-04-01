! Author: Simon Mason
Module scatter
!
! Modules
  Use data_gui_constants, Only: i_handle
#if GUI == 1
  Use graphics,           Only: image
#endif
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: regress
#if GUI == 1
  Public :: init_img_scatter, scatter_cv, scatter_ra
#endif
!
! Scalars
!
! Integer scalars
  Integer, Private :: ioutr0 ! - regression plot output unit number -
  Integer, Private :: ioutr1 ! - regression plot output unit number -
#if GUI == 1
  Integer, Private :: iouts0 ! - scatter plot output unit number -
#endif
!
  Integer(Kind=i_handle), Private :: ih_r1  ! - regression equations window handle -
!
#if GUI == 1
! Arrays
!
! Derived type arrays
  Type(image), Dimension(2), Private :: img_cv ! - cross-validated graphs -
  Type(image), Dimension(2), Private :: img_ra ! - retroactive graphs -
  Type(image), Dimension(2), Private :: img_rg ! - regression scatter plot -
!
#endif
Contains
!
!
#if GUI == 1
 Subroutine init_img_scatter (igsz1, igsz2)
!
! Initialises scatter and regression plots
!
! Modules
  Use graphics, Only: init_image
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: igsz1 ! - graphics area dimension -
  Integer, Intent(In) :: igsz2 ! - graphics area dimension -
!
! Executable Statements
!
! Initialise graph sizes
  Call init_image (img_cv(1), igsz1)
  Call init_image (img_cv(2), igsz1)
  Call init_image (img_ra(1), igsz1)
  Call init_image (img_ra(2), igsz1)
  Call init_image (img_rg(1), igsz2)
  Call init_image (img_rg(2), igsz2)
!
  Return
 End Subroutine init_img_scatter
!
!
!
 Function scatter_cv()
!
! Modules
  Use labels, Only: l_cvn
!
! Function type
  Integer :: scatter_cv
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Plot cross-validated graphs
  scatter_cv = scatters(Trim(l_cvn%c), img_cv(:), plot_cvs)
!
  Return
 End Function scatter_cv
!
!
!
 Function scatter_ra()
!
! Modules
  Use labels, Only: l_retro
!
! Function type
  Integer :: scatter_ra
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Plot cross-validated graphs
  scatter_ra = scatters(Trim(l_retro%c), img_ra(:), plot_ras)
!
  Return
 End Function scatter_ra
!
!
!
 Function scatters(ctype, img_rsp, plot_scatter_func)
!
! Plots a residuals plot and a scatter diagram
!
! Modules
  Use labels,      Only: l_axislims, l_axislimsy, l_residual, l_scatter
  Use windows_gui, Only: imgid, &
                         init_window_graphic
!
! Function type
  Integer :: scatters
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctype ! - validation type -
!
! Input/output arrays
  Type(image), Dimension(:), Intent(InOut) :: img_rsp ! - residuals/scatter plots -
!
! Procedure arguments
  Integer, External :: plot_scatter_func ! - plotting function -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Set image titles
  img_rsp(1)%title = Trim(ctype)//' Residuals'
  img_rsp(2)%title = 'Observations and '//Trim(ctype)//' Hindcasts'
!
! Plot graphs
  scatters = init_window_graphic(2,img_rsp(:), imgid, 'Graphs', (/l_residual, l_scatter/), &
             plot_scatter_func, f_title, f_export,                                         &
             cxy='Y', iloc=iouts0, limits_grph=f_limits, &
             llims=(/l_axislimsy, l_axislims/), f_help=help_scatter)
  scatters = 1
!
  Return
 End Function scatters
!
!
!
#endif
 Function regress()
!
! Displays regression results
!
#if GUI == 1
! Modules
  Use labels,      Only: l_lab, &
                         l_axislimsy, l_regress
  Use iofiles,     Only: mxa
  Use windows_gui, Only: imgid, &
                         init_window_graphic_0, init_window_graphic
!
#endif
! Function type
  Integer :: regress
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Set image title
  If (mxa == 1) Then
     img_rg(:)%title = 'Regression Plot'
!
! Plot graphs
     regress = init_window_graphic(1, img_rg(:), imgid, l_regress%c, (/l_lab/), &
               update_regr, f_title, f_export, &
               cxy='Y', iloc=ioutr0, add_win_v=show_regr, limits_grph=f_limits, llims=(/l_axislimsy/))
  Else
     regress = init_window_graphic_0(Trim(l_regress%c), update_regr, show_regr, &
               cxy='Y')
  End If
!
#else
! Print results
  regress = show_regr()
  regress = update_regr()
#endif
  regress = 1
!
  Return
 End Function regress
!
!
!
#if GUI == 1
 Function plot_cvs() &
          Result (plot)
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graphs
  plot = plot_res_cv()
  plot = plot_sct_cv()
!
  Return
 End Function plot_cvs
!
!
!
 Function plot_res_cv() &
          Result (plot)
!
! Modules
  Use data_time_constants, Only: iuseq
  Use time,                Only: Operator(+), &
                                 iseq
  Use labels,              Only: ca_seq
  Use arrays,              Only: res=>rwk2, y, yhat
  Use iofiles,             Only: yfile
  Use settings,            Only: iva, nt, nu
  Use missing,             Only: ktuse
  Use season,              Only: lb
  Use plotting,            Only: plot_line_graph
!
! Function type
  Integer :: plot
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Calculate residuals
  Allocate (res(nu,1))
  res(:,1) = yhat(iva,1:nu) - y(iva,1:nu,lb)
!
! Update graphs
! - added 0 because of FTN95 bug -
  plot = plot_line_graph(img_cv(1), Trim(ca_seq(iuseq(iseq))), 'Hindcast - Observation', .false., nt, 1, res(:,:), &
         d1=yfile%fdate+0, kuse=ktuse(:))
  Deallocate (res)
!
  Return
 End Function plot_res_cv
!
!
!
 Function plot_sct_cv() &
          Result (plot)
!
! Modules
  Use data_cpt_constants, Only: nts
  Use arrays,             Only: y, yhat
  Use settings,           Only: iva, nu
  Use categories,         Only: tobs
  Use season,             Only: lb
  Use plotting,           Only: plot_scatter
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graphs
  plot = plot_scatter(img_cv(2), 'Hindcast', 'Observed', iouts0, nu, nts, yhat(iva,:), y(iva,:,lb), .true., &
         tobs=tobs(iva,:))
!
  Return
 End Function plot_sct_cv
!
!
!
 Function plot_ras() &
          Result (plot)
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graphs
  plot = plot_res_ra()
  plot = plot_sct_ra()
!
  Return
 End Function plot_ras
!
!
!
 Function plot_res_ra() &
          Result (plot)
!
! Modules
  Use data_time_constants, Only: iuseq
  Use time,                Only: Operator(+), &
                                 iseq
  Use labels,              Only: ca_seq
  Use arrays,              Only: res=>rwk2, y, yret
  Use iofiles,             Only: yfile
  Use settings,            Only: iva, ntr, nt1, nu, nur, nu1
  Use missing,             Only: ktuse
  Use season,              Only: lb
  Use plotting,            Only: plot_line_graph
!
! Function type
  Integer :: plot
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Calculate residuals
  Allocate (res(nur,1))
  res(:,1) = yret(iva,1:nur) - y(iva,nu1+1:nu,lb)
!
! Update graphs
! - added 0 because of FTN95 bug -
  plot = plot_line_graph(img_ra(1), Trim(ca_seq(iuseq(iseq))), 'Hindcast - Observation', .false., ntr, 1, res(:,:), &
         d1=yfile%fdate+(nt1-1), kuse=ktuse(nt1:))
  Deallocate (res)
!
  Return
 End Function plot_res_ra
!
!
!
 Function plot_sct_ra() &
          Result (plot)
!
! Modules
  Use data_cpt_constants, Only: nts
  Use arrays,             Only: y, yret
  Use settings,           Only: iva, nu, nur, nu1
  Use categories,         Only: tobs
  Use season,             Only: lb
  Use plotting,           Only: plot_scatter
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graphs
  plot = plot_scatter(img_ra(2), 'Hindcast', 'Observed', iouts0, nur, nts, yret(iva,:), y(iva,nu1+1:nu,lb), .true., &
         tobs=tobs(iva,:))
!
  Return
 End Function plot_sct_ra
!
!
!
 Function plot_reg_scatter() &
          Result (plot)
!
! Modules
  Use data_cpt_constants, Only: nts
  Use arrays,             Only: b0, b, x, y
  Use settings,           Only: iva, nu
  Use season,             Only: lb
  Use axes,               Only: x_axis
  Use plotting,           Only: plot_scatter
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graph
  plot = plot_scatter(img_rg(1), 'Predictor', 'Predictand', ioutr0, nu, nts, x(1,:,1), y(iva,:,lb), .false., &
       b0=b0(iva,1), b1=b(1,iva), xmin=x_axis%amin, xmax=x_axis%amax, ymin=img_rg(1)%ymin, ymax=img_rg(1)%ymax)
!
  Return
 End Function plot_reg_scatter
!
!
!
 Function plot_regr() &
          Result (plot)
!
! Modules
  Use data_cpt_constants, Only: nts
  Use arrays,             Only: b0, b, x, y
  Use settings,           Only: igauss, iva, nu
  Use season,             Only: lb
  Use plotting,           Only: plot_scatter
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graphs
  Select Case (igauss)
   Case (0) ! - no transformaion -
     plot = plot_scatter(img_rg(1), 'Predictor', 'Predictand', ioutr0, nu, nts, x(1,:,1), y(iva,:,lb), .false., &
            b0=b0(iva,1), b1=b(1,iva))
   Case Default
     plot = plot_scatter(img_rg(1), 'Predictor', 'Predictand', ioutr0, nu, nts, x(1,:,1), y(iva,:,lb), .false., &
            b0=b0(iva,1))
  End Select
!
  Return
 End Function plot_regr
!
!
!
#endif
 Function show_regr() &
          Result (show)
!
! Modules
  Use labels,  Only: cg_regreqs
  Use gui,     Only: box_close, box_open, open_clearwin
#if GUI == 1
  Use gui,     Only: get_handle
#endif
  Use iofiles, Only: mxa, mya
!
! Function type
  Integer :: show
!
! Locals
!
! Local scalars
  Integer :: mlmx ! - maximum number of lines -
!
  Logical :: lhs ! - add horizontal scroll bar? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Max, Min, Trim
!
! Executable Statements
!
! Print regression equations
  CALL box_open (.true., .false., &
       title=Trim(cg_regreqs))
! - open Clearwin+ windows -
#if GUI == 1
  ih_r1 = get_handle()
#endif
  If (mxa == 1) Then
     lhs = .false.
     mlmx = 2
  Else
     lhs = .true.
     mlmx = Min(Max(2, mya), 41)
  End If
  If (mya <= mlmx) Then
     CALL open_clearwin (ioutr1, 88, mlmx, .true., &
          ih=ih_r1, lhs=lhs)
  Else
     CALL open_clearwin (ioutr1, 88, mlmx, .true., &
          ih=ih_r1, maxl=mya, lhs=lhs, lvs = .true.)
  End If
  CALL box_close (.true.)
  show = 1
!
  Return
 End Function show_regr
!
!
!
 Function update_regr() &
          Result (update)
!
! Updates Tools ~ Regression results and graphs
!
! Modules
  Use gui,          Only: print_text, set_cursor_waiting
  Use iofiles,      Only: mxa, mya
  Use fields,       Only: update_grid
  Use cpt_output,   Only: update_results
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Construct and add coordinate label for initial point
  CALL set_cursor_waiting (1)
  If (mxa == 1) Then
     CALL update_grid (ioutr0)
!
#if GUI == 1
! Plot graph
     update = plot_regr()
!
#endif
! Update regression equation
     CALL print_text (' ', &
          iout=ioutr1)
     update = print_regr(ioutr1)
  Else If (mya>1) Then
     update = update_results(ioutr1, ih_r1, 1, ' ', .false., print_regs)
  Else
     update = print_regr(ioutr1)
  End If
  CALL set_cursor_waiting (0)
  update = 1
!
  Return
 End Function update_regr
!
!
!
 Function print_regs(iout, ija, i2, istn, ifld) &
          Result (printr)
!
! Prints regression equation
!
! Modules
  Use data_numbers, Only: zero
  Use arrays,       Only: b, b0
  Use iofiles,      Only: mxa, yfile
  Use settings,     Only: iprecf
  Use space,        Only: coor, cstndx, cstndy
  Use fields,       Only: iffx
!
! Function type
  Integer :: printr
!
! Arguments
!
! Input scalaras
  Integer, Intent(In) :: iout ! - output unit number -
  Integer, Intent(In) :: ija  ! - location index -
  Integer, Intent(In) :: i2   ! - second (dummy) index -
  Integer, Intent(In) :: istn ! - station index -
  Integer, Intent(In) :: ifld ! - field index -
!
! Locals
!
! Local scalars
  Integer :: i ! - predictor index -
!
  Character(Len=32) :: cfmt  ! - format statement -
  Character(Len= 3) :: csign ! - plus/minus -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Abs
!
! Executable Statements
!
! Print regression equation
  Select Case (yfile%idstr)
   Case (1)
     Write (Unit=iout, Fmt='(A)', Advance='no') coor
   Case (2, 3)
     Write (Unit=iout, Fmt='(A)', Advance='no') cstndy(istn,ifld)
  End Select
  Write (Unit=iout, Fmt='(A)', Advance='no') ' = '
  Write (Unit=cfmt, Fmt='(A,I1,A)') '(F11.', iprecf, ')'
  Write (Unit=iout, Fmt=cfmt, Advance='no') b0(ija,i2)
  Write (Unit=cfmt, Fmt='(A,I1,A)') '(A,F11.', iprecf, ',2A)'
  Do i = 1, mxa
     If (b(i,ija) < zero) Then
        csign = ' - '
     Else
        csign = ' + '
     End If
     Write (Unit=iout, Fmt=cfmt, Advance='no') csign, Abs(b(i,ija)), ' * ', cstndx(i,iffx)
  End Do
  Write (Unit=iout, Fmt='()', Advance='yes')
  printr = 0
!
  Return
 End Function print_regs
!
!
!
 Function print_regr(iout) &
          Result (printr)
!
! Prints regression equation (RHS only)
!
! Modules
  Use data_numbers, Only: zero
  Use arrays,       Only: b, b0
  Use iofiles,      Only: mxa
  Use settings,     Only: iprecf, iva
  Use space,        Only: cstndx
  Use fields,       Only: iffx
!
! Function type
  Integer :: printr
!
! Arguments
!
! Input scalaras
  Integer, Intent(In) :: iout ! - output unit number -
!
! Locals
!
! Local scalars
  Integer :: i ! - predictor index -
!
  Character(Len=32) :: cfmt ! - format statement -
  Character(Len=32) :: ceqn ! - equation component -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Abs, AdjustL, Trim
!
! Executable Statements
!
! Print regression equation
  Write (Unit=cfmt, Fmt='(A,I1,A)') '(F11.', iprecf, ')'
  Write (Unit=ceqn, Fmt=cfmt) b0(iva,1)
  Write (Unit=iout, Fmt='(A)', Advance='no') Trim(AdjustL(ceqn))
  Do i = 1, mxa
     If (b(i,iva) < zero) Then
        Write (Unit=iout, Fmt='(A)', Advance='no') '  -  '
     Else
        Write (Unit=iout, Fmt='(A)', Advance='no') '  +  '
     End If
     Write (Unit=ceqn, Fmt=cfmt) Abs(b(i,iva))
     Write (Unit=iout, Fmt='(3A)', Advance='no') Trim(AdjustL(ceqn)), ' * ', Trim(AdjustL(cstndx(i,iffx)))
  End Do
  Write (Unit=iout, Fmt='()', Advance='yes')
  printr = 0
!
  Return
 End Function print_regr
!
!
!
#if GUI == 1
 Function f_title() &
          Result (title)
!
! Modules
  Use plotting, Only: get_title
!
! Function type
  Integer :: title
!
! Executable Statements
!
! Prompt for title
  title = replot_scatter(get_title)
!
  Return
 End Function f_title
!
!
!
 Function f_limits() &
          Result (limits)
!
! Modules
  Use plotting, Only: get_limits
!
! Function type
  Integer :: limits
!
! Executable Statements
!
! Prompt for limits
  limits = replot_scatter(get_limits)
!
  Return
 End Function f_limits
!
!
!
 Function replot_scatter(get_custom) &
          Result (replot)
!
! Modules
  Use plotting,    Only: replot_graph
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: replot
!
! Arguments
!
! Procedure arguments
  Interface
    Subroutine get_custom (img,lreset)
     Use graphics, Only: image
     Type(image), Intent(InOut) :: img ! - image -
     Logical, Intent(Out) :: lreset ! - reset? -
    End Subroutine get_custom
  End Interface
!
! Executable Statements
!
! Replot customised scatter diagrams
! - cross-validated residuals diagram -
  If (img_cv(1)%id == imgid) Then
     replot = replot_graph(img_cv(1), get_custom, plot_res_cv)
! - cross-validated scatter diagram -
  Else If (img_cv(2)%id == imgid) Then
     replot = replot_graph(img_cv(2), get_custom, plot_sct_cv)
! - retroactive residuals diagram -
  Else If (img_ra(1)%id == imgid) Then
     replot = replot_graph(img_ra(1), get_custom, plot_res_ra)
! - retroactive scatter diagram -
  Else If (img_ra(2)%id == imgid) Then
     replot = replot_graph(img_ra(2), get_custom, plot_sct_ra)
! - regression scatter diagram -
  Else If (img_rg(1)%id == imgid) Then
     replot = replot_graph(img_rg(1), get_custom, plot_reg_scatter)
  End If
!
  Return
 End Function replot_scatter
!
!
!
 Function f_export() &
          Result (export)
!
! Modules
  Use iofiles,     Only: ny
  Use settings,    Only: iv
  Use graphics,    Only: save_graphic
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: export
!
! Executable Statements
!
! Prompt for montage title
  If (imgid<0) Then
     export = f_export_montage(-imgid)
     Return     
  End If
!
! Save graph
! - cross-validated residuals diagram -
  If (img_cv(1)%id == imgid) Then
     export = save_graphic('Residuals_cv', img_cv(1), &
              imap=iv, mmap=ny)
! - cross-validated scatter diagram -
  Else If (img_cv(2)%id == imgid) Then
     export = save_graphic('Scatter_cv', img_cv(2), &
              imap=iv, mmap=ny)
! - retroactive residuals diagram -
  Else If (img_ra(1)%id == imgid) Then
     export = save_graphic('Residuals_retro', img_ra(1), &
              imap=iv, mmap=ny)
! - retroactive scatter diagram -
  Else If (img_ra(2)%id == imgid) Then
     export = save_graphic('Scatter_retro', img_ra(2), &
              imap=iv, mmap=ny)
! - regression scatter diagram -
  Else If (img_rg(1)%id == imgid) Then
     export = save_graphic('Regress_plot', img_rg(1), &
              imap=iv, mmap=ny)
  End If
!
  Return
 End Function f_export
!
!
!
 Function f_export_montage(imgid) &
          Result (export)
!
! Modules
  Use graphics, Only: export_montage
!
! Function type
  Integer :: export
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imgid ! - image ID -
!
! Executable Statements
!
! Save images
  If (img_cv(1)%id == imgid) Then
     export = export_montage('Residuals_cv', img_cv(1), img_cv(2))
  Else If (img_ra(1)%id == imgid) Then
     export = export_montage('Residuals_ra', img_ra(1), img_ra(2))
  End If
!
  Return
 End Function f_export_montage
!
!
!
 Function help_scatter()
!
! Modules
  Use help,        Only: help_file
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: help_scatter
!
! Executable Statements
!
! Open context-specific help pages
  If ((img_cv(1)%id == imgid) .or. (img_ra(1)%id == imgid)) Then
    Call help_file ('Tools_Validate_Scatter.html')
  End If
  help_scatter = 2
!
  Return
 End Function help_scatter
#endif
End Module scatter
