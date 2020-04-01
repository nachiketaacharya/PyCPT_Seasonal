! Author: Simon Mason
Module validate
!
! Modules
  Use data_numbers,       Only: rp
#if GUI == 1
  Use data_gui_constants, Only: i_handle
  Use graphics,           Only: image
#endif
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
#if GUI == 1
  Public :: init_img_validate
#else
  Public :: roc_cv, roc_ra, save_roc
#endif
  Public :: validate_cv, validate_ra, view_good
!
! Scalars
!
! Integer scalars
  Integer, Private :: ioutv0 ! - validation output unit number -
  Integer, Private :: ioutv1 ! - validation output unit number -
  Integer, Private :: ioutv2 ! - validation output unit number -
  Integer, Private :: ioutg1 ! - goodness index output unit number -
#if GUI == 1
!
  Integer(Kind=i_handle), Private :: ih_g ! - goodness index window handle -
!
! Arrays
!
! Derived type scalars
  Type(image), Dimension(2), Private :: img_cv ! - cross-validated hindcasts and ROC graphs -
  Type(image), Dimension(2), Private :: img_ra ! - retroactive hindcasts and ROC graph -
!
#endif
Contains
!
!
#if GUI == 1
 Subroutine init_img_validate (igsize)
!
! Initialises validation plots
!
! Modules
  Use graphics, Only: init_image
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: igsize ! - graphics area dimension -
!
! Executable Statements
!
! Initialise graph sizes
  Call init_image (img_cv(1), igsize)
  Call init_image (img_cv(2), igsize)
  Call init_image (img_ra(1), igsize)
  Call init_image (img_ra(2), igsize)
!
  Return
 End Subroutine init_img_validate
!
!
!
#endif
 Function validate_cv()
!
! Initialises for cross-validated predictions
!
! Function type
  Integer :: validate_cv
!
! Executable Statements
!
! Validate cross-validated forecasts
#if GUI == 1
  img_cv(1)%title = 'Observations and Cross-Validated Hindcasts'
  img_cv(2)%title = 'Relative Operating Characteristics'
  validate_cv = validates('Cross-validated scores', img_cv(:), plot_cv)
#else
  validate_cv = validates('Cross-validated scores', update_cv)
#endif
!
  Return
 End Function validate_cv
!
!
!
 Function validate_ra()
!
! Initialises for retroactive predictions
!
! Function type
  Integer :: validate_ra
!
! Executable Statements
!
! Validate retroactive forecasts
#if GUI == 1
  img_ra(1)%title='Observations and Retroactive Hindcasts'
  img_ra(2)%title = 'Relative Operating Characteristics'
  validate_ra = validates('Retroactive scores', img_ra(:), plot_ra)
#else
  validate_ra = validates('Retroactive scores', update_ra)
#endif
!
  Return
 End Function validate_ra
!
!
!
#if GUI == 1
 Function validates(ctitle, img_val, plot_validation)
!
! Initialises and creates a validation window
!
! Modules
  Use labels,      Only: l_axislimsy, l_hindplot, l_lab, l_rocdiag
  Use gui,         Only: set_cursor_waiting
  Use settings,    Only: nu, sw_elev
  Use scores,      Only: init_dscores
  Use windows_gui, Only: imgid, &
                         init_window_graphic
!
! Function type
  Integer :: validates
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctitle ! - window title -
!
! Input/output arrays
  Type(image), Dimension(:), Intent(InOut) :: img_val ! - hindcasts plot and ROC diagram -
!
! Procedure arguments
  Integer, External :: plot_validation ! - plotting function -
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Initialise memory
  validates = 1
  If (init_dscores(nu) /= 0) Then
     Call set_cursor_waiting (0)
     Return
  End If
!
! Create results window
  If (sw_elev%on) Then ! - expert mode -
     l_lab%c = ' '
     l_lab%iak = 0
     validates = init_window_graphic(2, img_val(:), imgid, ctitle, (/l_hindplot, l_rocdiag/), &
                 plot_validation, title_val, export_val,                                      &
                 cxy='Y', iloc=ioutv0, add_win_h=show_validations, limits_grph=limits_val, &
                 llims=(/l_axislimsy, l_lab/), f_help=help_valid)
  Else ! - beginner mode -
     validates = init_window_graphic(1, img_val(:), imgid, ctitle,(/l_hindplot/), &
                 plot_validation, title_val, export_val,                          &
                 cxy='Y', iloc=ioutv0, add_win_v=show_validations, limits_grph=limits_val, &
                 llims=(/l_axislimsy/), f_help=help_valid)
  End If
  validates = 1
!
  Return
 End Function validates
!
!
!
#else
 Function validates(title, validations)
!
! Initialises and prints validation scores
!
! Modules
  Use gui,        Only: get_output_unit, init_win
  Use settings,   Only: nu
  Use scores,     Only: init_dscores
  Use cpt_output, Only: print_results
!
! Function type
  Integer :: validates
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: title ! - window title -
!
! Procedure arguments
  Integer, External :: validations ! - validation function -
!
! Executable Statements
!
! Initialise memory
  validates = 1
  If (init_dscores(nu) /= 0) Return
!
! Print validation results
  Call init_win (title=title)
  ioutv0 = get_output_unit()
  validates = show_validations()
  validates = print_results(validations)
!
  Return
 End Function validates
!
!
!
#endif
 Function show_validations()
!
! Formats validation results
!
! Modules
  Use data_cpt_constants, Only: lbmscore, ndcat, ndcon
#if GUI == 1
  Use labels,             Only: l_saveroc
#endif
  Use screen,             Only: iwdth
  Use gui,                Only: boxes_open, open_clearwin
#if GUI == 1
  Use gui,                Only: box_close, print_text
  Use menus,              Only: add_menu
#endif
  Use settings,           Only: sw_elev
  Use scores,             Only: get_width_dscores
!
! Function type
  Integer :: show_validations
!
! Locals
!
! Local scalars
  Integer :: iwid ! - width of window -
  Integer :: ns1  ! - number of lines of scores (1 + number of scores) -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count, Max
!
! Executable Statements
!
#if GUI == 1
! Add extra menu item
  Call add_menu ('pm', 1, 3, l_saveroc, &
       cbf=save_roc)
!
#endif
! Identify number of scores to show
  If (sw_elev%on) Then ! - expert mode -
     ns1 = Max(ndcat, ndcon) + 1
  Else ! - beginner mode -
     ns1 = Max(Count(lbmscore(1:ndcon)), Count(lbmscore(ndcon+1:))) + 1
     Call boxes_open (2, 1, .true., .false.)
  End If
!
! Create validation statistics windows
  If (iwdth >= 768) Then
     iwid = get_width_dscores(itab=3)
  Else
     iwid = get_width_dscores(itab=0)
  End If
! - continuous forecast validation measures -
#if GUI == 1
  Call open_clearwin (ioutv1, iwid, ns1, .true., &
       cbox='Continuous measures')
  Call box_close (.false.) ! - close first row first column box -
  Call print_text (' ', &
       ladvance=.false.)
#else
  Call open_clearwin (ioutv1, iwid, ns1, .true.)
#endif
! - categorical forecast validation measures -
#if GUI == 1
  Call open_clearwin (ioutv2, iwid, ns1, .true., &
       cbox='Categorical measures')
  Call box_close (.false.) ! - close first row second column box -
#else
  Call open_clearwin (ioutv2, iwid, ns1, .true.)
#endif
  show_validations = 1
!
  Return
 End Function show_validations
!
!
!
#if GUI == 1
 Function plot_cv() &
          Result (plot)
!
! Plots graphs for cross-validated predictions
!
! Modules
  Use settings, Only: sw_elev
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Calculate new validation statistics
  plot = update_cv()
!
! Plot cross-validated graphs
  plot = plot_hin_cv()
  If (sw_elev%on) plot = plot_roc_cv()
!
  Return
 End Function plot_cv
!
!
!
 Function plot_hin_cv() &
          Result (plot)
!
! Plots cross-validated ROC graph
!
! Modules
  Use arrays,   Only: y, yhat
  Use iofiles,  Only: yfile
  Use settings, Only: iva, nt, nu
  Use season,   Only: lb
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Plot cross-validated graphs
  plot = plot_hin(img_cv(1), nt, nu, y(iva,1:nu,lb), yfile%fdate, nt, nu, yhat(iva,1:nu), yfile%fdate, 1)
!
  Return
 End Function plot_hin_cv
!
!
!
 Function plot_roc_cv() &
          Result (plot)
!
! Plots graphs for cross-validated predictions
!
! Modules
  Use data_cpt_constants, Only: nts
  Use settings,           Only: nu
  Use scores,             Only: far, hit, roca
  Use plotting,           Only: plot_roc
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Plot cross-validated graphs
  plot = plot_roc(img_cv(2), nu, nts, hit, far, roca)
  plot = 2
!
  Return
 End Function plot_roc_cv
!
!
!
 Function plot_ra() &
          Result (plot)
!
! Plots graphs for retroactive predictions
!
! Modules
  Use settings, Only: sw_elev
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Calculate new validation statistics
  plot = update_ra()
!
! Plot retroactive graphs
  plot = plot_hin_ra()
  If (sw_elev%on) plot = plot_roc_ra()
!
  Return
 End Function plot_ra
!
!
!
 Function plot_hin_ra() &
          Result (plot)
!
! Plots retroactive hindcasts
!
! Modules
  Use time,     Only: Operator(+)
  Use arrays,   Only: y, yret
  Use iofiles,  Only: yfile
  Use settings, Only: iva, nt, ntr, nt1, nu, nur
  Use season,   Only: lb
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Plot retroactive hindcasts
  plot = plot_hin(img_ra(1), nt, nu, y(iva,1:nu,lb), yfile%fdate, ntr, nur, yret(iva,1:nur), yfile%fdate+nt1, nt1+1)
!
  Return
 End Function plot_hin_ra
!
!
!
 Function plot_roc_ra() &
          Result (plot)
!
! Plots retroactive ROC graph
!
! Modules
  Use data_cpt_constants, Only: nts
  Use settings,           Only: nur
  Use scores,             Only: far, hit, roca
  Use plotting,           Only: plot_roc
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Plot retroactive ROC
  plot = plot_roc(img_ra(2), nur, nts, hit, far, roca)
  plot = 2
!
  Return
 End Function plot_roc_ra
!
!
!
 Function plot_hin(img, n1, nu1, y, date1, n2, nu2, y2, date2, i2) &
          Result (plot)
!
! Draws a hindcast plot
!
! Modules
  Use labels,   Only: cg_black_l, cg_green_l, cg_grey_l, cg_red_l
  Use time,     Only: pdate, &
                      nse
  Use gui,      Only: icbw, &
                      set_cursor_waiting
!  Use plotting, Only: plot_series  ! FTN95 BUG
!
! Function type
  Integer :: plot
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n1  ! - number of observations -
  Integer, Intent(In) :: nu1 ! - number of used observationss -
  Integer, Intent(In) :: n2  ! - number of hindcasts -
  Integer, Intent(In) :: nu2 ! - number of used hindcasts -
  Integer, Intent(In) :: i2  ! - offest for hindcasts -
!
  Type(pdate), Intent(In) :: date1 ! - first date of observations -
  Type(pdate), Intent(In) :: date2 ! - first date of hindcasts -
!
! Input/output scalars
  Type(image), Intent(InOut) :: img ! - hindcasts plot -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: y  ! - observed values -
  Real(Kind=rp), Dimension(:), Intent(In) :: y2 ! - forecast values -
!
! Locals
!
! Local scalars
  Character(Len=64) :: cylab ! - y-axis label -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Set appropriate y-axis label
  Call set_cursor_waiting (1)
  Select Case (icbw)
   Case (0)
     cylab = 'Observations ('//Trim(cg_red_l)//') / Hindcasts ('//Trim(cg_green_l)//')'
   Case (1)
     cylab = 'Observations ('//Trim(cg_black_l)//') / Hindcasts ('//Trim(cg_grey_l)//')'
  End Select
!
! Update hindcasts graph
  If (nse == 1) Then
!     plot = plot_series(img, Trim(cylab), n1, nu1, y(:), date1, n2, nu2, y2(:), date2, i2)
     plot = plot_series_local(img, Trim(cylab), n1, nu1, y(:), date1, n2, nu2, y2(:), date2, i2)
  Else
!     plot = plot_series(img, Trim(cylab), n1, nu1, y(:), date1, n2, nu2, y2(:), date2, i2, &
!            nse=nse)
     plot = plot_series_local(img, Trim(cylab), n1, nu1, y(:), date1, n2, nu2, y2(:), date2, i2, &
            nse=nse)
  End If
!
  Return
 End Function plot_hin
!
!
!
! FTN95 BUG - copied from plotting.F95
 Function plot_series_local(img, cylab, n1, nu1, y, d1, n2, nu2, y2, d2, i2, &
          df, nf, mf, yf, yfis, nse) &
          Result (plot)
!
! Draws a time series of observations and hindcasts, and optionally forecasts
!
! Modules
  Use data_numbers,        Only: rp
  Use data_cpt_constants,  Only: nts
  Use data_time_constants, Only: iuseq
  Use maths,               Only: interval
  Use time,                Only: pdate
  Use labels,              Only: ca_seq
  Use iofiles,             Only: yfile
  Use settings,            Only: iva
  Use missing,             Only: kfuse, ktuse
  Use categories,          Only: tobs
  Use plotting,            Only: plot_line_graph
!
! Function type
  Integer :: plot
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n1  ! - number of observations -
  Integer, Intent(In) :: nu1 ! - number of used observations -
  Integer, Intent(In) :: n2  ! - number of hindcasts -
  Integer, Intent(In) :: nu2 ! - number of used hindcasts -
  Integer, Intent(In) :: i2  ! - offset for hindcasts -
!
  Character(Len=*), Intent(In) :: cylab ! - y-axis label -
!
  Type(pdate), Intent(In) :: d1 ! - first date of observations -
  Type(pdate), Intent(In) :: d2 ! - first date of hindcasts -
!
! - optional input scalars
  Integer, Intent(In), Optional :: mf  ! - number of new data variables -
  Integer, Intent(In), Optional :: nf  ! - number of new data -
  Integer, Intent(In), Optional :: nse ! - number of subseasons -
!
  Type(pdate), Intent(In), Optional :: df ! - first date of forecasts -
!
! Input/output scalars
  Type(image), Intent(InOut) :: img ! - image -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: y  ! - observed values -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: y2 ! - hindcast values -
!
! - optional input arrays -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: yf ! - forecasts -
!
  Type(interval), Dimension(:), Intent(In), Optional :: yfis ! - prediction intervals -
!
! Locals
!
! Local scalars
  Character(Len=16) :: cxlab ! - x-axis label -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Reshape, Trim
!
! Executable Statements
!
! Determine appropriate x-axis label
  If (yfile%iseq < 0) Then
     cxlab = ca_seq(iuseq(yfile%iseq))
  Else
     cxlab = ca_seq(iuseq(0))
  End If
!
! Update graph
  If (Present(yf)) Then
     plot = plot_line_graph(img, Trim(cxlab), cylab, .false., n1, 1, Reshape(y(:), (/nu1, 1/)), &
            nse=nse, nts=nts, tobs=tobs(iva,:), d1=d1, kuse=ktuse(:),                           &
            d2=d2, n2=n2, m2=1, ilw2=2, k2use=ktuse(i2:), y2=Reshape(y2(:),(/nu2, 1/)), &
            df=df, nf=nf, mf=mf, icolf=3, kfuse=kfuse(:), yf=yf(:,:), yfis=yfis(:))
  Else
     plot = plot_line_graph(img, Trim(cxlab), cylab, .false., n1, 1, Reshape(y(:), (/nu1, 1/)), &
            nse=nse, nts=nts, tobs=tobs(iva,:), d1=d1, kuse=ktuse(:), &
            d2=d2, n2=n2, m2=1, ilw2=2, k2use=ktuse(i2:), y2=Reshape(y2(:), (/nu2, 1/)))
  End If
  plot = 1
!
  Return
 End Function plot_series_local
#endif
!
!
!
 Function update_cv() &
          Result (update)
!
! Prepares for cross-validation results 
!
! Modules
  Use arrays,     Only: y, yhat
  Use settings,   Only: iva, nu
  Use season,     Only: lb
  Use categories, Only: ifor, iobs, pobs
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Calculate validation statistics
  update = update_validations(nu, y(iva,:,lb), yhat(iva,:), iobs(iva,:), ifor(iva,:), pobs(iva,:))
!
  Return
 End Function update_cv
!
!
!
 Function update_ra() &
          Result (update)
!
! Prepares for retroactive results 
!
! Modules
  Use arrays,     Only: y, yret
  Use settings,   Only: iva, nu, nur, nu1
  Use season,     Only: lb
  Use categories, Only: irfor, irobs, pobs
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Calculate validation statistics
  update = update_validations(nur, y(iva,nu1+1:nu,lb), yret(iva,:), irobs(iva,:), irfor(iva,:), pobs(iva,:))
!
  Return
 End Function update_ra
!
!
!
 Function update_validations(n, y, yhat, iobs, ifor, cps) &
          Result (update)
!
! Calculates and prints results for current location
!
! Modules
  Use data_numbers,   Only: one
  Use progress_meter, Only: end_progress_meter, set_progress_increment, start_progress_meter
  Use fields,         Only: update_grid
  Use scores,         Only: score, &
                            calc_scores, print_dscores
!
! Function type
  Integer :: update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of cases -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iobs ! - observed categories -
  Integer, Dimension(:), Intent(In) :: ifor ! - forecast categories -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: y    ! - observed values -
  Real(Kind=rp), Dimension(:), Intent(In) :: yhat ! - forecast values -
  Real(Kind=rp), Dimension(:), Intent(In) :: cps  ! - climatological probabilities -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Initialise progress meter
  update = 1
#if GUI == 1
  Call start_progress_meter ('Calculating validation statistics', .true.)
#else
  Call start_progress_meter (' ', .true.)
#endif
  Call set_progress_increment (one, .false.)
!
! Construct and add coordinate label for current point
  Call update_grid (ioutv0)
!
! Calculate scores
  Call calc_scores (n, yhat, y, ifor, iobs, cps, score, ifail)
  If (ifail /= 0) Return
!
! Print scores
#if GUI == 1
  Call print_dscores (ioutv1, ioutv2, .false., score)
#else
  Call print_dscores (ioutv1, ioutv2, .true., score)
#endif
!
! Update progress meter
#if GUI == 1
  Call end_progress_meter (.false., .true.)
#else
  Call end_progress_meter (.false., .false.)
#endif
  update = 0
!
  Return
 End Function update_validations
!
!
!
#if GUI == 1
 Function title_val() &
          Result (title)
!
! Updates graphic title
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
  title = replot_val(get_title)
!
  Return
 End Function title_val
!
!
!
 Function limits_val() &
          Result (limits)
!
! Prompts for axis limits
!
! Modules
  Use plotting, Only: get_limits
!
! Function type
  Integer :: limits
!
! Executable Statements
!
! Prompt for y-axis limits
  limits = replot_val(get_limits)
!
  Return
 End Function limits_val
!
!
!
 Function replot_val(get_custom) &
          Result (replot)
!
! Updates axis limits
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
    Subroutine get_custom (img, lreset)
     Use graphics, Only: image
     Type(image), Intent(InOut) :: img ! - image -
     Logical, Intent(Out) :: lreset ! - reset? -
   End Subroutine get_custom
  End Interface
!
! Executable Statements
!
! Prompt for y-axis limits
  If (img_cv(1)%id == imgid) Then
     replot = replot_graph(img_cv(1), get_custom, plot_hin_cv)
  Else If (img_cv(2)%id == imgid) Then
     replot = replot_graph(img_cv(2), get_custom, plot_roc_cv)
  Else If (img_ra(1)%id == imgid) Then
     replot = replot_graph(img_ra(1), get_custom, plot_hin_ra)
  Else If (img_ra(2)%id == imgid) Then
     replot = replot_graph(img_ra(2), get_custom, plot_roc_ra)
  End If
!
  Return
 End Function replot_val
!
!
!
 Function export_val() &
          Result (export)
!
! Saves a plot
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
  If (imgid < 0) Then
     export = export_montage_val(-imgid)
     Return     
  End If
!
! Save graph
  If (img_cv(1)%id == imgid) Then
     export = save_graphic('Hind_cv', img_cv(1), &
            imap=iv, mmap=ny)
  Else If (img_cv(2)%id == imgid) Then
     export = save_graphic('ROC_cv', img_cv(2), &
            imap=iv, mmap=ny)
  Else If (img_ra(1)%id == imgid) Then
     export = save_graphic('Hind_retro', img_ra(1), &
            imap=iv, mmap=ny)
  Else If (img_ra(2)%id == imgid) Then
     export = save_graphic('ROC_retro', img_ra(2), &
            imap=iv, mmap=ny)
  End If
!
  Return
 End Function export_val
!
!
!
 Function export_montage_val(imgid) &
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
     export = export_montage('Cross_validated', img_cv(1), img_cv(2))
  Else If (img_ra(1)%id == imgid) Then
     export = export_montage('Retroactive', img_ra(1), img_ra(2))
  End If
!
  Return
 End Function export_montage_val
!
!
!
 Function help_valid()
!
! Modules
  Use help,        Only: help_file
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: help_valid
!
! Executable Statements
!
! Open context-specific help pages
  If (img_cv(1)%id == imgid .or. img_cv(2)%id == imgid .or. &
      img_ra(1)%id == imgid .or. img_ra(2)%id == imgid) Then
     Call help_file ('CPT_use06a.html')
  End If
  help_valid = 2
!
  Return
 End Function help_valid
#else
!
!
!
 Function roc_cv()
!
! Prepares to print cross-validated ROC
!
! Modules
  Use settings, Only: nu
!
! Function type
  Integer :: roc_cv
!
! Executable Statements
!
! Cross-validated ROC
  roc_cv = print_roc('Cross-validated ROC', nu, update_roc_cv)
!
  Return
 End Function roc_cv
!
!
!
 Function roc_ra()
!
! Prepares to print retroactive ROC
!
! Modules
  Use settings, Only: nur
!
! Function type
  Integer :: roc_ra
!
! Executable Statements
!
! Retroactive ROC
  roc_ra = print_roc('Retroactive ROC', nur, update_roc_ra)
!
  Return
 End Function roc_ra
!
!
!
 Function print_roc(title,n,update_roc)
!
! Prints ROC results
!
! Modules
  Use gui,        Only: init_win
#if GUI == 0
  Use gui,        Only: get_output_unit
#endif
  Use scores,     Only: close_dscores, init_dscores
  Use cpt_output, Only: print_results
!
! Function type
  Integer :: print_roc
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of cases -
!
  Character(Len=*), Intent(In) :: title ! - window title -
!
! Procedure arguments
  Integer, External :: update_roc ! - ROC function -
!
! Executable Statements
!
! Allocate workspace
  print_roc = 1
  If (init_dscores(n) /= 0) Return
#if GUI == 0
  ioutv0 = get_output_unit()
  ioutv1 = get_output_unit()
#endif
!
! Print ROC results
  Call init_win (title=title)
  print_roc = print_results(update_roc)
  print_roc = close_dscores()
!
  Return
 End Function print_roc
!
!
!
 Function update_roc_cv() &
          Result (update)
!
! Prepares for cross-validated ROC
!
! Modules
  Use arrays,     Only: yhat
  Use settings,   Only: iva, nu
  Use categories, Only: iobs
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Calculate ROC
  update = update_roc(nu, iobs(iva,:), yhat(iva,:))
!
  Return
 End Function update_roc_cv
!
!
!
 Function update_roc_ra() &
          Result (update)
!
! Prepares for retroactive ROC
!
! Modules
  Use arrays,     Only: yret
  Use settings,   Only: iva, nur
  Use categories, Only: irobs
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Calculate ROC
  update = update_roc(nur, irobs(iva,:), yret(iva,:))
!
  Return
 End Function update_roc_ra
! 
!
!
 Function update_roc(n, iobs, yhat) &
          Result (update)
!
! Calculates and prints ROC results for current location
!
! Modules 
  Use data_cpt_constants, Only: ng
  Use labels,             Only: ca_cat_a
  Use arrays,             Only: rnkf=>rwk, &
                                rank_data
  Use fields,             Only: update_grid
  Use scores,             Only: far, hit, roca, &
                                roc
!
! Function type
  Integer :: update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of cases -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: iobs ! - observed categories -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: yhat ! - forecasts -
!
! Locals
!
! Local scalars
  Integer :: k ! - case index -
!
! Executable Statements
! 
! Calculate ROC 
  Call rank_data (yhat(:), n, 'a', rnkf)
  Call roc (n, 1, iobs(:), rnkf(:), roca(1), hit(1:n,1), far(1:n,1))
  Call roc (n, ng, iobs(:), rnkf(:), roca(2), hit(1:n,2), far(1:n,2))
!
! Construct and add coordinate label for current point
  Call update_grid (ioutv0)
!
! Print ROC
  Write (Unit=ioutv1, Fmt='(A,F10.3)') ' ROC area ('//ca_cat_a( 1)//'):', roca(1)
  Write (Unit=ioutv1, Fmt='(A,F10.3)') ' ROC area ('//ca_cat_a(ng)//'):', roca(2)
  Write (Unit=ioutv1, Fmt=*)
  Write (Unit=ioutv1, Fmt='(4A)') '            Hit rates ('//ca_cat_a( 1)//')', '    False-alarm rates ('//ca_cat_a( 1)//')', &
                                 '            Hit rates ('//ca_cat_a(ng)//')', '    False-alarm rates ('//ca_cat_a(ng)//')'
  Do k = 1, n
     Write (Unit=ioutv1, Fmt='(4F25.3)') hit(k,1), far(k,1), hit(k,2), far(k,2)
  End Do
  Write (Unit=ioutv1, Fmt=*)
  update = 0
!
  Return
 End Function update_roc
#endif
!
!
!
 Function save_roc()
!
! Saves ROC results
!
! Modules
  Use iofiles,    Only: rofile
  Use cpt_output, Only: save_results
!
! Function type
  Integer :: save_roc
!
! Executable Statements
!
! Save requested results
  save_roc = save_results(rofile)
!
  Return
 End Function save_roc
!
!
!
 Function view_good()
!
! Creates results window for goodness index optimization results
!
! Modules
  Use labels,      Only: l_goodness
  Use errors,      Only: cpt_error
  Use analysis,    Only: lopt
  Use windows_gui, Only: init_window_graphic_0
!
! Function type
  Integer :: view_good
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Open results window
  If (Allocated(lopt)) Then
     view_good = init_window_graphic_0(l_goodness%c, update_good, show_good)
     Return
  End if
!
! Indicate if results are unavailable
  ifail = 1
  Call cpt_error ('view_good', ifail, .false.)
  view_good = 2
!
  Return
 End Function view_good
!
!
!
 Function show_good() &
          Result (show)
!
! Creates ClearWin+ window for goodness index results
!
! Modules
  Use gui, Only: open_clearwin
!
! Function type
  Integer :: show
!
! Executable Statements
!
! Create results window
#if GUI == 1
  Call open_clearwin (ioutg1, 88, 22, .true., &
       ih=ih_g, lvs=.true.)
#else
  Call open_clearwin (ioutg1, 88, 22, .true.)
#endif
  show = 1
!
  Return
 End Function show_good
!
!
!
 Function update_good() &
          Result (update)
!
! Prints goodness index optimization results
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr
  Use data_io_constants,  Only: lprd
  Use labels,             Only: cg_to_l, cg_tperiod
  Use time,               Only: Operator(+), &
                                get_cdate
  Use gui,                Only: print_text
  Use iofiles,            Only: yfile
  Use settings,           Only: iretro, nret, nretro, nt, nt1
  Use scores,             Only: gm
  Use pcs,                Only: mcc, mxe, mye, ncc, nxe, nye
  Use analysis,           Only: ianal, lopt
!
! Function type
  Integer :: update
!
! Locals
!
! Local scalars
  Integer :: ir   ! - retroactice index -
  Integer :: irl  ! - index of last retroactive training-period date -
  Integer :: ixe  ! - X EOF mode index -
  Integer :: iye  ! - Y EOF mode index -
  Integer :: icc  ! - current number of CCA modes -
  Integer :: nxo  ! - optimal number of X EOF modes -
  Integer :: nyo  ! - optimal number of Y EOF modes -
  Integer :: nco  ! - optimal number of CCA modes -
  Integer :: iopt ! - optimization step index -
  Integer :: jopt ! - index of optimized settings -
!
  Character(Len=lprd) :: cprd1 ! - period -
  Character(Len=lprd) :: cprd2 ! - period -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustR, Min, Trim
!
! Executable Statements
!
! Print headers for each retroactice step
  cprd1 = get_cdate(yfile%prd1 + (yfile%it1 - 1), 1, 2)
  irl = nt1
  Do ir = 1, nret*iretro+1
     If (ir == nret*iretro+1) Then
        cprd2 = get_cdate(yfile%prd1 + (yfile%it1 + nt - 2), 1, 2)
     Else
        cprd2 = get_cdate(yfile%prd1 + (yfile%it1 + irl - 2), 1, 2)
     End If
     Call print_text (Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2), &
          iout=ioutg1)
     Call print_text (' ', &
          iout=ioutg1)
     Call print_text ('               CURRENT                         OPTIMUM', &
          iout=ioutg1)
     Call print_text (' ', &
          iout=ioutg1)
     Call print_text ('      Number of Modes    Goodness      Number of Modes    Goodness', &
          iout=ioutg1)
     Select Case (ianal)
      Case (ia_cca)
        Call print_text ('         X    Y   CCA       Index         X    Y   CCA       Index', &
            iout=ioutg1)
      Case (ia_pcr)
        Call print_text ('         X                  Index         X                  Index', &
            iout=ioutg1)
     End Select
     Call print_text (' ', &
          iout=ioutg1)
     iopt = 0
     jopt = 1
!
! Print optimization
     Do ixe = nxe, mxe
        Do iye = nye, mye
           Do icc = ncc, Min(mcc, ixe, iye)
              iopt = iopt + 1
              If (lopt(iopt,ir)) Then
                 If (iopt == 1 .or. gm(iopt,ir) > gm(jopt,ir)) Then
                    jopt = iopt
                    nxo = ixe
                    nyo = iye
                    nco = icc
                 End If
                 Select Case (ianal)
                  Case (ia_cca)
                    Write (Unit=ioutg1, Fmt='(2(5X,2I5,I6,F12.3))') &
                                     ixe, iye, icc, gm(iopt,ir), nxo, nyo, nco, gm(jopt,ir)
                  Case (ia_pcr)
                    Write (Unit=ioutg1, Fmt='(2(5X,I5,11X,F12.3))') &
                                     ixe, gm(iopt,ir), nxo, gm(jopt,ir)
                 End Select
              Else
                 Select Case (ianal)
                  Case (ia_cca)
                    Write (Unit=ioutg1, Fmt='(5X,2I5,I6,A12,5X,2I5,I6,F12.3)') &
                                     ixe, iye, icc, AdjustR('N/A'), nxo, nyo, nco, gm(jopt,ir)
                  Case (ia_pcr)
                    Write (Unit=ioutg1, Fmt='(5X,I5,11X,A12,5X,I5,11X,F12.3)') &
                                     ixe, AdjustR('N/A'), nxo, gm(jopt,ir)
                 End Select
              End If
           End Do
        End Do
     End Do
     If (ir < nret*iretro+1) Then
        Call print_text (' ', &
             iout=ioutg1)
        Call print_text (' ', &
             iout=ioutg1)
     End If
     irl = irl + nretro
  End Do
  update = 1
!
  Return
 End Function update_good
End Module validate
