! Author: Simon Mason
Module maps
!
! Modules
  Use data_numbers, Only: rp
#if GUI == 1
  Use graphics,     Only: image
#endif
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: map_averages, map_correls, map_pscores, map_skill_cv, map_skill_ra, map_thresholds
#if GUI == 1
  Public :: init_img_maps
!
! Scalars
!
! Integer scalars
  Integer, Private :: ioutm ! - correlation maps output unit number -
!
! Arrays
!
! Derived-type arrays
  Type(image), Dimension(1), Private :: img_aves ! - averages -
  Type(image), Dimension(1), Private :: img_cvsk ! - cross-validated skill map -
  Type(image), Dimension(1), Private :: img_rask ! - retroactive skill map -
  Type(image), Dimension(1), Private :: img_prss ! - probabilistic scores map -
  Type(image), Dimension(1), Private :: img_cors ! - correlations -
!
  Type(image), Dimension(2), Private :: img_thrs ! - upper (1) and lower (2) thresholds -
#endif
!
Contains
!
!
#if GUI == 1
 Subroutine init_img_maps (igsize)
!
! Initialises map graphics
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
! Initialise images
  Call init_image (img_thrs(1), (igsize*3)/2)
  Call init_image (img_thrs(2), (igsize*3)/2)
  Call init_image (img_aves(1), (igsize*3)/2)
  Call init_image (img_cors(1), (igsize*3)/2)
  Call init_image (img_cvsk(1), (igsize*3)/2)
  Call init_image (img_rask(1), (igsize*3)/2)
  Call init_image (img_prss(1), (igsize*3)/2)
!
  Return
 End Subroutine init_img_maps
!
!
!
#endif
 Function map_averages() &
          Result (maps)
!
#if GUI == 1
! Creates window for a climatological average map
#else
! Saves climatological average map data
#endif
!
! Modules
#if GUI == 1
  Use labels,      Only: l_lab, &
                         l_averages, l_contourlims
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#else
  Use iofiles,     Only: avfile
  Use cpt_output,  Only: save_results
#endif
!
! Function type
  Integer :: maps
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Reset image
  img_aves(:)%title = Trim(l_averages%c)
!
! Open results window
  maps = init_window_graphic(1, img_aves(:), imgid, l_averages%c, (/l_lab/), plot_aves, title_maps, export_maps, &
         cxy='Y', limits_grph=limits_maps, llims=(/l_contourlims/), f_help=help_maps)
!
#else
! Prompt for output file
  maps = save_results(avfile)
!
#endif
  Return
 End Function map_averages
!
!
!
 Function map_thresholds() &
          Result (maps)
!
#if GUI == 1
! Creates window for a climatological thresholds map
#else
! Saves climatological thresholds map data
#endif
!
! Modules
#if GUI == 1
  Use labels,      Only: l_contourlims, l_threshs, &
                         ls_threshs
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#else
  Use iofiles,     Only: thfile
  Use cpt_output,  Only: save_results
#endif
!
! Function type
  Integer :: maps
!
! Executable Statements
!
#if GUI == 1
! Open results window
  maps = init_window_graphic(2, img_thrs(:), imgid, l_threshs%c, (/ls_threshs(2), ls_threshs(1)/), &
         plot_thresholds, title_maps, export_maps,                                                 &
         cxy='Y', limits_grph=limits_maps, llims=(/l_contourlims, l_contourlims/), f_help=help_maps)
!
#else
! Prompt for output file
  maps = save_results(thfile)
!
#endif
  Return
 End Function map_thresholds
!
!
!
 Function map_correls() &
          Result (maps)
!
! Creates window for a climatological correlations map
!
! Modules
#if GUI == 1
  Use labels,             Only: l_lab, &
                                l_correls
#else
  Use data_cpt_constants, Only: ncorr
  Use labels,             Only: ca_corrs, &
                                l_correls
  Use gui,                Only: box_close, ganged_rbs, init_win
#endif
  Use iofiles,            Only: mxa
#if GUI == 1
  Use scores,             Only: init_correls
  Use windows_gui,        Only: imgid, &
                                init_window_graphic
#else
  Use iofiles,            Only: cofile
  Use scores,             Only: icorrs, &
                                init_correls
  Use cpt_output,         Only: save_results
#endif
!
! Function type
  Integer :: maps
!
#if GUI == 0
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
#endif
! Executable Statements
!
! Allocate memory
  maps = 1
  If (init_correls(mxa) /= 0) Return
!
! Open results window
#if GUI == 1
  maps = init_window_graphic(1, img_cors(:), imgid, l_correls%c, (/l_lab/), update_correls, title_maps, export_maps, &
         cxy='Y', iloc=ioutm, add_win_h=add_correls, f_help=help_maps)
#else
  Call init_win (title=Trim(l_correls%c))
  icorrs(:) = 0
  Call ganged_rbs (ncorr, icorrs(:), ca_corrs)
! - close scores box -
  Call box_close (.true.)
!
! Plot correlation maps
  maps = update_correls()
!
! Construct and add coordinate label for current point
  maps = save_results(cofile)
#endif
!
  Return
 End Function map_correls
!
!
!
#if GUI == 1
 Function add_correls() &
          Result (adds)
!
! Adds correlations to results window
!
! Modules
  Use data_cpt_constants, Only: ncorr
  Use labels,             Only: ca_corrs, &
                                l_scores
  Use gui,                Only: box_open, box_close,ganged_rbs, print_text
  Use fields,             Only: prompt_fields
  Use scores,             Only: icorrs
!
! Function type
  Integer :: adds
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Prompt for score
  Call box_open (.true., .false., &
       title=Trim(l_scores%c))
  Call ganged_rbs (ncorr, icorrs(:), &
       crbs=ca_corrs, frb=update_correls)
! - close scores box -
  Call box_close (.true.)
!
! Prompt for X field and lagged fields
  Call print_text ('  ', &
       ladvance=.false.)
  Call prompt_fields ('X', update_correls)
  adds = 2
!
  Return
 End Function add_correls
!
!
!
#endif
 Function map_skill_cv() &
          Result (maps)
!
! Function type
  Integer :: maps
!
! Executable Statements
!
! Update maps
#if GUI == 1
  maps = map_skill('Cross-validated results', img_cvsk(:), update_skill_cv, prompt_skill_cv, select_skill_cv)
#else
  maps = map_skill('Cross-validated results', update_skill_cv, prompt_skill_cv)
#endif
!
  Return
 End Function map_skill_cv
!
!
!
 Function map_skill_ra() &
          Result (maps)
!
! Function type
  Integer :: maps
!
! Executable Statements
!
! Update maps 
#if GUI == 1
  maps = map_skill('Retroactive results', img_rask(:), update_skill_ra, prompt_skill_ra, select_skill_ra)
#else
  maps = map_skill('Retroactive results', update_skill_ra, prompt_skill_ra)
#endif
!
  Return
 End Function map_skill_ra
!
!
!
#if GUI == 1
 Function map_skill(title, img_sk, update_skill, prompt_skill, select_skill) &
          Result (maps)
#else
 Function map_skill(title, update_skill, prompt_skill) &
          Result (maps)
#endif
!
! Modules
  Use scores,    Only: init_skill
  Use bootstrap, Only: ipval
!
! Function type
  Integer :: maps
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: title ! - window title -
!
#if GUI == 1
! Input/output arrayss
  Type(image), Dimension(:), Intent(InOut) :: img_sk ! - skill map -
!
! Procedure arguments
  Integer, External :: select_skill
#endif
  Integer, External :: prompt_skill
  Integer, External :: update_skill
!
! Executable Statements
!
! Map skill scores
#if GUI == 1
  maps = map_scores(title, img_sk(:), ipval, &
         init_skill, update_skill, prompt_skill, save_skill, select_skill, &
         save_pvalues=save_pvalues)
#else
  maps = map_scores(title, ipval, init_skill, update_skill, prompt_skill)
#endif
!
  Return
 End Function map_skill
!
!
!
 Function map_pscores() &
          Result (maps)
!
! Plots probabilistic scores maps
!
! Modules
  Use labels,    Only: l_scoresp
  Use scores,    Only: init_pscores
  Use bootstrap, Only: ibcl
!
! Function type
  Integer :: maps
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Map probabilistic scores
#if GUI == 1
  maps = map_scores(Trim(l_scoresp%c), img_prss(:), ibcl, &
         init_pscores, update_pscores, prompt_pscores, save_pscores, select_pscores, &
         extra_prompt=prompt_cat)
#else
  maps = map_scores(Trim(l_scoresp%c), ibcl, init_pscores, update_pscores, prompt_pscores)
  maps = 1
#endif
!
  Return
 End Function map_pscores
!
!
!
#if GUI == 1
 Function map_scores(title, img, isig, init_scores, update_scores, prompt_scores, save_scores, select_scores, &
          save_pvalues, extra_prompt) &
          Result (maps)
#else
 Function map_scores(title, isig, init_scores, update_scores, prompt_scores) &
          Result (maps)
#endif
!
! Plots scores maps
!
! Modules
#if GUI == 1
  Use labels,      Only: l_lab
#endif
  Use iofiles,     Only: mya
  Use settings,    Only: nu, nur
#if GUI == 1
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#endif
!
! Function type
  Integer :: maps
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: title ! - window title -
!
! Input/output scalars
  Integer, Intent(InOut) :: isig ! - significance test -
!
#if GUI == 1
! Input/output arrays
  Type(image), Dimension(:), Intent(InOut) :: img ! - image -
!
#endif
! Procedure arguments
  Integer, External :: init_scores
  Integer, External :: prompt_scores
#if GUI == 1
  Integer, External :: save_scores
  Integer, External :: select_scores
#endif
  Integer, External :: update_scores
!
#if GUI == 1
! - optional procedure arguments -
  Integer, External, Optional :: save_pvalues
  Integer, External, Optional :: extra_prompt
!
#endif
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Max
!
! Executable Statements
!
! Allocate memory
  If (init_scores(mya, Max(nu, nur), isig) /= 0) Return
!
#if GUI == 1
! Open results window
  maps = init_window_graphic(1, img(:), imgid, title, (/l_lab/), update_scores, title_maps, export_maps, &
         select_grph=select_scores, save_rslt=save_scores, save_pvals=save_pvalues, ipv=isig, cxy='Y', &
         add_win_h=prompt_scores, extra_prompt=extra_prompt, f_help=help_maps)
!
#else
  maps = prompt_scores()
  maps = update_scores()
#endif
  Return
 End Function map_scores
!
!
!
 Function prompt_cat() &
          Result (prompt)
!
! Prompts for category
!
! Modules
#if GUI == 0
  Use data_cpt_constants, Only: ng
#endif
  Use gui,                Only: prompt_integer
  Use categories,         Only: icatc, &
                                change_cat
!
! Function type
  Integer :: prompt
!
! Executable Statements
!
! Prompt for mode
#if GUI == 1
  Call prompt_integer ('Category', icatc, &
       itab=17, fcb1=change_cat, fcb2=update_pscores)
#else
  Call prompt_integer ('Category', icatc, &
       ilow=1, ihgh=ng, itab=17, fcb1=change_cat)
#endif
  prompt = 0
!
  Return
 End Function prompt_cat
!
!
!
 Function update_correls() &
          Result (update)
!
! Updates a climatological correlations map for a new location
!
! Modules
#if GUI == 1
  Use data_numbers,       Only: one
#endif
  Use data_cpt_constants, Only: ncorr
  Use labels,             Only: ca_corrs, &
                                cg_done
#if GUI == 1
  Use labels,             Only: cg_field, cg_lagfield
#endif
  Use gui,                Only: set_cursor_waiting
  Use arrays,             Only: x, y
  Use iofiles,            Only: mxa
#if GUI == 1
  Use iofiles,            Only: yfile
#endif
  Use settings,           Only: iva, nu
#if GUI == 1
  Use fields,             Only: ilfy, ifdy, &
                                update_grid
#endif
  Use scores,             Only: corrs, icorr, icorrs, &
                                get_corrs
  Use season,             Only: lb
!
! Function type
  Integer :: update
!
#if GUI == 1
! Locals
!
! Local scalars
  Integer :: i ! - index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
#endif
  Intrinsic Trim
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Identify metric
  Do icorr = 1, ncorr
     If (icorrs(icorr) == 1) Exit
  End Do
!
! Calculate correlations
  Print *, 'Calculating '//Trim(ca_corrs(icorr))//' ...'
  Call get_corrs (icorr, nu, mxa, x(1:mxa,:,1), y(iva,:,lb), corrs(:))
  Print *, Trim(cg_done)//'!'
!
#if GUI == 1
! Construct and add coordinate label for current point
  Call update_grid (ioutm)
!
! Set contours
  img_cors(:)%lcontour = .true.
  img_cors(:)%ncs = 13
  Do i = 1, img_cors(1)%ncs
     img_cors(:)%contour(i) = -0.90_rp + 0.15_rp*Real(i-1, Kind=rp)
  End Do
  img_cors(:)%ymin = -one
  img_cors(:)%ymax = one
  If (yfile%nfs > 1) Then
     If (yfile%nls > 1) Then
        Write (img_cors(1)%title, Fmt='(4A,I0,3A,I0)') Trim(ca_corrs(icorr)), '; ', Trim(cg_lagfield), ': ', ifdy, ', ', &
                                                       Trim(cg_lagfield), ': ', ilfy
     Else
        Write (img_cors(1)%title, Fmt='(4A,I0)') Trim(ca_corrs(icorr)), '; ', Trim(cg_field), ': ', ifdy
     End If
  Else If (yfile%nls > 1) Then
     Write (img_cors(1)%title, Fmt='(4A,I0)') Trim(ca_corrs(icorr)), '; ', Trim(cg_lagfield), ': ', ilfy
  Else
     img_cors(1)%title = ca_corrs(icorr)
  End If
  img_cors(:)%scale = one
!
! Plot correlations map
  update = plot_correl()
#else
  update = 0
#endif
!
  Return
 End Function update_correls
!
!
!
 Function update_skill_cv() &
          Result (update)
!
! Modules
  Use arrays,     Only: y, yhat
  Use iofiles,    Only: mya
  Use settings,   Only: nu
  Use season,     Only: lb
  Use categories, Only: ifor, iobs, pobs
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Plot skill
#if GUI == 1
  update = update_skills(img_cvsk(1), mya, nu, yhat, y(:,:,lb), ifor, iobs, pobs, plot_skill_cv)
#else
  update = update_skills(mya, nu, yhat, y(:,:,lb), ifor, iobs, pobs)
#endif
!
  Return
 End Function update_skill_cv
!
!
!
 Function update_skill_ra() &
          Result (update)
!
! Modules
  Use arrays,     Only: y, yret
  Use iofiles,    Only: mya
  Use settings,   Only: nu, nur, nu1
  Use season,     Only: lb
  Use categories, Only: irfor, irobs, pobs
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Plot skill
#if GUI == 1
  update = update_skills(img_rask(1), mya, nur, yret, y(1:mya,nu1+1:nu,lb), irfor, irobs, pobs, plot_skill_ra)
#else
  update = update_skills(mya, nur, yret, y(1:mya,nu1+1:nu,lb), irfor, irobs, pobs)
#endif
!
  Return
 End Function update_skill_ra
!
!
!
#if GUI == 1
 Function update_skills(img_sk, m, n, yhat, y, ifor, iobs, pobs, plot_skill) &
          Result (update)
#else
 Function update_skills(m, n, yhat, y, ifor, iobs, pobs) &
          Result (update)
#endif
!
! Modules
#if GUI == 1
  Use data_numbers, Only: one
#endif
  Use labels,       Only: ca_skill_t, &
                          cg_done
  Use gui,          Only: set_cursor_waiting
  Use arrays,       Only: get_flag
  Use iofiles,      Only: mya
  Use settings,     Only: nu, nur
  Use scores,       Only: iskill, iskills, pvalues, skills, &
                          get_scores, init_skill
#if GUI == 1
  Use scores,       Only: set_dscore_range
#endif
  Use bootstrap,    Only: ipval, &
                          get_pvalues
!
! Function type
  Integer :: update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: m ! - number of variables -
  Integer, Intent(In) :: n ! - number of cases -
!
#if GUI == 1
! Input/output scalars
  Type(image), Intent(InOut) :: img_sk ! - skill map -
!
#endif
! Input arrays
  Integer, Dimension(:,:), Intent(In) :: ifor ! - estimated categories -
  Integer, Dimension(:,:), Intent(In) :: iobs ! - observed categories -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: y    ! - observed values -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: yhat ! - estimated values -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: pobs ! - climatological probabilities -
!
#if GUI == 1
! Procedure arguments
  Integer, External :: plot_skill ! - skill map -
!
#endif
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated, Max, Trim
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Identify metric
  iskill = get_flag(iskills)
!
! Calculate new skill scores
  Print *, 'Calculating ', Trim(ca_skill_t(iskill)), ' ...'
  Call get_scores (iskill, n, m, yhat(:,:), y(:,:), ifor(:,:), iobs(:,:), pobs(:,:), skills(:))
  If (ipval == 1) Then
     If (.not.Allocated(pvalues)) Then
        If (init_skill(mya, Max(nu, nur), ipval) /= 0) Return
     End If
     Call get_pvalues (iskill, n, m, yhat(:,:), y(:,:), ifor(:,:), iobs(:,:), pobs(:,:), skills(:), pvalues(:))
  End If
  Print *, Trim(cg_done)//'!'
!
#if GUI == 1
! Set contours
  img_sk%lcontour = .true.
  img_sk%scale = one
  Call set_dscore_range (iskill, img_sk%ymin, img_sk%ymax, img_sk%ncs, img_sk%contour)
!
! Plot scores
  Call init_score_plot (img_sk, ca_skill_t(iskill))
  update = plot_skill()
#else
  update = save_skill()
#endif
!
  Return
 End Function update_skills
!
!
!
 Function update_pscores() &
          Result (update)
!
! Modules
#if GUI == 1
  Use data_numbers,       Only: one
#endif
  Use data_cpt_constants, Only: npall
  Use labels,             Only: ca_pvs_t, &
                                cg_done
  Use gui,                Only: set_cursor_waiting
  Use arrays,             Only: get_flag
  Use iofiles,            Only: mya
  Use settings,           Only: ibs, nb, nur
  Use fields,             Only: yfield
  Use scores,             Only: ips, ipss, pscores, &
                                get_pscores
#if GUI == 1
  Use scores,             Only: set_pscore_range
#endif
  Use categories,         Only: icatc, irobs, pobs, rfps
!
! Function type
  Integer :: update
!
! Locals
!
! Local scalars
  Integer :: ipsp = 0
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Identify metric
  ips = get_flag(ipss)
  If (ips == ipsp .and. ips <= npall) Then
#if GUI == 1
     Call set_cursor_waiting (0)
     update = 2
     Return
  End If
  ipsp = ips
#endif
!
! Calculate new probabilistic scores
  Print *, 'Calculating ', Trim(ca_pvs_t(ips)), ' ...'
  Call get_pscores (ips, icatc, nur, mya, nb, ibs, irobs(:,:), rfps(:,:,:), pobs(:,:), pscores(:), yfield(1)%rmiss)
  Print *, Trim(cg_done)//'!'
!
#if GUI == 1
! Set contours
  img_prss%lcontour = .true.
  img_prss%scale = one
  Call set_pscore_range (ips, img_prss(1)%ymin, img_prss(1)%ymax, img_prss(1)%ncs, img_prss(1)%contour)
!
! Plot scores
  If (ips > npall) Then
     Call init_score_plot (img_prss(1), ca_pvs_t(ips), &
          icat=icatc)
  Else
     Call init_score_plot (img_prss(1), ca_pvs_t(ips))
  End If
  update = plot_pscore()
#else
  End If
  update = save_pscores()
  If (ips <= npall) Then
     ipsp = ips
     update = save_pscores()
  End If
#endif
!
  Return
 End Function update_pscores
!
!
!
#if GUI == 1
 Subroutine init_score_plot (img, ca_scores_t, &
            icat)
!
! Sets a default image title
!
! Modules
  Use labels,  Only: cg_field, cg_lagfield, &
                     ls_cats
  Use iofiles, Only: yfile
  Use fields,  Only: ilfy, ifdy
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ca_scores_t ! - score name -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: icat ! - category identifier -
!
! Input/output scalars
  Type(image), Intent(InOut) :: img ! - image -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Set title
  If (yfile%nfs > 1) Then
     If (yfile%nls > 1) Then
        Write (img%title, Fmt='(4A,I0,3A,I0)') Trim(ca_scores_t), '; ', Trim(cg_lagfield), ': ', ifdy, ', ', &
                                               Trim(cg_lagfield), ': ', ilfy
     Else
        Write (img%title, Fmt='(4A,I0)') Trim(ca_scores_t), '; ', Trim(cg_field), ':', ifdy
     End If
  Else If (yfile%nls > 1) Then
     Write (img%title, Fmt='(4A,I0)') Trim(ca_scores_t), '; ', Trim(cg_lagfield), ': ', ilfy
  Else
     img%title = ca_scores_t
  End If
!
! Add category
  If (Present(icat)) img%title = Trim(img%title)//' ('//Trim(ls_cats(icat)%c)//')'
!
  Return
 End Subroutine init_score_plot
!
!
!
 Function plot_aves() &
          Result (plot)
!
! Plots a climatological average map
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use labels,             Only: l_averages
  Use climate,            Only: ave
  Use season,             Only: lb
  Use plotting,           Only: plot_ydata
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
! Plot averages map
  plot = plot_ydata(img_aves(1), icol_sk, 0, 1, ave(:,lb), Trim(l_averages%c))
!
  Return
 End Function plot_aves
!
!
!
 Function plot_thresholds() &
          Result (plot)
!
! Plots a climatological thresholds map
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use labels,             Only: cg_threshslow, cg_threshsupp
  Use categories,         Only: tobs
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer :: plot
!
! Locals
!
! Local scalars
  Real(Kind=rp) :: vlow ! - minimum value -
  Real(Kind=rp) :: vhgh ! - maximum value -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal, MinVal, Trim
!
! Executable Statements
!
! Calculate data limits
  vlow = MinVal(tobs(:,1))
  vhgh = MaxVal(tobs(:,2))
!
! Set titles
  img_thrs(1)%title = Trim(cg_threshsupp)
  img_thrs(2)%title = Trim(cg_threshslow)
!
! Plot thresholds maps
  plot = plot_ydata(img_thrs(2), icol_sk, 0, 1, tobs(:,1), Trim(img_thrs(2)%title), &
         vlow=vlow, vhgh=vhgh)
  plot = plot_ydata(img_thrs(1), icol_sk, 0, 1, tobs(:,2), Trim(img_thrs(1)%title), &
         vlow=vlow, vhgh=vhgh)
!
  Return
 End Function plot_thresholds
!
!
!
 Function plot_tlow() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use labels,             Only: cg_threshslow
  Use categories,         Only: tobs
  Use plotting,           Only: plot_ydata
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
! Plot lower thresholds map
  plot = plot_ydata(img_thrs(2), icol_sk, 0, 1, tobs(:,1), Trim(cg_threshslow))
!
  Return
 End Function plot_tlow
!
!
!
 Function plot_tupp() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use labels,             Only: cg_threshsupp
  Use categories,         Only: tobs
  Use plotting,           Only: plot_ydata
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
! Replot upper thresholds map
  plot = plot_ydata(img_thrs(1), icol_sk, 0, 1, tobs(:,2), Trim(cg_threshsupp))
!
  Return
 End Function plot_tupp
!
!
!
 Function plot_skill_cv() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use scores,             Only: skills
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update skill map
  plot = plot_ydata(img_cvsk(1), icol_sk, 0, 1, skills, 'Score')
!
  Return
 End Function plot_skill_cv
!
!
!
 Function plot_skill_ra() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use scores,             Only: skills
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update skill map
  plot = plot_ydata(img_rask(1), icol_sk, 0, 1, skills, 'Score')
!
  Return
 End Function plot_skill_ra
!
!
!
 Function plot_pscore() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use scores,             Only: pscores
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update skill map
  plot = plot_ydata(img_prss(1), icol_sk, 0, 1, pscores, 'Score')
!
  Return
 End Function plot_pscore
!
!
!
 Function plot_correl() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use labels,             Only: cg_correl
  Use scores,             Only: corrs
  Use plotting,           Only: plot_xdata
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
! Update correlation map
  plot = plot_xdata(img_cors(1), icol_sk, 0, 1, corrs, Trim(cg_correl))
!
  Return
 End Function plot_correl
#endif
!
!
!
 Function save_skill() &
          Result (saves)
!
! Modules
  Use iofiles,    Only: skfile
  Use cpt_output, Only: save_results
!
! Function type
  Integer :: saves
!
! Executable Statements
!
! Save requested results
  saves = save_results(skfile)
!
  Return
 End Function save_skill
!
!
!
 Function save_pscores() &
          Result (saves)
!
! Modules
  Use iofiles,    Only: psfile
  Use cpt_output, Only: save_results
!
! Function type
  Integer :: saves
!
! Executable Statements
!
! Save requested results
  saves = save_results(psfile)
!
  Return
 End Function save_pscores
!
!
!
 Function save_pvalues() &
          Result (saves)
!
! Modules
  Use iofiles,    Only: pvfile
  Use cpt_output, Only: save_results
!
! Function type
  Integer :: saves
!
! Executable Statements
!
! Save requested results
  saves = save_results(pvfile)
!
  Return
 End Function save_pvalues
!
!
!
#if GUI == 1
 Function select_skill_cv() &
          Result (selects)
!
! Modules
  Use graphics, Only: select_location
  Use validate, Only: validate_cv
!
! Function type
  Integer :: selects
!
! Executable Statements
!
! Get mouse position
  selects = select_location(img_cvsk(1), validate_cv)
!
  Return
 End Function select_skill_cv
!
!
!
 Function select_skill_ra() &
          Result (selects)
!
! Modules
  Use graphics, Only: select_location
  Use validate, Only: validate_ra
!
! Function type
  Integer :: selects
!
! Executable Statements
!
! Get mouse position
  selects = select_location(img_rask(1), validate_ra)
!
  Return
 End Function select_skill_ra
!
!
!
 Function select_pscores() &
          Result (selects)
!
! Modules
  Use graphics, Only: select_location
  Use validate, Only: validate_ra
!
! Function type
  Integer :: selects
!
! Executable Statements
!
! Get mouse position
  selects = select_location(img_prss(1), validate_ra)
!
  Return
 End Function select_pscores
!
!
!
#endif
 Function prompt_skill_cv() &
          Result (prompt)
!
! Prompts for scores
!
! Function type
  Integer :: prompt
!
! Executable Statements
!
! Prompt for score
  prompt = prompt_skill(update_skill_cv)
!
  Return
 End Function prompt_skill_cv
!
!
!
 Function prompt_skill_ra() &
          Result (prompt)
!
! Prompts for scores
!
! Function type
  Integer :: prompt
!
! Executable Statements
!
! Prompt for score
  prompt = prompt_skill(update_skill_ra)
!
  Return
 End Function prompt_skill_ra
!
!
!
 Function prompt_skill(update_skill) &
          Result (prompt)
!
! Prompts for scores
!
! Modules
  Use data_cpt_constants, Only: nskill
  Use labels,             Only: ca_skill_t, &
                                l_scores
  Use gui,                Only: prompt_bullet
  Use scores,             Only: iskills
!
! Function type
  Integer :: prompt
!
! Arguments
!
! Procedure arguments
  Integer, External :: update_skill
!
! Executable Statements
!
! Prompt for score
  prompt = prompt_bullet(l_scores%c, nskill, iskills, ca_skill_t, update_skill)
!
  Return
 End Function prompt_skill
!
!
!
 Function prompt_pscores() &
          Result (prompt)
!
! Prompts for scores
!
! Modules
  Use data_cpt_constants, Only: ipsi, npall, npcat
#if GUI == 1
  Use data_cpt_constants, Only: nps
#endif
  Use labels,             Only: ca_pavscores_t, ca_pcvscores_t
  Use gui,                Only: prompt_bullet
#if GUI == 1
  Use gui,                Only: box_close, boxes_open, ganged_rbs, print_text
#else
  Use gui,                Only: prompt_integer
#endif
  Use scores,             Only: ipss
!
! Function type
  Integer :: prompt
!
#if GUI == 0
! Locals
!
! Local scalars
  Integer :: i     ! - selected score -
  Integer :: ifail ! - error indicator -
#endif
! Executable Statements
!
#if GUI == 1
! Format
  Call boxes_open (1, 2, .false., .false.)
!
#endif
! Prompt for scores
! - prompt for all-category scores -
  prompt = prompt_bullet('All-category scores', npall, ipss(1:npall), ca_pavscores_t, update_pscores, &
           lprompt=.false., ns=ipsi(1:npall))
#if GUI == 1
  Call box_close (.false.)
  Call print_text (' ')
#endif
! - prompt for category-specific scores -
  prompt = prompt_bullet('Category-specific scores', npcat, ipss(npall+1:), ca_pcvscores_t, update_pscores, &
           lprompt=.false., ns=ipsi(npall+1:))
#if GUI == 0
  Call prompt_integer ('Score', i, ifail, &
       ns=ipsi(:), lzero=.true.)
  If (i > 0) ipss(i) = 1
! - prompt for category -
  If (i > npall) prompt = prompt_cat()
#else
  Call box_close (.false.)
!
! Gang the two sets of scores
  Call ganged_rbs (nps, ipss(:))
#endif
!
  Return
 End Function prompt_pscores
!
!
!
#if GUI == 1
 Function title_maps() &
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
  title = replot_maps(get_title)
!
  Return
 End Function title_maps
!
!
!
 Function limits_maps() &
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
  limits = replot_maps(get_limits)
!
  Return
 End Function limits_maps
!
!
!
 Function replot_maps(get_custom) &
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
    Subroutine get_custom (img, lreset)
     Use graphics, Only: image
     Type(image), Intent(InOut) :: img ! - image -
     Logical, Intent(Out) :: lreset ! - reset? -
    End Subroutine get_custom
  End Interface
!
! Executable Statements
!
! Replot customised maps
! - climatological maps -
  If (img_aves(1)%id == imgid) Then
     replot = replot_graph(img_aves(1), get_custom, plot_aves)
  Else If (img_thrs(1)%id == imgid) Then
     replot = replot_graph(img_thrs(1), get_custom, plot_tupp)
  Else If (img_thrs(2)%id == imgid) Then
     replot = replot_graph(img_thrs(2), get_custom, plot_tlow)
  Else If (img_cors(1)%id == imgid) Then
     replot = replot_graph(img_cors(1), get_custom, plot_correl)
! - skill maps -
  Else If (img_cvsk(1)%id == imgid) Then
     replot = replot_graph(img_cvsk(1), get_custom, plot_skill_cv)
  Else If (img_rask(1)%id == imgid) Then
     replot = replot_graph(img_rask(1), get_custom, plot_skill_ra)
  Else If (img_prss(1)%id == imgid) Then
     replot = replot_graph(img_prss(1), get_custom, plot_pscore)
  End If
!
  Return
 End Function replot_maps
!
!
!
 Function export_maps() &
          Result (export)
!
! Modules
  Use labels,      Only: l_averages, l_correls
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
     export = export_montage_maps(-imgid)
     Return     
  End If
!
! Save graph
! - climatological maps -
  If (img_aves(1)%id == imgid) Then
     export = save_graphic(l_averages%c, img_aves(1), plot_aves)
  Else If (img_thrs(1)%id == imgid) Then
     export = save_graphic('Upper_thresholds', img_thrs(1))
  Else If (img_thrs(2)%id == imgid) Then
     export = save_graphic('Lower_thresholds', img_thrs(2))
  Else If (img_cors(1)%id == imgid) Then
     export = save_graphic(l_correls%c, img_cors(1))
! - skill maps -
  Else If (img_cvsk(1)%id == imgid) Then
     export = save_graphic('Skill_cv', img_cvsk(1))
  Else If (img_rask(1)%id == imgid) Then
     export = save_graphic('Skill_retro', img_rask(1))
  Else If (img_prss(1)%id == imgid) Then
     export = save_graphic('Skill_retro', img_prss(1))
  End If
!
  Return
 End Function export_maps
!
!
!
 Function export_montage_maps(imgid) &
          Result (export)
!
! Modules
  Use labels,   Only: l_threshs
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
! - climatological maps -
  If (img_thrs(1)%id == imgid) Then
     export = export_montage(l_threshs%c, img_thrs(2), img_thrs(1))
  End If
!
  Return
 End Function export_montage_maps
!
!
!
 Function help_maps()
!
! Modules
  Use help,        Only: help_file
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: help_maps
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any
!
! Executable Statements
!
! Open context-specific help pages
! - climatological maps -
  If (img_aves(1)%id == imgid) Then
     Call help_file ('Tools_Climate_Averages.html')
  Else If (Any(img_thrs(:)%id == imgid)) Then
     Call help_file ('Tools_Climate_Thresholds.html')
  Else If (img_cors(1)%id == imgid) Then
     Call help_file ('Tools_Climate_Correlations.html')
! - skill maps -
  Else If (img_cvsk(1)%id == imgid .or. img_rask(1)%id == imgid) Then
     Call help_file ('Tools_Validate_Map.html')
  Else If (img_prss(1)%id == imgid) Then
     Call help_file ('Tools_Verify_PSkillMaps.html')
  End If
  help_maps = 2
!
  Return
 End Function help_maps
#endif
End Module maps
