! Author: Simon Mason
!
! Climate Predictability Tool for MOS Downscaling
!
! Written by Dr Simon J. Mason, Dr Michael K. Tippett, and Lulin Song
!
#if GUI == 1 && FTN95 == 1
Winapp
!
#endif
Program cpt
!
! Modules
  Use projects,      Only: get_project_startup
  Use cpt_interface, Only: cpt_window, title_page
!
! Implicit declarations
  Implicit None
!
! Executable Statements
!
! Initialise CPT
  Call init_cpt1 ()
!
! Check for project from comand line
  Call get_project_startup ()
!
! Print out title page
  Call title_page ()
!
! Complete initialisation
  Call init_cpt2 ()
!
! Prompt for program settings
  Call cpt_window ()
!
Contains
!
!
 Subroutine init_cpt1 ()
!
! Initialises CPT settings required before CPT can start up
!
! Modules
  Use data_numbers,      Only: sp, oneh=>oneh_sp
  Use data_io_constants, Only: ddir
  Use data_version,      Only: cpt_name, cpt_title, cpt_ver, ver, ver_major, ver_minor, ver_patch
  Use maths,             Only: init_numbers
  Use labels,            Only: init_chars, init_labels
  Use time,              Only: init_time
  Use screen,            Only: init_screen
#if GUI == 1
  Use screen,            Only: init_fonts
  Use gui,               Only: init_cols
#endif
  Use errors,            Only: cpt_error, init_errorl
  Use iofiles,           Only: init_dirs
#if GUI == 1
  Use help,              Only: init_help
#endif
  Use verif,             Only: init_verif
  Use user,              Only: get_user
!
! Implicit declarations
  Implicit None
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
  Character(Len=11) :: cproc ! - procedure -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Trim
!
! External routines
  External init_sets
!
! Executable Statements
!
! Initialise version
  ver = Real(ver_major, Kind=sp) + Real(ver_minor, Kind=sp)/oneh
  Write (Unit=cpt_title, Fmt='(A,I2)') cpt_name//' ', ver_major
  Write (Unit=cpt_ver, Fmt='(I0,A,I0,A,I0)') ver_major, '.', ver_minor, '.', ver_patch
!
! Initialise time
  Call init_time ()
!
! Initialise numbers
  Call init_numbers ()
!
! Initialise screen size settings
  Call init_screen ()
!
#if GUI == 1
! Initialise colours
  Call init_cols ()
!
! Initialise GUI settings
  Call init_fonts ()
!
#endif
! Initialise characters
  Call init_chars ()
!
! Initialise directories
  Call init_dirs (ifail)
  If (ifail /= 0) Then
     cproc = 'init_dirs'
     GoTo 1
  End If
!
! Initialise labels
  Call init_labels (1,Trim(ddir),ifail)
  If (ifail /= 0) Then
     cproc = 'init_labels'
     GoTo 1
  End If
!
! Initialise verification
  Call init_verif ()
!
! Identify user
  Call get_user ()
!
! Initialise settings
  Call init_sets ()
!
#if GUI == 1
! Identify help file
  Call init_help ()
!
! Initialise font settings
  Call init_fonts () ! - re-initialise for GUI-version in case of text resizing -
#endif
!
!
  Return
!
! Errors
1 Call init_errorl ()
  Call cpt_error (cproc, ifail, .false.)
  Stop
 End Subroutine init_cpt1
!
!
!
 Subroutine init_cpt2 ()
!
! Initialises CPT settings required before CPT can prompt for options
!
! Modules
#if GUI == 1
  Use gui,      Only: mcw, &
                      init_cw
  Use errors,   Only: cpt_error
#endif
  Use iofiles,  Only: init_iofiles
#if GUI == 1
  Use settings, Only: init_elevel
  Use pcs,      Only: lcde
#endif
!
! Implicit declarations
  Implicit None
#if GUI == 1
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! External routines
  Interface
   Subroutine init_image_sizes ()
   End Subroutine init_image_sizes
  End Interface
!
! Executable Statements
!
! Re-initialise expertise level
  Call init_elevel (lcde)
!
! Initialise images
  Call init_image_sizes ()
!
! Initialise child windows
  Call init_cw (ifail)
  If (ifail /= 0) Then
     ifail = 1
     Call cpt_error ('init_cw', ifail, .false., &
          i_arg1=mcw)
     Stop
  End If
#endif
!
! Initialise structures and output file descriptions
  Call init_iofiles ()
!
  Return
 End Subroutine init_cpt2
End Program cpt
!
!

#if GUI == 1
 Subroutine init_image_sizes ()
!
! Initialises CPT settings required before CPT can prompt for options
!
! Modules
  Use screen,   Only: igsize, ihsize
  Use validate, Only: init_img_validate
  Use maps,     Only: init_img_maps
  Use forecast, Only: init_img_fcst
  Use eofs,     Only: init_img_eofs
  Use scatter,  Only: init_img_scatter
  Use verif,    Only: init_img_verify
!
! Implicit declarations
  Implicit None
!
! Executable Statements
!
! Initialise validation plots
  Call init_img_validate (igsize)
!
! Initialise maps
  Call init_img_maps (igsize)
!
! Initialise forecast images
  Call init_img_fcst (igsize, ihsize)
!
! Initialise scree plots
  Call init_img_eofs (igsize, ihsize)
!
! Initialise scatter plots
  Call init_img_scatter (igsize, ihsize)
!
! Initialise verification plots
  Call init_img_verify (igsize, ihsize)
!
  Return
 End Subroutine init_image_sizes
#endif
!
!
!
Subroutine init_sets ()
!
! Initialises CPT default settings
!
! Modules
  Use time,       Only: reset_time
  Use screen,     Only: init_gfonts, window_update
  Use gui,        Only: set_win_title
#if GUI == 1
  Use gui,        Only: reset_jpg
#endif
  Use iofiles,    Only: xfile, yfile, zfile, &
                        reset_iofiles
  Use settings,   Only: igcms, istd, isynch, nf, nt, sw_elev, &
                        reset_limits, reset_settings, training_refresh
  Use space,      Only: reset_space
  Use fields,     Only: nmfd, &
                        reset_fields
  Use scores,     Only: reset_scores
  Use season,     Only: icsv, &
                        reset_season
  Use pcs,        Only: init_pcs
  Use missing,    Only: xmiss,ymiss
  Use categories, Only: init_refdates, init_icat
  Use climate,    Only: ncy, &
                        reset_climate
  Use analysis,   Only: ianal, icalc, &
                        reset_analysis
  Use markers,    Only: init_marker
  Use cpt_view,   Only: set_menuflags
  Use bootstrap,  Only: reset_bootstrap
  Use forecast,   Only: reset_forecast
  Use verif,      Only: reset_verif
#if GUI == 1
  Use cpt_file,   Only: files_refresh
#endif
  Use user,       Only: cpt_ini, &
                        read_ini
  Use downloads,  Only: init_download
!
! Implicit declarations
  Implicit None
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
! Read customised defaults
  Call read_ini (Trim(cpt_ini))
!
! Update window
  Call training_refresh (ncy)
!
! Initialise settings
  Call reset_settings (ianal)
!
! Reset input file settings
  Call reset_iofiles ()
!
! Reset field settings
  Call reset_fields (xmiss, ymiss, ifail)
  If (ifail /= 0) Return
  nmfd = 0
!
! Reset space settings
  Call reset_space ()
!
#if GUI == 1
! Update window
  Call files_refresh ()
!
#endif
! Reset graphics font scaling
  Call init_gfonts ()
!
! Initialise analysis settings
  Call reset_analysis ()
!
! Current EOFs
  Call init_pcs ()
!
! Scores selections
  Call reset_scores ()
!
! Switch off climatological period
  Call reset_climate ()
!
! Initialise category indicator
  Call init_icat ()
!
! Switch off reference years
  Call init_refdates ()
!
! Reset marker settings
  Call init_marker ()
!
! Reset verification settings
  Call reset_verif ()
!
! Reset forecast settings
  Call reset_forecast ()
!
! Reset time settings
  Call reset_time (.true.)
!
! Initialise seasonal values settings
  Call reset_season (.false., .true.)
!
! Initialise download settings
  Call init_download ()
!
! Initialise resampling settings
  Call reset_bootstrap ()
!
#if GUI == 1
! Initialise JPEG quality
  Call reset_jpg ()
#endif
!
! Set view flags
  Call set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, sw_elev%on, &
       ixdstr=xfile%idstr, iydstr=yfile%idstr, izdstr=zfile%idstr)
!
! Reset window caption and update training period settings
  Call set_win_title (ianal)
  If (ianal > 0) Then
     Call reset_limits (0, .false., 0, (/0/), 0, 0)
     Call window_update (nt)
     Call window_update (nf)
  End If
!
  Return
End Subroutine init_sets
!
!
!
#if GUI == 1
Subroutine draw_domain (isize, mlim, alim, ctitle, ifail)
!
! Modules
  Use space,    Only: area
  Use graphics, Only: draw_simple_map
!
! Implicit declarations
  Implicit None
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: isize ! - image size -
!
  Character(Len=*), Intent(In) :: ctitle ! - map title -
!
  Type(area), Intent(In) :: mlim ! - map domain -
  Type(area), Intent(In) :: alim ! - selected domain -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Draw map
  Call draw_simple_map (isize, mlim, alim, ctitle, ifail)
!
  Return
End Subroutine draw_domain
#endif
