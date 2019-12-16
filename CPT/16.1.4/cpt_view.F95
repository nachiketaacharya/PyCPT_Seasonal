! Author: Simon Mason
Module cpt_view
!
! Modules
  Use menus, Only: switch
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: change_view, set_view, set_menuflags
!
! Scalars
!
! Integer scalars
! - grey menu flag complements -
  Integer, Public :: jm_dia   ! - complement of diagnostics grey menu flag -
  Integer, Public :: jm_for   ! - complement of forecast grey menu flag -
!
! Derived type scalars
! - analyses -
  Type(switch), Public :: sw_cca ! - CCA menu flag -
  Type(switch), Public :: sw_gcm ! - GCM menu flag -
  Type(switch), Public :: sw_mlr ! - MLR menu flag -
  Type(switch), Public :: sw_pcr ! - PCR menu flag -
  Type(switch), Public :: sw_pfv ! - PFV menu flag -
!
Contains
!
!
 Function change_view()
!
! Changes analysis option
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_dia, ia_for, ia_pcr, ia_mlr, ia_gcm, ia_pfv, mnt
  Use labels,             Only: cls_l
  Use screen,             Only: window_update
  Use gui,                Only: set_win_title
  Use iofiles,            Only: xfile, yfile, zfile
  Use settings,           Only: iensc, iensc_bk, igauss, igauss_bk, igcms, iretro, istd, istd_bk, isynch, mnu, sw_elev
  Use season,             Only: icsv
  Use pcs,                Only: ieofx, ieofx_bk, im_yeof, npx, npy
  Use categories,         Only: ithr, ithr_bk
  Use climate,            Only: ieclim, ieclim_bk
  Use analysis,           Only: ianal, ianaln, icalc
!
! Function type
  Integer :: change_view
!
! Executable Statements
!
! Convert from beginner mode
  change_view = 2
  Select Case (ianaln)
   Case (ia_dia, ia_for)
     If (ianal == ia_pfv) Then
        ianaln = ia_cca
     Else
        ianaln = ianal
     End If
     If (ianaln == ia_cca) Then
        im_yeof = 1
     Else
        im_yeof = 0
     End If
     If (ianal == ia_dia) Then
        isynch = 1
     Else
        isynch = 0
     End If
!     jdia=1-isynch
!     jfor=isynch
  End Select
!
! Return If no change in analysis
  If (ianaln == ianal) Return
!
! Update CPT window
  Select Case (ianaln)
   Case (ia_cca, ia_pcr) ! - CCA and PCR -
     ieofx = ieofx_bk
     If (ianal == ia_pfv) Then
        igauss = igauss_bk
        ieclim = ieclim_bk
        istd = istd_bk
        ithr = ithr_bk
     End If
     If (ianal == ia_gcm) iensc_bk = iensc
     iensc = 1
     mnu = mnt
   Case (ia_mlr) ! - MLR -
     If (ianal == ia_cca) ieofx_bk = ieofx
     ieofx = 2
     If (ianal == ia_pfv) Then
        igauss = igauss_bk
        ieclim = ieclim_bk
        istd = istd_bk
        ithr = ithr_bk
     End If
     If (ianal == ia_gcm) iensc_bk = iensc
     iensc = 1
     mnu = mnt
   Case (ia_gcm) ! - GCM -
     If (ianal == ia_cca) ieofx_bk = ieofx
     ieofx = 2
     npx = 0
     npy = 0
     If (ianal == ia_pfv) Then
        ieclim = ieclim_bk
        istd = istd_bk
        ithr = ithr_bk
     End If
     iensc = iensc_bk
     mnu = mnt
   Case (ia_pfv) ! - PFV -
     If (ianal == ia_cca) ieofx_bk = ieofx
     ieofx = 2
     igauss_bk = igauss
     igauss = 0
     npx = 0
     npy = 0
     ieclim_bk = ieclim
     ieclim = 1
     If (ianal == ia_gcm) iensc_bk = iensc
     iensc = 1
     istd_bk = istd
     istd = 0
     ithr_bk = ithr
     ithr = 1
     mnu = 1
  End Select
  ianal = ianaln
  Call set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, sw_elev%on, &
                      ixdstr=xfile%idstr, iydstr=yfile%idstr, izdstr=zfile%idstr, iretro=iretro)
  Call set_win_title (ianal)
!
! Update lags/seasons label
  Call window_update (cls_l)
!
  Return
 End Function change_view
!
!
!
 Function set_view()
!
! Sets analysis option
!
! Modules
  Use data_cpt_constants, Only: ia_mlr, ia_gcm, ia_pfv, mnt
  Use settings,           Only: igauss, istd, mnu
  Use pcs,                Only: ieofx
  Use climate,            Only: ieclim
  Use analysis,           Only: ianal
!
! Function type
  Integer :: set_view
!
! Executable Statements
!
! Set X EOF option for MLR, GCM, and PFV
  set_view = 0
  If (ianal == ia_mlr .or. ianal == ia_gcm .or. ianal == ia_pfv) Then
     ieofx = 2
!
! Switch on extended climatology, and switch off standardization and transformation flags, for PFV
     If (ianal == ia_pfv) Then 
        ieclim = 1
        igauss = 0
        istd = 0
        mnu = 1
     End If
  End If
  If (ianal /= ia_pfv) mnu = mnt
!
  Return
 End Function set_view
!
!
!
 Subroutine set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, lelev, &
            ixdstr, iydstr, izdstr, iretro, f_scree)
!
! Sets which menu items to display or grey
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use labels,             Only: cg_lags_l, cg_seasons_l, &
                                cls_l
  Use menus,              Only: im_edit, im_ov, im_seas, im_trns
  Use iofiles,            Only: im_sat, im_saz, &
                                set_grey_files
  Use space,              Only: igeog
  Use pcs,                Only: set_grey_modes
  Use analysis,           Only: iaction, jm_gauss
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal  ! - analysis flag -
  Integer, Intent(In) :: icalc  ! - calculation completed flag -
  Integer, Intent(In) :: istd   ! - standardization option -
  Integer, Intent(In) :: igcms  ! - GCM standardization option -
  Integer, Intent(In) :: icsv   ! - calculate seasonal values option -
  Integer, Intent(In) :: isynch ! - allow synchronous predictots flag -
!
  Logical, Intent(In) :: lelev  ! - expert level? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: ixdstr ! - X file structure -
  Integer, Intent(In), Optional :: iydstr ! - Y file structure -
  Integer, Intent(In), Optional :: izdstr ! - Z file structure -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: iretro ! - retroactive analysis flag -
!
! - optional arguments -
  Integer, External, Optional :: f_scree ! - scree plot settings Callback function -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Set grey View menu items
  If (Present(ixdstr) .and. Present(iydstr)) Then
     Call set_grey_view (ianal, isynch, ixdstr, lelev)
!
! Set grey Modes menu items
    Call set_grey_modes (ianal, ixdstr, &
         f_scree=f_scree)
!
! Set grey menu items
     Select Case (ianal)
      Case (ia_cca) ! - CCA -
        im_seas = icsv
        jm_gauss = 1
        cls_l = cg_lags_l
      Case (ia_pcr) ! - PCR -
        im_seas = icsv
        jm_gauss = 1
        cls_l = cg_lags_l
      Case (ia_mlr) ! - MLR -
        im_seas = icsv
        jm_gauss = 1
        cls_l = cg_lags_l
      Case (ia_gcm) ! - GCM -
        im_seas = icsv
        If (igcms == 3) Then
           jm_gauss = 1
        Else
           jm_gauss = 0
        End If
        cls_l = cg_lags_l
      Case (ia_pfv) ! - PFV -
        im_seas = 0
        jm_gauss = 0
        cls_l = cg_seasons_l
        im_ov = 1
        If (Present(iretro)) iretro = 1
     End Select
!
! Set file open and edit menu flag
     Call set_grey_files (ixdstr, iydstr, iaction, igeog)
!
! Set options flags
! - verification options flag -
     If (ianal /= ia_pfv) Then
        If (Present(iretro)) Then
           im_ov = iretro
        Else
           im_ov = 0
        End If
!
! - transformation option flag -
        Select Case (istd)
         Case Default
           im_trns = 1
         Case (3)
           im_trns = 0
        End Select
     Else
        im_trns = 0
     End If
  End If
!
! Set forecast file flags
  If (Present(izdstr)) Then
     If (izdstr > 0) Then
        im_edit = 1
        im_sat = 1
        im_saz = 1
     Else
        If (icalc == 1) Then
           im_edit = 0
           im_sat = 0
        End If
        im_saz = 0
     End If
  End If
!
  Return
 End Subroutine set_menuflags
!
!
!
 Subroutine set_grey_view (ianal, isynch, ixdstr, lelev)
!
! Sets which View menu items to display or grey
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use menus,              Only: Assignment(=)
  Use analysis,           Only: ihind
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal  ! - analysis flag -
  Integer, Intent(In) :: isynch ! - allow synchronous predictots flag -
  Integer, Intent(In) :: ixdstr ! - X file structure -
!
  Logical, Intent(In) :: lelev ! - expert level? -
!
! Executable Statements
!
! Set grey View menu items
  Select Case (ianal)
   Case (ia_cca) ! - CCA -
     ihind = 1
     sw_cca = 1
     sw_pcr = 0
     sw_mlr = 0
     sw_gcm = 0
   Case (ia_pcr) ! - PCR -
     ihind = 1
     sw_cca = 0
     sw_pcr = 1
     sw_mlr = 0
     sw_gcm = 0
   Case (ia_mlr) ! - MLR -
     ihind = 1
     sw_cca = 0
     sw_pcr = 0
     sw_mlr = 1
     sw_gcm = 0
   Case (ia_gcm) ! - GCM -
     ihind = 1
     sw_cca = 0
     sw_pcr = 0
     sw_mlr = 0
     sw_gcm = 1
   Case (ia_pfv) ! - PFV -
     ihind = 0
     sw_cca = 0
     sw_pcr = 0
     sw_mlr = 0
  End Select
!
! Set view flags depending on whether input files are open
  If (ianal /= ia_pfv) Then
     sw_pfv%ion = 0
     Select Case (ixdstr)
      Case (0)
        If (ianal /= ia_mlr) sw_mlr%iacc = 1
        sw_pfv%iacc = 1
      Case (3)
        If (ianal /= ia_mlr) sw_mlr%iacc = 1
        sw_pfv%iacc = 0
      Case Default
        If (ianal /= ia_mlr) sw_mlr%iacc = 0
        sw_pfv%iacc = 0
     End Select
     If (.not.lelev) Then
        jm_dia = 1 - isynch
        jm_for = isynch
     End If
  Else ! - switch off all but PFV if a probabilistic X input file is open -
     sw_pfv%ion = 1
     sw_pfv%iacc = 0
     Select Case (ixdstr)
      Case (1:)
        sw_cca%iacc = 0
        sw_pcr%iacc = 0
        sw_mlr%iacc = 0
        sw_gcm%iacc = 0
      Case Default
        sw_cca%iacc = 1
        sw_pcr%iacc = 1
        sw_mlr%iacc = 1
        sw_gcm%iacc = 1
     End Select
     Select Case (ixdstr)
      Case (1:)
        jm_dia = 0
        jm_for = 0
      Case Default
        jm_dia = 1
        jm_for = 1
     End Select
  End If
  sw_pfv%ioff=1-sw_pfv%ion
!
  Return
 End Subroutine set_grey_view
End Module cpt_view
