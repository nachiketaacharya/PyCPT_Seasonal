! Author: Simon Mason
Module menus
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: add_menu
!
! Derived-type definitions
!
! - grey -
  Public :: switch
  Type switch
     Sequence
     Integer :: ion  ! - on -
     Integer :: ioff ! - off -
     Integer :: iacc ! - accessible -
     Logical :: on   ! - on -
  End Type switch
!
#if GUI == 0
! Parameters
!
! Integer parameters
  Integer, Parameter, Private :: lmax = 9 ! - maximum number of menu levels -
!
#endif
! Scalars
!
! Integer scalars
! - grey menu flags - 
  Integer, Public :: im_edit ! - edit grey menu flag -
  Integer, Public :: im_ov   ! - verification options grey menu flag -
  Integer, Public :: im_seas ! - target season grey menu flag -
  Integer, Public :: im_trns ! - transformation option grey menu flag -
! - miscellaneous flags -
  Integer, Public :: inew = 0   ! - new (i.e., no open files) grey menu flag -
#if GUI == 1
  Integer, Public :: iexit      ! - exit program flag -
#else
  Integer, Public :: imenus = 1 ! - print menus (default is on)? -
!
! Logical arrays
  Logical, Dimension(lmax), Public :: lgrey ! - greyed menu items? -
#endif
!
! Interfaces
!
! Interface assignments
  Public :: Assignment(=)
  Interface Assignment(=)
     Module Procedure init_switch
  End Interface
!
! Generic interfaces
  Interface add_menu
   Module Procedure add_menu_lab
   Module Procedure add_menu_txt1
#if FTN95 == 0 && GFORTRAN == 0
   Module Procedure add_menu_txt3
#endif
  End Interface add_menu
!
Contains
!
!
 Subroutine init_switch (s, i)
!
! Initialises a switch
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i ! - assignment value -
!
! Output scalars
  Type(switch), Intent(Out) :: s ! - switch -
!
! Executable Statements
!
! Assign switch values
  Select Case (i)
   Case (0)
     s%ion = 0
   Case Default
     s%ion = 1
  End Select
  s%ioff = 1 - s%ion
  s%iacc = s%ioff
  s%on = (s%ion == 1)
!
  Return
 End Subroutine init_switch
!
!
!
 Subroutine add_menu_lab (cmn, lmenu, imenu, lab, &
            c_suf, lsep, itab, igrey, icheck, itoggle, iset1, ival1, iset2, ival2, iset3, ival3, scbf, cbf, acc_key)
!
! Adds a menu item
!
! Modules
#if GUI == 0
  Use data_text, Only: c_tab
#endif
  Use labels,    Only: label, &
                       c_lab, &
                       make_label
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imenu ! - menu identifier -
  Integer, Intent(In) :: lmenu ! - level of menu -
!
  Character(Len=*), Intent(In) :: cmn ! - menu type -
!
  Type(label), Intent(In) :: lab ! - label -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: igrey ! - toggle for greyed menu -
  Integer, Intent(In), Optional :: itab  ! - tab point -
  Integer, Intent(In), Optional :: ival1 ! - value to set to iset1 -
  Integer, Intent(In), Optional :: ival2 ! - additional value to set to iset2 -
  Integer, Intent(In), Optional :: ival3 ! - additional value to set to iset3 -
!
  Character(Len=*), Intent(In), Optional :: acc_key ! - accelerator key -
  Character(Len=*), Intent(In), Optional :: c_suf   ! - label suffix -
  Character(Len=*), Intent(In), Optional :: scbf    ! - Silverfrost call-back function -
!
  Logical, Intent(In), Optional :: lsep ! - include separator? -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: icheck  ! - check mark -
  Integer, Intent(InOut), Optional :: itoggle ! - value to toggle -
  Integer, Intent(InOut), Optional :: iset1   ! - variable to set -
  Integer, Intent(InOut), Optional :: iset2   ! - additional variable to set -
  Integer, Intent(InOut), Optional :: iset3   ! - additional variable to set -
!
! Procedure arguments
! - optional procedure arguments -
  Integer, External, Optional :: cbf ! - call-back function -
!
! Locals
!
! Local scalars
#if GUI == 0
  Integer :: ltrim ! - trimmed length of label without the suffix -
!
#endif
  Character(Len=1) :: cmnfmt ! - menu format -
!
  Logical :: llak ! - include accelerator key? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 0
  Intrinsic Index, Repeat
#endif
!
! Executable Statements
!
! Construct label
#if GUI == 1
  cmnfmt = 'T'
  llak = .true.
#else
  Select Case (lmenu)
   Case (1)
     cmnfmt = 'U'
   Case Default
     cmnfmt = 'T'
  End Select
  llak = .false.
#endif
  c_lab = make_label(cmnfmt, lab, &
          lak=llak)
  If (Present(c_suf)) Then
     c_lab = Trim(c_lab)//c_suf
#if GUI == 0
     If (Index(c_suf, c_tab) == 1 .and. Present(itab)) Then
        ltrim = Index(c_lab, c_tab) - 1
        c_lab = c_lab(1:ltrim)//Repeat(' ', itab-ltrim)//Trim(c_lab(ltrim+2:))
     End If
#endif
  End If
!
! Add menu item
  Call add_menu_txt1 (cmn, lmenu, imenu, c_lab, &
       lsep=lsep, igrey=igrey, icheck=icheck, itoggle=itoggle, &
       iset1=iset1, ival1=ival1, iset2=iset2, ival2=ival2, iset3=iset3, ival3=ival3, &
       scbf=scbf, cbf=cbf, acc_key=acc_key)
!
  Return
 End Subroutine add_menu_lab
!
!
!
 Subroutine add_menu_txt1 (cmn, lmenu, imenu, cmenu, &
            lsep, igrey, icheck, itoggle, iset1, ival1, iset2, ival2, iset3, ival3, scbf, cbf, acc_key)
!
! Adds a menu item
!
! Modules
#if GUI == 1
  Use clrwin$,   Only: winio$
  Use data_text, Only: c_tab
  Use gui,       Only: iw
#else
  Use labels,    Only: cg_turnoff_l, cg_turnon_l
  Use gui,       Only: print_text
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imenu ! - menu identifier -
  Integer, Intent(In) :: lmenu ! - level of menu -
!
  Character(Len=*), Intent(In) :: cmn   ! - menu type -
  Character(Len=*), Intent(In) :: cmenu ! - menu label -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: igrey ! - toggle for greyed menu -
  Integer, Intent(In), Optional :: ival1 ! - value to set to iset1 -
  Integer, Intent(In), Optional :: ival2 ! - additional value to set to iset2 -
  Integer, Intent(In), Optional :: ival3 ! - additional value to set to iset3 -
!
  Character(Len=*), Intent(In), Optional :: acc_key ! - accelerator key -
  Character(Len=*), Intent(In), Optional :: scbf    ! - Silverfrost call-back function -
!
  Logical, Intent(In), Optional :: lsep ! - include separator? -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: icheck  ! - check mark -
  Integer, Intent(InOut), Optional :: itoggle ! - value to toggle -
  Integer, Intent(InOut), Optional :: iset1   ! - variable to set -
  Integer, Intent(InOut), Optional :: iset2   ! - additional variable to set -
  Integer, Intent(InOut), Optional :: iset3   ! - additional variable to set -
!
! Procedure arguments
  Integer, External, Optional :: cbf ! - call-back function -
!
! Locals
!
! Local scalars
#if GUI == 1
  Character(Len=256) :: cfmt  ! - menu format text -
#else
  Character(Len=  8) :: cmnum ! - menu number -
#endif
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 1
  Intrinsic Repeat
#else
  Intrinsic All
#endif
!
! Executable Statements
!
! Construct menu item
#if GUI == 1
  cfmt = '%'//cmn//Repeat('[', lmenu)                   ! - open levels -
  If (Present(lsep)) Then
     If (lsep) cfmt = Trim(cfmt)//'|, '                 ! - add separator -
  End If
  If (Present(icheck)) Then
     cfmt = Trim(cfmt)//'#'                             ! - add check mark -
  End If
  If (Present(igrey)) Then
     cfmt = Trim(cfmt)//'~'                             ! - add greying -
  End If
  cfmt = Trim(cfmt)//Trim(cmenu)                        ! - add menu item -
  If (Present(acc_key)) Then
     cfmt = Trim(cfmt)//c_tab//'('//Trim(acc_key)//')'  ! - add accelerator key indicator -
  End If
  If (.not.(Present(cbf) .or. Present(scbf) .or. Present(icheck) .or. Present(itoggle) .or. &
           (Present(iset1) .and. Present(ival1)))) Then
     cfmt = Trim(cfmt)//'[]'                            ! - add empty level if there is no callback -
  End If
  cfmt = Trim(cfmt)//Repeat(']', lmenu)//'&'            ! - close levels -
!
! Add menu item
  If (Present(icheck)) Then
     If (Present(igrey)) Then
        If (Present(itoggle)) Then
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                           'SET', iset2, ival2, &
                                                                           'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, &
                                                                 'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                            'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', 'TOGGLE', itoggle, 'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'TOGGLE', itoggle, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', 'TOGGLE', itoggle, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', 'TOGGLE', itoggle, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               'TOGGLE', itoggle)
                       End If
                    End If
           End If
        Else ! - no toggle -
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', '+', 'SET', iset1, ival1, &
                                                   'SET', iset2, ival2, &
                                                   'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, &
                                         'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', 'SET', iset1, ival1, &
                                    'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', '+', 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               '+', cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey, &
                               scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, igrey)
                       End If
                    End If
           End If
        End If
     Else ! - no grey -
        If (Present(itoggle)) Then
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                           'SET', iset2, ival2, &
                                                                           'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, &
                                                                 'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                            'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', 'TOGGLE', itoggle, 'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'TOGGLE', itoggle, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', 'TOGGLE', itoggle, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', 'TOGGLE', itoggle, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               'TOGGLE', itoggle)
                       End If
                    End If
           End If
        Else ! - no toggle -
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', '+', 'SET', iset1, ival1, &
                                                   'SET', iset2, ival2, & 
                                                   'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, &
                                         'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', 'SET', iset1, ival1, &
                                    'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', '+', 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               '+', cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               icheck, &
                               scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               icheck)
                       End If
                    End If
           End If
        End If
     End If
  Else ! - no check -
     If (Present(igrey)) Then
        If (Present(itoggle)) Then
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                           'SET', iset2, ival2, &
                                                                           'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, &
                                                                 'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                            'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', 'TOGGLE', itoggle, 'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'TOGGLE', itoggle, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', 'TOGGLE', itoggle, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', 'TOGGLE', itoggle, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               'TOGGLE', itoggle)
                       End If
                    End If
           End If
        Else ! - no toggle -
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', '+', 'SET', iset1, ival1, &
                                                   'SET', iset2, ival2, &
                                                   'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, &
                                         'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', 'SET', iset1, ival1, &
                                    'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', '+', 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               '+', cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               igrey, &
                               scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               igrey)
                       End If
                    End If
           End If
        End If
     Else ! - no grey -
        If (Present(itoggle)) Then
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                           'SET', iset2, ival2, &
                                                                           'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, &
                                                                      'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, &
                                                                 'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                      'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                                 'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, &
                                                            'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'TOGGLE', itoggle, 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', 'TOGGLE', itoggle, 'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'TOGGLE', itoggle, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', 'TOGGLE', itoggle, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', 'TOGGLE', itoggle, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               'TOGGLE', itoggle)
                       End If
                    End If
           End If
        Else ! - no toggle -
           If (Present(iset1) .and. Present(ival1)) Then
              If (Present(iset2) .and. Present(ival2)) Then
                 If (Present(iset3) .and. Present(ival3)) Then
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', '+', 'SET', iset1, ival1, &
                                                   'SET', iset2, ival2, &
                                                   'SET', iset3, ival3, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, &
                                              'SET', iset3, ival3, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, &
                                         'SET', iset3, ival3)
                       End If
                    End If
                 Else ! - no iset3 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', '+', 'SET', iset1, ival1, &
                                              'SET', iset2, ival2, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'SET', iset1, ival1, &
                                         'SET', iset2, ival2, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                             '+', 'SET', iset1, ival1, &
                                  'SET', iset2, ival2)
                       End If
                    End If
                 End If
              Else ! - no iset2 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', '+', 'SET', iset1, ival1, cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               '+', 'SET', iset1, ival1, cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', 'SET', iset1, ival1, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               'SET', iset1, ival1)
                       End If
                    End If
              End If
           Else ! - no iset1 -
                    If (Present(cbf)) Then
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               '+', cbf, scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt), &
                               cbf)
                       End If
                    Else ! - no user callback -
                       If (Present(scbf)) Then
                          iw = winio$(Trim(cfmt), &
                               scbf)
                       Else ! - no Silverfrost callback -
                          iw = winio$(Trim(cfmt))
                       End If
                    End If
           End If
        End If
     End If
  End If
!
! Set accelerator key
  If (Present(acc_key)) Then
     If (Present(iset1) .and. Present(ival1)) Then
        If (Present(iset2) .and. Present(ival2)) Then
           If (Present(iset3) .and. Present(ival3)) Then
              If (Present(cbf)) Then
                 If (Present(scbf)) Then
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', '+', '+', 'SET', iset1, ival1, &
                                             'SET', iset2, ival2, &
                                             'SET', iset3, ival3, scbf, cbf)
                 Else ! - no Silverfrost callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', '+', 'SET', iset1, ival1, &
                                        'SET', iset2, ival2, &
                                        'SET', iset3, ival3, cbf)
                 End If
              Else ! - no user callback -
                 If (Present(scbf)) Then  ! - no user callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', '+', 'SET', iset1, ival1, &
                                        'SET', iset2, ival2, &
                                        'SET', iset3, ival3, scbf)
                 Else
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', 'SET', iset1, ival1, &
                                   'SET', iset2, ival2, &
                                   'SET', iset3, ival3)
                 End If
              End If
           Else ! - no iset3 -
              If (Present(cbf)) Then
                 If (Present(scbf)) Then
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', '+', 'SET', iset1, ival1, &
                                        'SET', iset2, ival2, scbf, cbf)
                 Else ! - no Silverfrost callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', 'SET', iset1, ival1, &
                                   'SET', iset2, ival2, cbf)
                 End If
              Else ! - no user callback -
                 If (Present(scbf)) Then  ! - no user callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', 'SET', iset1, ival1, &
                                   'SET', iset2, ival2, scbf)
                 Else
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', 'SET', iset1, ival1, &
                              'SET', iset2, ival2)
                 End If
              End If
           End If
        Else ! - no iset2 -
              If (Present(cbf)) Then
                 If (Present(scbf)) Then
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', '+', 'SET', iset1, ival1, scbf, cbf)
                 Else ! - no Silverfrost callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', 'SET', iset1, ival1, cbf)
                 End If
              Else ! - no user callback -
                 If (Present(scbf)) Then  ! - no user callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', 'SET', iset1, ival1, scbf)
                 Else
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         'SET', iset1, ival1)
                 End If
              End If
        End If
     Else ! - no iset1 -
              If (Present(cbf)) Then
                 If (Present(scbf)) Then
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         '+', scbf, cbf)
                 Else ! - no Silverfrost callback -
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         cbf)
                 End If
              Else ! - no user callback -
                 If (Present(scbf)) Then
                    iw = winio$('%ac['//Trim(acc_key)//']&', &
                         scbf)
                 End If
              End If
     End If
  End If
#else
! Set greying
  If (Present(igrey)) Then
     lgrey(lmenu) = (igrey > 0)
  Else
     lgrey(lmenu) = .true.
  End If
!
! Print menu
  If (All(lgrey(1:lmenu))) Then
     Select Case (lmenu)
      Case (1)
        Call print_text (Trim(cmenu))
      Case Default
        Write (Unit=cmnum, Fmt='(I5,A)') imenu,'.  '
        If (Present(icheck)) Then
           If (icheck == 0) Then
              Call print_text (cmnum//Trim(cmenu)//' ('//Trim(cg_turnon_l)//')', &
                   itab=5*(lmenu-1))
           Else
              Call print_text (cmnum//Trim(cmenu)//' ('//Trim(cg_turnoff_l)//')', &
                   itab=5*(lmenu-1))
           End If
        Else
           Call print_text (cmnum//Trim(cmenu), &
                itab=5*(lmenu-1))
        End If
     End Select
  End If
#endif
!
  Return
 End Subroutine add_menu_txt1
#if FTN95 == 0 && GFORTRAN == 0
!
!
!
 Subroutine add_menu_txt3 (cmn, lmenu, imenu, cmenu, &
            lsep, igrey, icheck, itoggle, iset1, ival1, iset2, ival2, iset3, ival3, scbf, cbf, acc_key)
!
! Adds a menu item
!
! Modules
  Use data_text, Only: i_cyrillic
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imenu ! - menu identifier -
  Integer, Intent(In) :: lmenu ! - level of menu -
!
  Character(Len=*), Intent(In) :: cmn   ! - menu type -
!
  Character(Len=*, Kind=i_cyrillic), Intent(In) :: cmenu ! - menu label -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: igrey ! - toggle for greyed menu -
  Integer, Intent(In), Optional :: ival1 ! - value to set to iset1 -
  Integer, Intent(In), Optional :: ival2 ! - additional value to set to iset2 -
  Integer, Intent(In), Optional :: ival3 ! - additional value to set to iset3 -
!
  Character(Len=*), Intent(In), Optional :: acc_key ! - accelerator key -
  Character(Len=*), Intent(In), Optional :: scbf    ! - Silverfrost call-back function -
!
  Logical, Intent(In), Optional :: lsep ! - include separator? -
!
! Input/output scalars
! - optional input/output scalars -
  Integer, Intent(InOut), Optional :: icheck  ! - check mark -
  Integer, Intent(InOut), Optional :: itoggle ! - value to toggle -
  Integer, Intent(InOut), Optional :: iset1   ! - variable to set -
  Integer, Intent(InOut), Optional :: iset2   ! - additional variable to set -
  Integer, Intent(InOut), Optional :: iset3   ! - additional variable to set -
!
! Procedure arguments
  Integer, External, Optional :: cbf ! - call-back function -
!
! Locals
!
! Local scalars
  Integer :: i ! - index -
!
  Character(Len=128) :: c1menu ! - menu label -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Char, Ichar
!
! Executable Statements
!
! Add menu item
  Do i = 1, Len_Trim(cmenu)
     c1menu(i:i) = Char(Ichar(cmenu(i:i)))
  End Do
  Call add_menu_txt1 (cmn, lmenu, imenu, c1menu, &
       lsep=lsep, igrey=igrey, icheck=icheck, itoggle=itoggle, &
       iset1=iset1, ival1=ival1, iset2=iset2, ival2=ival2, iset3=iset3, ival3=ival3, &
       scbf=scbf, cbf=cbf, acc_key=acc_key)
!
  Return
 End Subroutine add_menu_txt3
#endif
End Module menus
