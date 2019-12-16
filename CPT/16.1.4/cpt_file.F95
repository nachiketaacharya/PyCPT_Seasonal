! Author: Simon Mason
Module cpt_file
!
! Modules
  Use space, Only: domain
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: close_file, files_refresh, merge_files, input_files, update_file
#if GUI == 0
  Public :: get_file
#endif
!
! Scalars
!
! Integer scalars
  Integer, Public :: i_file ! - current input file identifier -
  Integer, Public :: ixyz   ! - file indicator (1 = X; 2 = Y; 3 = Z) -
!
! Arrays
!
! Derived-type arrays
  Type(domain), Dimension(:), Allocatable, Private :: domain_bk ! - back-up domain limits -
!
Contains
!
!
 Function get_file()
!
! Opens an input file
!
! Modules
  Use data_cpt_constants, Only: c_file
#if GUI == 1
  Use screen,             Only: window_update
  Use gui,                Only: set_colour_direction
#endif
  Use gui,                Only: set_cursor_waiting
#if GUI == 1
  Use iofiles,            Only: nz, xfile, zfile
#endif
  Use settings,           Only: sw_elev
#if GUI == 1
  Use fields,             Only: yfield
#endif
  Use pcs,                Only: get_eofx_settings, get_eofy_settings
  Use forecast,           Only: set_fcst
!
! Function type
  Integer :: get_file
!
! Locals
!
! Executable Statements
!
! Open input file
  Select Case (c_file(i_file))
! - X file -
   Case ('X')
     If (sw_elev%on) Then ! - expert mode -
        get_file = get_infile('X', .true., get_xfile, &
                   get_eof_settings=get_eofx_settings)
     Else ! - beginner mode -
        get_file = get_infile('X', .true., get_xfile)
     End If
! - Y file -
   Case ('Y')
     If (sw_elev%on) Then ! - expert mode -
        get_file = get_infile('Y', .true., get_yfile, &
                   get_eof_settings=get_eofy_settings)
     Else ! - beginner mode -
        get_file = get_infile('Y', .true., get_yfile)
     End If
#if GUI == 1
     If (get_file == 0) Call set_colour_direction (yfield(1)%var) ! - set default colour direction -
#endif
! - Z file -
   Case ('Z')
     get_file = get_infile('Z', .false., get_zfile)
     get_file = set_fcst(get_file)
#if GUI == 1
     If (get_file == 0) Call file_refresh (zfile, &
                             nv=nz)
#endif
! - first file -
   Case ('1')
     get_file = get_infile('1', .true., get_xfile)
#if GUI == 1
     If (get_file == 2) Call window_update (xfile%idstr)
#endif
! - second file -
   Case ('2')
     get_file = get_infile('2', .true., get_yfile)
! - update file -
   Case ('3')
     get_file = get_infile('3', .true., get_zfile)
  End Select
! - switch off cursor waiting -
  get_file = 2
  Call set_cursor_waiting (0)
!
  Return
 End Function get_file
!
!
!
 Function get_infile(cxy, lseason, get_afile, &
          get_eof_settings)
!
! Opens an input file, and prompts for numbers of modes if appropriate
!
! Modules
  Use data_cpt_constants,  Only: ia_pfv
  Use data_time_constants, Only: isq_mn, isq_sn, nmn
  Use time,                Only: ims, ims_bk, lmsq, lssq, mpsn, mpsn_bk
  Use screen,              Only: window_update
#if GUI == 1
  Use gui,                 Only: nclr
#endif
  Use gui,                 Only: set_cursor_waiting
  Use menus,               Only: inew
  Use iofiles,             Only: bkfile, bk2file=>pfile, mfile, nv, nx, ny, nz, xfile, yfile, zfile, &
                                 init_ifile
  Use settings,            Only: igcms, iretro, istd, isynch, nt, ntm, sw_elev, &
                                 get_maxnt, record_change
  Use space,               Only: cstndx, cstnx, idomx, rlatdx, rlatx, rlngdx, rlngx, &
                                 cstndy, cstny, idomy, rlatdy, rlaty, rlngdy, rlngy, &
                                 cstndz, cstnz, idomz, rlatdz, rlatz, rlngdz, rlngz
  Use fields,              Only: xfield, yfield, zfield, mfield,&
                                 get_domains, set_mfields
  Use season,              Only: icsv, lb, lensn, lprst, lsn, nsn, &
                                 get_cssn, get_seasons
  Use pcs,                 Only: get_eofs_settings
  Use missing,             Only: kavx, kavy, kavz, xmiss, ymiss
  Use analysis,            Only: ianal, icalc
  Use cpt_view,            Only: set_menuflags
  Use cpt_actions,         Only: reset
!
! Function type
  Integer :: get_infile
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cxy ! - X/Y flag -
!
  Logical, Intent(In) :: lseason ! - prompt for season? -
!
! Procedure arguments
  Integer, External :: get_afile
!
! - optional procedure arguments -
  Integer, External, Optional :: get_eof_settings
!
! Locals
!
! Local scalars
  Integer :: nv_old ! - backup number of variables -
  Integer :: ifail  ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Check for completed calculations
  get_infile = 2
#if GUI == 1
  nclr = 3
#endif
  Select Case (cxy)
   Case ('X', 'x', 'Y', 'y')
     If (reset('Opening a new '//cxy//' file', .true.) == 1) Return
  End Select
!
! Open file
  Select Case (cxy)
   Case ('X', 'x', '1')
     nv_old = nx
     bkfile = xfile
   Case ('Y', 'y', '2')
     nv_old = ny
     bkfile = yfile
     bk2file = mfile
     Call init_ifile (mfile)
     If (xfile%iseq == isq_sn) Then ! - backup season settings -
        ims_bk(:) = ims(:)
        mpsn_bk(:) = mpsn(:)
     End If
   Case ('Z', 'z', '3')
     nv_old = nz
     bkfile = zfile
  End Select
  lssq = (ianal == ia_pfv)
  If (get_afile() /= 0) GoTo 1
!
! Prompt for season
  Select Case (cxy)
   Case ('Y', 'y', '2')
     If (xfile%iseq == isq_sn) Then ! - backup season settings -
        ims(:) = ims_bk(:)
        mpsn(:) = mpsn_bk(:)
     End If
  End Select
  If (lseason) Then
     Call set_cursor_waiting (0)
     mfield(1)%tprd = yfield(1)%tprd
     Call get_seasons (ianal, xfield, yfield, ifail)
     If (ifail /= 0) GoTo 1
     Call set_cursor_waiting (1)
  End If
!
! Set default length of training period
  Select Case (cxy)
   Case ('X', 'x', 'Y', 'y', '1', '2')
     Call get_maxnt (icsv, lprst, lensn, lsn(lb:), lb, nsn, nt, ntm)
     Call window_update (nt)
  End Select
!
! Reformat monthly files if merge or update are used
  Select Case (cxy)
   Case ('1')
     If (xfile%iseq == isq_mn) Call reformat_monthly (xfile%nt, xfile%ntm, xfile%nfl, kavx, kavy, ifail)
   Case ('2')
     If (yfile%iseq == isq_mn) Call reformat_monthly (yfile%nt, yfile%ntm, yfile%nfl, kavy, kavz, ifail)
   Case ('3')
     If (zfile%iseq == isq_mn) Call reformat_monthly (zfile%nt, zfile%ntm, zfile%nfl, kavz, kavy, ifail)
     zfile%cssn = get_cssn(zfile%prd1, zfile%iseq, nsn)
  End Select
!
! Prompt for domain
  Select Case (cxy)
   Case ('X', 'x', '1')
     Select Case (xfile%idstr)
      Case (1, 2)
        lmsq = .false.
        lssq = (ianal == ia_pfv)
     End Select
     Call get_domains (cxy, xfile, xfield, rlatx, rlngx, cstnx, nv, rlatdx, rlngdx, cstndx, idomx, ifail)
   Case ('Y', 'y', '2')
     Select Case (yfile%idstr)
       Case (1, 2)
        lmsq = .true.
        lssq = (ianal == ia_pfv)
     End Select
     Call get_domains (cxy, yfile, yfield, rlaty, rlngy, cstny, nv, rlatdy, rlngdy, cstndy, idomy, ifail)
     Call set_mfields ()
   Case ('Z', 'z', '3')
     lmsq = .false.
     lssq = (ianal == ia_pfv)
     ifail = 0 ! - ignore forecast file domain, which must be the same as the X file -
  End Select
  If (ifail /= 0) GoTo 1
!
! Prompt for mode settings
! - prompt for EOF mode and CCA mode settings -
  If (Present(get_eof_settings)) Then
     If (xfile%cfile%lset .and. yfile%cfile%lset) Then
        get_infile = get_eofs_settings(cxy, ianal, .true., get_eof_settings)
! - prompt for mode settings -
     Else
        get_infile = get_eofs_settings(cxy, ianal, .false., get_eof_settings)
     End If
  End If
!
! Update availability of menu items depending on which files are open -
  Select Case (cxy)
   Case ('X', 'x')
     Call reset_zfile ()
     If (xfile%cfile%ffile /= bkfile%cfile%ffile .or. xfile%mtime /= bkfile%mtime) ifail = record_change()
     Call set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, sw_elev%on, &
          ixdstr=xfile%idstr, iydstr=yfile%idstr, izdstr=zfile%idstr, iretro=iretro)
   Case ('Y', 'y')
     If (yfile%cfile%ffile /= bkfile%cfile%ffile .or. yfile%mtime /= bkfile%mtime) ifail = record_change()
     Call set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, sw_elev%on, &
          ixdstr=xfile%idstr, iydstr=yfile%idstr, izdstr=zfile%idstr, iretro=iretro)
   Case ('Z', 'z')
     If (zfile%cfile%ffile /= bkfile%cfile%ffile .or. zfile%mtime /= bkfile%mtime) ifail = record_change()
     Call set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, sw_elev%on, &
          izdstr=zfile%idstr, iretro=iretro)
  End Select
!
  get_infile = 0
  GoTo 2
!
! Errors
! - reset old file -
1 Select Case (cxy)
   Case ('X', 'x', '1')
     Call reset_oldfile (lmsq, lssq, yfile%cfile%lset, xfile, bkfile, xfield, xmiss, &
          rlatx, rlngx, cstnx, rlatdx, rlngdx, cstndx, idomx, kavx)
     nx = nv_old
   Case ('Y', 'y', '2')
     Call reset_oldfile (lmsq, lssq, xfile%cfile%lset, yfile, bkfile, yfield, ymiss, &
          rlaty, rlngy, cstny, rlatdy, rlngdy, cstndy, idomy, kavy)
     mfile = bk2file
     ny = nv_old
   Case ('Z', 'z')
     Call reset_oldfile (lmsq, lssq, .true., zfile, bkfile, zfield, xmiss, &
          rlatz, rlngz, cstnz, rlatdz, rlngdz, cstndz, idomz, kavz)
     nz = nv_old
   Case ('3')
     zfile = bkfile
  End Select
  If (.not.xfile%cfile%lset .and. .not.yfile%cfile%lset) inew = 1
  get_infile = 1
!
! Refresh files
2 Select Case (cxy)
   Case ('X', 'x')
     Call file_refresh (xfile, &
          nv=nx)
     If (yfile%cfile%lset) Call file_refresh (yfile) ! - update Y-file in case of modified start date -
   Case ('Y', 'y')
     Call file_refresh (yfile, &
          nv=ny)
     Call file_refresh (mfile)
     If (xfile%cfile%lset) Call file_refresh (xfile) ! - update X-file in case of modified start date -
   Case ('Z', 'z')
     Call file_refresh (zfile, &
          nv=nz)
   Case ('1')
     Call file_refresh (xfile)
   Case ('2')
     Call file_refresh (yfile)
   Case ('3')
     Call file_refresh (zfile)
  End Select
  Deallocate (domain_bk)
! - switch off cursor waiting -
  Call set_cursor_waiting (0)
!
  Return
!
 Contains
!
!
  Subroutine reformat_monthly (nt, ntm, nfl, kava, kavb, ifail)
!
! Reformats monthly files if merge or update are used
!
! Modules
  Use missing, Only: init_kav
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl ! - total number of fields and lagged fields -
  Integer, Intent(In) :: nt  ! - number of time steps -
  Integer, Intent(In) :: ntm ! - number of time steps (monthly) -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Logical, Dimension(:,:), Pointer :: kava ! - cases available? -
  Logical, Dimension(:,:), Pointer :: kavb ! - cases available? (workspace) -
!
! Locals
!
! Local scalars
  Integer :: k  ! - year index -
  Integer :: l  ! - month index -
  Integer :: kl ! - year and month index -
!
! Executable Statements
!
! Reformats monthly files if merge or update are used
  Call init_kav (ntm, nfl, 1, kavb, ifail)
  kl = 0
  Do k = 1, nt
     Do l = 1, nmn
        kl = kl + 1
        If (kl > nt) Exit
        kavb(kl,1) = kava(k,l)
     End Do
  End Do
  Call init_kav (ntm, nfl, 1, kava, ifail)
  kava(:,:) = kavb(:,:)
!
  Return
  End Subroutine reformat_monthly
 End Function get_infile
!
!
!
 Function get_xfile() &
          Result (get_afile)
!
! Accesses an X file, and determines its structure
!
! Error indicators:
!    get_afile =  0 Successful
!    get_afile =  1 Problem accessing file
!
! Modules
  Use data_time_constants, Only: isq_yr
  Use time,                Only: Operator(+), Operator(<), &
                                 iseq, lmsq, lssq, nss
  Use iofiles,             Only: mfile, nv, nx, xfile, yfile, &
                                 set_mfile
  Use settings,            Only: xfd_old
  Use space,               Only: cstndx, cstnx, idomx, rlatdx, rlatx, rlngdx, rlngx
  Use fields,              Only: mfield, xfield, yfield, &
                                 init_field
  Use season,              Only: icsv
  Use missing,             Only: kavx, xmiss
!
! Function type
  Integer :: get_afile
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Ubound
!
! Executable Statements
!
! Get X input file
  get_afile = 1
  lmsq = .false.
  nv => nx
  Call get_inputfile ('X', xfile, mfile%cfile%lset, xfield, nv, xmiss, check_xfile, &
       rlatx, rlngx, rlatdx, rlngdx, cstnx, cstndx, idomx, kavx, ifail)
  If (ifail /= 0) Return
  xfd_old = xfile%fdate%iyr
!
! Check for calculation of seasonal values
  If (icsv > 0) Then
     If (xfile%iseq == 0) Then ! - assign X-file sequencing if there is only one case -
        Select Case (yfile%iseq)
         Case (:-1)
           xfile%iseq = isq_yr
         Case Default
            xfile%iseq = yfile%iseq
        End Select
     Else
        yfile%iseq = xfile%iseq
     End If
  Else
     iseq = xfile%iseq
  End If
  Call set_mfile (icsv, .false.)
!
! Reset Y-fields if Y file is set
  If (yfile%cfile%lset) Then
     yfile%bdate = mfile%bdate
     yfile%cprdb = mfile%cprdb
     If (lssq .and. Ubound(yfield, Dim=1) /= nss) Then
        Call init_field (yfield, nss, mfield(1:1)%rmiss, ifail)
        If (ifail /= 0) Return  
        yfield(:) = mfield(1)
        If (nss > 0) yfield(2:)%nv = 0
     End If
  End If
!
! Set start dates if Y-file is open -
  If (mfile%cfile%lset) Call set_fdates ('X', xfile, mfile)
!
  get_afile = 0
!
  Return
 End Function get_xfile
!
!
!
 Function get_yfile() &
          Result (get_afile)
!
! Accesses a Y file, and determines its structure
!
! Error indicators:
!    get_afile =  0 Successful
!    get_afile =  1 Problem accessing file
!
! Modules
  Use data_time_constants, Only: isq_yr
  Use labels,              Only: ca_dstruct
  Use time,                Only: Operator(>=), Operator(<=), &
                                 iseq, lmsq
  Use iofiles,             Only: nv, ny, xfile, yfile, &
                                 set_mfile
  Use settings,            Only: izero, yfd_old
  Use space,               Only: cstndy, cstny, idomy, rlatdy, rlaty, rlngdy, rlngy
  Use fields,              Only: dsdy, yfield, &
                                 set_yfields
  Use season,              Only: icsv, lb
  Use missing,             Only: kavy, ymiss
  Use climate,             Only: iclim, climate_bak, climate_per, &
                                 init_climatep
!
! Function type
  Integer :: get_afile
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
  Intrinsic All
!
! Executable Statements
!
! Get Y input file
  get_afile = 1
  lmsq = .true.
  nv => ny
  Call get_inputfile ('Y', yfile, xfile%cfile%lset, yfield, nv, ymiss, check_yfile, &
       rlaty, rlngy, rlatdy, rlngdy, cstny, cstndy, idomy, kavy, ifail)
  If (ifail /= 0) Return
  yfd_old = yfile%fdate%iyr
!
! Check for calculation of seasonal values
  Call set_mfile (icsv, .true.)
  If (icsv > 0) Then
     yfile%iseq = xfile%iseq
  Else
     If (yfile%iseq /= 0) Then
        iseq = yfile%iseq
     Else
        iseq = xfile%iseq
     End If
  End If
!
! Set start dates if X-file is open -
  If (xfile%cfile%lset) Call set_fdates ('Y', yfile, xfile)
!
! Set default zero-bound setting
  Call set_yfields (lb, yfile%nse)
!
! Get data structure description
  dsdy = ca_dstruct(yfile%idstr)
!
! Set a guess at zero-bound
  If (yfile%nfs == 1) Then
     Select Case (yfield(1)%var)
      Case ('precipitation', 'prcp', 'prec', 'pre', 'rain', 'rainfall')
        izero = 1
      Case Default
        izero = 0
     End Select
  Else If (All(yfield(2:)%var == yfield(1)%var)) Then
     Select Case (yfield(1)%var)
      Case ('precipitation', 'prcp', 'prec', 'pre', 'rain', 'rainfall')
        izero = 1
      Case Default
        izero = 0
     End Select
  Else
     izero = 0
  End If
!
! Reset climatology
  climate_bak = climate_per
  If (climate_bak%d1 >= yfile%prd1%sdate .and. climate_bak%d2 <= yfile%prdn%sdate) Then
! - set unused month and day for yearly sequencing -
     If (yfile%iseq == isq_yr) Then
        climate_per%d1%imn = yfile%prd1%sdate%imn
        climate_per%d1%idy = 0
        climate_per%d2%imn = yfile%prd1%sdate%imn
        climate_per%d2%idy = 0
     End If
  Else
     Call init_climatep ()
  End If
  iclim = 1
!
  get_afile = 0
!
  Return
 End Function get_yfile
!
!
!
 Function get_zfile() &
          Result (get_afile)
!
! Accesses a forecast / update file, and determines its structure
!
! Error indicators:
!    get_afile =  0 Successful
!    get_afile =  1 Problem accessing file
!
! Modules
  Use time,     Only: Operator(+), Operator(<), &
                      lmsq, &
                      date_diff
  Use iofiles,  Only: mza, nv, nz, xfile, zfile
  Use settings, Only: nf
  Use space,    Only: cstndz, cstnz, idomz, rlatdz, rlatz, rlngdz, rlngz, &
                      nullify_grids
  Use fields,   Only: xfield, zfield
  Use missing,  Only: kavz, xmiss
!
! Function type
  Integer :: get_afile
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Get forecast input file
  nv => nz
  lmsq = .false.
  Call get_inputfile ('Z', zfile, xfile%cfile%lset, zfield, nv, xmiss, check_zfile, &
       rlatz, rlngz, rlatdz, rlngdz, cstnz, cstndz, idomz, kavz, ifail)
  get_afile = ifail
  Select Case (ifail)
   Case (0) ! - no errors -
     If (zfile%iseq == 0) zfile%iseq = xfile%iseq
     If (xfile%cfile%lset) Then
        Call set_fdates ('Z', zfile, xfile) ! - set start-at dates -
        If (zfile%cfile%ffile == xfile%cfile%ffile) Then
           nf = date_diff(zfile%fdate, zfile%prdn%sdate, zfile%iseq) + 1
        Else
           nf = zfile%nt
        End If
     End If
     zfield(:)%region = xfield(:)%region ! - copy domain from X -
   Case (2) ! - information unset -
     Continue
   Case (1) ! - errors -
     Return
  End Select
  mza = 0
! - free unnecessary array memory -
  Call nullify_grids (rlatz, rlngz, cstnz, rlatdz, rlngdz, cstndz, idomz)
!
  Return
 End Function get_zfile
!
!
!
 Subroutine get_inputfile (cxy, afile, lbfile, afield, nv, amiss, check_file, &
            rlat, rlng, rlatd, rlngd, cstn, cstnd, idom, kav, ifail)
!
! Accesses an input file, and determines its structure
!
! Error indicators:
!    ifail =  0 Successful
!    ifail =  1 Problem getting file
!    ifail =  2 Successful, but location information not set
!
! Modules
  Use data_numbers,      Only: rp
  Use data_io_constants, Only: lstn
  Use time,              Only: Operator(<), Operator(>)
  Use iofiles,           Only: ifile, &
                               bkfile
  Use fields,            Only: field
  Use season,            Only: lprst
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In) :: amiss ! - default missing value flag -
!
  Logical, Intent(In) :: lbfile ! - other file open? -
!
  Character(Len=1), Intent(In) :: cxy ! - X / Y variables flag -
!
! Pointer scalars
  Integer, Pointer :: nv ! - number of variables -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - input file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Integer, Dimension(:,:), Pointer :: idom ! - variable is within domain -
!
  Real(Kind=rp), Dimension(:,:), Pointer :: rlat  ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlng  ! - longitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlatd ! - used latitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlngd ! - realigned/used longitudes -
!
  Character(Len=lstn), Dimension(:,:), Pointer :: cstn  ! - names of stations -
  Character(Len=lstn), Dimension(:,:), Pointer :: cstnd ! - names of stations within domains -
!
  Logical, Dimension(:,:), Pointer :: kav ! - cases available? -
!
  Type(field), Dimension(:), Pointer :: afield ! - field information -
!
! Procedure arguments
  Integer, External :: check_file ! - file checking function -
!
! Locals
!
! Local scalars
  Logical :: lfile ! - prompt for file? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim, Size
!
! Executable Statements
!
! Backup fields
  Allocate (domain_bk(Size(afield)))
  domain_bk(:) = afield(:)%region
!
! Get X input file
  If (Len_Trim(afile%cfile%desc) == 0) Then
     Select Case (cxy)
      Case ('X', 'x', 'Y', 'y')
        afile%cfile%desc = cxy//' Input File'
      Case ('Z', 'z')
        afile%cfile%desc = 'Forecast Input File'
     End Select
  End If
  lfile = .true.
1 Call get_ifile (lfile, lprst, afile, lbfile, bkfile, afield, amiss, rlat, rlng, rlatd, rlngd, cstn, cstnd, idom, kav, ifail)
  Select Case (ifail)
   Case (0) ! - no errors -
     ifail = check_file() ! - check for valid input file formats -
     Select Case (ifail)
      Case (0)
        Continue
      Case (-1)
        lfile = .false.
        GoTo 1 ! - re-read the file -
      Case Default
        GoTo 1 ! - check for valid input file formats -
     End Select
   Case (1) ! - no new file -
     ifail = 0
   Case (2) ! - errors -
     GoTo 2
  End Select
!
! Retain prior start date unless it is invalid
  If (afile%fdate < afile%prd1%sdate .or. afile%fdate > afile%prdn%sdate) afile%fdate = afile%prd1%sdate
!
! Update screen
  Call file_refresh (afile, &
       nv=nv)
!
  Return
!
! Errors
2 ifail = 1
!
  Return
 End Subroutine get_inputfile
!
!
!
 Subroutine get_ifile (lfile, lprsst, afile, lbfile, bkfile, afield, amiss, rlat, rlng, rlatd, rlngd, cstn, cstnd, idom, kav, ifail)
!
! Accesses an input file, and determines its structure
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Successful but file is unchanged
!    ifail =  2 Unsuccessful
!
! Modules
  Use data_numbers,        Only: rp
  Use data_io_constants,   Only: lstn
  Use data_time_constants, Only: isq_mn, isq_sn, isq_so, isq_yr
  Use time,                Only: pprd, &
                                 Operator(-), &
                                 iseq, lmsq, lssq, &
                                 add_to_month, get_cdate, reset_time
  Use gui,                 Only: set_cursor_waiting
  Use iofiles,             Only: ifile, &
                                 file_mtime, get_old_file
  Use fields,              Only: field
  Use season,              Only: icsv, lsn, &
                                 get_cssn
  Use file_structure,      Only: get_structure
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In) :: amiss ! - default missing value flag -
!
  Logical, Intent(In) :: lfile  ! - prompt for file? -
  Logical, Intent(In) :: lbfile ! - other file open? -
  Logical, Intent(In) :: lprsst ! - persistence component? -
!
  Type(ifile), Intent(In) :: bkfile ! - backup file -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - input file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Integer, Dimension(:,:), Pointer :: idom ! - variable is within domain -
!
  Real(Kind=rp), Dimension(:,:), Pointer :: rlat  ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlng  ! - longitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlatd ! - used latitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlngd ! - realigned/used longitudes -
!
  Character(Len=lstn), Dimension(:,:), Pointer :: cstn  ! - names of stations -
  Character(Len=lstn), Dimension(:,:), Pointer :: cstnd ! - names of stations within domains -
!
  Logical, Dimension(:,:), Pointer :: kav ! - cases available? -
!
  Type(field), Dimension(:), Pointer :: afield ! - field information -
!
! Locals
!
! Local scalars
  Integer :: l ! - field index -
!
  Type(pprd) :: apprd ! - period -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Initialise time
  If (.not.lbfile) Call reset_time (.true.)
!
! Select input file
  If (lfile) Then
     Call get_old_file (Trim(afile%cfile%desc), afile%cfile%ffile, afile%cfile%fdir, (/afile%cfile%desc/), (/'*.*'/), 1, &
          afile%cfile%fname, ifail)
     If (ifail /= 0) GoTo 1
     Call set_cursor_waiting (1)
!
! Check modification date
     afile%mtime = file_mtime(afile%cfile%ffile)
!
! Check whether file is new
     If (afile%cfile%ffile == bkfile%cfile%ffile .and. afile%mtime == bkfile%mtime .and. afile%mtime /= 0) Then
        If (.not.lbfile) Then
           iseq = afile%iseq
           icsv = 0
        End If
        ifail = 1
        Return
     End If
  End If
!
! Determine dataset structure
  Call get_structure (lmsq, lssq, lbfile, afile, afield, amiss, rlat, rlng, cstn, rlatd, rlngd, cstnd, idom, kav, ifail)
  If (ifail /= 0) GoTo 1
  afile%cfile%lset = .true.
!
! Set time sequence
  If (.not.lbfile) Then
     iseq = afile%iseq
     icsv = 0
  End If
!
! Update file settings
  afile%cfile%lset = .true.
  If (icsv == 0) lsn(1) = 1
  Select Case (iseq)
   Case Default
     If (lprsst) Then
        apprd = afile%prd1
        apprd%sdate%imn = add_to_month(apprd%sdate%imn, -lsn(0), &
                          iyr=apprd%sdate%iyr)
        afile%cprd1 = get_cdate(apprd, lsn(1), 2)
     Else
        afile%cprd1 = get_cdate(afile%prd1, lsn(1), 2)
     End If
     afile%cprdn = get_cdate(afile%prdn, lsn(1), 2)
   Case (isq_sn, isq_so)
     afile%cprd1 = get_cdate(afield(1)%tprd, 1, 2)
     afile%cprdn = get_cdate(afile%prdn, 1, 2)
  End Select
! - set season indicator if file version is CPT10+ -
  If (afile%cfile%ffmt%iver >= 10) Then
     Do l = 1, afile%nfl
        afield(l)%cssn = get_cssn(afield(l)%tprd, afile%iseq, 1) ! - REDUNDANT? -
     End Do
  Else
     afield(:)%cssn = 'N/A'
  End If
! - delete any old first date extraneous information -
  Select Case (afile%iseq)
   Case (isq_yr)
     afile%fdate%idy = afile%prd1%sdate%idy
     afile%fdate%imn = afile%prd1%sdate%imn
   Case (isq_mn, isq_sn, isq_so)
     afile%fdate%idy = afile%prd1%sdate%idy
  End Select
!
  ifail = 0
  Return
!
! Errors
1 ifail = 2
!
  Return
 End Subroutine get_ifile
!
!
!
 Function check_xfile() &
          Result (check)
!
! Checks for a valid X input file given current analysis settings
!
! On exit:
!    check =  0 Successful
!    check = -1 Probabilistic file needs to be re-read
!    check =  4 Multi-field or multi-lagged input files cannot be used for GCM
!    check =  5 Insufficient number of grids / stations / indices in file
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_gcm, ia_mlr, ia_pfv
  Use data_io_constants,  Only: lfil
  Use labels,             Only: ca_atypes_a
  Use time,               Only: lssq
  Use errors,             Only: cpt_error, query
  Use iofiles,            Only: mfile, xfile
  Use settings,           Only: mnu, sw_elev
  Use analysis,           Only: ianal, ianaln
  Use cpt_view,           Only: change_view
!
! Function type
  Integer :: check
!
! Locals
!
! Local parameters
  Character(Len=*), Parameter :: cproc = 'check_xfile' ! - procedure -
!
! Local scalars
  Integer :: ifail ! - error indicator -
  Integer :: iarg1 ! - integer error indicator -
!
  Character(Len=lfil) :: carg1 ! - first character error message argument -
  Character(Len=   8) :: carg2 ! - second character error message argument -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Disable gridded and station X input files for MLR
  If (xfile%idstr <= 2 .and. ianal == ia_mlr) Then
     If (sw_elev%on) Then
        carg1 = ca_atypes_a(ianal)
        carg2 = ca_atypes_a(ia_cca)
        ifail = 1
        If (query(cproc, ifail, .true., c_arg1=Trim(carg1), c_arg2=Trim(carg2)) /= 1) GoTo 2
     End If
     ianaln = ia_cca
     check = change_view()
!
! Force probabilistic X files for PFV
  Else If (xfile%ngs == 1 .and. ianal == ia_pfv) Then
     If (sw_elev%on) Then
        carg1 = ca_atypes_a(ianal)
        carg2 = ca_atypes_a(ia_cca)
        ifail = 2
        If (query(cproc, ifail, .true., c_arg1=Trim(carg1), c_arg2=Trim(carg2)) /= 1) GoTo 2
     End If
     ianaln = ia_cca ! - change to CCA if file is not probabilistic -
     check = change_view()
     lssq = .false.
     check = -1 ! - re-read the file -
     Return
!
! Enable probabilistic X files for PFV only
  Else If (xfile%ngs > 1 .and. ianal /= ia_pfv) Then
     If (sw_elev%on) Then
        carg1 = ca_atypes_a(ianal)
        carg2 = ca_atypes_a(ia_pfv)
        ifail = 3
        If (query(cproc, ifail, .true., c_arg1=Trim(carg1), c_arg2=Trim(carg2)) /= 1) GoTo 2
     End If
     ianaln = ia_pfv ! - change to PFV if file is probabilistic -
     check = change_view()
     lssq = .true.
     check = -1 ! - re-read the file -
     Return
  End If
!
! Disable multi-field input files for GCM
  If (xfile%nfl /= xfile%nms .and. ianal == ia_gcm) Then 
     carg1 = ca_atypes_a(ianal)
     ifail = 4
     GoTo 1
  End If
!
! Check for sufficient data
  If (xfile%nt < mnu) Then
     carg1 = xfile%cfile%ffile
     iarg1 = xfile%iseq
     ifail = 5
     GoTo 1
  End If
!
! Perform checks for X and Y files
  check = check_xyfiles('X', xfile, mfile)
  Return
!
! Errors
1 Call cpt_error (cproc, ifail, .false., &
       i_arg1=iarg1, c_arg1=Trim(carg1), c_arg2=Trim(carg2))
2 check = 1
!
  Return
 End Function check_xfile
!
!
!
 Function check_yfile() &
          Result (check)
!
! Checks for a valid Y input file given current analysis settings
!
! On exit:
!    check =  0 Successful
!    check =  1 Multi-model Y input files cannot be used for GCM
!    check =  2 Probabilistic Y input files are not permitted
!    check =  4 Multi-field or multi-lagged input files cannot be used for GCM
!    check =  5 Insufficient number of cases in file
!    check =  6 Insufficient number of cases in Y file to calculate a climatology
!
! Modules
  Use data_cpt_constants, Only: ia_gcm
  Use labels,             Only: ca_atypes_a
  Use errors,             Only: cpt_error
  Use iofiles,            Only: xfile, yfile
  Use settings,           Only: mnu
  Use analysis,           Only: ianal
!
! Function type
  Integer :: check
!
! Locals
!
! Local parameters
  Character(Len=*), Parameter :: cproc = 'check_yfile'
!
! Local scalars
  Integer :: ifail ! - error indicator -
  Integer :: iarg1 ! - integer error indicator -
!
  Character(Len=8) :: carg ! - character error message argument -
!
! Executable Statements
!
! Check for valid input file formats
  If (ianal == ia_gcm) Then
     If (yfile%nms > 1) Then ! - disable multi-model Y input files for GCM -
        carg = ca_atypes_a(ianal)
        ifail = 1
        GoTo 1
     End If
     If (yfile%nfl /= yfile%nms) Then ! - disable multi-field input files for GCM -
        carg = ca_atypes_a(ianal)
        ifail = 4
        GoTo 1
     End If
  End If
  If (yfile%ngs > 1) Then ! - disable probabilistic Y files -
     carg = ' '
     ifail = 2
     GoTo 1
  End If
!
! Check for sufficient data
  If (yfile%nt < mnu) Then
     iarg1 = yfile%iseq
     If (yfile%nt == 1 .and. mnu == 2) Then
        ifail = 6
     Else
        ifail = 5
     End If
     carg = ' '
     GoTo 1
  End If
!
! Perform checks for X and Y files
  check = check_xyfiles('Y', xfile, yfile)
  Return
!
! Errors
1 Call cpt_error (cproc, ifail, .false., &
       i_arg1=iarg1, c_arg1=Trim(carg))
  check = ifail
!
  Return
 End Function check_yfile
!
!
!
 Function check_xyfiles(cxy, xfile, yfile) &
          Result (check)
!
! Checks for consistency between X and Y input files given current analysis settings
!
! Modules
  Use data_cpt_constants,  Only: ia_pfv
  Use data_io_constants,   Only: lfil
  Use data_time_constants, Only: isq_mn
  Use time,                Only: iseq
  Use errors,              Only: cpt_error
  Use iofiles,             Only: ifile
  Use fields,              Only: yfield
  Use season,              Only: icsv, &
                                 set_icsv
  Use analysis,            Only: ianal
!
! Function type
  Integer :: check
!
! Arguments
!
! Input scalars
  Integer :: ifail ! - error indicator -
!
  Character(Len=1), Intent(In) :: cxy ! - X / Y variables flag -
!
  Type(ifile), Intent(In) :: xfile ! - X file -
  Type(ifile), Intent(In) :: yfile ! - Y file -
!
! Locals
!
! Local parameters
  Character(Len=*), Parameter :: cproc = 'check_xyfiles'
!
! Local scalars
  Integer :: iarg1 ! - integer error indicator -
  Integer :: iarg2 ! - integer error indicator -
!
  Character(Len=lfil) :: carg3 ! - character error message argument -
  Character(Len=lfil) :: carg4 ! - character error message argument -
!
! Executable Statements
!
! Initialise routine
  carg3 = ' '
  carg4 = ' '
  iarg1 = 0
  iarg2 = 0
!
! Check for compatibility of X and Y files for PFV
  If (xfile%cfile%lset .and. yfile%cfile%lset) Then
     If (ianal == ia_pfv) Then
        ifail = 0
        Select Case (xfile%idstr)
         Case (1) ! - gridded forecasts -
           If (yfile%idstr == 3) ifail = 1 ! - gridded or station observations only -
         Case (2) ! - station forecasts -
           If (yfile%idstr /= 2) ifail = 1 ! - station observations only -
         Case (3) ! - unreferenced forecasts -
           If (yfile%idstr /= 3) ifail = 1 ! - unreferenced observations only -
        End Select
        If (ifail /= 0) GoTo 1
     End If
!
! Check for consistent file sequences
! - permit monthly sequenced Y-file -
     If (xfile%iseq /= yfile%iseq .and. (xfile%iseq /= 0 .and. yfile%iseq /= 0)) Then
        If (yfile%iseq == isq_mn) Call set_icsv (icsv, yfield(1)%unit)
        Select Case (cxy) 
         Case ('X', 'x')
           If (icsv /= 0) Then
              If (xfile%iseq /= 0) Then
                 iseq = xfile%iseq
              Else
                 iseq = yfile%iseq
              End If
           End If
           If (icsv == 0) Then
              carg3 = xfile%cfile%ffile
              carg4 = yfile%cfile%ffile
              iarg1 = xfile%iseq
              iarg2 = yfile%iseq
              ifail = 2
              GoTo 1
           End If
         Case ('Y', 'y')
           If (xfile%iseq /= 0) Then
              iseq = xfile%iseq
           Else
              iseq = yfile%iseq
           End If
           If (icsv == 0) Then
              carg3 = yfile%cfile%ffile
              carg4 = xfile%cfile%ffile
              iarg1 = yfile%iseq
              iarg2 = xfile%iseq
              ifail = 2
              GoTo 1
           End If
        End Select
     Else If (xfile%iseq == 0 .and. yfile%iseq /= 0) Then
        iseq = yfile%iseq
     Else If (xfile%iseq /= 0 .and. yfile%iseq == 0) Then
        iseq = xfile%iseq
     Else
        If (yfile%iseq == isq_mn) Then
           Call set_icsv (icsv, yfield(1)%unit)
        Else
           icsv = 0
        End If
     End If
  Else
     icsv = 0
  End If
!
  check = 0
  Return
!
! Errors
1 Select Case (cxy)
   Case ('X', 'x')
     Call cpt_error (cproc, ifail, .false., &
          c_arg1='X', c_arg2='Y', c_arg3=carg3, c_arg4=carg4, i_arg1=iarg1, i_arg2=iarg2)
   Case ('Y', 'y')
     Call cpt_error (cproc, ifail, .false., &
          c_arg1='Y', c_arg2='X', c_arg3=carg3, c_arg4=carg4, i_arg1=iarg1, i_arg2=iarg2)
  End Select
  check = 1
!
  Return
 End Function check_xyfiles
!
!
!
 Function check_zfile() &
          Result (check)
!
! Checks for consistency between forecast and X files
!
! Error indicators:
!    ifail =  0 Successful
!    ifail =  1 Problem allocating memory
!    ifail =  2 File version is inconsistent with corresponding X file
!    ifail =  3 File structure is inconsistent with corresponding X file
!    4  <=  ifail  <=  8 Fields are inconsistent with corresponding X file
!
! Modules
  Use data_io_constants,   Only: lfil
  Use data_time_constants, Only: isq_yr
  Use time,                Only: Operator(+), Operator(<)
  Use errors,              Only: cpt_error
  Use iofiles,             Only: xfile, zfile
  Use space,               Only: rlatx, rlatz, rlngx, rlngz
  Use fields,              Only: xfield, zfield
!
! Function type
  Integer :: check
!
! Locals
!
! Local parameters
  Character(Len=*), Parameter :: cproc = 'check_zfile' ! - procedure -
!
! Local scalars
  Integer :: l     ! - field / lagged-field index -
  Integer :: iarg1 ! - error optional argument -
  Integer :: iarg2 ! - error optional argument -
  Integer :: ifd   ! - current field -
  Integer :: ilf   ! - current lagged-field -
  Integer :: ifail ! - error indicator -
!
  Character(Len=lfil) :: carg2 ! - error optional argument -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any, Trim
!
! Executable Statements
!
! Initialise routine
  carg2 = ' '
  iarg1 = 1
  iarg2 = 1
!
! Check file versions
  If (zfile%cfile%ffmt%iver /= xfile%cfile%ffmt%iver) Then
     ifail = 2
     carg2 = xfile%cfile%ffile
     GoTo 1
  End If
!
! Check for consistency with xfile
  If (zfile%idstr /= xfile%idstr) Then
     ifail = 3
     iarg1 = zfile%idstr
     iarg2 = xfile%idstr
     GoTo 1
  End If
  iarg1 = 0
  iarg2 = 0
!
! Check for consistency of field settings
  Do ifd = 1, zfile%nfs
     If (zfile%nfs > 1) iarg1 = ifd
     Do ilf = 1, zfile%nls
        If (zfile%nls > 1) iarg2 = ilf
        l = (ifd - 1)*zfile%nls + ilf
! - permit inconsistent precipitation abbreviations -
        If (zfield(l)%var /= xfield(l)%var) Then
           If (Trim(zfield(l)%var) == 'prec' .and. Trim(xfield(l)%var) == 'prcp') Then
              zfield(l)%var = 'prcp'
           Else If (Trim(zfield(l)%var) == 'prcp' .and. Trim(xfield(l)%var) == 'prec') Then
              zfield(l)%var = 'prec'
           Else
              ifail = 4
              carg2 = 'field'
              GoTo 1
           End If
        End If
! - check for consistent units -
        If (zfield(l)%unit /= xfield(l)%unit) Then
           ifail = 4
           carg2 = 'unit'
           GoTo 1
        End If
! - check for consistent domains -
        If (zfield(l)%nlt /= xfield(l)%nlt) Then
           ifail = 5
           carg2 = 'latitudes'
           GoTo 1
        End If
        If (zfield(l)%nlg /= xfield(l)%nlg) Then
           ifail = 5
           carg2 = 'longitudes'
           GoTo 1
        End If
        If (zfield(l)%ln2s .neqv. xfield(l)%ln2s) Then
           ifail = 6
           GoTo 1
        End If
        If (xfile%igeog == 0) Cycle
        If (Any(rlatz(:,l) /= rlatx(:,l))) Then
           ifail = 7
           carg2 = 'Latitudes'
           GoTo 1
        End If
        If (Any(rlngz(:,l) /= rlngx(:,l))) Then
           ifail = 7
           carg2 = 'Longitudes'
           GoTo 1
        End If
     End Do
  End Do
! - check for consistent seasons -
  If (zfile%cssn /= xfile%cssn) Then
     Select Case (zfile%iseq)
      Case (isq_yr)
        ifail = 8
        GoTo 1
     End Select
  End If
!
  check = 0
!
  Return
!
! Errors
1 Call cpt_error (cproc, ifail, .false., &
       c_arg1=Trim(zfile%cfile%ffile), c_arg2=Trim(carg2), i_arg1=iarg1, i_arg2=iarg2)
  check = ifail
!
  Return
 End Function check_zfile
!
!
!
 Subroutine reset_zfile ()
!
! Resets forecast file after opening a new X file
!
! Modules
  Use time,     Only: Operator(+), Operator(<), &
                      nss, &
                      date_diff, set_bdate
  Use screen,   Only: window_update
  Use iofiles,  Only: bkfile, idir, nx, nz, xfile, zfile, &
                      init_ifile
  Use settings, Only: nf, nt, xfd_old
  Use fields,   Only: larea, xfield, zfield, &
                      init_field
  Use season,   Only: lsn, nsn
  Use missing,  Only: kavx, kavz, xmiss, &
                      init_kav
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Reset forecast file
  ifail = 0  
  If (xfile%cfile%ffile /= bkfile%cfile%ffile) Then
     xfd_old = xfile%fdate%iyr
     If (.not.larea) Then
        ifail = 0
        Return
     End If
     zfile = xfile
     zfile%cfile%desc = ' '
! - set default first forecast date to after training period - 
     If (nt > 0 .and. xfile%fdate+(nt-1) < xfile%prdn%sdate) Then
        zfile%fdate = xfile%fdate + nt
        nf = date_diff(zfile%fdate, zfile%prdn%sdate, zfile%iseq) + 1
     Else
        zfile%fdate = zfile%prdn%sdate
        nf = 1
     End If
     Call set_bdate (1, lsn(1:), nsn, xfile%prd1%sdate, zfile%fdate, zfile%bdate, &
          cprdb=zfile%cprdb)
! - reset Z fields -
     Call init_field (zfield, zfile%nfl*zfile%nse, zfield(:)%rmiss, ifail)
     If (ifail == 0) Then
        zfield(:) = xfield(:)
        nz = nx
     Else
        ifail = 1
        Call init_ifile (zfile, &
             dir=idir)
        Call init_field (zfield, 1, (/xmiss/), ifail)
        nz = 0
     End If
  End If
  If (ifail == 0) Then
! - update forecast file settings in case of modified xfile contents -
     If (zfile%cfile%ffile == xfile%cfile%ffile) Then
        zfile%cprd1 = xfile%cprd1
        zfile%cprdn = xfile%cprdn
        zfile%cssn = xfile%cssn
        zfile%prd1 = xfile%prd1
        zfile%prdn = xfile%prdn
        Call init_kav (xfile%nt, xfile%nfl, nss, kavz, ifail)
        If (ifail == 0) Then
           kavz(:,:) = kavx(:,:)
        Else
           ifail = 1
           Call init_ifile (zfile, &
                dir=idir)
           Call init_field (zfield, 1, (/xmiss/), ifail)
           nz = 0
        End If
        ifail = 0
     End If
  End If
  Call file_refresh (zfile, &
       nv=nz)
  Call window_update (nf)
!
  Return
 End Subroutine reset_zfile
!
!
!
 Subroutine files_refresh ()
!
! Refresh input file settings
!
! Modules
  Use iofiles, Only: mxa, mya, mza, nx, ny, nz, xfile, yfile, zfile, mfile
!
! Executable Statements
!
! Refresh file information
  Call file_refresh (xfile, &
       nv=nx, nva=mxa)
  Call file_refresh (yfile, &
       nv=ny, nva=mya)
  Call file_refresh (zfile, &
       nv=nz, nva=mza)
  Call file_refresh (mfile, &
       nv=ny)
!
  Return
 End Subroutine files_refresh
!
!
!
 Subroutine file_refresh (afile, &
            nv, nva)
!
! Refresh file settings
!
! Modules
  Use screen,  Only: window_update
  Use iofiles, Only: ifile
!
! Arguments
!
! Input scalars
  Type(ifile), Intent(In) :: afile ! - file -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: nv  ! - number of variables -
  Integer, Intent(In), Optional :: nva ! - number of used variables -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Refresh file information
  Call window_update (afile%cfile%fname)
  Call window_update (afile%cgss)
  Call window_update (afile%cprd1)
  Call window_update (afile%cprdn)
  Call window_update (afile%cprdb)
  Call window_update (afile%nfs)
  Call window_update (afile%nls)
  If (Present(nv))  Call window_update (nv)
  If (Present(nva)) Call window_update (nva)
!
  Return
 End Subroutine file_refresh
!
!
!
 Subroutine reset_oldfile (lmsq, lssq, lbfile, afile, bkfile, afield, amiss, rlat, rlng, cstn, rlatd, rlngd, cstnd, idom, kav)
!
! Resets a backup file
!
! Arguments
!
! Modules
  Use data_numbers,      Only: rp
  Use data_io_constants, Only: lstn
  Use iofiles,           Only: ifile
  Use space,             Only: nullify_grids
  Use fields,            Only: field, &
                               init_field
  Use file_structure,    Only: get_structure
!
! Input scalars
  Real(Kind=rp), Intent(In) :: amiss ! - missing value flag -
!
  Logical, Intent(In) :: lbfile ! - other file open? -
  Logical, Intent(In) :: lmsq   ! - permit monthly-sequenced files? -
  Logical, Intent(In) :: lssq   ! - permit seasonally-sequenced files? -
!
  Type(ifile), Intent(In) :: bkfile ! - backup file -
!
! Output scalars
  Type(ifile), Intent(Out) :: afile ! - input file -
!
! Pointer arrays
  Integer, Dimension(:,:), Pointer :: idom ! - variables within domain -
!
  Real(Kind=rp), Dimension(:,:), Pointer :: rlat  ! - latitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlng  ! - longitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlatd ! - realigned used latitudes -
  Real(Kind=rp), Dimension(:,:), Pointer :: rlngd ! - realigned used longitudes -
!
  Character(Len=lstn), Dimension(:,:), Pointer :: cstn  ! - names of stations -
  Character(Len=lstn), Dimension(:,:), Pointer :: cstnd ! - names of stations within domain -
!
  Logical, Dimension(:,:), Pointer :: kav ! - cases available? -
!
  Type(field), Dimension(:), Pointer :: afield ! - field -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Associated
!
! Executable Statements
!
! Resets a backup file
  afile = bkfile
  If (afile%cfile%lset) Then
     Call get_structure (lmsq, lssq, lbfile, afile, afield, amiss, rlat, rlng, cstn, rlatd, rlngd, cstnd, idom, kav, ifail)
     If (ifail == 0) Then
        afield(:)%region = domain_bk(:) ! - restore domain settings -
     Else
        Deallocate (afield)
     End If
  Else
     Call nullify_grids (rlat, rlng, cstn, rlatd, rlngd, cstnd, idom)
  End If
  If (.not.Associated(afield)) Call init_field (afield, 1, (/amiss/), ifail) ! - reset fields if file is opened unsuccessfully -
!
  Return
 End Subroutine reset_oldfile
!
!
!
 Subroutine set_fdates (cxy, afile, bfile)
!
! Sets first dates of interest
!
! Modules
  Use data_cpt_constants,  Only: ia_pfv
  Use data_time_constants, Only: isq_mn, isq_sn
  Use time,                Only: pdate, &
                                 Operator(<), Operator(<=), Operator(>), Operator(>=), Operator(==), Operator(+), Operator(-), &
                                 iafter, iseq, &
                                 date_diff, get_cdate, reset_iseq, set_fdate
  Use iofiles,             Only: ifile, &
                                 mfile, yfile
  Use settings,            Only: nt
  Use season,              Only: lprst, lsn
  Use analysis,            Only: ianal
!
! Arguments
!
! Input scalars
  Character(Len=1), Intent(In) :: cxy ! - X / Y variables flag -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - X input file -
  Type(ifile), Intent(InOut) :: bfile ! - Y/M input file -
!
! Locals
!
! Local scalars
  Integer :: lag    ! - lag -
  Integer :: lagspi ! - lag-adjustment for persistence component of SPI -
  Integer :: ifail  ! - error indicator -
!
  Type(pdate) :: xbdate ! - X-file beginning date -
  Type(pdate) :: xfdate ! - X-file start date -
  Type(pdate) :: ybdate ! - Y-file beginning date -
  Type(pdate) :: yfdate ! - Y-file start date -
!
! Executable Statements
!
! Set default first date of interest for X file
  Select Case (cxy)
   Case ('X', 'x', '1')
     Call init_bfdates (afile, bfile) ! - identify beginning and start dates -
! - Y's start date is within X's data range -
     If (date_diff(afile%prd1%sdate, yfdate, iseq) >= 0 .and. yfdate <= afile%prdn%sdate) Then
        Call set_fdate (-iafter, afile%iseq, yfdate, afile%prd1%sdate, xbdate, xfdate, ifail)
        If (ifail /= 0) Call set_fdate (iafter, afile%iseq, xfdate, mfile%prd1%sdate, ybdate, yfdate, ifail)
! - set start dates to earliest commonly available if there is overlap -
     Else If (afile%prd1%sdate <= bfile%prdn%sdate .and. afile%prdn%sdate >= bfile%prd1%sdate) Then
        If (iseq < 0 .and. iseq /= isq_sn) Then
           Call reset_iseq (isq_new=isq_mn)
        Else
            Call reset_iseq (isq_new=iseq)
        End If
        lag = date_diff(bfile%prd1%sdate, afile%prd1%sdate, iseq)
        Call reset_iseq ()
! - Y file starts before X file -
        If (lag >= 0) Then
           Call set_fdate (iafter, afile%iseq, xfdate, bfile%prd1%sdate, ybdate, yfdate, ifail)
! - X file starts before Y file -
        Else
           yfdate = bfile%prd1%sdate
           ybdate = yfdate
           If (lagspi > 0) Then
              If (iseq < 0 .and. iseq /= isq_sn) Then
                 Call reset_iseq (isq_new=isq_mn)
              Else
                 Call reset_iseq (isq_new=iseq)
              End If
              ybdate = yfdate - lagspi
              Call reset_iseq ()
           End If
           Call set_fdate (-iafter, afile%iseq, yfdate, afile%prd1%sdate, xbdate, xfdate, ifail)
        End If
! - no overlap -
     Else
        xfdate = afile%prd1%sdate
        xbdate = xfdate
     End If
!
! Set default first date of interest for Y file
   Case ('Y', 'y', '2')
     Call init_bfdates (bfile, afile) ! - identify beginning and start dates -
! - X start date is within Y's data range -
     If (date_diff(afile%prd1%sdate, xfdate, iseq) >= -1 .and. xfdate <= afile%prdn%sdate) Then
        Call set_fdate (iafter, bfile%iseq, xfdate, mfile%prd1%sdate, ybdate, yfdate, ifail)
        If (ifail /= 0) Call set_fdate (-iafter, bfile%iseq, yfdate, bfile%prd1%sdate, xbdate, xfdate, ifail)
! - set start dates to earliest commonly available if there is overlap -
     Else If (afile%prd1%sdate <= bfile%prdn%sdate .and. afile%prdn%sdate >= bfile%prd1%sdate) Then
        If (iseq<0 .and. iseq /= isq_sn) Then
           Call reset_iseq (isq_new=isq_mn)
        Else
           Call reset_iseq (isq_new=iseq)
        End If
        lag = date_diff(bfile%prd1%sdate, afile%prd1%sdate, iseq)
        Call reset_iseq ()
! - X file starts before Y file -
        If (lag > 0) Then
           Call set_fdate (-iafter, bfile%iseq, yfdate, bfile%prd1%sdate, xbdate, xfdate, ifail)
! - X file starts before Y file -
        Else
           xfdate = bfile%prd1%sdate
           xbdate = xfdate
           Call set_fdate (iafter, bfile%iseq, xfdate, bfile%prd1%sdate, ybdate, yfdate, ifail)
        End If
! - no overlap -
     Else
        yfdate = afile%prd1%sdate
        ybdate = yfdate
        If (lagspi > 0) Then
           If (iseq < 0 .and. iseq /= isq_sn) Then
              Call reset_iseq (isq_new=isq_mn)
           Else
              Call reset_iseq (isq_new=iseq)
           End If
           ybdate = yfdate - lagspi
           Call reset_iseq ()
        End If
     End If
  End Select
  Select Case (cxy)
   Case ('X', 'x', '1')
     afile%bdate = xbdate
     afile%fdate = xfdate
     bfile%bdate = ybdate
     bfile%fdate = yfdate
   Case ('Y', 'y', '2')
     bfile%bdate = xbdate
     bfile%fdate = xfdate
     afile%bdate = ybdate
     afile%fdate = yfdate
   Case ('Z', 'z', '3')
     If (afile%cfile%ffile == bfile%cfile%ffile) Then
        If (bfile%fdate+(nt-1) < bfile%prdn%sdate) Then
           afile%fdate = bfile%fdate + nt
        Else
           afile%fdate = afile%prdn%sdate
        End If
     Else
        afile%fdate = afile%prd1%sdate
        afile%cprdb = get_cdate(afile%fdate, 2)
     End If
  End Select
! - make sure any changes to Y file start date are applied -
  Select Case (cxy)
   Case ('X', 'x', '1') ! - if current file is X-file, bfile is M-file -
     If (iseq > 0) Then
        yfile%fdate = mfile%fdate
        yfile%bdate = mfile%bdate
     Else
        yfile%fdate%iyr = mfile%fdate%iyr
        yfile%bdate%iyr = mfile%bdate%iyr
        If (ianal == ia_pfv) Then
           yfile%fdate%imn = mfile%fdate%imn
           yfile%bdate%imn = mfile%bdate%imn
        End If
     End If
   Case ('Y', 'y', '2') ! - copy to M-file if current file is Y-file -
     mfile%fdate = yfile%fdate
     mfile%bdate = yfile%bdate
  End Select
!
  Return
!
 Contains
!
!
  Subroutine init_bfdates (filex, filey)
!
! Sets initial estimates for beginning and start at dates
!
! Arguments
!
! Input scalars
  Type(ifile), Intent(In) :: filex ! - X input file -
  Type(ifile), Intent(In) :: filey ! - Y/M input file -
!
! Executable Statements
!
! Set initial estimates for beginning and start at dates
! - set appropriate lag -
  xfdate = filex%fdate
  xbdate = xfdate
  yfdate = filey%fdate
  If (lprst) Then ! - calculate adjustment for persistence component -
     If (iseq < 0 .and. iseq /= isq_sn) Then
        Call reset_iseq (isq_new=isq_mn)
     Else
        Call reset_iseq (isq_new=iseq)
     End If
     lagspi = lsn(0)
     ybdate = yfdate - lagspi
     Call reset_iseq ()
     Do
        If (ybdate >= mfile%prd1%sdate) Exit
        ybdate = ybdate + 1
     End Do
     If (iseq < 0 .and. iseq /= isq_sn) Then
        Call reset_iseq (isq_new=isq_mn)
     Else
        Call reset_iseq (isq_new=iseq)
     End If
     yfdate = ybdate + lagspi
     Call reset_iseq ()
  Else ! - no persistence component -
     lagspi = 0
     ybdate = yfdate
  End If
  End Subroutine init_bfdates
 End Subroutine set_fdates
!
!
!
 Function merge_files()
!
! Merges two input files to produce a multi-field file
!
! Modules
  Use data_time_constants, Only: isq_mn, isq_yr
  Use time,                Only: iafter
  Use screen,              Only: window_update
  Use progress_meter,      Only: lcalc, lstop, &
                                 end_progress_meter, start_progress_meter
  Use errors,              Only: cpt_error
  Use arrays,              Only: iusex, iusey
  Use iofiles,             Only: mfile, nx, ny, xfile, x2file, yfile, &
                                 get_new_file
  Use settings,            Only: isynch, mnu, ntm, nt, &
                                 get_maxnt, set_options
  Use fields,              Only: xfield, yfield
  Use season,              Only: nsn  
  Use missing,             Only: kuse=>kcuse, &
                                 init_missing
  Use cpt_input,           Only: get_data
  Use cpt_output,          Only: write_results
  Use cpt_actions,         Only: reset_mem
!
! Function type
  Integer :: merge_files
!
! Locals
!
! Local scalars
  Integer :: i     ! - index -
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Backup settings
  ifail = set_options(i1=nt, i2=isynch, i3=mnu) ! - backup settings -
  isynch = 1                                    ! - enable synchronous files -
  iafter = 0
  mnu = 1
!
! Prompt for files to merge
  xfile%cfile%desc = 'First file'
  yfile%cfile%desc = 'Second file'
1 Call get_2files ('Merge files', 4, xfile, '1', 5, yfile, '2', ifail)
  If (ifail /= 0) GoTo 2
!
! Check for consistency
! - data structure -
  If (xfile%idstr /= yfile%idstr) Then
     ifail = 1
     Call cpt_error ('merge_files', ifail, .false., &
          c_arg1=xfile%cgss, c_arg2=yfile%cgss)
     GoTo 1
  End If
! - sequencing -
  If (xfile%iseq /= yfile%iseq) Then
     If (yfile%iseq == 0) Then
        yfile%iseq = xfile%iseq
     Else If (xfile%iseq == 0) Then
        xfile%iseq = yfile%iseq
     Else
        ifail = 2
        Call cpt_error ('merge_files', ifail, .false., &
             i_arg1=xfile%iseq, i_arg2=yfile%iseq)
        GoTo 1
     End If
  End If
!
! Reset sequencing if data are monthly
  If (xfile%iseq == isq_mn) Then
     xfile%iseq = isq_yr
     xfile%nt = xfile%ntm
     mfile%iseq = isq_yr
     mfile%nt = mfile%ntm
     yfile%iseq = isq_yr
     yfile%nt = yfile%ntm
  End If
!
! Read data
  Call start_progress_meter ('Reading data', .false.)
  Call get_maxnt (0, .false., 0, (/1/), 1, nsn, nt, ntm) ! - check number of cases -
  lcalc = .false.                                        ! - switch off calculation flag -
  ifail = get_data(0, 0, .false.)                        ! - read data -
  Call end_progress_meter (.false., .false.)
  If (ifail /= 0 .or. lstop) GoTo 2
!
! Prompt for output file
  x2file%cfile%fdir = xfile%cfile%fdir
  x2file%cfile%fname = ' '
  x2file%cfile%fext = xfile%cfile%fext
  x2file%cfile%ffmt = xfile%cfile%ffmt
  Call get_new_file ('Combined file', .true., x2file%cfile%ffile, x2file%cfile%fdir, x2file%cfile%fext, &
       (/'Data Files'/), (/'*.*'/), 1, x2file%cfile%fname, ifail)
  If (ifail /= 0) Then
     merge_files = reset_mem()
     GoTo 1
  End If
  x2file%cfile%lset = .true.
!
! Reset lagged fields to fields
  xfile%nls = 1
  xfile%nfs = xfile%nfl
  yfile%nls = 1
  yfile%nfs = yfile%nfl
!
! Save data
  merge_files = init_missing(.false., ntm, nx, ny)
  xfield(:)%nva = xfield(:)%nv
  yfield(:)%nva = yfield(:)%nv
  iusex(:) = (/(i,i=1,nx)/)
  iusey(:) = (/(i,i=1,ny)/)
  kuse(:) = .true.
  Call write_results ()
!
! Restore settings
2 merge_files = reset_mem()
  ixyz = 1
  merge_files = close_file()
  ixyz = 2
  merge_files = close_file()
  merge_files = set_options(iw=0, i1=nt, i2=isynch, i3=mnu)
  Call window_update (nt)
!
  Return
 End Function merge_files
!
!
!
 Function update_file()
!
! Updates a file using data in a second file
!
! Modules
  Use data_numbers,   Only: rp, one
  Use data_time_constants, Only: isq_mn, isq_yr
  Use time,           Only: Operator(<), Operator(>), &
                            iafter, &
                            date_diff
  Use screen,         Only: window_update
  Use progress_meter, Only: lstop, &
                            end_progress_meter, set_progress_increment, start_progress_meter
  Use errors,         Only: cpt_error
  Use arrays,         Only: iusex, x
  Use iofiles,        Only: nx, ny, xfile, xofile, zfile, &
                            get_new_file
  Use settings,       Only: mnu, ntm, nf, nt, &
                            set_it1, set_options
  Use space,          Only: idomx
  Use fields,         Only: larea, xfield, zfield
  Use missing,        Only: kavx, kavz, kuse=>kcuse, &
                            init_missing
  Use cpt_input,      Only: init_read, num_read, read_data
  Use cpt_output,     Only: write_results
  Use cpt_actions,    Only: reset_mem
!
! Function type
  Integer :: update_file
!
! Locals
!
! Local scalars
  Integer :: i     ! - index -
  Integer :: ifail ! - error indicator -
  Integer :: lgap  ! - gap between primary and update file -
  Integer :: nact  ! - number of lines to read -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Sum
!
! Executable Statements
!
! Backup settings
  ifail = set_options(i1=nt, i2=nf, i3=mnu) ! - backup settings -
  iafter = 0
  mnu = 1
!
! Prompt for files to update
  larea = .false.
  xfile%cfile%desc = 'Primary file'
  zfile%cfile%desc = 'Update file'
1 Call get_2files ('Update files', 4, xfile, '1', 6, zfile, '3', ifail, &
       iopen=xfile%idstr)
  If (ifail /= 0) GoTo 3
!
! Read data
! - check start-at dates -
  xfile%it1 = set_it1(xfile%iseq, xfile%prd1%sdate, xfile%fdate)
  zfile%it1 = set_it1(zfile%iseq, zfile%prd1%sdate, zfile%fdate)
  If (zfile%fdate < xfile%fdate) Then
     ifail = 1
     Call cpt_error ('update_file', ifail, .false.)
     GoTo 1
  End If
! - determine number of cases -
  lgap = date_diff(xfile%prdn%sdate, zfile%fdate, xfile%iseq) ! - identify any gaps -
  Select Case(lgap)
   Case (1) ! - update file immediately follows primary file -
     nt = date_diff(xfile%fdate, xfile%prdn%sdate, xfile%iseq) + 1
     ntm = nt + nf
     lgap = 0
   Case (:0) ! - overlap -
     nt = date_diff(xfile%fdate, zfile%fdate, xfile%iseq)
     ntm = nt + nf
     lgap = 0
   Case Default ! - gap -
     nt = date_diff(xfile%fdate, xfile%prdn%sdate, xfile%iseq) + 1
     ntm = nt + lgap + nf - 1
     x(:,nt+1:nt+lgap,:) = xfield(1)%rmiss
  End Select
  nf = date_diff(zfile%fdate, zfile%prdn%sdate, zfile%iseq) + 1
  nact=num_read(xfile%idstr, xfile%lstack, xfile%nt, nt, Sum(xfield(:)%nlt), xfile%nfs, xfile%nls, xfile%ngs, xfile%it1) + &
       num_read(zfile%idstr, zfile%lstack, zfile%nt, nf, Sum(zfield(:)%nlt), zfile%nfs, zfile%nls, zfile%ngs, zfile%it1)
  Call set_progress_increment (one/Real(nact+1, Kind=rp), .false.)
! - read primary data -
  Call start_progress_meter ('Reading data', .false.)
  If (init_read(nx, 0, 1, ntm, 0, 0, 0) /= 0) Then
     ifail = 2
     Call cpt_error ('update_file', ifail, .false.)
     GoTo 3
  End If
  If (xfile%iseq == isq_mn) xfile%iseq = isq_yr
  Call read_data (xfile, xfield, nt, 1, idomx, kavx, x, ifail)
  If (ifail>0) Then
     ifail = 3
     Call cpt_error ('update_file', ifail, .false., &
          c_arg1 = Trim(xfile%cfile%ffile))
     GoTo 3
  End If
! - read update data -
  nt = nt + lgap
  If (zfile%iseq == isq_mn) zfile%iseq = isq_yr
  zfield(:)%region = xfield(:)%region
  Call read_data (zfile, zfield, nf, 1, idomx, kavz, x(:,nt+1:nt+nf,:), ifail)
  If (ifail > 0) Then
     ifail = 3
     Call cpt_error ('update_file', ifail, .false., &
          c_arg1=Trim(zfile%cfile%ffile))
     GoTo 3
  End If
  Call end_progress_meter (.false., .false.)
  If (ifail /= 0 .or. lstop) GoTo 3
!
! Update X file
  If (zfile%prdn%sdate > xfile%prdn%sdate) Then
     xfile%prdn = zfile%prdn
     xfile%nt = ntm
  End If
!
! Prompt for output file
  xofile%cfile%ffile = xfile%cfile%ffile
  xofile%cfile%fdir = xfile%cfile%fdir
  xofile%cfile%fname = xfile%cfile%fname
  xofile%cfile%fext = xfile%cfile%fext
  xofile%cfile%ffmt = xfile%cfile%ffmt
  nt = ntm
2 Call get_new_file ('Updated file', .true., xofile%cfile%ffile, xofile%cfile%fdir, xofile%cfile%fext, &
       (/'Data Files'/), (/'*.*'/), 1, xofile%cfile%fname,ifail)
  If (ifail /= 0) GoTo 2
  xofile%cfile%lset = .true.
!
! Save data
  update_file=init_missing(.false., ntm, nx, ny)
  iusex(:) = (/(i,i=1,nx)/)
  kuse(:) = .true.
  Call write_results ()
!
! Restore settings
3 update_file = reset_mem()
  ixyz = 1
  update_file = close_file()
  ixyz = 3
  update_file = close_file()
  update_file = set_options(iw=0, i1=nt, i2=nf, i3=mnu)
  Call window_update (nt)
  Call window_update (nf)
  larea = .true.
!
  Return
 End Function update_file
!
!
!
 Subroutine get_2files (ctitle, id1, file1, cfile1, id2, file2, cfile2, ifail, &
            iopen)
!
! Prompts for two input files
!
! Modules
  Use labels,  Only: cg_infiles_t
#if GUI == 1
  Use gui,     Only: iw, &
                     win_prompt
#endif
  Use gui,     Only: box_close, box_open, init_win
  Use iofiles, Only: ifile
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: id1 ! - identifier for first file -
  Integer, Intent(In) :: id2 ! - identifier for second file -
!
  Character(Len=*), Intent(In) :: ctitle ! - title -
  Character(Len=*), Intent(In) :: cfile1 ! - first file -
  Character(Len=*), Intent(In) :: cfile2 ! - second file -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: iopen ! - ready for opening flag -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: file1 ! - first file -
  Type(ifile), Intent(InOut) :: file2 ! - second file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for input files
1 Call init_win (title=ctitle)


  Call box_open (.true., .false., &
       title=Trim(cg_infiles_t), c_opts='bottom_exit')
! - first input file -
  Call box_open (.false., .false.)
  Call input_files (cfile1, Trim(file1%cfile%desc), file1, &
       id_file=id1)
  Call box_close (.false.)
! - second input file -
  Call box_open (.false., .false.)
  Call input_files (cfile2, Trim(file2%cfile%desc), file2, &
       id_file=id2, iopen=iopen)
  Call box_close (.false.)
! - close input box -
  Call box_close (.true.)
! - prompt -
#if GUI == 1
  iw = win_prompt(.true.)
  If (iw == 3) Then
     If (.not.(file1%cfile%lset) .or. .not.(file2%cfile%lset)) GoTo 1
#endif
     file2%bdate=file2%fdate
     ifail = 0
#if GUI == 1
  Else
     ifail = 1
  End If 
#endif
!
  Return
 End Subroutine get_2files
!
!
!
 Subroutine input_files (cxyz, cdesc, afile, &
            iopen, id_file, nv, nva, cprd1, cprdn)
!
! Displays CPT input file information
!
! Modules
#if GUI == 1
  Use data_numbers,       Only: one
#endif
  Use data_cpt_constants, Only: c_file, m_file
  Use labels,             Only: cg_data1, cg_datan, &
                                l_startat
#if GUI == 1
  Use labels,             Only: cg_browse_l, cg_fields_l, cg_filename, cg_nof, cg_nused, &
                                cls_l
  Use screen,             Only: dts_adj, iwdth
  Use gui,                Only: add_button, adjust_y, box_open, gui_creturn, print_ivalue
#endif
  Use gui,                Only: box_close, print_cvalue
  Use iofiles,            Only: ifile
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cxyz  ! - X, Y, or Z file / first or second input file -
  Character(Len=*), Intent(In) :: cdesc ! - file decription -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: id_file ! - file ID -
  Integer, Intent(In), Optional :: iopen   ! - ready for opening flag -
  Integer, Intent(In), Optional :: nv      ! - total number of variables -
  Integer, Intent(In), Optional :: nva     ! - number of used variables -
!
  Character(Len=*), Intent(In), Optional :: cprd1 ! - first period available -
  Character(Len=*), Intent(In), Optional :: cprdn ! - last period available -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - input file -
!
! Locals
!
! Local scalars
  Integer :: itab  ! - tab position -
  Integer :: i_fid ! - input file ID -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
#if GUI == 1
  Intrinsic Min, Nint
#endif
!
! Executable Statements
!
! Identify input file type
  Do i_fid = 1, m_file
     If (cxyz == c_file(i_fid)) Exit
  End Do
!
#if GUI == 1
! Prompt for input file settings
  Call box_open (.true., .false., &
       title=cdesc)
! - prompt for file -
  If (dts_adj > one) Call adjust_y (dts_adj-one)
  Call add_button ('tt', 1, Trim(cg_browse_l), &
       iwdth=5, igrey=iopen, iset1=i_file, ival1=i_fid, cbf=get_file) 
  Call gui_creturn (.false.)
  itab=12
  Call print_cvalue (Trim(cg_filename), afile%cfile%fname, &
       itab=itab, iwid=20, ldepth=.false.)
#else
  If (Present(id_file)) i_file = id_file
  itab = get_file()
  itab = 20
#endif
! - print date of first data -
  If (Present(cprd1)) Then
     Call print_cvalue (Trim(cg_data1), cprd1, &
          itab=itab, iwid=20)
  Else
     Call print_cvalue (Trim(cg_data1), afile%cprd1, &
          itab=itab, iwid=20)
  End If
! - print date of last data -
  If (Present(cprdn)) Then
     Call print_cvalue (Trim(cg_datan), cprdn, &
          itab=itab, iwid=20)
  Else
     Call print_cvalue (Trim(cg_datan), afile%cprdn, &
          itab=itab, iwid=20)
  End If
! - print start date -
  Select Case (cxyz)
   Case ('X', 'x', 'Y', 'y', 'Z', 'z')
     Call print_cvalue (Trim(l_startat%c), afile%cprdb, &
          itab=itab, iwid=20)
  End Select
#if GUI == 1
! - print numbers of fields, lagged fields, variables, and used variables -
  itab = Min(Nint(0.03*Real(iwdth)), 35)
  Call print_ivalue (Trim(cg_nof)//' '//Trim(cg_fields_l), afile%nfs, &
       itab=itab)
  Call print_ivalue (Trim(cg_nof), afile%nls, &
       cqual=cls_l, itab=itab)
  If (Present(nv)) Call print_ivalue (Trim(cg_nof), nv, &
                        cqual=afile%cgss, itab=itab)
  If (Present(nva)) Call print_ivalue (Trim(cg_nused), nva, &
                         itab=itab)
#endif
  Call box_close (.true.)
!
  Return
 End Subroutine input_files
!
!
!
 Function close_file()
!
! Closes an input file
!
! Modules
  Use time,        Only: iseq, &
                         reset_time
  Use screen,      Only: window_update
  Use iofiles,     Only: mxa, mya, mza, nx, ny, nz, &
                         mfile, xfile, yfile, zfile
  Use settings,    Only: igcms, istd, isynch, nf, sw_elev
  Use fields,      Only: nmfd
  Use season,      Only: icsv, &
                         close_season
  Use analysis,    Only: ianal, icalc
  Use cpt_view,    Only: set_menuflags
  Use forecast,    Only: ifc, ifcast, &
                         close_fcast
  Use cpt_actions, Only: reset
!
! Function type
  Integer :: close_file
!
! Locals
!
! Local scalars
  Integer :: ireset ! - reset flag -
!
! Executable Statements
!
! Close input file
  close_file = 1
1 Select Case (ixyz)
! - X file -
   Case (1)
     If (icalc == 1) Then
        ireset = reset('Closing the X file', .true.)
        If (ireset == 1) Return
     End If
     Call reset_time (.false.)
     Call close_season ()
     Call close_ifile (xfile, nx, &
          nva=mxa)
     If (yfile%cfile%lset) Then ! - restore monthly file settings -
        yfile%cprd1 = mfile%cprd1
        yfile%cprdn = mfile%cprdn
        yfile%cssn = mfile%cssn
        yfile%iseq = mfile%iseq
        yfile%prd1 = mfile%prd1
        yfile%prdn = mfile%prdn
        iseq = mfile%iseq
     End If
     Call window_update (nx)
     Call window_update (mxa)
     ixyz = 3
     GoTo 1
! - Y file -
   Case (2)
     If (icalc == 1) Then
        ireset = reset('Closing the Y file', .true.)
        If (ireset == 1) Return
     End If
     Call close_season ()
     Call close_ifile (yfile, ny, &
          nva=mya)
     Call close_ifile (mfile, ny)
     If (xfile%nt == 1) xfile%iseq = 0
     Call window_update (ny)
     Call window_update (mya)
     icsv = 0
     nmfd = 1
! - forecast file -
   Case (3)
     ifc = close_fcast()
     ifcast = 0
     Call close_ifile (zfile, nz, &
          nva=mza)
     nf = 0
     Call window_update (nf)
     Call window_update (mza)
  End Select
  ixyz = 0
!
! Enable analyses
  Call set_menuflags (ianal, icalc, istd, igcms, icsv, isynch, sw_elev%on, &
       ixdstr=xfile%idstr, iydstr=yfile%idstr, izdstr=zfile%idstr)
!
  Return
 End Function close_file
!
!
!
 Subroutine close_ifile (afile, nv, &
            nva)
!
! Closes an input file
!
! Modules
  Use iofiles, Only: ifile, &
                     init_ifile
!
! Arguments
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - input file -
!
! Output scalars
  Integer, Intent(Out) :: nv ! - number of variables -
!
! - optional output scalars -
  Integer, Intent(Out), Optional :: nva ! - number of used variables -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Reset variables
  nv = 0
  If (Present(nva)) nva = 0
!
! Close input file
  Call init_ifile (afile)
  Call file_refresh (afile, &
       nv=nv, nva=nva)
!
  Return
 End Subroutine close_ifile
End Module cpt_file
