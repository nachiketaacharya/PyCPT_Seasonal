! Author: Simon Mason
Module pcs
!
! Modules
  Use data_numbers,       Only: rp
  Use data_gui_constants, Only: i_handle
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: advanced_cca, advanced_eof, check_pcs, close_svds, get_cca_settings, get_eofs_settings, get_eofx_settings, &
            get_eofy_settings, init_ccamaps, init_pc, init_pcs, init_svds, init_svdwk, make_eof_title, max_cca_modes,  &
            max_eof_modes, proj_read_pcs, proj_write_pcs, set_grey_modes, set_scree
#if GUI == 1
  Public :: change_cca, change_eofx, change_eofy, prompt_eofs
#endif
!
! Scalars
!
! Integer scalars
  Integer, Target, Public :: ieofx ! - X EOF option -
  Integer, Target, Public :: ieofy ! - Y EOF option -
!
  Integer, Pointer, Public :: ieofp ! - EOF option -
!
  Integer, Public :: icco     ! - perfect canonical correlation option flag (0 = continue; 1 = avoid) -
  Integer, Public :: iec      ! - current CCA mode -
  Integer, Public :: ieofx_bk ! - backup X EOF option -
  Integer, Public :: iex      ! - current X EOF mode -
  Integer, Public :: iey      ! - current Y EOF mode -
  Integer, Public :: lc1      ! - first dimension of ce -
  Integer, Public :: lr1      ! - first dimension of r -
  Integer, Public :: lxt      ! - minimum of nxa and nu -
  Integer, Public :: lxyt     ! - Min(Max(nxa,nya),nu) -
  Integer, Public :: lyt      ! - minimum of nya and nu -
  Integer, Public :: mcc      ! - maximum number of CCA modes -
  Integer, Public :: mxe      ! - maximum number of X EOF modes -
  Integer, Public :: mye      ! - maximum number of Y EOF modes -
  Integer, Public :: ncc      ! - number of CCA modes -
  Integer, Public :: nco      ! - optimal number of CCA modes -
  Integer, Public :: nccu     ! - number of canonical modes to calculate -
  Integer, Public :: npx      ! - number of non-zero X eigenvalues -
  Integer, Public :: npy      ! - number of non-zero Y eigenvalues -
  Integer, Public :: nxe      ! - number of X EOF modes -
  Integer, Public :: nxo      ! - optimal number of X EOF modes -
  Integer, Public :: nye      ! - number of Y EOF modes -
  Integer, Public :: nyo      ! - optimal number of Y EOF modes -
! - grey menu flags -
  Integer, Public :: im_modes ! - modes grey menu flag -
  Integer, Public :: im_ogcm  ! - GCM options grey menu flag -
  Integer, Public :: im_xeof  ! - X EOF grey menu flag -
  Integer, Public :: im_yeof  ! - Y EOF grey menu flag -
  Integer, Public :: im_scree ! - scree plot grey menu flag -
!
#if GUI == 1
  Integer, Private :: iecp ! - previous CCA mode -
  Integer, Private :: iexp ! - previous X EOF mode -
  Integer, Private :: ieyp ! - previous Y EOF mode -
!
  Integer, Parameter, Private :: ilceofs = 8 ! - size of ceofs -
!
#endif
  Integer(i_handle), Public :: ih_eof   ! - handle for EOF options window -
!
! Real scalars
  Real(Kind=rp), Public :: cancor ! - canonical correlation -
!
! Logical scalars
  Logical, Public :: lcde = .false. ! - calculate EOF options? -
!
! Arrays
!
! Real arrays
  Real(Kind=rp), Dimension(:), Allocatable, Public :: mu   ! - canonical correlations -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: prjc ! - projections onto CCA modes -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: svx  ! - singular values of x -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: svx2 ! - singular values of x -
  Real(Kind=rp), Dimension(:), Allocatable, Public :: svy  ! - singular values of y -
!
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: bz     ! - principal component regression coefficients -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: bz2    ! - principal component regression coefficients -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: ce     ! - EOF cross-correlations -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: eofx   ! - x EOF patterns -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: eofx2  ! - x EOF patterns -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: eofy   ! - y EOF patterns -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: fsx    ! - forecast time scores of X EOFs -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: fsxc   ! - forecast time scores of X CCA modes -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: hx_map ! - X homogeneous covariance maps -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: hx_ser ! - X homogeneous covariance map time series -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: hy_map ! - Y homogeneous covariance maps -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: hy_ser ! - Y homogeneous covariance map time series -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: r      ! - canonical Y EOF weights -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: s      ! - canonical X EOF weights (transposed) -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: tsx    ! - time-series of x EOFs (transposed) -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: tsx2   ! - time-series of x EOFs (transposed) -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: tsy    ! - time-series of y EOFs (transposed) -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: xwk    ! - workspace -
!
#if GUI == 1
! Character arrays
  Character(Len=ilceofs), Dimension(:), Allocatable, Private :: ceofs ! - list of modes -
!
#endif
Contains
!
!
 Subroutine init_pcs ()
!
! Initialises PC settings
!
! Modules
  Use settings, Only: init_elevel
!
! Executable Statements
!
! Reset calculation options
  ieofx = 1
  ieofy = 1
  ieofx_bk = 1
  icco = 0
!
! Reset numbers of positive modes
  npx = 0
  npy = 0
!
! Reset current modes
  iex = 1
  iey = 1
  iec = 1
#if GUI == 1
  iexp = iex
  ieyp = iey
  iecp = iec
#endif
!
! Reset workspace settings
  lc1 = 0
  lr1 = 0
  nco = 0
  nccu = 0
  nxo = 0
  nyo = 0
!
! Set whether to calculate EOF mode settings
  Call init_elevel (lcde)
!
  Return
 End Subroutine init_pcs
!
!
!
 Function init_ccamaps(nco) &
          Result (init)
!
! Initialises memory and settings to calculate CCA maps
!
! Modules
  Use iofiles,  Only: mxa, mya
  Use settings, Only: nu
!
! Function type
  Integer :: init
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nco ! - number of modes -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Allocate memory
! - X homogeneous maps -
  Allocate (hx_map(mxa,nco), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - Y homogeneous maps -
  Allocate (hy_map(mya,nco), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - X homogeneous maps time series -
  Allocate (hx_ser(nco,nu), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - Y homogeneous maps time series -
  Allocate (hy_ser(nco,nu), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
  init = 0
  Return
!
! Errors
1 init = 1
!
  Return
 End Function init_ccamaps
!
!
!
 Function init_pc(j, k, l, sv, eof, ts) &
          Result (init)
!
! Initialises memory to calculate principal components
!
! Function type
  Integer :: init
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: j ! - number of locations / variables -
  Integer, Intent(In) :: k ! - number of time-steps -
!
! Output scalars
  Integer, Intent(Out) :: l ! - Min(j, k) -
!
! Output arrays
  Real(Kind=rp), Dimension(:), Allocatable :: sv ! - singular values -
!
  Real(Kind=rp), Dimension(:,:), Allocatable :: eof ! - EOF patterns -
  Real(Kind=rp), Dimension(:,:), Allocatable :: ts  ! - time-series of EOFs (transposed) -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Min
!
! Executable Statements
!
! Calculate minimum dimension
  l = Min(j, k)
!
! Initialise memory to calculate principal components
  Allocate (sv(l), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - independent variables EOF patterns -
  Allocate (eof(j,l), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - independent variables EOF scores -
  Allocate (ts(l,k), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
1 init = ifail
!
  Return
 End Function init_pc
!
!
!
 Function init_svds(ianalc) &
          Result (init)
!
! Initialises memory and settings to perform SVDs
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr
  Use iofiles,            Only: mxa, mya
  Use settings,           Only: nu
!
! Function type
  Integer :: init
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianalc ! - analysis identifier -
!
! Locals
!
! Local scalars
#if GUI == 1
  Integer :: i     ! - counter -
  Integer :: mxy   ! - maximum number of modes -
#endif
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic routines
#if GUI == 1
  Intrinsic AdjustL
#endif
  Intrinsic Max, Min
!
! Executable Statements
!
! Allocate memory
! - independent variables SVD -
  If (ianalc==ia_cca .or. ianalc==ia_pcr .or. ianalc==ia_mlr) Then
     If (init_pc(mxa, nu, lxt, svx, eofx, tsx) /= 0) GoTo 1
  End If
!
  Select Case (ianalc)
   Case (ia_cca)
! - dependent variables SVD -
     If (init_pc(mya, nu, lyt, svy, eofy, tsy) /= 0) GoTo 1
     lxyt = Min(Max(mxa, mya), nu)
#if GUI == 1
     mxy = Max(lxt, lyt)
#endif
! - canonical correlations -
     Allocate (mu(lxyt), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - canonical Y EOF weights -
     lr1 = mye
     Allocate (r(lr1,lxyt), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - canonical X EOF weights -
     Allocate (s(lxyt,mxe), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - projections onto CCA modes -
     Allocate (prjc(mcc), Stat=ifail)
     If (ifail /= 0) GoTo 1
   Case (ia_pcr, ia_mlr)
! - principal component regression coefficients -
#if GUI == 1
     mxy = lxt
#endif
     Allocate (bz(lxt,mya), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End Select
#if GUI == 1
! - list of modes -
  Allocate (ceofs(mxy), Stat=ifail)
  If (ifail /= 0) GoTo 1
  Do i = 1, mxy
     Write (Unit=ceofs(i), Fmt='(I8)') i
     ceofs(i) = AdjustL(ceofs(i))
  End Do
#endif
!
  init = 0
  Return
!
! Errors
1 init = 1
!
  Return
 End Function init_svds
!
!
!
 Function init_svdwk() &
          Result (init)
!
! Initialises workspace and settings to perform SVDs
!
! Function type
  Integer :: init
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Allocate memory
! - EOF cross-correlations -
  lc1 = mye
  Allocate (ce(lc1,mxe), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
  init = 0
  Return
!
! Errors
1 init = 1
!
  Return
 End Function init_svdwk
!
!
!
 Function set_scree()
!
! Sets scree plot flag
!
! Function type
  Integer :: set_scree
!
! Executable Statements
!
! Identify whether scree plot can be drawn
  If (npx > 1 .or. npy > 1) Then
     im_scree = 1
  Else
     im_scree = 0
  End If
  set_scree = 0
!
  Return
 End Function set_scree
!
!
!
 Function get_eofs_settings(cxy, ianal, lbf, feof) &
          Result (eofs)
!
! Sets EOF settings
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr
  Use gui,                Only: set_cursor_waiting
!
! Function type
  Integer :: eofs
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis indicator -
!
  Character(Len=*), Intent(In) :: cxy ! - X/Y indicator -
!
  Logical, Intent(In) :: lbf ! - both files set? -
!
! Procedure arguments
  Integer, External :: feof ! - EOF function -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator
!
! Executable Statements
!
! Prompt for numbers of modes if relevant
  Select Case (ianal)
   Case (ia_cca) ! - CCA -
     Call set_cursor_waiting (0)
     ifail = feof()
     If (ifail /= 1 .and. lbf) ifail = get_cca_settings()
   Case (ia_pcr) ! - PCR -
     Select Case (cxy)
      Case ('X', 'x')
        Call set_cursor_waiting (0)
        ifail = feof()
      Case Default
        Continue
     End Select
   Case Default
     Continue
  End Select
  eofs = 2
!
  Return
 End Function get_eofs_settings
!
!
!
 Function get_eofx_settings() &
          Result (modes)
!
! Calls function to set X EOF options
!
! Modules
  Use labels,  Only: l_modeoptsx
  Use iofiles, Only: nx
!
! Function type
  Integer :: modes
!
! Executable Statements
!
! Prompt for EOF options
  modes = get_eof_settings(l_modeoptsx%c, 'X', nx, nxe, mxe, ieofx)
  ieofx_bk = ieofx
!
  Return
 End Function get_eofx_settings
!
!
!
 Function get_eofy_settings() &
          Result (modes)
!
! Calls function to set Y EOF options
!
! Modules
  Use labels,  Only: l_modeoptsy
  Use iofiles, Only: ny
!
! Function type
  Integer :: modes
!
! Executable Statements
!
! Prompt for EOF options
  modes = get_eof_settings(l_modeoptsy%c, 'Y', ny, nye, mye, ieofy)
!
  Return
 End Function get_eofy_settings
!
!
!
 Function get_eof_settings(cmodes, cxyc, nv, ne, me, ieof) &
          Result (modes)
!
! Calls function to set EOF options
!
! Function type
  Integer :: modes
!
! Input scalars
  Integer, Intent(In) :: nv ! - number of variables -
!
  Character(Len=*), Intent(In) :: cmodes ! - options title -
  Character(Len=*), Intent(In) :: cxyc   ! - X/Y variable flag -
!
! Input/output scalars
  Integer, Intent(InOut) :: ne ! - minimum number of modes -
  Integer, Intent(InOut) :: me ! - maximum number of modes -
!
  Integer, Intent(InOut), Target :: ieof ! - EOF options -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for EOF options
  modes = mode_settings(Trim(cmodes), cxyc, nv, ne, me, ieof, advanced_eof)
!
  Return
 End Function get_eof_settings
!
!
!
 Function get_cca_settings() &
          Result (modes)
!
! Sets CCA options
!
! Modules
  Use labels, Only: l_ccaopts
!
! Function type
  Integer :: modes
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Calculate maximum number of modes
  modes = 2
  modes = mode_settings(Trim(l_ccaopts%c), 'CCA', 0, ncc, mcc, iec, advanced_cca)
!
  Return
 End Function get_cca_settings
!
!
!
 Function mode_settings(ctitle, cxyc, nv, ne, me, ieof, fadvanced) &
          Result (modes)
!
! Sets EOF options
!
! Modules
  Use iofiles,  Only: cxy
  Use settings, Only: nt, lcw
!
! Function type
  Integer :: modes
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nv ! - number of variables -
!
  Character(Len=*), Intent(In) :: ctitle ! - title -
  Character(Len=*), Intent(In) :: cxyc   ! - X/Y variable flag -
!
! Input/output scalars
  Integer, Intent(InOut) :: ne ! - minimum number of modes -
  Integer, Intent(InOut) :: me ! - maximum number of modes -
!
  Integer, Intent(InOut), Target :: ieof ! - EOF options -
!
! Procedure arguments
  Integer, External :: fadvanced ! - advanced callback function -
!
! Locals
!
! Local scalars
  Integer :: mm ! - maximum number of modes -
!
! Executable Statements
!
! Calculate maximum number of modes
  modes = 2
  Select Case (cxyc)
   Case ('X', 'x', 'Y', 'y')
     mm = max_eof_modes(lcw, nt, nv, ne, me)
   Case ('CCA')
     mm = max_cca_modes()
  End Select
  If (mm == 1) Return
!
! Prompt for modes options
  If (.not.lcde) Then
     cxy = cxyc
     ieofp => ieof
     modes = get_mode_settings(ctitle, cxyc, ne, me, mm, fadvanced)
  Else
!
! Set CCA options
     me = mm
     Select Case (cxyc)
      Case ('X', 'x', 'Y', 'y')
        Call calc_default_nmodes (ne, me)
     End Select
  End If
!
  Return
 End Function mode_settings
!
!
!
 Function get_mode_settings(ctitle, cxyc, ne, me, mm, fadvanced) &
          Result (modes)
!
! Sets EOF options
!
! Modules
#if GUI == 1
  Use labels,   Only: l_advanced
#endif
  Use labels,   Only: cg_maxnof, cg_minnof, cg_nof, &
                      c_lab=>c_lab, &
                      l_modes, &
                      make_label
  Use gui,      Only: iw, &
                      box_close, box_open, init_win, prompt_integer, win_prompt
#if GUI == 1
  Use errors,   Only: cpt_error
#endif
  Use settings, Only: set_options
!
! Function type
  Integer :: modes
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ctitle ! - title -
  Character(Len=*), Intent(In) :: cxyc   ! - X/Y/CCA mode -
!
! Input/output scalars
  Integer, Intent(InOut) :: ne ! - minimum number of modes -
  Integer, Intent(InOut) :: me ! - maximum number of modes -
  Integer, Intent(InOut) :: mm ! - maximum number of modes -
!
! Procedure arguments
  Integer, External :: fadvanced ! - advanced callback function -
!
#if GUI == 1
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
#endif
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for mode options
  modes = set_options(i1=me, i2=ne)
#if GUI == 1
  Do
     ifail = 0
#endif
     Call init_win (title=ctitle, &
          ihandle=ih_eof)
     c_lab = make_label('L', l_modes)
     Call box_open (.false., .false., &
          title=Trim(cg_nof)//' '//Trim(c_lab)//' ('//Trim(cxyc)//')')
#if GUI == 1
     Call prompt_integer (Trim(cg_maxnof)//' '//Trim(c_lab), me, &
          itab=37, ilow=1, ihgh=mm)
#endif
     Call prompt_integer (Trim(cg_minnof)//' '//Trim(c_lab), ne, &
          itab=37, ilow=1, ihgh=mm)
#if GUI == 0
     Call prompt_integer (Trim(cg_maxnof)//' '//Trim(c_lab), me, &
          ilow=ne, ihgh=mm)
#endif
     Call box_close (.false.)
#if GUI == 1
     c_lab = make_label('T', l_advanced, &
             lak=.true.)
     iw = win_prompt(.true., &
          cextra=Trim(c_lab), fextra=fadvanced)
     If (me < ne) Then
        ifail = 1
        Call cpt_error ('get_mode_settings', ifail, .false., &
             c_arg1=cxyc)
        iw = 0
     End If
#else
     iw = win_prompt(.true.)
#endif
!
! Set options
     modes = set_options(iw=iw, i1=me, i2=ne)
     If (iw == 1) lcde = .false.
#if GUI == 1
     If (ifail == 0) Exit
  End Do
#endif
  modes = 2
!
  Return
 End Function get_mode_settings
!
!
!
 Function advanced_eof()
!
! Sets EOF options
!
! Modules
  Use data_text,          Only: ilang
  Use data_cpt_constants, Only: npo
  Use labels,             Only: ca_pccos, &
                                cg_modecopts_t
  Use gui,                Only: iw, &
                                box_close, box_open, ganged_rbs, gui_creturn, init_win, win_prompt
  Use arrays,             Only: get_flag
  Use iofiles,            Only: cxy
  Use settings,           Only: set_options
!
! Function type
  Integer :: advanced_eof
!
! Locals
!
! Local scalars
  Character(Len=64) :: ctitle ! - title -
!
! Local arrays
  Integer, Dimension(npo) :: ipo ! - principal components calculation option flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Identify current option
  ipo(:) = 0
  ipo(ieofp) = 1
  advanced_eof = set_options(i1=ieofp)
!
! Prompt for EOF calculation options
  Select Case (ilang)
   Case (1, 3, 4)
     ctitle = cxy//' '//Trim(cg_modecopts_t)
   Case Default
     ctitle = Trim(cg_modecopts_t)//' '//cxy
  End Select
  Call init_win (title=Trim(ctitle), &
       idisable=ih_eof) ! - disable host window -
  Call gui_creturn (.false.) 
  Call box_open (.true., .false., &
       title=Trim(ctitle))
  Call ganged_rbs (npo, ipo, &
       crbs=ca_pccos)
  Call box_close (.true.)
  iw = win_prompt(.true.)
!
! Set option
  If (iw == 1) ieofp = get_flag(ipo)
  advanced_eof = set_options(iw=iw, i1=ieofp)
  advanced_eof = 2
!
  Return
 End Function advanced_eof
!
!
!
 Function advanced_cca()
!
! Sets CCA options for too many combined numbers of modes
!
! Modules
  Use labels,   Only: cg_ccacopts_t
  Use gui,      Only: iw, &
                      box_close, box_open, ganged_rbs, init_win, win_prompt
  Use arrays,   Only: get_flag
  Use settings, Only: set_options
!
! Function type
  Integer :: advanced_cca
!
! Locals
!
! Local arrays
  Integer, Dimension(2) :: iccs ! - canonical correlation calculation option flags -
!
! Executable Statements
!
! Backup current option
  advanced_cca = set_options(i1=icco)
  iccs(:) = 0
  iccs(2-icco) = 1
!
! Prompt for CCA calculation options
  Call init_win (title=Trim(cg_ccacopts_t), &
       idisable=ih_eof) ! - disable host window -
  Call box_open (.false., .false., &
       title=Trim(cg_ccacopts_t))
  Call ganged_rbs (2, iccs, (/'Ignore modes with perfect correlation', 'Use modes with perfect correlation   '/))
  Call box_close (.false.)
  iw = win_prompt(.true.)
!
! Set option
  If (iw == 1) icco = get_flag(iccs)
  advanced_cca = set_options(iw=iw, i1=icco)
  advanced_cca = 2
!
  Return
 End Function advanced_cca
!
!
!
 Function max_eof_modes(lcw, nt, nv, ne, me) &
          Result (mn)
!
! Sets maximum number of EOF modes
!
! Modules
  Use data_numbers, Only: ihuge
!
! Function type
  Integer :: mn
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nv  ! - number of variables -
  Integer, Intent(In) :: nt  ! - number of time steps -
  Integer, Intent(In) :: lcw ! - length of cross-validation window -
!
! Input/output scalars
  Integer, Intent(InOut) :: ne ! - minimum number of modes -
  Integer, Intent(InOut) :: me ! - maximum number of modes -
!
! Locals
!
! Local scalars
  Integer :: n ! - nt-lcw-1 -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Min
!
! Executable Statements
!
! Calculate maximum number of modes
  n = nt - lcw - 1
  If (nv > 0) Then
     If (n > 1) Then
        mn = Min(n, nv)
     Else
        mn = nv
     End If
  Else
     If (n > 1) Then
        mn = n
     Else
        mn = ihuge
     End If
  End If
  If (mn == 1) Then
     ne = 1
     me = 1
     Return
  End If
!
  Return
 End Function max_eof_modes
!
!
!
 Function max_cca_modes() &
          Result (mn)
!
! Sets maximum number of CCA modes
!
! Function type
  Integer :: mn
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Min
!
! Executable Statements
!
! Calculate maximum number of modes
  mn = Min(mxe, mye)
  If (mn == 1) Then
     ncc = 1
     mcc = 1
  End If
!
  Return
 End Function max_cca_modes
!
!
!
 Subroutine calc_default_nmodes (ne, me)
!
! Sets default maximum number of EOF modes
!
! Modules
  Use data_cpt_constants, Only: mre
!
! Arguments
!
! Input/output scalars
  Integer, Intent(InOut) :: ne ! - minimum number of modes -
  Integer, Intent(InOut) :: me ! - maximum number of modes -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Max, Min, Nint, Real, Sqrt
!
! Executable Statements
!
! Calculate reasonable number of modes
  ne = 1
  me = Nint(Sqrt(Real(Max(1, me/2))))
  me = Min(me, mre)
!
  Return
 End Subroutine calc_default_nmodes
!
!
!
 Function check_pcs(ianalc) &
          Result (check)
!
! Checks PC settings
!
! On exit:
!    check =  0 All checks passed
!    check =  1 nxe > nt-lcw-1, nxe reset to nt-lcw-1
!    check =  2 nxe > mxa, nxe reset to mxa
!    check =  3 mxe > nxe, mxe reset to nxe
!    check =  4 mxe > nt-lcw-1, mxe reset to nt-lcw-1
!    check =  5 mxe > mxa, mxe reset to mxa
!    check =  6 nye > nt-lcw-1, nye reset to nt-lcw-1
!    check =  7 nye > mya, nye reset to mya
!    check =  8 mye > nye, mye reset to nye
!    check =  9 mye > nt-lcw-1, mye reset to nt-lcw-1
!    check = 10 mye > mya, mye reset to mya
!    check = 11 ncc > nxe, ncc reset to nxe
!    check = 12 ncc > nye, ncc reset to nye
!    check = 13 mcc < ncc, mcc reset to ncc
!    check = 14 mcc > mxe, mcc reset to mxe
!    check = 15 mcc > mye, mcc reset to mye
!    check = 16 mxe+mye > nt
!    check = -1 abort
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr
  Use errors,             Only: cpt_error
  Use iofiles,            Only: mxa, mya
  Use settings,           Only: nt
!
! Function type
  Integer :: check
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianalc ! - analysis indicator -
!
! Locals
!
! Local scalars
  Integer :: ierr  ! - error indicator -
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Check that number of modes is not greater than number of variables or number of cases
  check = -1
  ifail = 0
  If (ianalc == ia_cca .or. ianalc == ia_pcr) Then
     Call check_ne (mxa, nt, nxe, mxe, ifail)
     If (ifail /= 0) Then
        ierr = ifail
        Call cpt_error ('check_pcs', ierr, .false.)
        If (ierr > 0) Return
     End If
     If (ianalc == ia_cca) Then
        Call check_ne (mya, nt, nye, mye, ifail)
        If (ifail /= 0) Then
           ifail = ifail + 5
           ierr = ifail
           Call cpt_error ('check_pcs', ierr, .false.)
           If (ierr > 0) Return
        End If
!
! Check that number of CCA modes is not greater than number of EOF modes
        Call check_nc (nxe, mxe, nye, mye, ncc, mcc, ifail)
        If (ifail /= 0) Then
           ifail = ifail + 10
           ierr = ifail
           Call cpt_error ('check_pcs', ierr, .false.)
           If (ierr > 0) Return
        End If
!
! Check that combined number of modes is not greater than number of cases
        Call check_nes (nt, mxe, mye, ifail)
        If (ifail /= 0) Then
           ifail = 16
           ierr = ifail
           Call cpt_error ('check_pcs', ierr, .false.)
           If (ierr > 0) Return
        End If
     End If
  End If
!
  check = ifail
  Return
 End Function check_pcs
!
!
!
 Subroutine check_ne (nv, nt, ne, me, ifail)
!
! Checks that number of EOF modes is not greater than number of cases-1 or variables
! and that maximum number of modes is not less than minimum number.
!
! On exit:
!    ifail = 0 All checks passed
!    ifail = 1 ne > nt-lcw-1, ne reset to nt-lcw-1
!    ifail = 2 ne > nv, ne reset to nv
!    ifail = 3 me > ne, me reset to ne
!    ifail = 4 me > nt-lcw-1, me reset to nt-lcw-1
!    ifail = 5 me > nv, me reset to nv
!
! Modules
  Use settings, Only: lcw
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nv ! - number of variables -
  Integer, Intent(In) :: nt ! - number of cases -
!
! Input/output scalars
  Integer, Intent(InOut) :: ne ! - minimum number of modes -
  Integer, Intent(InOut) :: me ! - maximum number of modes -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: lvt ! - minimum of nv and nt -
  Integer :: n   ! - number of cross-validation cases - 1  -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Min
!
! Executable Statements
!
! Set maximum number of modes
  n = nt - lcw - 1
  lvt = Min(nv, n)
!
! Force minimum number of modes to be no more than number of variables or number of cases
  If (ne > lvt) Then
     ne = lvt
     If (lvt == n) Then
        ifail = 1
     Else
        ifail = 2
     End If
     Return
  End If
!
! Force maximum number of modes to be at least the minimum number
  If (me < ne) Then
     me = ne
     ifail = 3
     Return
  End If
!
! Force maximum number of modes to be no more than number of variables or number of cases
  If (me > lvt) Then
     me = lvt
     If (lvt == n) Then
        ifail = 4
     Else
        ifail = 5
     End If
     Return
  End If
  ifail = 0
!
  Return
 End Subroutine check_ne
!
!
!
 Subroutine check_nes (nt, mxe, mye, ifail)
!
! Checks that combined number of EOF modes is not greater than number of cases.
!
! On exit:
!    ifail = 0 All checks passed
!    ifail = 1 mxe+mye > nt
!
! Modules
  Use settings, Only: lcw
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: mxe ! - minimum number of X modes -
  Integer, Intent(In) :: mye ! - maximum number of Y modes -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: n ! - number of cross-validation cases  -
  Integer :: m ! - combined number of EOF modes  -
!
! Executable Statements
!
! Identify maximum number of combined modes
  n = nt - lcw
  m = mxe + mye
!
! Warn if number of modes is too large
  If (m > n .and. icco == 1) Then
     ifail = 1
  Else
     ifail = 0
  End If
!
  Return
 End Subroutine check_nes
!
!
!
 Subroutine check_nc (nxe, mxe, nye, mye, ncc, mcc, ifail)
!
! Checks that number of CCA modes is not greater than number of EOF modes
! and that maximum number of modes is not less than minimum number
!
! On exit:
!    ifail = 0 All checks passed
!    ifail = 1 nxe > nxe, ncc reset to nxe
!    ifail = 2 ncc > nye, ncc reset to nye
!    ifail = 3 mcc < ncc, mcc reset to ncc
!    ifail = 4 mcc > mxe, mcc reset to mxe
!    ifail = 5 mcc > mye, mcc reset to mye
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nxe ! - minimum number of X EOF modes -
  Integer, Intent(In) :: nye ! - minimum number of Y EOF modes -
  Integer, Intent(In) :: mxe ! - maximum number of X EOF modes -
  Integer, Intent(In) :: mye ! - maximum number of Y EOF modes -
!
! Input/output scalars
  Integer, Intent(InOut) :: ncc ! - number of CCA modes -
  Integer, Intent(InOut) :: mcc ! - maximum number of CCA modes -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: lxy ! - minimum of nxe and nye -
!
! Executable Statements
!
! Force minimum number of CCA modes to be no more than minimum numbers of EOF modes
  If (nxe <= nye) Then
     lxy = nxe
     ifail = 1
  Else
     lxy = nye
     ifail = 2
  End If
  If (ncc > lxy) Then
     ncc = lxy
     Return
  End If
!
! Force maximum number of CCA modes to be at least the minimum number
  If (mcc < ncc) Then
     mcc = ncc
     ifail = 3
     Return
  End If
!
! Force maximum number of CCA modes to be no more than smallest of the maximum numbers of EOF modes
  If (mxe <= mye) Then
     lxy = mxe
     ifail = 4
  Else
     lxy = mye
     ifail = 5
  End If
  If (mcc > lxy) Then
     mcc = lxy
     Return
  End If
  ifail = 0
!
  Return
 End Subroutine check_nc
!
!
!
#if GUI == 1
 Function change_eofx() &
          Result (change)
!
! Change X mode
!
! Modules
  Use settings, Only: change_setting
!
! Function type
  Integer :: change
!
! Executable Statements
!
! Change number of X modes
  change = change_setting(iex, iexp, 1, npx)
!
  Return
 End Function change_eofx
!
!
!
 Function change_eofy() &
          Result (change)
!
! Change Y mode
!
! Modules
  Use settings, Only: change_setting
!
! Function type
  Integer :: change
!
! Executable Statements
!
! Change number of Y modes
  change = change_setting(iey, ieyp, 1, npy)
!
  Return
 End Function change_eofy
!
!
!
 Function change_cca() &
          Result (change)
!
! Change CCA mode
!
! Modules
  Use screen,   Only: window_update
  Use settings, Only: change_setting
!
! Function type
  Integer :: change
!
! Executable Statements
!
! Change number of CCA modes
  change = change_setting(iec, iecp, 1, nccu)
!
! Update canonical correlation
  cancor = mu(iec)
  Call window_update (cancor)
!
  Return
 End Function change_cca
!
!
!
 Function prompt_eofs(ieof, ne, change_eof, update_eof) &
          Result (prompt)
!
! Prompts for mode number
!
! Modules
  Use labels, Only: cg_mode
  Use gui,    Only: prompt_list
!
! Function type
  Integer :: prompt
!
! Arguments
!
! Input scalars
  Integer, Intent(InOut) :: ieof ! - current mode -
  Integer, Intent(InOut) :: ne   ! - total number of modes -
!
! Procedure arguments
  Integer, External :: change_eof
  Integer, External :: update_eof
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for mode
!  Call prompt_integer (Trim(cg_mode), ieof, &
!       itab=17, fcb1=change_eof, fcb2=update_eof)
  Call prompt_list (Trim(cg_mode), .true., ceofs, ne, ieof, &
       itab=17, iwdth=ilceofs, cbf1=change_eof, cbf2=update_eof)
  prompt = 0
!
  Return
 End Function prompt_eofs
!
!
!
#endif
 Subroutine make_eof_title (nfs, nes, ie, afield, ceof, cmode, img_title)
!
! Constructs title for EOF map
!
! Modules
  Use labels, Only: cg_field
  Use fields, Only: field
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfs ! - number of fields -
  Integer, Intent(In) :: nes ! - number of EOF dimensions -
  Integer, Intent(In) :: ie  ! - current EOF -
!
  Character(Len=1), Intent(In) :: ceof  ! - EOF flag -
  Character(Len=*), Intent(In) :: cmode ! - mode flag -
!
  Type(field), Intent(In) :: afield ! - field -
!
! Output scalars
  Character(Len=*), Intent(Out) :: img_title ! - image title -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Construct default title
  If (nfs > 1 .and. nes > 1) Then
     Write (img_title, Fmt='(A,I0,A)') &
        ceof//' Spatial Loadings ('//cmode, ie, '); '//Trim(cg_field)//': '//Trim(afield%var)//'; Season: '//Trim(afield%cssn)
  Else If (nfs > 1) Then
     Write (img_title, Fmt='(A,I0,A)') &
        ceof//' Spatial Loadings ('//cmode, ie, '); '//Trim(cg_field)//': '//Trim(afield%var)
  Else If (nes > 1) Then
     Write (img_title, Fmt='(A,I0,A)') &
        ceof//' Spatial Loadings ('//cmode, ie, '); Season: '//Trim(afield%cssn)
  Else
     Write (img_title, Fmt='(A,I0,A)') &
        ceof//' Spatial Loadings ('//cmode, ie, ')'
  End If
!
  Return
 End Subroutine make_eof_title
!
!
!
 Subroutine set_grey_modes (ianal, ixdstr, &
            f_scree)
!
! Sets which menu items to display or grey
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal  ! - analysis flag -
  Integer, Intent(In) :: ixdstr ! - X file structure -
!
! Procedure arguments
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
! Set grey Modes menu items
  Select Case (ianal)
   Case (ia_cca) ! - CCA -
     im_modes = 1
     If (Present(f_scree)) im_xeof = f_scree() ! - assign function result to dummy variable -
     im_xeof = 1
     im_yeof = 1
     im_ogcm = 0
   Case (ia_pcr) ! - PCR -
     im_modes = 1
     im_ogcm = 0
     If (Present(f_scree)) im_xeof = f_scree() ! - assign function result to dummy variable -
     im_xeof = 1
     im_yeof = 0
   Case (ia_mlr) ! - MLR -
     im_modes = 0
     im_ogcm = 0
     im_scree = 0
     im_xeof = 0
     im_yeof = 0
   Case (ia_gcm) ! - GCM -
     im_modes = 0
     im_ogcm = 1
     im_scree = 0
     im_xeof = 0
     im_yeof = 0
   Case (ia_pfv) ! - PFV -
     im_modes = 0
     im_ogcm = 1
     im_scree = 0
     im_xeof = 0
     im_yeof = 0
  End Select
!
! Set view flags depending on whether input files are open
  If (ianal /= ia_pfv) Then
     Select Case (ixdstr)
      Case (0)
        im_scree = 0
     End Select
  Else
     im_scree = 0
  End If
!
  Return
 End Subroutine set_grey_modes
!
!
!
 Subroutine proj_read_pcs (iin, ianal, icalc, lres, ifail)
!
! Reads PC settings from project file
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin   ! - input file unit number -
  Integer, Intent(In) :: ianal ! - analysis type -
  Integer, Intent(In) :: icalc ! - results available flag -
!
  Logical, Intent(In) :: lres ! - read results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Read PC settings
  If (.not.lres) Then
     Read (Unit=iin, Err=1, End=2) icco, ieofx, ieofy, mcc, mxe, mye, ncc, nxe, nye
     If (icalc==1) Read (Unit=iin, Err=1, End=2) lc1, lr1, nco, nccu, npx, npy, nxo, nyo
!
! Read PC results
  Else
     Select Case (ianal)
      Case (ia_cca) ! - CCA -
        Read (Unit=iin, Err=1, End=2) eofx
        Read (Unit=iin, Err=1, End=2) eofy
        Read (Unit=iin, Err=1, End=2) hx_map
        Read (Unit=iin, Err=1, End=2) hx_ser
        Read (Unit=iin, Err=1, End=2) hy_map
        Read (Unit=iin, Err=1, End=2) hy_ser
        Read (Unit=iin, Err=1, End=2) mu
        Read (Unit=iin, Err=1, End=2) r
        Read (Unit=iin, Err=1, End=2) s
        Read (Unit=iin, Err=1, End=2) svx
        Read (Unit=iin, Err=1, End=2) svy
        Read (Unit=iin, Err=1, End=2) tsx
        Read (Unit=iin, Err=1, End=2) tsy
      Case (ia_pcr, ia_mlr) ! - PCR, MLR -
        Read (Unit=iin, Err=1, End=2)
        Read (Unit=iin, Err=1, End=2) eofx
        Read (Unit=iin, Err=1, End=2) svx
        Read (Unit=iin, Err=1, End=2) tsx
     End Select
  End If
!
  ifail = 0
  Return
!
! Errors in project file
! - problem reading file -
1 ifail = 1
  Return
!
! - end of file -
2 ifail = 2
  Return
!
 End Subroutine proj_read_pcs
!
!
!
 Subroutine proj_write_pcs (iout, ianal, icalc, lres, ifail)
!
! Writes PC settings to project file
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout  ! - output file unit number -
  Integer, Intent(In) :: ianal ! - analysis type -
  Integer, Intent(In) :: icalc ! - results available flag -
!
  Logical, Intent(In) :: lres ! - read results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write PC settings
  If (.not.lres) Then
     Write (Unit=iout, Err=1) icco, ieofx, ieofy, mcc, mxe, mye, ncc, nxe, nye
     If (icalc == 1) Write (Unit=iout, Err=1) lc1, lr1, nco, nccu, npx, npy, nxo, nyo
!
! Read PC results
  Else
    Select Case (ianal)
      Case (ia_cca) ! - CCA -
        Write (Unit=iout, Err=1) eofx
        Write (Unit=iout, Err=1) eofy
        Write (Unit=iout, Err=1) hx_map
        Write (Unit=iout, Err=1) hx_ser
        Write (Unit=iout, Err=1) hy_map
        Write (Unit=iout, Err=1) hy_ser
        Write (Unit=iout, Err=1) mu
        Write (Unit=iout, Err=1) r
        Write (Unit=iout, Err=1) s
        Write (Unit=iout, Err=1) svx
        Write (Unit=iout, Err=1) svy
        Write (Unit=iout, Err=1) tsx
        Write (Unit=iout, Err=1) tsy
      Case (ia_pcr, ia_mlr) ! - PCR, MLR -
        Write (Unit=iout, Err=1)
        Write (Unit=iout, Err=1) eofx
        Write (Unit=iout, Err=1) svx
        Write (Unit=iout, Err=1) tsx
     End Select
  End If
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
!
 End Subroutine proj_write_pcs
!
!
!
 Function close_svds()
!
! Frees memory allocated for SVDs
!
! Function type
  Integer :: close_svds
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory
  If (Allocated(hy_ser)) Deallocate (hy_ser)
  If (Allocated(hx_ser)) Deallocate (hx_ser)
  If (Allocated(hy_map)) Deallocate (hy_map)
  If (Allocated(hx_map)) Deallocate (hx_map)
  If (Allocated(ce))     Deallocate (ce)
#if GUI == 1
  If (Allocated(ceofs))  Deallocate (ceofs)
#endif
  If (Allocated(prjc))   Deallocate (prjc)
  If (Allocated(s))      Deallocate (s)
  If (Allocated(r))      Deallocate (r)
  If (Allocated(mu))     Deallocate (mu)
  If (Allocated(tsy))    Deallocate (tsy)
  If (Allocated(svy))    Deallocate (svy)
  If (Allocated(eofy))   Deallocate (eofy)
  If (Allocated(bz))     Deallocate (bz)
  If (Allocated(tsx))    Deallocate (tsx)
  If (Allocated(eofx))   Deallocate (eofx)
  If (Allocated(svx))    Deallocate (svx)
!
! Reset workspace settings
  lxt = 0
  lxyt = 0
  lyt = 0
!
  close_svds = 2
!
  Return
 End Function close_svds
End Module pcs
