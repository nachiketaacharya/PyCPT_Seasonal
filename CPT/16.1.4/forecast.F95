! Author: Simon Mason
Module forecast
!
! Modules
  Use data_cpt_constants, Only: ng
  Use data_gui_constants, Only: i_handle
  Use data_io_constants,  Only: lprd
  Use time,               Only: pprd
#if FTN95 == 1
  Use iofiles, Only: zfile ! - for FTN95 compiler bug -
#endif
  Use graphics,           Only: image
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: change_fcast, close_fcast, fcst_fens, fcst_fps, fcst_fser, fcst_fval, fcst_odds, fcst_pexc, &
            proj_read_forecasts, proj_write_forecasts, reset_forecast, set_fcst, update_fprds
#if GUI == 1
  Public :: fcst_hps, init_img_fcst
#else
  Public :: get_nf
#endif
!
! Scalars
!
! Integer scalars
  Integer, Public :: ifcast ! - operational forecast flag -
  Integer, Public :: ifc    ! - recalculate forecast flag -
!                               (0=allocate and forecast;
!                                1=forecast;
!                                2=do nothing;
!                               -1=deallocate first)
!
  Integer, Private :: ifm    ! - current forecast map -
  Integer, Private :: ifmy   ! - year of current forecast map -
#if GUI == 1
  Integer, Private :: ifps    ! - forecast probabilities diagram index -
  Integer, Private :: ifvl    ! - forecast values diagram index -
  Integer, Private :: iodd    ! - forecast odds diagram index -
#endif
  Integer, Private :: ioutf0  ! - forecast output unit number -
  Integer, Private :: ioutf1  ! - forecast output unit number -
  Integer, Private :: ioutf2  ! - forecast output unit number -
  Integer, Private :: ioutf3  ! - forecast output unit number -
  Integer, Private :: ioutf4  ! - forecast output unit number -
  Integer, Private :: ioutf5  ! - forecast output unit number -
  Integer, Private :: ioutf6  ! - forecast output unit number -
  Integer, Private :: ioutf7  ! - forecast output unit number -
  Integer, Private :: ioutf8  ! - forecast output unit number -
  Integer, Private :: ioutf10 ! - forecast output unit number -
  Integer, Private :: itpo    ! - target period offset -
!
  Integer(Kind=i_handle), Private :: ih_f1  ! - forecasts window handle -
  Integer(Kind=i_handle), Private :: ih_f2  ! - forecast ranges window handle -
  Integer(Kind=i_handle), Private :: ih_f3  ! - map forecasts window handle -
  Integer(Kind=i_handle), Private :: ih_f4  ! - map forecasts window handle -
  Integer(Kind=i_handle), Private :: ih_f5  ! - map forecasts window handle -
  Integer(Kind=i_handle), Private :: ih_f8  ! - forecast ensembles window handle -
  Integer(Kind=i_handle), Private :: ih_f10 ! - forecast odds window handle -
!
! Character scalars
  Character(Len=32), Private :: cfmt ! - format statement -
!
! Derived type scalars
  Type(image), Dimension(1), Private :: img_pexc ! - probability of exceedance graph -
!
#if GUI == 1
! Derived type arrays
  Type(image), Dimension( 1), Private :: img_fser ! - forecasts graph -
  Type(image), Dimension( 1), Private :: img_fens ! - forecast ensemble graph -
  Type(image), Dimension( 3), Private :: img_fval ! - forecasts map -
  Type(image), Dimension( 3), Private :: img_hps  ! - probabilistic forecast map -
  Type(image), Dimension(ng), Private :: img_fps  ! - forecast probabilities maps -
  Type(image), Dimension(ng), Private :: img_odds ! - forecast odds map -
!
#endif
  Type(pprd), Dimension(:), Allocatable, Private :: fprd ! - target periods -
!
! Arrays
!
! Character arrays
  Character(Len=lprd), Dimension(:), Allocatable, Public :: cdates ! - list of dates -
!
Contains
!
!
 Subroutine reset_forecast ()
!
! Resets forecast settings
!
! Executable Statements
!
! Reset forecast settings
  ifcast = 0
  ifc = 0
!
  Return
 End Subroutine reset_forecast
!
!
!
#if GUI == 1
 Subroutine init_img_fcst (igsize,ihsize)
!
! Initialises forecast images
!
! Modules
  Use graphics, Only: init_image
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: igsize ! - graphics area dimension -
  Integer, Intent(In) :: ihsize ! - reduced graphics area dimension -
!
! Locals
!
! Local scalars
  Integer :: i ! - diagram index -
!
! Executable Statements
  Call init_image (img_fser(1), (igsize*4)/3)
  Call init_image (img_fens(1), (igsize*4)/3)
  Call init_image (img_pexc(1), (igsize*3)/2)
  Call init_image (img_hps(1), (igsize*3)/2)
  Do i = 1, 3
     Call init_image (img_fval(i), ihsize)
  End Do
  Do i = 1, ng
     Call init_image (img_fps(i), ihsize)
  End Do
  Do i = 1, ng
     Call init_image (img_odds(i), ihsize)
  End Do
!
  Return
 End Subroutine init_img_fcst
!
!
!
 Function reset_img_pexc() &
          Result (reset)
!
! Resets probability of exceedance images
!
! Modules
  Use data_numbers, Only: zero, one
!
! Function type
  Integer :: reset
!
! Executable Statements
!
! Set y-axis limits
  img_pexc(:)%lcontour = .false.
  img_pexc(:)%ymin = zero
  img_pexc(:)%ymax = one
  img_pexc(:)%scale = one
  reset = 0
!
  Return
 End Function reset_img_pexc
!
!
!
 Function reset_img_fps() &
          Result (reset)
!
! Resets forecast probability images
!
! Modules
  Use data_numbers, Only: rp, zero, one, oneh
  Use iofiles,      Only: yfile
!
! Function type
  Integer :: reset
!
! Locals
!
! Local scalars
  Integer :: i ! - index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
!
! Executable Statements
!
! Set forecast probability contours
  img_fps(:)%lcontour= .true.
  Select Case (yfile%idstr)
   Case (1, 2)
     img_fps(:)%ncs = 9
     Do i = 1, img_fps(1)%ncs
        img_fps(:)%contour(i) = Real(i*10, Kind=rp)
     End Do
  End Select
  img_fps(:)%ymin = zero
  img_fps(:)%ymax = oneh
  img_fps(:)%scale = one
  reset = 0
!
  Return
 End Function reset_img_fps
!
!
!
 Function reset_img_hps() &
          Result (reset)
!
! Modules
  Use data_numbers, Only: rp, one
!
! Function type
  Integer :: reset
!
! Locals
!
! Local scalars
  Integer :: i ! - index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
!
! Executable Statements
!
! Set forecast probability contours
  img_hps(:)%lcontour = .true.
  img_hps(:)%ncs = 21
  Do i = 1, img_hps(1)%ncs/2-2
     img_hps(:)%contour(img_hps(1)%ncs/2-1-i) = Real(-35-i*5, Kind=rp)
  End Do
  img_hps(:)%contour(img_hps(1)%ncs/2-1) = Real(-33, Kind=rp)
  img_hps(:)%contour(img_hps(1)%ncs/2) = Real(40, Kind=rp)
  img_hps(:)%contour(img_hps(1)%ncs/2+1) = Real(45, Kind=rp)
  img_hps(:)%contour(img_hps(1)%ncs/2+2) = Real(50, Kind=rp)
  img_hps(:)%contour(img_hps(1)%ncs/2+3) = Real(133, Kind=rp)
  Do i  =1, img_hps(1)%ncs/2-2
     img_hps(:)%contour(img_hps(1)%ncs/2+3+i) = Real(135+i*5, Kind=rp)
  End Do
  img_hps(:)%scale = one
  reset = 0
!
  Return
 End Function reset_img_hps
!
!
!
#endif
 Function init_fcasts(ianalc) &
          Result (init)
!
! Initialises memory for forecasting
!
! Modules
  Use data_cpt_constants,  Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use data_time_constants, Only: nmn
  Use arrays,              Only: fcast, fcastt, fpls, pev, xc, xiny, xvp, yfit, z, ziny, zyp
  Use iofiles,             Only: mxa, mya, nx, ny, zfile
  Use settings,            Only: igauss, nenf, nmf, nf, nu
  Use season,              Only: lprst
  Use pcs,                 Only: fsx, fsxc, mxe, nco
  Use missing,             Only: kfuse
  Use categories,          Only: fps, odds
  Use analysis,            Only: ihind
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
  Integer :: ifail ! - memory allocation flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Lbound, Max, Ubound
!
! Executable Statements
!
! Allocate workspace
! - forecast data -
  Allocate (z(nx,nf,zfile%ngs), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - copy of x -
  Allocate (xc(mxa,Max(nu,nf),1), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - fitted values -
  If (ihind == 1) Then
     Allocate (yfit(mya,nu), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - forecasts -
     Allocate (fcast(mya,nf,0:nenf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - transformed forecasts -
     If (igauss > 0) Then
        Allocate (fcastt(mya,nf), Stat=ifail)
        If (ifail /= 0) GoTo 1
     End If
! - prediction error variance -
     Allocate (pev(mya,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - forecast prediction limits -
     Allocate (fpls(mya,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
! - forecast probabilities -
  Allocate (fps(mya,nf,ng), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - odds -
  Allocate (odds(mya,nf,ng,2), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - predictors over forecast period -
  Select Case (ianalc)
   Case (ia_cca)
     Allocate (xvp(1,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - projections onto X EOFs -
     Allocate (fsx(mxe,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - projections onto X EOFs -
     Allocate (fsxc(nco,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - predictors over forecast period -
   Case (ia_pcr, ia_mlr)
     Allocate (xvp(1,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - projections onto X EOFs -
     Allocate (fsx(mxe,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - predictors over forecast period -
   Case (ia_gcm)
     Allocate (xvp(mya,nf), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - forecast data -
     Allocate (ziny(mya,nf,Lbound(xiny,Dim=3):Ubound(xiny,Dim=3)), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - forecast data -
   Case (ia_pfv)
     Allocate (ziny(mya,nf,zfile%ngs), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End Select
! - used forecast flag -
  Allocate (kfuse(nf), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - persistence component -
  If (lprst) Then
     Allocate (zyp(ny,nmf,nmn), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
!
1 init = ifail
!
  Return
 End Function init_fcasts
!
!
!
 Function init_exc() &
          Result (init)
!
! Initialises memory for exceedence probabilities
!
! Modules
  Use data_cpt_constants, Only: nep
  Use arrays,             Only: c, e, f, tex
!
! Function type
  Integer :: init
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
! Set forecast probability contours
  init = 1
  If (.not.Allocated(f)) Allocate (f(0:nep), Stat=ifail)
  If (ifail /= 0) Return
  If (.not.Allocated(e)) Allocate (e(0:nep), Stat=ifail)
  If (ifail /= 0) Return
  If (.not.Allocated(c)) Allocate (c(0:nep), Stat=ifail)
  If (ifail /= 0) Return
  If (.not.Allocated(tex)) Allocate (tex(2), Stat=ifail)
  If (ifail /= 0) Return
  init = 0
!
  Return
 End Function init_exc
!
!
!
#if GUI == 1
 Function init_hps() &
          Result (init)
!
! Initialises memory for highest probability mapping
!
! Modules
  Use categories, Only: hps
  Use iofiles,    Only: mya
  Use settings,   Only: nf
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
! Set forecast probability contours
  init = 1
  Allocate (hps(mya,nf), Stat=ifail)
  If (ifail /= 0) Return
  init = 0
!
  Return
 End Function init_hps
!
!
!
#else
 Function get_nf ()
!
! Prompts for number of forecasts
!
! Modules
  Use gui,      Only: prompt_integer
  Use settings, Only: nf
!
! Function type
  Integer :: get_nf
!
! Executable Statements
!
! Get training data settings
  Call prompt_integer ('Number of forecasts', nf, &
       ilow=1)
  get_nf = change_fcast()
!
  Return
 End Function get_nf
!
!
!
#endif
 Function set_fcst(iopen)
!
! Sets forecasting to ready state
!
! Modules
  Use screen,   Only: window_update
  Use iofiles,  Only: bkfile, mza, nz, zfile, &
                      init_ifile
  Use settings, Only: nf
  Use fields,   Only: tfield, xfield, zfield
!
! Function type
  Integer :: set_fcst
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iopen ! - forecast file is open (0 = false; 1 = true) -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Associated, Size
!
! Executable Statements
!
! Prompt for forecast file
  set_fcst = iopen
  Select Case (set_fcst)
   Case (0)
     If (zfile%cfile%lset) Then
        If (zfile%cfile%ffile == bkfile%cfile%ffile) Then
           ifc = -1
        Else
           ifc = 0
        End If
! - reset if cancelled -
     Else
        zfield(:) = tfield(:)
        ifc = 2
     End If
   Case Default
     set_fcst = close_fcast()
     Call init_ifile (zfile)
     ifc = 0
     If (set_fcst == 1) Then
        If (Associated(tfield)) Nullify (tfield)
        Allocate (zfield(Size(xfield)))
        zfield(:) = xfield(:)
     End If
  End Select
  If (ifc == 0) ifc = close_fcast()
! - set Tools ~ Forecast menu item grey state -
  If (zfile%cfile%lset) Then
     ifcast = 1
  Else
     ifcast = 0
  End If
!
! Update window
  If (set_fcst == 0) Then
     Call window_update (nf)
     Call window_update (nz)
     Call window_update (mza)
  End If
  Call window_update (ifcast)
  set_fcst = 2
!
  Return
 End Function set_fcst
!
!
!
 Function read_fcast()
!
! Reads data for forecasting
!
! Modules
  Use data_numbers,        Only: rp, one
  Use data_time_constants, Only: isq_yr, nmn
  Use time,                Only: Operator(+), &
                                 iseq, lead, &
                                 date_diff
  Use progress_meter,      Only: lstop, &
                                 end_progress_meter, set_progress_increment, start_progress_meter
  Use arrays,              Only: z, zyp, &
                                 convert_units
  Use iofiles,             Only: mfile, pfile, zfile
  Use settings,            Only: nf, nmf, &
                                 set_it1
  Use space,               Only: idomx, idomy
  Use fields,              Only: mfield, yfield, zfield
  Use season,              Only: lsn, lprst, nsn, &
                                 seasonal_zyp
  Use missing,             Only: kavy, kavz, &
                                 replace_missing_z, replace_missing_zyp
  Use climate,             Only: ave, sdev
  Use analysis,            Only: ianal, lconvert
  Use cpt_input,           Only: num_read, read_data
!
! Function type
  Integer :: read_fcast
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real, Sum
!
! Executable Statements
!
! Check first date of interest
  ifail = 3
  zfile%bdate = zfile%fdate
  zfile%it1 = set_it1(zfile%iseq, zfield(1)%tprd%sdate, zfile%bdate)
  If (lprst) Then
     ifail = 4
     pfile%it1=set_it1(pfile%iseq, mfile%prd1%sdate, pfile%bdate)
!
! Determine number of months to read
     Select Case (iseq)
      Case (isq_yr)
        nmf = (nf - 1)*nmn + lsn(0)
      Case Default
        nmf = date_diff(pfile%bdate, pfile%fdate+(nf-1), pfile%iseq) + 1
     End Select
  Else
     nmf = nf
  End If
!
! Initialise progress meter
  Call start_progress_meter (' ', .false.)
  Call set_progress_increment &
     (one/Real(num_read(zfile%idstr, zfile%lstack, zfile%nt, nmf, Sum(zfield(:)%nlt), &
                        zfile%nfs, zfile%nls, zfile%ngs, zfile%it1), Kind=rp), .false.)
!
! Allocate workspace
  If (init_fcasts(ianal) > 0) Then
     ifail = 1
     GoTo 1
  End If
!
! Read forecast data
  Call read_data (zfile, zfield(:), nf, 1, idomx(:,:), kavz(:,:), z(:,:,:), ifail)
  If (ifail > 0) Then
     ifail = 2
     GoTo 1
  End If
!
! Read persistence data
  If (lprst) Then
     Call read_data (pfile, mfield(:), nf, 1, idomy(:,:), kavy(:,:), zyp(:,:,:), ifail)
     If (ifail > 0) Then
        ifail = 2
        GoTo 1
     End If
!
! Calculate seasonal aggregates
     Call seasonal_zyp (yfield(0:)%tprd%sdate%imn, yfield(0:)%rmiss)
  End If
!
! Replace missing values and check for missing variables
  ifail = replace_missing_z(ianal, nsn)
  If (ifail /= 0) Then
     ifail = 5 - ifail
     GoTo 1
  End If
!
  If (lprst) Then
     ifail = replace_missing_zyp(ave(:,0:), sdev(:,0:))
     If (ifail /= 0) Then
        ifail = 5
        GoTo 1
     End If
  End If
!
! Check for consistent units
  If (lconvert) Then
     Call convert_units (zfield(1)%unitc, yfield(1)%unit, z(:,:,1), ifail)
     If (ifail /= 0) ifail = 2
     If (lprst) Then
        Call convert_units (zfield(1)%unitc, yfield(1)%unit, zyp(:,:,0), ifail)
        If (ifail /= 0) ifail = 2
     End If
  End If
!
! Set initial values
  ifm = 1
  ifmy = zfile%fdate%iyr + lead
!
1 If (lstop) ifail = -1
  read_fcast = ifail
!
  Call end_progress_meter (.false., .false.)
!
  Return
 End Function read_fcast
!
!
!
 Function calc_fcast() &
          Result (calc)
!
! Calculate forecasts
!
! Modules
  Use data_numbers,       Only: rp, one
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use progress_meter,     Only: lstop, &
                                end_progress_meter, set_progress_increment, start_progress_meter
  Use arrays,             Only: b, b0, fcast, fcastt, fpls, iusex, iusey, pev, x, xvp, y, yfit, yhat, yhatt, yret, yrett, yt, &
                                z, ziny, zyp
  Use iofiles,            Only: mxa, mya, yfile, zfile
  Use settings,           Only: clf, iev, igauss, ig0, intp, istd, istdo, izero, lzero, nenf, nf, nu, nur, nu1
  Use space,              Only: get_interpolated
  Use fields,             Only: imfs, yfield, zfield
  Use distribs,           Only: restdize, zero_bound
  Use season,             Only: lb, lprst
  Use pcs,                Only: fsx, fsxc, ieofx, ieofy, nco, nxo, nyo
  Use categories,         Only: fps, odds, pobs, tobst, &
                                calc_odds, calc_reg_probs
  Use climate,            Only: ave, sdev
  Use ensemble,           Only: get_combined
  Use cca,                Only: cca_predict
  Use pcr,                Only: pcr_predict
  Use gcm,                Only: gcm_predict
  Use analysis,           Only: ianal
!
! Function type
  Integer :: calc
!
! Locals
!
! Local scalars
  Real(Kind=rp) :: dprog ! - progess increment -
  Real(Kind=rp) :: r     ! - missing values -
!
  Logical :: lstd ! - apply zero bound if given restandardization? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
!
! Executable Statements
!
! Initialise progress meter
  calc = -1
  Call start_progress_meter ('Calculating forecasts', .false.)
  Select Case (igauss)
   Case (0)
     dprog = Real(2*(nf+izero), Kind=rp)
   Case Default
     dprog = Real(2*(nf+1+izero), Kind=rp)
  End Select
  If (iev == 3) dprog = dprog + Real(nu+izero, Kind=rp)
  Call set_progress_increment (one/dprog, .false.)
  lstd = (istd /= istdo)
!
! Calculate fitted values
  If (iev == 3) Then
     Select Case (ianal)
      Case (ia_cca) ! - CCA -
        Call cca_predict (nu, mxa, mya, nxo, nyo, nco, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
             x(1:mxa,1:nu,1), ieofy, yfit)
      Case (ia_pcr, ia_mlr) ! - PCR/MLR -
        Call pcr_predict (nu, mxa, mya, ieofx, (lstd).and.(lzero).and.(igauss /= 2), &
             x(1:mxa,1:nu,1), b, yfit)
     End Select
  End If
!
! Predict anomaly
  Select Case (ianal)
   Case (ia_cca) ! - CCA -
     Select Case (igauss)
      Case (0)
        If (.not.lprst) Then
           Call cca_predict (nf, mxa, mya, nxo, nyo, nco, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), ieofy, fcast(:,:,0), &
                nt=nu, xvp=xvp(1,:), fsx=fsx, fsxc=fsxc)
        Else
           Call cca_predict (nf, mxa, mya, nxo, nyo, nco, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1),ieofy,fcast(:,:,0), &
                nt=nu, xvp=xvp(1,:), fsx=fsx, fsxc=fsxc, yprst=zyp(1:mya,:,0))
        End If
      Case Default
        If (.not.lprst) Then
           Call cca_predict (nf, mxa, mya, nxo, nyo, nco, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), ieofy, fcast(:,:,0), &
                nt=nu, xvp=xvp(1,:), fsx=fsx, fsxc=fsxc, nc=nu, tclim1=y(1:mya,1:nu,1), fcastt=fcastt)
        Else
           Call cca_predict (nf, mxa, mya, nxo, nyo, nco, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), ieofy, fcast(:,:,0), &
                nt=nu, xvp=xvp(1,:), fsx=fsx, fsxc=fsxc, nc=nu, tclim1=y(1:mya,1:nu,1), tclim0=y(1:mya,1:nu,0), &
                fcastt=fcastt, yprst=zyp(1:mya,:,0))
        End If
     End Select
   Case (ia_pcr, ia_mlr) ! - PCR/MLR -
     Select Case (igauss)
      Case (0)
        If (.not.lprst) Then
           Call pcr_predict (nf, mxa, mya, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), b, fcast(:,:,0), &
                nt=nu, nxe=nxo, xvp=xvp(1,:), fsx=fsx)
        Else
           Call pcr_predict (nf, mxa, mya, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), b, fcast(:,:,0), &
                nt=nu, nxe=nxo, xvp=xvp(1,:), fsx=fsx, yprst=zyp(1:mya,:,0))
        End If
      Case Default
        If (.not.lprst) Then
           Call pcr_predict (nf, mxa, mya, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), b, fcast(:,:,0), &
                nt=nu, nxe=nxo, xvp=xvp(1,:), fsx=fsx, nc=nu, tclim1=y(1:mya,1:nu,1), fcastt=fcastt)
        Else
           Call pcr_predict (nf, mxa, mya, ieofx, (lstd).and.(lzero).and.(igauss/=2), &
                z(1:mxa,:,1), b, fcast(:,:,0), &
                nt=nu, nxe=nxo, xvp=xvp(1,:), fsx=fsx, nc=nu, tclim1=y(1:mya,1:nu,1), tclim0=y(1:mya,1:nu,0), &
                fcastt=fcastt, yprst=zyp(1:mya,:,0))
        End If
     End Select
   Case (ia_gcm) ! - GCM -
     Call get_interpolated (ianal, zfile%idstr, yfile%idstr, intp, mya, iusex, iusey, zfile%nms, yfile%nfl, nf, zfile%ngs, &
          zfield(:)%nv, yfield(1:)%nva, z, yfield(1:)%rmiss, zfield(:)%region, yfield(1:)%region, ziny(:,:,1:))
     If (zfile%nms > 1) Call get_combined (mya, zfile%nms, b0(1:,:), b(1:,:), ziny(:,:,1:), ziny(:,:,0))
     Select Case (igauss)
      Case (0)
        If (.not.lprst) Then
           Call gcm_predict (nf, mya, (lstd).and.(lzero).and.(igauss/=2), &
                ziny(:,:,ig0), b0(ig0,:), b(ig0,:), fcast(:,:,0), &
                nt=nu, xvp=xvp)
        Else
           Call gcm_predict (nf, mya, (lstd).and.(lzero).and.(igauss/=2), &
                ziny(:,:,ig0), b0(ig0,:), b(ig0,:), fcast(:,:,0), &
                nt=nu, xvp=xvp, yprst=zyp(1:mya,:,0))
        End If
      Case Default
        If (.not.lprst) Then
           Call gcm_predict (nf, mya, (lstd).and.(lzero).and.(igauss/=2), &
                ziny(:,:,ig0), b0(ig0,:), b(ig0,:), fcast(:,:,0), &
                nt=nu, xvp=xvp, tclim1=y(1:mya,1:nu,1), fcastt=fcastt)
        Else
           Call gcm_predict (nf, mya, (lstd).and.(lzero).and.(igauss/=2), &
                ziny(:,:,ig0), b0(ig0,:), b(ig0,:), fcast(:,:,0), &
                nt=nu, xvp=xvp, tclim1=y(1:mya,1:nu,1), tclim0=y(1:mya,1:nu,0), fcastt=fcastt, yprst=zyp(1:mya,:,0))
        End If
     End Select
   Case (ia_pfv) ! - PFV -
     Call get_interpolated (ianal, zfile%idstr, yfile%idstr, 0, mya, iusex, iusey, zfile%nfl, yfile%nfl, nf, zfile%ngs, &
          zfield(:)%nv, yfield(1:)%nva, z, zfield(:)%rmiss, zfield(:)%region, yfield(1:)%region, fps(:,:,:),            &
          imfs=imfs)
     r = zfield(1)%rmiss
     Where (fps(:,:,:) == r) fps(:,:,:)=-one
     Call calc_odds (nf,ng,mya,pobs,fps,odds)
     GoTo 1
  End Select
!
! Restandardize if necessary
  If (lstd) Then
     If ((.not.lprst) .and. (igauss == 0)) Then
        If (iev == 3) Then
           Call restdize (istdo, istd, igauss, 1, mya, nu, ave(:,1), sdev(:,1), yfit(:,:))
           Call zero_bound (istd, mya, nu, ave(:,1), sdev(:,1), yfit(:,:))
        End If
        Call restdize (istdo, istd, igauss, 1, mya, nf, ave(:,1), sdev(:,1), fcast(:,:,0))
        Call zero_bound (istd, mya, nf, ave(:,1), sdev(:,1), fcast(:,:,0))
     End If
  End If
!
! Calculate probabilities
  Select Case (igauss)
   Case (0)
     Select Case (iev)
      Case (1)
        Call calc_reg_probs (ianal, igauss, nu, nf, ng, nxo, mya, lb,                                      &
             ave(:,lb), sdev(:,lb), y(:,:,lb), yhat, xvp, fcast(:,:,0), tobst, pobs, clf, fps, odds, fpls, &
             nenf=nenf, fens=fcast(:,:,1:), pev=pev)
      Case (2)
        Call calc_reg_probs (ianal, igauss, nur, nf, ng, nxo, mya, lb,                                            &
             ave(:,lb), sdev(:,lb), y(:,nu1+1:nu,lb), yret, xvp, fcast(:,:,0), tobst, pobs, clf, fps, odds, fpls, &
             nenf=nenf, fens=fcast(:,:,1:), pev=pev)
      Case (3)
        Call calc_reg_probs (ianal, igauss, nu, nf, ng, nxo, mya, lb,                                   &
             ave(:,1), sdev(:,1), y(:,:,1), yfit, xvp, fcast(:,:,0), tobst, pobs, clf, fps, odds, fpls, &
             nenf=nenf, fens=fcast(:,:,1:), pev=pev)
     End Select
   Case Default
     Select Case (iev)
      Case (1)
        Call calc_reg_probs (ianal, igauss, nu, nf, ng, nxo, mya, lb,                          &
             ave(:,lb), sdev(:,lb), yt, yhatt, xvp, fcastt, tobst, pobs, clf, fps, odds, fpls, &
             nenf=nenf, fens=fcast(:,:,1:), nc=nu, clim=y(1:mya,1:nu,lb), pev=pev)
      Case (2)
        Call calc_reg_probs (ianal, igauss, nur, nf, ng, nxo, mya, lb,                                     &
             ave(:,lb), sdev(:,lb), yt(:,nu1+1:nu), yrett, xvp, fcastt, tobst, pobs, clf, fps, odds, fpls, &
             nenf=nenf, fens=fcast(:,:,1:), nc=nu, clim=y(1:mya,1:nu,lb), pev=pev)
      Case (3)
        Call calc_reg_probs (ianal, igauss, nu, nf, ng, nxo, mya, lb,                         &
             ave(:,lb), sdev(:,lb), yt, yfit, xvp, fcastt, tobst, pobs, clf, fps, odds, fpls, &
             nenf=nenf, fens=fcast(:,:,1:), nc=nu, clim=y(1:mya,1:nu,lb), pev=pev)
     End Select
  End Select
  If (lstop) Return
!
1 Call end_progress_meter (.false., .true.)
  calc = 0
!
  Return
 End Function calc_fcast
!
!
!
 Function calc_pexc() &
          Result (calc)
!
! Calculate exceedence probabilities
!
! Modules
  Use data_numbers,       Only: rp, zero, one, oneh
  Use data_cpt_constants, Only: nep
  Use gui,                Only: set_cursor_waiting
  Use progress_meter,     Only: set_progress_increment
  Use arrays,             Only: c, rwk, e, f, fcast, fcastt, pev, tex, y, &
                                insertion_sort
  Use settings,           Only: dofr, igauss, istd, iva, nu
  Use distribs,           Only: alpha, beta, pcen, &
                                cdf_gamma, cdf_empirical, cdf_normal, q_normal, s_student_t, ustdize
  Use season,             Only: lb
  Use climate,            Only: ave, clim, nc, sdev
  Use graphics,           Only: init_exceed
!
! Function type
  Integer :: calc
!
! Locals
!
! Local scalars
  Integer :: i ! - probability of exceedance index -
!
  Real(Kind=rp) :: x    ! - x-value -
  Real(Kind=rp) :: xinc ! - x-value increment -
  Real(Kind=rp) :: xmin ! - x minimum -
  Real(Kind=rp) :: df   ! - number of cases -
  Real(Kind=rp) :: vbar ! - mean -
  Real(Kind=rp) :: vsd  ! - standard deviation -
!
! Local arrays
  Real(Kind=rp), Dimension(1,1) :: xx ! - x-value -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count, Real
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Determine mean and variance
  Select Case (istd)
   Case (0) ! - no standardization -
     vbar = ave(iva,1)
     vsd = sdev(iva,1)
   Case (1) ! - anomalies -
     vbar = zero
     vsd = sdev(iva,1)
   Case (2) ! - standardization -
     vbar = zero
     vsd = one
   Case (3) ! - % of average -
     vbar = oneh
     vsd = oneh*sdev(iva,1)/ave(iva,1)
  End Select
!
! Calculate probabilities of exceedance
  Call set_progress_increment (zero, .false.)
#if GUI == 1
  If (.not.img_pexc(1)%lcontour) Then
#endif
     Call init_exceed (img_pexc(1), nc, clim(iva,1:nc,lb))
#if GUI == 1
  End If
#endif
  df = Real(nc, Kind=rp)
  xinc = (img_pexc(1)%ymax - img_pexc(1)%ymin)/(Real(nep, Kind=rp)*img_pexc(1)%scale)
  xmin = img_pexc(1)%ymin/img_pexc(1)%scale
  tex(1) = xmin
  tex(2) = img_pexc(1)%ymax/img_pexc(1)%scale
  x = xmin
  e(0) = Real(100*Count(clim(iva,1:nc,lb) > x), Kind=rp)/df
  Select Case (igauss)
! - no transformaion -
   Case (0)
     f(0) = oneh*s_student_t((x-fcast(iva,ifm,0))/pev(iva,ifm), dofr)
     c(0) = oneh - oneh*cdf_normal((x - vbar)/vsd)
     Do i = 1, nep
        Call set_cursor_waiting (1)
        x = (xmin + Real(i, Kind=rp)*xinc)
        c(i) = oneh - oneh*cdf_normal((x - vbar)/vsd)                    ! - climatological probability -
        e(i) = Real(100*Count(clim(iva,1:nc,lb) > x), Kind=rp)/df        ! - empirical probability -
        f(i) = oneh*s_student_t((x-fcast(iva,ifm,0))/pev(iva,ifm), dofr) ! - forecast probability -
        Call set_cursor_waiting (0)
     End Do
     Call set_cursor_waiting (1)
! - empirical distribution -
   Case (1)
     If (nc < nu) Then
        rwk(1:nu) = y(iva,1:nu,lb)
        Call insertion_sort (nu, rwk, 'a')
     Else
        rwk(1:nu) = clim(iva,1:nu,lb)
     End If
     x = cdf_empirical(rwk(1:nu), nu, x)
     c(0) = oneh - oneh*x
     x = q_normal(x)
     f(0) = oneh*s_student_t((x-fcastt(iva,ifm))/pev(iva,ifm), dofr)
     Do i = 1, nep
        Call set_cursor_waiting (1)
        x = xmin + Real(i, Kind=rp)*xinc
        e(i) = Real(100*Count(clim(iva,1:nc,lb) > x), Kind=rp)/df       ! - empirical probability -
        x = cdf_empirical(clim(iva,1:nc,lb), nc, x)
        c(i) = oneh - oneh*x                                            ! - climatological probability -
        x = q_normal(x)
        f(i) = oneh*s_student_t((x-fcastt(iva,ifm))/pev(iva,ifm), dofr) ! - forecast probability -
        Call set_cursor_waiting (0)
     End Do
     Call set_cursor_waiting (1)
! - gamma distribution -
   Case (2)
     xx(1,1) = x
     Call ustdize (istd, igauss, lb, 1, 1, ave(iva:iva,lb), sdev(iva:iva,lb), xx)
     x = cdf_gamma(xx(1,1)/beta(iva,lb), alpha(iva,lb), pcen(iva,lb))
     c(0) = oneh - oneh*x
     x = q_normal(x)
     f(0) = oneh*s_student_t((x-fcast(iva,ifm,0))/pev(iva,ifm), dofr)
     Do i = 1, nep
        Call set_cursor_waiting (1)
        x = xmin + Real(i, Kind=rp)*xinc
        e(i) = Real(100*Count(clim(iva,1:nc,lb) > x), Kind=rp)/df        ! - empirical probability -
        xx(1,1) = x
        Call ustdize (istd, igauss, lb, 1, 1, ave(iva:iva,lb), sdev(iva:iva,lb), xx)
        x = cdf_gamma(xx(1,1)/beta(iva,lb), alpha(iva,lb), pcen(iva,lb))
        c(i) = oneh - oneh*x                                                         ! - climatological probability -
        x = q_normal(x)
        f(i) = oneh*s_student_t((x-fcastt(iva,ifm))/pev(iva,ifm), dofr)              ! - forecast probability -
        Call set_cursor_waiting (0)
     End Do
     Call set_cursor_waiting (1)
  End Select
  calc = 1
!
  Return
 End Function calc_pexc
!
!
!
 Function get_forecast()
!
! Checks whether forecasts need to be (re-)calculated
!
! Modules
  Use data_time_constants, Only: isq_mn
  Use time,                Only: Operator(+), Operator(-), Operator(>), Operator(<), &
                                 iseq, lead, &
                                 date_diff, reset_iseq
  Use gui,                 Only: set_cursor_waiting
  Use progress_meter,      Only: lcalc,lstop
  Use errors,              Only: cpt_error
  Use iofiles,             Only: mfile, pfile, xfile, yfile, zfile
  Use settings,            Only: nf
  Use season,              Only: lprst, lsn
  Use categories,          Only: fps
  Use analysis,            Only: analysis_flags, close_analysis
!
! Function type
  Integer :: get_forecast
!
! Locals
!
! Local scalars
  Integer :: ierr  ! - error indicator -
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated, Trim
!
! Executable Statements
!
! Initialise flags
  Call set_cursor_waiting (1)
! -  calculation flags -
  lcalc = .true.
! -  stop calculation flag -
  lstop = .false.
!
! Initialise forecasts if necessary
  If (.not.Allocated(fps)) ifc = 0
!
! Calculate forecasts if necessary
! - check for invalid number of forecasts -
  If (ifc <= 0) Then
     If (zfile%fdate+(nf-1) > zfile%prdn%sdate) Then
        ifail = 3
        GoTo 1
     End If
! - check for availability of persistence component -
     If (lprst) Then
        pfile = mfile
        pfile%fdate = zfile%fdate
        Select Case (iseq)
         Case (:-1)
           Call reset_iseq (isq_new=isq_mn)
         Case Default
           Call reset_iseq (isq_new=iseq)
        End Select
        pfile%bdate = pfile%fdate + date_diff(xfile%bdate, yfile%bdate, iseq)
        pfile%fdate = pfile%bdate + (lsn(0) - 1)
        Call reset_iseq ()
        If (pfile%fdate + (nf-1) > pfile%prdn%sdate) Then
           ifail = 4
           GoTo 1
        End If
        Select Case (iseq)
         Case (:-1)
           Call reset_iseq (isq_new=isq_mn)
         Case Default
           Call reset_iseq (isq_new=iseq)
        End Select
        pfile%fdate = pfile%fdate + 1
        Call reset_iseq ()
        If (pfile%bdate < pfile%prd1%sdate) Then
           ifail = 5
           GoTo 1
        End If
     End If
! - clear memory -
     If (ifc == -1) Then
        lcalc = .false.
        ifail = close_fcast()
        lcalc = .true.
     End If
! - initialise -
     ifail = read_fcast()
     Select Case (ifail)
      Case (0)
        Continue
      Case (3, 4, 5)
        Call cpt_error ('read_fcast', ifail, .true., &
             i_arg1=zfile%iseq)
      Case (6:)
        ierr = ifail - 5
        ifail = 6
        Call cpt_error ('read_fcast', ifail, .true., &
             c_arg1=Trim(zfile%cfile%ffile), i_arg2=ierr)
        GoTo 1
      Case Default
        GoTo 1
     End Select
  End If
  If (ifm == 0) Then
     ifm = 1
     ifmy = zfile%fdate%iyr + lead
  End If
! - calculate forecasts -
  If (ifc <= 1) Then
     ifail = calc_fcast()
     If (ifail /= 0) GoTo 1
     Call get_fprds (nf, 4)
  End If
  ifc = 2
  lcalc = .false.
  ifail = 0
!
! Errors
1 Call set_cursor_waiting (0)
  Select Case (ifail)
   Case (0)
     get_forecast = 0
     Return
   Case (-1)
     Call analysis_flags ('off')
     get_forecast = close_fcast()
     get_forecast = close_analysis()
     get_forecast = 2
     Return
   Case (1)
     Call cpt_error ('get_forecast', ifail, .false.)
     lcalc = .false.
     get_forecast = close_fcast()
   Case (2)
     lcalc = .false.
     get_forecast = close_fcast()
   Case (3, 4, 5)
     Call cpt_error ('get_forecast', ifail, .false.)
     lcalc = .false.
   Case Default
     lcalc = .false.
     get_forecast = close_fcast()
  End Select
  get_forecast = 1
!
  Return
 End Function get_forecast
!
!
!
 Subroutine get_fprds (nf, ifmt)
!
! Calculates target periods
!
! Modules
  Use time,     Only: pdate, &
                      Operator(+), &
                      iseq, &
                      date_diff, get_cdate, set_bdate
  Use iofiles,  Only: pfile, xfile, yfile, zfile
  Use settings, Only: set_it1
  Use season,   Only: lb, lprst, lsn, nsn
!
! Arguments
!
! Input scalars
  Integer :: nf   ! - number of forecasts -
  Integer :: ifmt ! - format -
!
! Locals
!
! Local scalars
  Integer :: k ! - forecast index -
!
! Local scalars
  Type(pdate) :: sdate ! - start date -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated, Size
!
! Executable Statements
!
! Allocate memory
  If (Allocated(fprd)) Then
     If (Size(fprd, Dim=1) /= nf) Then
        Deallocate (fprd)
        Allocate (fprd(nf))
     End If
  Else
     Allocate (fprd(nf))
  End If
  If (Allocated(cdates)) Then
     If (Size(cdates, Dim=1) /= nf) Then
        Deallocate (cdates)
        Allocate (cdates(nf))
     End If
  Else
     Allocate (cdates(nf))
  End If
!
! Calculate target periods
  If (yfile%it1 > 0) Then
     fprd(1) = yfile%prd1 + (yfile%it1 + date_diff(xfile%fdate, zfile%fdate, iseq) - 1)
     If (lprst) Then ! - adjust for persistence component -
        itpo = fprd(1)%sdate%iyr - pfile%bdate%iyr
        fprd(1)%sdate = pfile%bdate
     Else
        itpo = 0
     End if
     cdates(1) = get_cdate(fprd(1), 0, ifmt)
  Else
     fprd(1) = yfile%prd1 + (set_it1(yfile%iseq, yfile%prd1%sdate, yfile%fdate) + date_diff(xfile%fdate, zfile%fdate, iseq) - 1)
     If (lprst) Then ! - adjust for persistence component -
        sdate = fprd(1)%sdate
        Call set_bdate (lb, lsn, nsn, yfile%prd1%sdate, sdate, fprd(1)%sdate, &
                        cprdb=cdates(1))
     Else
        cdates(1) = get_cdate(fprd(1), 0, ifmt)
     End if
  End If
! - repeat for all forecasts -
  If (nf > 1) Then
     Do k = 2, nf
        fprd(k) = fprd(k-1) + 1
        cdates(k) = get_cdate(fprd(k), 0, ifmt)
     End Do
  End If
!
  Return
 End Subroutine get_fprds
!
!
!
 Function update_fprds() &
          Result (update)
!
! Updates forecast periods
!
! Modules
  Use screen,  Only: window_update
  Use iofiles, Only: yfile
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Update season selection indicator
  If (yfile%idstr > 0) Then
     Call get_fprds(1, 2)
     Call window_update (cdates(1))
  End If
  update = 2
!
  Return
 End Function update_fprds
!
!
!
 Function fcst_fser() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Series
!
! Modules
#if GUI == 1
  Use labels,      Only: l_axislimsy, l_fcasts, l_lab
#endif
  Use gui,         Only: set_cursor_waiting
#if GUI == 1
  Use settings,    Only: istd
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#else
  Use gui,         Only: get_output_unit
  Use cpt_output,  Only: print_results
#endif
!
! Function type
  Integer :: fcst
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
#if GUI == 1
! Set image title
  Select Case (istd)
   Case (0)
     img_fser(:)%title = 'Forecasts and Cross-Validated Hindcasts'
   Case (1)
     img_fser(:)%title = 'Forecast and Cross-Validated Hindcast Anomalies'
   Case (2)
     img_fser(:)%title = 'Standardized Forecasts and Cross-Validated Hindcasts'
   Case (3)
     img_fser(:)%title = 'Forecast and Cross-Validated Hindcast SPIs'
   Case (4)
     img_fser(:)%title = 'Forecast and Cross-Validated Hindcast %s of Average'
  End Select
!
! Create window
  fcst = init_window_graphic(1, img_fser(:), imgid, l_fcasts%c, (/l_lab/), &
        update_fser, f_title, f_export,                                    &
        cxy='Y', iloc=ioutf0, limits_grph=f_limits, llims=(/l_axislimsy/), add_win_h=show_fser, f_help=help_fcast)
  fcst = 1
#else
!
! Print forecasts
  ioutf0 = get_output_unit()
  fcst = show_fser()
  fcst = print_results(update_fser)
#endif
!
  Return
 End Function fcst_fser
!
!
!
 Function fcst_fens() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Ensemble
!
! Modules
#if GUI == 1
  Use labels,      Only: l_axislimsy, l_lab
#endif
  Use gui,         Only: set_cursor_waiting
#if GUI == 1
  Use settings,    Only: istd
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#else
  Use gui,         Only: get_output_unit
  Use cpt_output,  Only: print_results
#endif
!
! Function type
  Integer :: fcst
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
#if GUI == 1
! Set image title
  Select Case (istd)
   Case (0)
     img_fens(:)%title='Ensemble Forecasts and Cross-Validated Hindcasts'
   Case (1)
     img_fens(:)%title='Ensemble Forecasts and Cross-Validated Hindcast Anomalies'
   Case (2)
     img_fens(:)%title='Standardized Ensemble Forecasts and Cross-Validated Hindcasts'
   Case (3)
     img_fens(:)%title='Ensemble Forecasts and Cross-Validated Hindcast SPIs'
   Case (4)
     img_fens(:)%title='Ensemble Forecasts and Cross-Validated Hindcast %s of Average'
  End Select
!
! Create window
  fcst = init_window_graphic(1, img_fens(:), imgid, 'Forecast Ensembles', (/l_lab/), &
         update_fens, f_title, f_export,                                             &
         cxy='Y', iloc=ioutf7, limits_grph=f_limits, llims=(/l_axislimsy/), extra_prompt=prompt_forecast, &
         add_win_h=show_fens, f_help=help_fcast)
#else
!
! Print forecasts
  ioutf7 = get_output_unit()
  fcst = show_fens()
  fcst = print_results(update_fens)
#endif
  fcst = 1
!
  Return
 End Function fcst_fens
!
!
!
 Function fcst_pexc() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Exceedence Probabilities
!
! Modules
#if GUI == 1
  Use labels,      Only: l_axislimsx, l_lab
#endif
  Use gui,         Only: set_cursor_waiting
#if GUI == 1
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#else
  Use gui,         Only: get_output_unit
  Use cpt_output,  Only: print_results
#endif
!
! Function type
  Integer :: fcst
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
!
! Initialise memory
  If (init_exc() /= 0) Return
  Call set_cursor_waiting (1)
!
#if GUI == 1
! Initialise graphic
  img_pexc(:)%title = 'Exceedance probabilities'
!
! Create window
  fcst = init_window_graphic(1, img_pexc(:), imgid, 'Exceedance probabilities', (/l_lab/), &
         update_pexc, f_title, f_export,                                                   &
         reset_img=reset_img_pexc, save_rslt=save_pexc, close_grph=close_pexc, cxy='Y', iloc=ioutf6, extra_prompt=prompt_forecast, &
         limits_grph=f_limits, llims=(/l_axislimsx/), f_help=help_fcast)
!
#else
! Print forecasts
  ioutf6 = get_output_unit()
  fcst = print_results(update_pexc)
  fcst = close_pexc()
#endif
  fcst = 1
!
  Return
 End Function fcst_pexc
!
!
!
 Function fcst_fval() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Maps ~ Values
!
! Modules
#if GUI == 1
  Use labels,      Only: l_contourlims, l_fcasts, &
                         ls_threshs
#endif
  Use gui,         Only: set_cursor_waiting
#if GUI == 1
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#endif
!
! Function type
  Integer :: fcst
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
!
! Create window
#if GUI == 1
  fcst = init_window_graphic(3, img_fval(:), imgid, l_fcasts%c, (/ls_threshs(1), l_fcasts, ls_threshs(2)/), &
         update_fval, f_title, f_export, &
         cxy='Y', select_grph=select_fval, extra_prompt=prompt_forecast, add_win_h=show_fval, &
         limits_grph=f_limits, llims=(/l_contourlims, l_contourlims, l_contourlims/), f_help=help_fcast)
#else
!
! Print results
  fcst = show_fval()
  fcst = update_fval()
#endif
  fcst = 1
!
  Return
 End Function fcst_fval
!
!
 Function fcst_fps() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Maps ~ Probabilities
!
! Modules
#if GUI == 1
  Use labels,      Only: l_contourlims, &
                         ls_cats
#endif
  Use gui,         Only: set_cursor_waiting
#if GUI == 1
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#endif
!
! Function type
  Integer :: fcst
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
!
#if GUI == 1
! Create window
  fcst = init_window_graphic(3, img_fps(:), imgid, 'Forecast probabilities maps', ls_cats(1:), update_fps, f_title, f_export, &
         cxy='Y', select_grph=select_fps, extra_prompt=prompt_forecast, add_win_h=show_fps, &
         reset_img=reset_img_fps, limits_grph=f_limits, llims=(/l_contourlims, l_contourlims, l_contourlims/), &
         f_help=help_fcast)
#else
!
! Print results
  fcst = show_fps()
  fcst = update_fps()
#endif
  fcst = 1
!
  Return
 End Function fcst_fps
!
!
!
 Function fcst_odds() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Maps ~ Odds
!
! Modules
#if GUI == 1
  Use labels,      Only: l_contourlims, &
                         ls_cats
#endif
  Use gui,         Only: set_cursor_waiting
#if GUI == 1
  Use settings,    Only: iodds
  Use windows_gui, Only: imgid, &
                         init_window_graphic
#endif
!
! Function type
  Integer :: fcst
!
#if GUI == 1
! Locals
!
! Local scalars
  Character(Len=27) :: codds ! - odds / relative odds -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
#endif
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
!
#if GUI == 1
! Plot odds
  Select Case (iodds)
   Case (0)
     codds = 'Forecast odds maps'
   Case Default
     codds = 'Forecast relative odds maps'
  End Select
!
! Create window
  fcst = init_window_graphic(3, img_odds(:), imgid, Trim(codds), ls_cats(1:), update_odds, f_title, f_export, &
         cxy='Y', select_grph=select_odds, extra_prompt=prompt_forecast, add_win_h=show_odds, &
         limits_grph=f_limits, llims=(/l_contourlims, l_contourlims, l_contourlims/), f_help=help_fcast)
#else
!
! Print results
  fcst = show_odds()
  fcst = update_odds()
#endif
  fcst = 1
!
  Return
 End Function fcst_odds
!
!
!
#if GUI == 1
 Function fcst_hps() &
          Result (fcst)
!
! Create window for Tools ~ Forecast ~ Maps ~ Probabilistic
!
! Modules
  Use data_numbers, Only: one, oneh
  Use labels,       Only: l_lab
  Use gui,          Only: set_cursor_waiting
  Use fields,       Only: yfield
  Use categories,   Only: fps, hps
  Use windows_gui,  Only: imgid, &
                          init_window_graphic
!
! Function type
  Integer :: fcst
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal
!
! Executable Statements
!
! Set cursor
  Call set_cursor_waiting (1)
!
! Check whether forecasts require initialization
  fcst = 1
  If (get_forecast() /= 0) Return
!
! Initialise memory
  If (init_hps() /= 0) Return
!
! Identify highest probabilities
  hps(:,:) = MaxVal(fps(:,:,:), Dim=3)
  Where (fps(:,:,1) == -one) ! - missing values -
     hps(:,:) = yfield(1)%rmiss
  Else Where (hps(:,:) == fps(:,:,1)) ! - below-normal -
     hps(:,:) = -hps(:,:)
  Else Where (hps(:,:) == fps(:,:,3)) ! - above-normal -
     hps(:,:) = hps(:,:) + oneh
  End Where
!
! Set title
  img_hps(:)%title = 'Probabilistic forecasts'
!
! Plot probabilistic forecasts
  fcst = init_window_graphic(1, img_hps(:), imgid, 'Probabilistic forecast map', (/l_lab/), &
         update_hps, f_title, f_export, &
         reset_img=reset_img_hps, close_grph=close_hps, cxy='Y', select_grph=select_hps, extra_prompt=prompt_forecast, &
         f_help=help_fcast)
  fcst = 1
!
  Return
 End Function fcst_hps
!
!
!
#endif
 Function show_fser() &
          Result (show)
!
! Creates output window for Forecast ~ Series results
!
! Modules
#if GUI == 1
  Use labels,   Only: l_fcasts, l_pintervals, l_probs, l_threshs
#endif
  Use gui,      Only: open_clearwin
#if GUI == 1
  Use gui,      Only: box_close, box_open, gui_creturn, print_text
  Use menus,    Only: add_menu
#endif
  Use settings, Only: istd, nf
#if GUI == 1
  Use climate,  Only: print_climatology
  Use graphics, Only: ipis
#endif
!
! Function type
  Integer :: show
!
! Locals
!
! Local scalars
  Integer :: maxc = 8 + ng*4 + 2*(1 + ng*6) ! - maximum number of columns -
  Integer :: maxl = 8                       ! - maximum number of lines -
  Integer :: nl                             ! - number of lines -
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Add pop-up menu
  Call add_menu ('pm', 1, -1, l_pintervals, &
       icheck=ipis, itoggle=ipis, cbf=plot_fser)
!
! Print thresholds
  Call box_open (.true., .false., &
       title=Trim(l_threshs%c))
  Call print_climatology (.true.)
  Call print_text (Trim(l_threshs%c)//':', &
       lcolour=.true.)
  Call gui_creturn (.false.)
#endif
  If (istd /= 3) Then
     Call open_clearwin (ioutf1, maxc, 7, .true., &
          ih=ih_f1)
  Else
     Call open_clearwin (ioutf1, maxc, 11, .true., &
          ih=ih_f1)
  End If
#if GUI == 1
  Call box_close (.true.)
!
! Print forecast information
  Call gui_creturn (.false.)
  Call gui_creturn (.true.)
  Call box_open (.true., .false., &
       title=Trim(l_fcasts%c))
  Call print_text (Trim(l_probs%c)//':', &
       lcolour=.true.)
! - open Clearwin+ windows -
  Call gui_creturn (.false.)
#endif
  nl = nf + 3
  If (nl <= maxl) Then
     Call open_clearwin (ioutf2, maxc, maxl, .true., &
          ih=ih_f2)
  Else
     Call open_clearwin (ioutf2, maxc, maxl, .true., &
          ih=ih_f2, maxl=nl, lvs=.true.)
  End If
#if GUI == 1
!
! Print forecast ranges
  Call gui_creturn (.false.)
  Call gui_creturn (.true.)
  Call print_text ('Forecast ranges:', &
       lcolour=.true.)
! - open Clearwin+ windows -
  Call gui_creturn (.false.)
#endif
  If (nl<=maxl) Then
     Call open_clearwin (ioutf3, maxc, maxl, .true., &
          ih=ih_f3)
  Else
     Call open_clearwin (ioutf3, maxc, maxl, .true., &
          ih=ih_f3, maxl=nl, lvs=.true.)
  End If
#if GUI == 1
  Call box_close (.true.)
#endif
  show = 1
!
  Return
 End Function show_fser
!
!
!
 Function show_fens() &
          Result (show)
!
! Creates output window for Forecast ~ Ensemble results
!
! Modules
#if GUI == 1
  Use labels,   Only: l_fcasts
#endif
  Use gui,      Only: open_clearwin
#if GUI == 1
  Use gui,      Only: box_close, box_open, gui_creturn, print_text
#endif
  Use settings, Only: nenf
!
! Function type
  Integer :: show
!
! Locals
!
! Local scalars
  Integer :: maxc = 28 ! - maximum number of columns -
  Integer :: maxl = 14 ! - maximum number of lines -
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Print forecast information
  Call box_open (.true., .false., &
       title=Trim(l_fcasts%c))
  Call print_text ('Ensemble members:', &
       lcolour=.true.)
! - open Clearwin+ windows -
  Call gui_creturn (.false.)
#endif
  If (nenf+3 <= 14) Then
     Call open_clearwin (ioutf8, maxc, maxl, .true., &
          ih=ih_f8)
  Else
     Call open_clearwin (ioutf8, maxc, maxl, .true., &
          ih=ih_f8, maxl=nenf+3, lvs=.true.)
  End If
#if GUI == 1
  Call box_close (.true.)
#endif
  show = 1
!
  Return
 End Function show_fens
!
!
!
 Function show_fval() &
          Result (show)
!
! Creates output window for Forecast ~ Values results
!
! Modules
  Use labels,      Only: l_fcasts
  Use fields,      Only: yfield
  Use windows_gui, Only: show_results
!
! Function type
  Integer :: show
!
! Locals
!
! Local scalars
  Integer :: maxc = 18 + 3*11 ! - maximum number of columns -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal, Trim
!
! Executable Statements
!
! Print forecast information
  show = show_results(Trim(l_fcasts%c), maxc, MaxVal(yfield(:)%nva)+2, ioutf4, ih_f4)
!
  Return
 End Function show_fval
!
!
!
 Function show_fps() &
          Result (show)
!
! Creates output window for Forecast ~ Probabilities results
!
! Modules
  Use fields,      Only: yfield
  Use windows_gui, Only: show_results
!
! Function type
  Integer :: show
!
! Locals
!
! Local scalars
  Integer :: maxc = 18 + ng*7 ! - maximum number of columns -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal
!
! Executable Statements
!
! Print forecast information
  show = show_results('Forecast probabilities', maxc, MaxVal(yfield(:)%nva)+2, ioutf5, ih_f5)
!
  Return
 End Function show_fps
!
!
!
 Function show_odds() &
          Result (show)
!
! Creates output window for Forecast ~ Odds results
!
! Modules
  Use labels,      Only: l_odds, l_oddsrel
  Use settings,    Only: iodds
  Use fields,      Only: yfield
  Use windows_gui, Only: show_results
!
! Function type
  Integer :: show
!
! Locals
!
! Local scalars
  Integer :: maxc = 18 + ng*7 ! - maximum number of columns -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal, Trim
!
! Executable Statements
!
! Print forecast information
  Select Case (iodds)
   Case (0)
     show = show_results(Trim(l_odds%c), maxc, MaxVal(yfield(:)%nva)+2, ioutf10, ih_f10)
   Case Default
     show = show_results(Trim(l_oddsrel%c), maxc, MaxVal(yfield(:)%nva)+2, ioutf10, ih_f10)
  End Select
!
  Return
 End Function show_odds
!
!
!
 Function update_fser() &
          Result (update)
!
! Updates Forecast ~ Series results and graphs
!
! Modules
#if GUI == 1
  Use clrwin$,             Only: clear_window$
#endif
  Use data_numbers,        Only: rp, zero, eps, one, oneh
  Use data_cpt_constants,  Only: nts
  Use data_io_constants,   Only: lprd
  Use data_time_constants, Only: isq_yr, iuseq
  Use maths,               Only: prob_rounding
  Use labels,              Only: ca_cat_a, ca_seq, &
                                 l_fcast, l_odds, l_oddsrel, l_probs, &
                                 ls_cats, ls_threshs, &
                                 make_label
#if GUI == 0
  Use labels,              Only: l_fcasts, l_threshs
#endif
  Use time,                Only: Operator(+), &
                                 iseq, &
                                 get_cdate
#if GUI == 0
  Use gui,                 Only: print_subtitle
#endif
  Use gui,                 Only: print_text, set_cursor_waiting
  Use arrays,              Only: fcast, fpls, y
  Use settings,            Only: ifpr, iprec, istd, iva, nf
  Use fields,              Only: update_grid
  Use distribs,            Only: cdf_normal
  Use season,              Only: lb
  Use categories,          Only: fps, irefs, ithr, odds, pobs, refs, tobs
#if GUI == 0
  Use climate,             Only: print_climatology
#endif
!
! Function type
  Integer :: update
!
! Locals
!
! Local scalars
  Integer :: i    ! - index -
  Integer :: j    ! - threshold/category index -
  Integer :: k    ! - forecast index -
  Integer :: jj   ! - threshold index -
  Integer :: jmin ! - minimum threshold index -
!
  Real(Kind=rp) :: terr ! - threshold difference -
  Real(Kind=rp) :: tmin ! - minimum threshold difference -
!
  Character(Len=   8) :: c_cat ! - category -
  Character(Len=lprd) :: cdate ! - date -
!
! Local arrays
  Real(Kind=rp), Dimension(ng) :: rps   ! - rounded probabilities -
  Real(Kind=rp), Dimension(ng) :: tpobs ! - theoretical probabilities -
!
  Real(Kind=rp), Dimension(ng,2) :: codds ! - climatological odds -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Abs, AdjustR, All, Huge, Len_Trim, Nint, Trim
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
!
! Construct and add coordinate label for initial point
  Call set_cursor_waiting (1)
  Call update_grid (ioutf0)
!
! Update thresholds
#if GUI == 1
  Call clear_window$ (ih_f1)
#else
  Call print_subtitle (Trim(l_fcasts%c), &
       iout=ioutf0)
  Call print_text (' ', ioutf0)
  Call print_climatology (.true., &
       iout=ioutf0)
  Call print_subtitle (Trim(l_threshs%c), &
       iout=ioutf0)
#endif
  If (ithr /= 3) Then
     Write (Unit=cfmt, Fmt='(A,I1,A)') '(4X,A8,F12.', iprec, ')'
     Do j = nts, 1, -1
        Write (Unit=ioutf1, Fmt=cfmt) ls_threshs(j)%c, tobs(iva,j)
     End Do
  Else
     Write (Unit=cfmt, Fmt='(A,I1,A)') '(4X,A8,F12.', iprec, ',A,I0,A)'
     Do j = nts, 1, -1
        jmin = 0
        tmin = Huge(tmin)
        Do jj = 1, nts
           terr = Abs(tobs(iva,j)-y(iva,irefs(jj), lb))
           If (terr < tmin) Then
              jmin = jj
              tmin = terr
           End If
        End Do
        Write (Unit=ioutf1, Fmt=cfmt) ls_threshs(j)%c, tobs(iva,j), '  (', refs(jmin)%iyr, ')'
     End Do
  End If
  Write (Unit=ioutf1, Fmt=*)
!
! Calculate climatological odds
  If (All(pobs(iva,:) > eps)) Then
     codds(:,1 ) =pobs(iva,:)/(one-pobs(iva,:))
  Else
     codds(:,1) = -one
  End If
  If (istd == 3) Then
     tpobs(ng) = one
     Do j = 1, ng-1 
        tpobs(j) = cdf_normal(tobs(iva,j))
     End Do
     Do j = ng, 2, -1
        tpobs(j) = tpobs(j) - tpobs(j-1)
     End Do
     If (All(tpobs(:) > eps)) Then
        codds(:,2) = tpobs(:)/(one - tpobs(:))
     Else
        codds(:,2) = -one
     End If
  End If
!
! Update climatological probabilities and odds
#if GUI == 0
  Call print_subtitle ('Climatological probabilities and odds',iout=ioutf1)
#endif
  If (istd == 3) Call print_text ('Empirical', &
                      iout=ioutf1)
  If (codds(1,1) >= zero) Then
     Do j = ng, 1, -1
        c_cat = Trim(make_label('T', ls_cats(j)))
        Write (Unit=ioutf1, Fmt='(4X,A,I5,A,F12.2)') AdjustR(c_cat), Nint(pobs(iva,j)*oneh), '%', codds(j,1)
     End Do
  Else
     Do j = ng, 1, -1
        c_cat = Trim(make_label('T', ls_cats(j)))
        Write (Unit=ioutf1, Fmt='(4X,A,I5,A,9X,A)') AdjustR(c_cat), Nint(pobs(iva,j)*oneh), '%', 'N/A'
     End Do
  End If
  If (istd == 3) Then ! - print theoretical probabilities and odds -
     Call print_text ('Theoretical', &
          iout=ioutf1)
     If (codds(1,2) >= zero) Then
        Do j = ng, 1, -1
           c_cat = Trim(make_label('T', ls_cats(j)))
           Write (Unit=ioutf1, Fmt='(4X,A,I5,A,F12.2)') AdjustR(c_cat), Nint(tpobs(j)*oneh), '%',codds(j,2)
        End Do
     Else
        Do j = ng, 1, -1
           c_cat = Trim(make_label('T', ls_cats(j)))
           Write (Unit=ioutf1, Fmt='(4X,A,I5,A,9X,A)') AdjustR(c_cat), Nint(tpobs(j)*oneh), '%', 'N/A'
        End Do
     End If
  End If
!
! Update forecast probabilities
#if GUI == 1
  Call clear_window$ (ih_f2)
  Write (cfmt, '(3(A,I2.2),A)') &
        '(A,T', 20-Len_Trim(l_probs%c), ',A,T', 39-Len_Trim(l_odds%c), ',A,T', 58-Len_Trim(l_oddsrel%c), ',A)'
#else
  Write (Unit=ioutf2, Fmt=*)
  Write (cfmt, '(3(A,I2.2),A)') &
        '(1X,A,T', 37-Len_Trim(l_probs%c), ',A,T', 57-Len_Trim(l_odds%c), ',A,T', 77-Len_Trim(l_oddsrel%c), ',A)'
#endif
  Write (Unit=ioutf2, Fmt=cfmt) ca_seq(iuseq(iseq)), Trim(l_probs%c), Trim(l_odds%c), Trim(l_oddsrel%c)
#if GUI == 1
  Write (Unit=ioutf2, Fmt='(7X,3(3X,A),2(1X,3(5X,A)))') ((ca_cat_a(j), j = 1, ng), k = 1, 3)
  cfmt='(A4,3X,3I4,2(1X,3F6.2))'
#else
  Write (Unit=ioutf2, Fmt='(T25,3(3X,A),2(2X,3(5X,A)))') ((ca_cat_a(j), j = 1, ng), k = 1, 3)
  cfmt='(1X,A,T25,3I4,2(2X,3F6.2))'
#endif
  Write (Unit=ioutf2, Fmt=*)
  Do k = 1, nf
     cdate = Trim(get_cdate(fprd(k)+itpo, 1, 1))
     rps(:) = fps(iva,k,:)
     Call prob_rounding (ifpr, ng, pobs(iva,:), rps(:))
     Write (Unit=ioutf2, Fmt=cfmt) &
        Trim(cdate), (Nint(rps(j)), j=1,ng), ((odds(iva,k,j,i), j=1,ng), i=1,2)
  End Do
!
! Update forecasts
#if GUI == 1
  Call clear_window$ (ih_f3)
#else
  Write (Unit=ioutf3, Fmt=*)
#endif
  Write (cfmt, '(3(A,I2.2),A)') &
        '(1X,A,T', 17-Len_Trim(l_fcast%c), ',A,T', 28-Len_Trim(ls_threshs(1)%c), ',A,T', 39-Len_Trim(ls_threshs(2)%c), ',A)'
  Write (Unit=ioutf3, Fmt=cfmt) ca_seq(iuseq(iseq)), Trim(l_fcast%c), Trim(ls_threshs(1)%c), Trim(ls_threshs(2)%c)
  Write (Unit=ioutf3, Fmt=*)
  Select Case (iseq)
   Case (isq_yr)
     Write (Unit=cfmt, Fmt='(A,I1,A)') '(1X,I4,3F11.', iprec, ')'
     Do k = 1, nf
        Write (Unit=ioutf3, Fmt=cfmt) fprd(k)%sdate%iyr+itpo, fcast(iva,k,0), fpls(iva,k)%lower, fpls(iva,k)%upper
     End Do
   Case Default
     Write (Unit=cfmt, Fmt='(A,I1,A)') '(1X,A4,3F11.', iprec, ')'
     Do k = 1, nf
        cdate = get_cdate(fprd(k)%sdate, 1)
        Write (Unit=ioutf3, Fmt=cfmt) Trim(cdate), fcast(iva,k,0), fpls(iva,k)%lower, fpls(iva,k)%upper
     End Do
  End Select
!
! Update graph
#if GUI == 1
  update = plot_fser()
#else
  Write (Unit=ioutf3, Fmt=*)
  update = 0
#endif
  Call set_cursor_waiting (0)
!
  Return
 End Function update_fser
!
!
!
 Function update_fens() &
          Result (update)
!
! Updates Forecast ~ Ensemble results and graphs
!
! Modules
#if GUI == 1
  Use clrwin$,      Only: clear_window$
#endif
  Use data_numbers, Only: rp
#if GUI == 0
  Use labels,       Only: cg_ensemble, &
                          l_fcasts
  Use gui,          Only: print_subtitle
#endif
  Use gui,          Only: print_text, set_cursor_waiting
  Use arrays,       Only: fcast
  Use settings,     Only: iprec, iva, nenf
#if GUI == 0
  Use settings,     Only: nf
#endif
  Use fields,       Only: update_grid
#if GUI == 0
  Use climate,      Only: print_climatology
#endif
!
! Function type
  Integer :: update
!
! Locals
!
! Local scalars
  Integer :: i ! - ensemble member -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
#if GUI == 0
  Intrinsic Trim
#endif
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
!
! Construct and add coordinate label for initial point
  Call set_cursor_waiting (1)
  Call update_grid (ioutf7)
!
! Update forecasts
  Write (Unit=cfmt, Fmt='(A,I1,A)') '(1X,I6,F10.3,F11.', iprec, ')'
#if GUI == 1
  Call clear_window$ (ih_f8)
#else
  Call print_subtitle (Trim(l_fcasts%c), &
       iout=ioutf7)
  Call print_text (' ', &
       iout=ioutf7)
  Call print_climatology (.true., &
       iout=ioutf7)
  Call print_subtitle (Trim(cg_ensemble), &
       iout=ioutf7)
  Do ifm = 1, nf
     update = change_ifm()
#endif
     Call print_text ('Member  Quantile   Forecast', &
          iout=ioutf8)
     Call print_text (' ', &
          iout=ioutf8)
     Do i = 1, nenf
        Write (Unit=ioutf8, Fmt=cfmt) i, Real(i, Kind=rp)/Real(nenf+1, Kind=rp), fcast(iva,ifm,i)
     End Do
!
#if GUI == 1
! Update graph
     update = plot_fens()
#else
     Call print_text (' ', &
          iout=ioutf8)
  End Do
!
  update = 0
#endif
  Call set_cursor_waiting (0)
!
  Return
 End Function update_fens
!
!
!
 Function update_pexc() &
          Result (update)
!
! Updates Forecast ~ Probability of exceedance results and graphs
!
! Modules
#if GUI == 0
  Use data_cpt_constants, Only: nep
  Use labels,             Only: ca_stds_t
  Use gui,                Only: print_text
#endif
  Use gui,                Only: set_cursor_waiting
#if GUI == 0
  Use arrays,             Only: c, e, f
  Use settings,           Only: istd, nf
#endif
  Use fields,             Only: update_grid
!
! Function type
  Integer :: update
!
#if GUI == 0
! Locals
!
! Local scalars
  Integer :: i ! - probability of exceedance index -
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
#endif
! Construct and add coordinate label for initial point
  Call update_grid (ioutf6)
!
#if GUI == 1
! Update graph
  update = plot_pexc()
!
#else
! Repeat for each forecast
  Do ifm = 1, nf
     update = change_ifm()
!
! Calculate probabilities of exceedance
     update = calc_pexc()
!
! Print probabilities of exceedance
     Call print_text (Trim(ca_stds_t(istd+1)))
     Call print_text ('Probabilities of exceedance:')
     Call print_text ('     Given forecast   Empirical   Climatological')
     Do i = 1, nep
        Write (Unit=*, Fmt=*) i, f(i), e(i), c(i)
     End Do
     Call print_text (' ')
  End Do
#endif
  Call set_cursor_waiting (0)
  update = 1
!
  Return
 End Function update_pexc
!
!
!
 Function update_fval() &
          Result (update)
!
! Updates Forecast ~ Values results and maps
!
! Modules
  Use labels,     Only: ca_dstruct
#if GUI == 1
  Use labels,     Only: ca_stds_t, &
                        ls_threshs
#endif
  Use gui,        Only: set_cursor_waiting
  Use iofiles,    Only: yfile
#if GUI == 1
  Use settings,   Only: clf, istd
#else
  Use settings,   Only: nf
#endif
  Use settings,   Only: iprec
  Use cpt_output, Only: update_results
!
! Function type
  Integer :: update
!
  Character(Len=64) :: chead ! - header line -
#if GUI == 1
  Character(Len= 6) :: cclf  ! - forecast confidence level -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
! Set titles
  Write (Unit=cclf, Fmt='(F6.3)') clf
  If (cclf(1:1) /= '*') Then
     img_fval(3)%title = Trim(ls_threshs(2)%c)//' ('//cclf//'%)'
     img_fval(1)%title = Trim(ls_threshs(1)%c)//' ('//cclf//'%)'
  Else
     img_fval(3)%title = Trim(ls_threshs(2)%c)//' (>99.999%)'
     img_fval(1)%title = Trim(ls_threshs(1)%c)//' (>99.999%)'
  End If
  img_fval(2)%title = ca_stds_t(istd+1)
!
#else
  Do ifm = 1, nf
     update = change_ifm()
#endif
! Print forecasts
  Write (Unit=cfmt, Fmt='(A,I1,A)') '(A18,3F11.', iprec, ')'
  Select Case (yfile%idstr)
   Case (1)
     Write (Unit=chead, Fmt='(A)') '    Lat.    Long.   Forecast      Lower      Upper'
   Case (2, 3)
     Write (Unit=chead, Fmt='(1X,2A)') ca_dstruct(yfile%idstr), '    Forecast      Lower      Upper'
  End Select
  update = update_results(ioutf4, ih_f4, ifm, Trim(chead), .true., print_fval)
!
! Plot maps
#if GUI == 1
  update = plot_fvals()
#else
  End Do
!
  update = 0
#endif
  Call set_cursor_waiting (0)
!
  Return
 End Function update_fval
!
!
!
 Function update_fps() &
          Result (update)
!
! Updates Forecast ~ Probabilities results and maps
!
! Modules
#if GUI == 1
  Use clrwin$,    Only: clear_window$
  Use labels,     Only: ls_cats
#endif
  Use labels,     Only: ca_dstruct
  Use gui,        Only: set_cursor_waiting
  Use iofiles,    Only: yfile
#if GUI == 0
  Use settings,   Only: nf
#endif
  Use cpt_output, Only: update_results
!
! Function type
  Integer :: update
!
! Locals
!
! Local scalars
#if GUI == 1
  Integer :: i ! - image index -
!
#endif
  Character(Len=40) :: chead ! - header line -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustR, Trim
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
#if GUI == 1
! Clear window
  Call clear_window$ (ih_f5)
!
! Set titles
  Do i = 1, ng
     img_fps(i)%title = Trim(ls_cats(i)%c)
  End Do
!
#else
  Do ifm = 1, nf
     update = change_ifm()
#endif
! Print forecasts
  Select Case (yfile%idstr)
   Case (1)
     cfmt = '(A18,3(I6,A))'
     Write (Unit=chead, Fmt='( A)') '    Lat.    Long.  Below Normal  Above'
   Case (2, 3)
     cfmt = '(A,T18,3(I6,A))'
     Write (Unit=chead, Fmt='(2A)') AdjustR(ca_dstruct(yfile%idstr)), '   Below Normal  Above'
  End Select
  update = update_results(ioutf5, ih_f5, ifm, Trim(chead), .true., print_fps)
!
! Plot maps
#if GUI == 1
  update = plot_fpss()
#else
  End Do
  update = 0
#endif
  Call set_cursor_waiting (0)
!
  Return
 End Function update_fps
!
!
!
 Function update_odds() &
          Result (update)
!
! Updates Forecast ~ Odds results and maps
!
! Modules
#if GUI == 1
  Use clrwin$,    Only: clear_window$
  Use labels,     Only: ls_cats
#endif
  Use labels,     Only: ca_dstruct
  Use gui,        Only: set_cursor_waiting
  Use iofiles,    Only: yfile
#if GUI == 0
  Use settings,   Only: nf
#endif
  Use cpt_output, Only: update_results
!
! Function type
  Integer :: update
!
! Locals
!
! Local scalars
#if GUI == 1
  Integer :: i ! - image index -
!
#endif
  Character(Len=40) :: chead ! - header line -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustR, Trim
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
#if GUI == 1
! Clear window
  Call clear_window$ (ih_f10)
!
! Set titles
  Do i = 1, ng
     img_odds(i)%title = Trim(ls_cats(i)%c)
  End Do
!
#else
  Do ifm = 1, nf
     update = change_ifm()
#endif
! Print forecast odds
  Select Case (yfile%idstr)
   Case (1)
     cfmt = '(A18,3F7.2)'
     Write (Unit=chead, Fmt='( A)') '    Lat.    Long.  Below Normal  Above'
   Case (2,3)
     cfmt = '(A,T18,3F7.2)'
     Write (Unit=chead, Fmt='(2A)') AdjustR(ca_dstruct(yfile%idstr)), '   Below Normal  Above'
  End Select
  update = update_results(ioutf10, ih_f10, ifm, Trim(chead), .true., print_odds)
!
! Plot maps
#if GUI == 1
  update = plot_oddss()
#else
  End Do
  update = 0
#endif
  Call set_cursor_waiting (0)
!
  Return
 End Function update_odds
!
!
!
#if GUI == 1
 Function update_hps() &
          Result (update)
!
! Updates Forecast ~ Probabilistic results and maps
!
! Modules
  Use gui, Only: set_cursor_waiting
!
! Function type
  Integer :: update
!
! Executable Statements
!
! Check whether forecasts require initialization
  update = 1
  If (get_forecast() /= 0) Return
  Call set_cursor_waiting (1)
!
! Plot maps
  update = plot_hps()
  Call set_cursor_waiting (0)
!
  Return
 End Function update_hps
!
!
!
#endif
 Function print_fval(ioutf, ija, ifm, istn, ifld) &
          Result (printf)
!
! Prints Forecast ~ Values
!
! Modules
  Use arrays,  Only: fcast, fpls
  Use iofiles, Only: yfile
  Use space,   Only: coor, cstndy
!
! Function type
  Integer :: printf
!
! Arguments
!
! Input scalaras
  Integer, Intent(In) :: ioutf ! - output unit number -
  Integer, Intent(In) :: ija   ! - location index -
  Integer, Intent(In) :: ifm   ! - forecast index -
  Integer, Intent(In) :: istn  ! - station index -
  Integer, Intent(In) :: ifld  ! - field index -
!
! Executable Statements
!
! Print forecast values
  Select Case (yfile%idstr)
   Case (1)
     Write (Unit=ioutf, Fmt=cfmt) coor, fcast(ija,ifm,0), fpls(ija,ifm)%lower, fpls(ija,ifm)%upper
   Case (2, 3)
     Write (Unit=ioutf, Fmt=cfmt) cstndy(istn,ifld), fcast(ija,ifm,0), fpls(ija,ifm)%lower, fpls(ija,ifm)%upper
  End Select
  printf = 0
!
  Return
 End Function print_fval
!
!
!
 Function print_fps(ioutf, ija, ifm, istn, ifld) &
          Result (printf)
!
! Prints Forecast ~ Probabilities
!
! Modules
  Use iofiles,    Only: yfile
  Use space,      Only: coor, cstndy
  Use categories, Only: fps
!
! Function type
  Integer :: printf
!
! Arguments
!
! Input scalaras
  Integer, Intent(In) :: ioutf ! - output unit number -
  Integer, Intent(In) :: ija   ! - location index -
  Integer, Intent(In) :: ifm   ! - forecast index -
  Integer, Intent(In) :: istn  ! - station index -
  Integer, Intent(In) :: ifld  ! - field index -
!
! Locals
!
! Local scalars
  Integer :: j ! - forecast index -
!
! Functions and Subroutines

! Intrinsic functions
  Intrinsic Nint
!
! Executable Statements
!
! Print forecast probabilities
  Select Case (yfile%idstr)
   Case (1)
     Write (Unit=ioutf, Fmt=cfmt) coor, (Nint(fps(ija,ifm,j)), '%', j=1,ng)
   Case (2, 3)
     Write (Unit=ioutf, Fmt=cfmt) cstndy(istn,ifld), (Nint(fps(ija,ifm,j)), '%', j=1,ng)
  End Select
  printf = 0
!
  Return
 End Function print_fps
!
!
!
 Function print_odds(ioutf,ija,ifm,istn,ifld) &
          Result (printf)
!
! Prints Forecast ~ Odds
!
! Modules
  Use iofiles,    Only: yfile
  Use settings,   Only: iodds
  Use space,      Only: coor, cstndy
  Use categories, Only: odds
!
! Function type
  Integer :: printf
!
! Arguments
!
! Input scalaras
  Integer, Intent(In) :: ioutf ! - output unit number -
  Integer, Intent(In) :: ija   ! - location index -
  Integer, Intent(In) :: ifm   ! - forecast index -
  Integer, Intent(In) :: istn  ! - station index -
  Integer, Intent(In) :: ifld  ! - field index -
!
! Locals
!
! Local scalars
  Integer :: j ! - forecast index -
!
! Functions and Subroutines

!
! Executable Statements
!
! Print forecast odds
  Select Case (yfile%idstr)
   Case (1)
     Write (Unit=ioutf, Fmt='(A18,3F7.2)') coor, (odds(ija,ifm,j,iodds+1), j=1,ng)
   Case (2, 3)
     Write (Unit=ioutf, Fmt='(A,T18,3F7.2)') cstndy(istn,ifld), (odds(ija,ifm,j,iodds+1), j=1,ng)
  End Select
  printf = 0
!
  Return
 End Function print_odds
!
!
!
 Function change_fcast() &
          Result (change)
!
! Determines whether a change in the forecast date requires forecast recalculation
!
! Modules
  Use iofiles,  Only: zfile
  Use settings, Only: record_change
!
! Function type
  Integer :: change
!
! Executable Statements
!
! Check for changes
  If (zfile%cfile%lset) Then
     ifc = record_change()
     ifc = -1
  End If
  change = 2
!
  Return
 End Function change_fcast
!
!
!
 Function change_ifm() &
          Result (change)
!
! Changes current forecast date
!
! Modules
#if GUI == 1
  Use time,    Only: lead
  Use iofiles, Only: zfile
#else
  Use gui,     Only: print_text
#endif
!
! Function type
  Integer :: change
!
#if GUI == 0
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
#endif
!
! Executable Statements
!
#if GUI == 1
! Update year
  ifmy = ifm + (zfile%fdate%iyr+lead) - 1
  change = 2
#else
  Call print_text (' ')
  Call print_text ('Date: '//Trim(cdates(ifm)))
  change = 0
#endif
!
  Return
 End Function change_ifm
!
!
!
#if GUI == 1
 Function prompt_forecast(fcallback) &
          Result (prompt)
!
! Prompts for forecast date
!
! Modules
  Use labels,   Only: l_fcast
  Use gui,      Only: print_ivalue, prompt_list
  Use settings, Only: nf
!
! Function type
  Integer :: prompt
!
! Arguments
!
! Procedure arguments
  Integer, External :: fcallback ! - callback function -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for forecast
  If (nf > 1) Then
!       Call prompt_list (Trim(l_fcast%c), .true., cdates(:), nf, ifm, & ! - FTN95 bug -
       Call prompt_list (Trim(l_fcast%c), .true., cdates, nf, ifm, &
          itab=15, iwdth=lprd, cbf1=change_ifm, cbf2=fcallback)
  Else
     Call print_ivalue (Trim(l_fcast%c), ifmy, &
          itab=15)
  End If
  prompt = 1
!
  Return
 End Function prompt_forecast
!
!
!
 Function save_pexc()
!
! Modules
  Use iofiles,    Only: exfile
  Use cpt_output, Only: save_results
!
! Function type
  Integer :: save_pexc
!
! Executable Statements
!
! Save requested results
  save_pexc = save_results(exfile)
!
  Return
 End Function save_pexc
!
!
!
 Function plot_fser() &
          Result (plot)
!
! Plots forecast series graph
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Plot graph
  plot = plot_fcasts(img_fser(1), 0, 1)
!
  Return
 End Function plot_fser
!
!
!
 Function plot_fens() &
          Result (plot)
!
! Plots forecast ensembles graph
!
! Modules
  Use settings, Only: nenf
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Update graph
  plot = plot_fcasts(img_fens(1), 1, nenf)
!
  Return
 End Function plot_fens
!
!
!
 Function plot_fcasts(img, i1, mf) &
          Result (plot)
!
! Plots forecast series graph
!
! Modules
  Use time,     Only: Operator(+), &
                      nse
  Use labels,   Only: cg_black_l, cg_green_l, cg_grey_l, cg_red_l, &
                      l_fcasts
  Use gui,      Only: icbw
  Use arrays,   Only: fcast, fpls, y, yhat
  Use iofiles,  Only: yfile
  Use settings, Only: iva, nf, nt, nu
  Use season,   Only: lb
  Use plotting, Only: plot_series
!
! Function type
  Integer :: plot
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: i1 ! - index of first ensemble member -
  Integer, Intent(In) :: mf ! - number of ensemble members -
!
! Input/output scalars
  Type(image), Intent(InOut) :: img ! - image -
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
  Select Case (icbw)
   Case (0)
     cylab = 'Observations ('//Trim(cg_red_l)//') / '//Trim(l_fcasts%c)//' ('//Trim(cg_green_l)//')'
   Case (1)
     cylab = 'Observations ('//Trim(cg_black_l)//') / '//Trim(l_fcasts%c)//' ('//Trim(cg_grey_l)//')'
  End Select
!
! Plot graph
  If (nse == 1) Then
     plot = plot_series(img, Trim(cylab), nt, nu, y(iva,1:nu,lb), yfile%fdate, nt, nu, yhat(iva,1:nu), yfile%fdate, 1, &
            df=fprd(1)%sdate+itpo, nf=nf, mf=mf, yf=fcast(iva,1:nf,i1:i1+mf-1), yfis=fpls(iva,1:nf))
  Else
     plot = plot_series(img, Trim(cylab), nt, nu, y(iva,1:nu,lb), yfile%fdate, nt, nu, yhat(iva,1:nu), yfile%fdate, 1, &
            df=fprd(1)%sdate+itpo, nf=nf, mf=mf, yf=fcast(iva,1:nf,i1:i1+mf-1), yfis=fpls(iva,1:nf), nse=nse)
  End If
!
  Return
 End Function plot_fcasts
!
!
!
 Function replot_pexc() &
          Result (plot)
!
! Replots probability of exceedance graph without resetting y-axis
!
! Modules
  Use data_cpt_constants, Only: nep, nts
  Use labels,             Only: ca_stds_t
  Use arrays,             Only: c, e, f
  Use settings,           Only: istd, iva
  Use categories,         Only: tobs
  Use climate,            Only: clim, ncu
  Use plotting,           Only: plot_exceed
!
! Function type
  Integer :: plot
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for modifications
  plot = plot_exceed(img_pexc(1), Trim(ca_stds_t(istd+1)), nep, f(:), e(:), c(:), ncu(1), clim(iva,1:ncu(1),1), &
         nts=nts, tobs=tobs(iva,:))
  plot = 2
!
  Return
 End Function replot_pexc
!
!
!
 Function plot_pexc() &
          Result (plot)
!
! Plots probability of exceedance graph initializing y-axis
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Calculate probabilities of exceedence
  plot = calc_pexc()
!
! Plot probabilities of exceedance
  plot = replot_pexc()
  plot = 1
!
  Return
 End Function plot_pexc
!
!
!
 Function plot_fvals() &
          Result (plot)
!
! Plots Forecast ~ Maps ~ Values results
!
! Modules
  Use data_numbers, Only: rp
  Use arrays,       Only: fcast
!
! Function type
  Integer :: plot
!
! Locals
!
! Local scalars
  Integer :: i ! - image index -
!
  Real(Kind=rp) :: vlow ! - minimum value -
  Real(Kind=rp) :: vhgh ! - maximum value -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal, MinVal
!
! Executable Statements
!
! Calculate data limits
  vlow = MinVal(fcast(:,:,0))
  vhgh = MaxVal(fcast(:,:,0))
!
! Plot forecast maps
  Do i = 1, 3
     plot = plot_fval(i, &
            vlow=vlow, vhgh=vhgh)
  End Do
  plot = 1
!
  Return
 End Function plot_fvals
!
!
!
 Function plot_fval(img_n, &
          vlow, vhgh, vmin, vmax) &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use data_numbers,       Only: rp
  Use labels,             Only: l_fcasts
  Use arrays,             Only: fcast, fpls
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer :: plot
!
! Arguments
!
! Input arguments
  Integer, Intent(In) :: img_n ! - image number -
!
! - optional input scalars -
  Real(Kind=rp), Intent(In), Optional :: vlow ! - lowest value -
  Real(Kind=rp), Intent(In), Optional :: vhgh ! - highest value -
  Real(Kind=rp), Intent(In), Optional :: vmin ! - absolute minimum -
  Real(Kind=rp), Intent(In), Optional :: vmax ! - absolute maximum -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Plot forecast maps
  Select Case (img_n)
   Case (1)
     plot = plot_ydata(img_fval(img_n), icol_sk, 0, 1, fpls(:,ifm)%lower, 'Lower prediction limits', &
            vlow=vlow, vhgh=vhgh, vmin=vmin, vmax=vmax)
   Case (2)
     plot = plot_ydata(img_fval(img_n), icol_sk, 0, 1, fcast(:,ifm,0), Trim(l_fcasts%c), &
            vlow=vlow, vhgh=vhgh, vmin=vmin, vmax=vmax)
   Case (3)
     plot = plot_ydata(img_fval(img_n), icol_sk, 0, 1, fpls(:,ifm)%upper, 'Upper prediction limits', &
            vlow=vlow, vhgh=vhgh, vmin=vmin, vmax=vmax)
  End Select
  plot = 1
!
  Return
 End Function plot_fval
!
!
!
 Function replot_fval() &
          Result (plot)
!
! Modules
  Use arrays, Only: fcast
!
! Function type
  Integer :: plot
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal, MinVal
!
! Executable Statements
!
! Prompt for forecasts map title
  img_fval(ifvl)%lcontour = .false.
  If (img_fval(ifvl)%ymin == img_odds(ifvl)%ymax) Then
     plot = plot_fval(ifvl, &
            vlow=MinVal(fcast(:,:,0)), vhgh=MaxVal(fcast(:,:,0)))
  Else
     plot = plot_fval(ifvl, &
            vmin=img_fval(ifvl)%ymin, vmax=img_fval(ifvl)%ymax)
  End If
!
  Return
 End Function replot_fval
!
!
!
 Function plot_fpss() &
          Result (plot)
!
! Function type
  Integer plot
!
! Locals
!
! Local scalars
  Integer :: i ! - image index -
!
! Executable Statements
!
! Plot forecast maps
  Do i = 1, ng
     plot = plot_fps(i)
  End Do
  plot = 1
!
  Return
 End Function plot_fpss
!
!
!
 Function plot_fps(img_n, &
          vlow, vhgh, vmin, vmax) &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_p2
  Use data_numbers,       Only: rp, one
  Use labels,             Only: ls_cats
  Use fields,             Only: iffy, yfield
  Use categories,         Only: fps
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer plot
!
! Arguments
!
! Input arguments
  Integer, Intent(In) :: img_n ! - image number -
!
! - optional input scalars -
  Real(Kind=rp), Intent(In), Optional :: vlow ! - lowest value -
  Real(Kind=rp), Intent(In), Optional :: vhgh ! - highest value -
  Real(Kind=rp), Intent(In), Optional :: vmin ! - absolute minimum -
  Real(Kind=rp), Intent(In), Optional :: vmax ! - absolute maximum -
!
! Locals
!
! Local scalars
  Real(Kind=rp) :: tmiss ! - missing value -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Reset missing values
  tmiss = yfield(iffy)%rmiss
  yfield(iffy)%rmiss = -one
!
! Plot forecast maps
  plot = plot_ydata(img_fps(img_n), icol_p2, 2, 1, fps(:,ifm,img_n), 'Probabilities of '//Trim(ls_cats(img_n)%c), &
         vlow=vlow, vhgh=vhgh, vmin=vmin, vmax=vmax)
  yfield(iffy)%rmiss = tmiss
!
  Return
 End Function plot_fps
!
!
!
 Function replot_fps() &
          Result (plot)
!
! Function type
  Integer :: plot
!
! Executable Statements
!
! Prompt for forecasts map title
  img_fps(ifps)%lcontour = .false.
  If (img_fps(ifps)%ymin == img_fps(ifps)%ymax) Then
     plot = reset_img_fps()
  End If
  plot = plot_fps(ifps, &
         vmin=img_fps(ifps)%ymin, vmax=img_fps(ifps)%ymax)
!
  Return
 End Function replot_fps
!
!
!
 Function plot_oddss() &
          Result (plot)
!
! Modules
  Use data_numbers, Only: rp
  Use settings,     Only: iodds
  Use categories,   Only: odds
!
! Function type
  Integer plot
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
  Intrinsic MaxVal, MinVal
!
! Executable Statements
!
! Calculate data limits
  vlow = MinVal(odds(:,:,:,iodds+1))
  vhgh = MaxVal(odds(:,:,:,iodds+1))
!
! Plot forecast odds maps
  Do iodd = 1, ng
     plot = plot_odds(vlow=vlow,vhgh=vhgh)
  End Do
  plot = 1
!
  Return
 End Function plot_oddss
!
!
!
 Function plot_odds(vlow,vhgh,vmin,vmax) &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_sk
  Use data_numbers,       Only: rp, one
  Use settings,           Only: iodds
  Use fields,             Only: iffy, yfield
  Use categories,         Only: odds
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer plot
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In), Optional :: vlow ! - lowest value -
  Real(Kind=rp), Intent(In), Optional :: vhgh ! - highest value -
  Real(Kind=rp), Intent(In), Optional :: vmin ! - minimum value -
  Real(Kind=rp), Intent(In), Optional :: vmax ! - maximum value -
!
! Locals
!
! Local scalars
  Real(Kind=rp) :: tmiss ! - backup missing value -
!
  Character(Len=26) :: codds ! - odds / relative odds -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Reset missing values
  tmiss = yfield(iffy)%rmiss
  yfield(iffy)%rmiss = -one
!
! Plot forecast odds maps
  Select Case (iodds)
   Case (0)
     codds = 'Forecast odds'
   Case Default
     codds = 'Forecast relative odds'
  End Select
  plot = plot_ydata(img_odds(iodd), icol_sk, 0, 1, odds(:,ifm,iodd,iodds+1), Trim(codds), &
         vlow=vlow, vhgh=vhgh, vmin=vmin, vmax=vmax)
  yfield(iffy)%rmiss = tmiss
  plot = 1
!
  Return
 End Function plot_odds
!
!
!
 Function replot_odds() &
          Result (plot)
!
! Modules
  Use settings,   Only: iodds
  Use categories, Only: odds
!
! Function type
  Integer :: plot
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MaxVal, MinVal
!
! Executable Statements
!
! Plot map
  img_odds(iodd)%lcontour = .false.
  If (img_odds(iodd)%ymin == img_odds(iodd)%ymax) Then
     plot = plot_odds(vlow=MinVal(odds(:,:,:,iodds+1)), vhgh=MaxVal(odds(:,:,:,iodds+1)))
  Else
     plot = plot_odds(vmin=img_odds(iodd)%ymin, vmax=img_odds(iodd)%ymax)
  End If
!
  Return
 End Function replot_odds
!
!
!
 Function plot_hps() &
          Result (plot)
!
! Modules
  Use data_gui_constants, Only: icol_p3
  Use categories,         Only: hps
  Use plotting,           Only: plot_ydata
!
! Function type
  Integer plot
!
! Executable Statements
!
! Plot forecast maps
  plot = plot_ydata(img_hps(1), icol_p3, 0, 3, hps(:,ifm), ' ')
!
  Return
 End Function plot_hps
!
!
!
 Function select_fval() &
          Result (selectf)
!
! Opens Forecast ~ Series for a selected point on a forecast values map
!
! Modules
  Use graphics,    Only: select_location
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: selectf
!
! Executable Statements
!
! Select location
  selectf = select_location(img_fval(imgid+1-img_fval(1)%id), fcst_fser)
!
  Return
 End Function select_fval
!
!
!
 Function select_fps() &
          Result (selectf)
!
! Opens Forecast ~ Series for a selected point on a forecast probabilities map
!
! Modules
  Use graphics,    Only: select_location
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: selectf
!
! Executable Statements
!
! Select location
  selectf=select_location(img_fps(imgid+1-img_fps(1)%id),fcst_fser)
!
  Return
 End Function select_fps
!
!
!
 Function select_odds() &
          Result (selectf)
!
! Opens Forecast ~ Series for a selected point on a forecast odds map
!
! Modules
  Use graphics,    Only: select_location
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: selectf
!
! Executable Statements
!
! Select location
  selectf = select_location(img_odds(imgid+1-img_odds(1)%id), fcst_fser)
!
  Return
 End Function select_odds
!
!
!
 Function select_hps() &
          Result (selectf)
!
! Opens Forecast ~ Series for a selected point on a probabilistic forecast map
!
! Modules
  Use graphics, Only: select_location
!
! Function type
  Integer :: selectf
!
! Executable Statements
!
! Select location
  selectf = select_location(img_hps(1), fcst_fser)
!
  Return
 End Function select_hps
!
!
!
 Function f_title() &
          Result (title)
!
! Prompts for new title
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
  title = f_replot(get_title)
!
  Return
 End Function f_title
!
!
!
 Function f_limits() &
          Result (limits)
!
! Prompts for contour limits
!
! Modules
  Use plotting, Only: get_limits
!
! Function type
  Integer :: limits
!
! Executable Statements
!
! Prompt for contour limits
  limits = f_replot(get_limits)
!
  Return
 End Function f_limits
!
!
!
 Function f_replot(get_custom) &
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
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any
!
! Executable Statements
!
! Replot customised forecast diagrams
  If (img_fser(1)%id == imgid) Then ! - Forecast ~ Series -
     replot = replot_graph(img_fser(1), get_custom, plot_fser)
  Else If (img_fens(1)%id == imgid) Then ! - Forecast ~ Ensembles -
     replot = replot_graph(img_fens(1), get_custom, plot_fens)
  Else If (img_pexc(1)%id == imgid) Then ! - Forecast ~ Exceedance -
     replot = replot_graph(img_pexc(1), get_custom, replot_pexc)
  Else If (Any(img_fval(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Values -
     ifvl = imgid + 1 - img_fval(1)%id
     replot = replot_graph(img_fval(ifvl), get_custom, replot_fval)
  Else If (Any(img_fps(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Probabilities -
     ifps = imgid + 1 - img_fps(1)%id
     replot = replot_graph(img_fps(ifps), get_custom, replot_fps)
  Else If (img_hps(1)%id == imgid) Then ! - Forecast ~ Maps ~ Probabilistic -
     replot = replot_graph(img_hps(1), get_custom, plot_hps)
  Else If (Any(img_odds(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Odds -
     iodd = imgid + 1 - img_odds(1)%id
     replot = replot_graph(img_odds(iodd), get_custom, replot_odds)
  End If
!
  Return
 End Function f_replot
!
!
!
 Function f_export() &
          Result (export)
!
! Saves images
!
! Modules
  Use labels,      Only: cg_ensemble, &
                         l_fcast, &
                         ls_cats
  Use iofiles,     Only: ny
  Use settings,    Only: iv
  Use graphics,    Only: save_graphic
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: export
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any, Trim
!
! Executable Statements
!
! Prompt for montage title
  If (imgid < 0) Then
     export = f_export_montage(-imgid)
     Return     
  End If
!
! Save graph
  If (imgid == img_fser(1)%id) Then ! - Forecast ~ Series -
     export = save_graphic(Trim(l_fcast%c), img_fser(1), &
              imap=iv, mmap=ny)
  Else If (imgid == img_pexc(1)%id) Then ! - Forecast ~ Ensembles -
     export = save_graphic('Exceedance_probability', img_pexc(1), &
              imap=ifmy)
  Else If (imgid == img_fens(1)%id) Then ! - Forecast ~ Ensembles -
     export = save_graphic(Trim(cg_ensemble), img_fens(1), &
              imap=iv, mmap=ny)
  Else If (imgid == img_fval(1)%id) Then ! - Forecast ~ Maps ~ Values -
     ifvl = 1
     export = save_graphic('Lower_map', img_fval(ifvl), &
              imap=ifmy)
  Else If (imgid == img_fval(2)%id) Then
     ifvl = 2
     export = save_graphic('Forecast_map', img_fval(ifvl), &
              imap=ifmy)
  Else If (imgid == img_fval(3)%id) Then
     ifvl = 3
     export = save_graphic('Upper_map', img_fval(ifvl), &
              imap=ifmy)
  Else If (Any(img_fps(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Probabilities -
     ifps = imgid + 1 - img_fps(1)%id
     export = save_graphic(Trim(ls_cats(ifps)%c)//'_map', img_fps(ifps), &
              imap=ifmy)
  Else If (imgid == img_hps(1)%id) Then ! - Forecast ~ Maps ~ Probabilistic -
     export = save_graphic('Probabilistic_map', img_hps(1), &
              imap=ifmy)
  Else If (Any(img_odds(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Odds -
     iodd = imgid + 1 - img_odds(1)%id
     export = save_graphic(Trim(ls_cats(iodd)%c)//'_map', img_odds(iodd), &
              imap=ifmy)
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
  Use labels,   Only: l_fcasts, l_odds, l_probs
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
  If (imgid == img_fval(1)%id) Then ! - Forecast ~ Maps ~ Values -
     export = export_montage(l_fcasts%c, img_fval(1), img_fval(2), &
              lsquare=.true., img_3=img_fval(3))
  Else If (img_fps(1)%id == imgid) Then ! - Forecast ~ Maps ~ Probabilities -
     export = export_montage(l_probs%c, img_fps(1), img_fps(2), &
              lsquare=.true., img_3=img_fps(3))
  Else If (img_odds(1)%id == imgid) Then ! - Forecast ~ Maps ~ Odds -
     export = export_montage(l_odds%c, img_odds(1), img_odds(2), &
              lsquare=.true., img_3=img_odds(3))
  End If
!
  Return
 End Function f_export_montage
!
!
!
 Function help_fcast()
!
! Modules
  Use help,        Only: help_file
  Use windows_gui, Only: imgid
!
! Function type
  Integer :: help_fcast
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any
!
! Executable Statements
!
! Open context-specific help pages
  If (img_fser(1)%id == imgid) Then ! - Forecast ~ Series -
     Call help_file ('Tools_Forecast_Series.html')
  Else If (img_fens(1)%id == imgid) Then ! - Forecast ~ Ensembles -
     Call help_file ('Tools_Forecast_Ensemble.html')
  Else If (img_pexc(1)%id == imgid) Then ! - Forecast ~ Exceedance -
     Call help_file ('Tools_Forecast_Exceed.html')
  Else If (Any(img_fval(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Values -
     Call help_file ('Tools_Forecast_Values.html')
  Else If (Any(img_fps(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Probabilities -
     Call help_file ('Tools_Forecast_Probs.html')
  Else If (img_hps(1)%id == imgid) Then ! - Forecast ~ Maps ~ Probabilistic -
     Call help_file ('Tools_Forecast_Prob.html')
  Else If (Any(img_odds(:)%id == imgid)) Then ! - Forecast ~ Maps ~ Odds -
     Call help_file ('Tools_Forecast_Odds.html')
  End If
  help_fcast = 2
!
  Return
 End Function help_fcast
!
!
!
#endif
 Subroutine proj_read_forecasts (iin, ianal, igauss, ihind, pver, ifail)
!
! Reads forecast settings from project file
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Problem reading file
!    ifail =  2 Premature end of file
!    ifail =  4 Problem allocating memory
!
! Modules
  Use data_numbers,       Only: sp, tolv
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use arrays,             Only: fcast, fcastt, fpls, pev, xvp, yfit, z, ziny
  Use settings,           Only: iev, nf
  Use pcs,                Only: fsx, fsxc
  Use missing,            Only: kfuse
  Use categories,         Only: fps, odds
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin    ! - input file unit number -
  Integer, Intent(In) :: ianal  ! - analysis identifier -
  Integer, Intent(In) :: igauss ! - transform Y data flag -
  Integer, Intent(In) :: ihind  ! - hindcasts flag -
!
  Real(Kind=sp), Intent(In) :: pver ! - project version -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Read forecast settings
  Read (Unit=iin, Err=1, End=2) ifc, ifcast
  ifm = 1
!
! Read forecast results
  If (ifc > 0) Then
     If (pver < 15.06_sp-tolv) Then
        ifc = 0
        GoTo 3
     End if
     If (init_fcasts(ianal) /= 0) GoTo 4
     Read (Unit=iin, Err=1, End=2) z
     If (ihind == 1) Then
        If (iev == 3) Read (Unit=iin, Err=1, End=2) yfit
        Read (Unit=iin, Err=1, End=2) fcast
        If (igauss > 0) Read (Unit=iin, Err=1, End=2) fcastt
        Read (Unit=iin, Err=1, End=2) pev
        Read (Unit=iin, Err=1, End=2) fpls
     End If
     Read (Unit=iin, Err=1, End=2) fps
     Read (Unit=iin, Err=1, End=2) odds(:,:,:,2)
     Read (Unit=iin, Err=1, End=2) odds(:,:,:,1)
     Select Case (ianal)
      Case (ia_cca)
        Read (Unit=iin, Err=1, End=2) xvp
        Read (Unit=iin, Err=1, End=2) fsx
        Read (Unit=iin, Err=1, End=2) fsxc
      Case (ia_pcr, ia_mlr)
        Read (Unit=iin, Err=1, End=2) xvp
        Read (Unit=iin, Err=1, End=2) fsx
      Case (ia_gcm)
        Read (Unit=iin, Err=1, End=2) xvp
        Read (Unit=iin, Err=1, End=2) ziny
      Case (ia_pfv)
        Read (Unit=iin, Err=1, End=2) ziny
     End Select
     Read (Unit=iin, Err=1, End=2) kfuse
     Call get_fprds (nf, 2)
  End If
!
3 ifail = 0
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
! - problem allocating memory -
4 ifail = 4
  Return
!
 End Subroutine proj_read_forecasts
!
!
!
 Subroutine proj_write_forecasts (iout, ianal, igauss, ihind, ifail)
!
! Writes forecast settings to project file
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Problem writing to file
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
  Use arrays,             Only: fcast, fcastt, fpls, pev, xvp, yfit, z, ziny
  Use settings,           Only: iev
  Use pcs,                Only: fsx, fsxc
  Use missing,            Only: kfuse
  Use categories,         Only: fps, odds
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout   ! - output file unit number -
  Integer, Intent(In) :: ianal  ! - analysis identifier -
  Integer, Intent(In) :: igauss ! - transform Y data flag -
  Integer, Intent(In) :: ihind  ! - hindcasts flag -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write forecast settings
  Write (Unit=iout, Err=1) ifc, ifcast
  If (ifc > 0) Then
     Write (Unit=iout, Err=1) z
     If (ihind == 1) Then
        If (iev == 3) Write (Unit=iout, Err=1) yfit
        Write (Unit=iout, Err=1) fcast
        If (igauss > 0) Write (Unit=iout, Err=1) fcastt
        Write (Unit=iout, Err=1) pev
        Write (Unit=iout, Err=1) fpls
     End If
     Write (Unit=iout, Err=1) fps
     Write (Unit=iout, Err=1) odds(:,:,:,2)
     Write (Unit=iout, Err=1) odds(:,:,:,1)
     Select Case (ianal)
      Case (ia_cca)
        Write (Unit=iout, Err=1) xvp
        Write (Unit=iout, Err=1) fsx
        Write (Unit=iout, Err=1) fsxc
      Case (ia_pcr, ia_mlr)
        Write (Unit=iout, Err=1) xvp
        Write (Unit=iout, Err=1) fsx
      Case (ia_gcm)
        Write (Unit=iout, Err=1) xvp
        Write (Unit=iout, Err=1) ziny
      Case (ia_pfv)
        Write (Unit=iout, Err=1) ziny
     End Select
     Write (Unit=iout, Err=1) kfuse
  End If
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = close_fcast()
  ifail = 1
  Return
!
 End Subroutine proj_write_forecasts
!
!
!
 Function close_fcast()
!
! Frees memory allocated for forecasts
!
! Modules
  Use arrays,     Only: fcast, fcastt, fpls, pev, xc, xvp, yfit, z, ziny, zyp
  Use pcs,        Only: fsx, fsxc
  Use missing,    Only: kfuse
  Use categories, Only: fps, odds
!
! Function type
  Integer :: close_fcast
!
! Executable Statements

! Free memory
  close_fcast=close_pexc()
  close_fcast=close_hps()
  If (Allocated(cdates)) Deallocate (cdates)
  If (Allocated(fprd))   Deallocate (fprd)
  If (Allocated(zyp))    Deallocate (zyp)
  If (Allocated(kfuse))  Deallocate (kfuse)
  If (Allocated(ziny))   Deallocate (ziny)
  If (Allocated(fsxc))   Deallocate (fsxc)
  If (Allocated(fsx))    Deallocate (fsx)
  If (Allocated(xvp))    Deallocate (xvp)
  If (Allocated(odds))   Deallocate (odds)
  If (Allocated(fps))    Deallocate (fps)
  If (Allocated(fpls))   Deallocate (fpls)
  If (Allocated(pev))    Deallocate (pev)
  If (Allocated(fcastt)) Deallocate (fcastt)
  If (Allocated(fcast))  Deallocate (fcast)
  If (Allocated(yfit))   Deallocate (yfit)
  If (Allocated(xc))     Deallocate (xc)
  If (Allocated(z))      Deallocate (z)
  close_fcast = 0
!
  Return
 End Function close_fcast
!
!
!
 Function close_pexc()
!
! Frees memory allocated for probabilities of exceedence
!
! Modules
  Use arrays, Only: c, e, f, tex
!
! Function type
  Integer :: close_pexc
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory
  If (Allocated(tex)) Deallocate (tex)
  If (Allocated(c))   Deallocate (c)
  If (Allocated(e))   Deallocate (e)
  If (Allocated(f))   Deallocate (f)
  close_pexc = 0
!
  Return
 End Function close_pexc
!
!
!
 Function close_hps()
!
! Frees memory allocated for highest probability mapping
!
! Modules
  Use categories, Only: hps
!
! Function type
  Integer :: close_hps
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory
  If (Allocated(hps)) Deallocate (hps)
  close_hps = 0
!
  Return
 End Function close_hps
End Module forecast
