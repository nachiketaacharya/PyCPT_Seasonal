! Author: Simon Mason
Module analysis
!
! Declarations

! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: analysis_flags, close_analysis, init_analysis1, init_analysis2, init_results, &
            num_calcs, perform_analysis, proj_read_analysis, proj_write_analysis, reset_analysis
!
! Scalars
!
! Integer scalars
  Integer, Public :: ianal = 0 ! - analysis identifier (1=CCA, 2=PCR, 3=MLR, 4=GCM, 5=PFV) -
  Integer, Public :: ianaln    ! - intended new analysis identifier (1=CCA, 2=PCR, 3=MLR, 4=GCM, 5=PFV) -
  Integer, Public :: iaction   ! - ready for calculation flag -
  Integer, Public :: icalc     ! - calculation completed flag -
  Integer, Public :: ihind     ! - hindcasts flag -
  Integer, Public :: jcalc     ! - calculation incomplete flag -
  Integer, Public :: nopt      ! - total number of optimization steps -
! - grey menu flag complements -
  Integer, Public :: jm_gauss ! - transform Y data grey menu flag -
!
! Logical scalars
  Logical, Public :: lconvert ! - convert units? -
!
! Arrays
!
! Logical arrays
  Logical, Dimension(:,:), Allocatable, Public :: lopt ! - optimized? -
!
Contains
!
!
 Subroutine reset_analysis ()
!
! Resets analysis flags
!
! Executable Statements
!
! Initialise analysis settings
  Call analysis_flags ('off')
  iaction = 0
  lconvert = .false.
!
  Return
 End Subroutine reset_analysis
!
!
!
 Subroutine analysis_flags (onoff)
!
! Switches calculation flags on or off
!
! Modules
  Use data_numbers,       Only: zero, one
  Use data_cpt_constants, Only: ia_pfv
  Use progress_meter,     Only: lcalc, lstop, &
                                close_progress_meter, set_progress_meter
  Use settings,           Only: iretro
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: onoff
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Random_Seed
!
! Executable Statements
!
! Set flags on
  Select Case (onoff)
   Case ('on', 'ON', 'On')
     icalc = 0                       ! - calculation completed flag -
     lcalc = .true.                  ! - calculation active flag -
     lstop = .false.                 ! - stop calculation flag -
     Call set_progress_meter (zero)  ! - progress meter -
!
! - initialise random number generator -
     Call Random_Seed ()
!
! Set flags to indicate end of calculation
   Case ('end', 'END', 'End')
     icalc = 1                       ! - calculation completed flag -
     jcalc = 0                       ! - calculation incomplete flag -
     lcalc = .false.                 ! - calculation active flag -
     Call set_progress_meter (one)   ! - progress meter -
!
! Set all flags off
   Case ('off', 'OFF', 'Off')
     icalc = 0                       ! - calculation completed flag -
     jcalc = 1                       ! - calculation incomplete flag -
     lcalc = .false.                 ! - calculation active flag -
     lstop = .false.                 ! - stop calculation flag -
     If (ianal /= ia_pfv) iretro = 0 ! - retroactive calculation flag -
     nopt = 0                        ! - number of optimizations -
     Call close_progress_meter ()
  End Select
!
  Return
 End Subroutine analysis_flags
!
!
!
 Subroutine init_analysis1 (ianalc, lproj, ng, nret, nopt, nlt, nlg, ndom, ifail)
!
! Initialises memory and settings to perform EOF prefiltering and model fitting
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm
  Use arrays,             Only: b, b0, rwk, iwk, liwk, lrwk, xc, xiny, xm, xsd, xvp, yc, yhat, yhatt, ym, ysd, yopt, yret, yrett, &
                                yrpls
  Use iofiles,            Only: mxa, mya
  Use settings,           Only: iensc, igauss, ig0, imc, iretro, hcw, lcw, nu, nur
  Use space,              Only: init_space
  Use fields,             Only: first_field
  Use scores,             Only: gm
  Use pcs,                Only: lxt, lxyt, mcc, mxe, mye, ncc, nxe, nye, &
                                init_svds, init_svdwk
  Use categories,         Only: efps, rfps, rodds
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianalc ! - analysis identifier -
  Integer, Intent(In) :: ng     ! - number of categories -
  Integer, Intent(In) :: nret   ! - number of retroactive iterations -
  Integer, Intent(In) :: nlt    ! - number of latitudes -
  Integer, Intent(In) :: nlg    ! - number of longitudes -
  Integer, Intent(In) :: ndom   ! - number of domains -
!
  Logical, Intent(In) :: lproj ! - restore from project? -
!
! Output scalars
  Integer, Intent(Out) :: nopt  ! - number of optimization steps -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: lnn  ! - minimum of ndom, and nu -
  Integer :: mxt  ! - maximum of mxa, and nu -
  Integer :: mxyt ! - maximum of mxa, mya, and nu -
  Integer :: mopt ! - number of optimization steps -
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Allocated, Max, Min
!
! Executable Statements
!
! Free existing memory
  If (Allocated(rwk)) Deallocate (rwk)
!
! Allocate additional memory
! - means -
  Select Case (ianalc)
   Case (ia_cca, ia_pcr, ia_mlr)
     Allocate (xm(mxa), Stat=ifail)
   Case (ia_gcm)
     Allocate (xm(mya), Stat=ifail)
  End Select
  If (ifail /= 0) GoTo 1
  Allocate (ym(mya), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - standard deviations -
  Select Case (ianalc)
   Case (ia_cca, ia_pcr, ia_mlr)
     Allocate (xsd(mxa), Stat=ifail)
   Case (ia_gcm)
     Allocate (xsd(mya), Stat=ifail)
  End Select
  If (ifail /= 0) GoTo 1
  Allocate (ysd(mya), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - cross-validated predictions -
  Allocate (yhat(mya,nu), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - cross-validated transformed predictions -
  If (igauss > 0) Then
     Allocate (yhatt(mya,nu), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
! - optimization results -
  If (nret > 0) Then
     Select Case (ianalc)
      Case (ia_cca)
        nopt = get_nopt(nxe, mxe, &
               nye=nye, mye=mye, ncc=ncc, mcc=mcc)
      Case (ia_pcr)
        nopt = get_nopt(nxe, mxe)
      Case (ia_mlr)
        nopt = 1
      Case Default
        nopt = 0
     End Select
     If (nopt > 0) Then
        Allocate (lopt(nopt,nret*iretro+1), Stat=ifail)
        If (ifail /= 0) GoTo 1
        Allocate (gm(nopt,nret*iretro+1), Stat=ifail)
        If (ifail /= 0) GoTo 1
     End If
  Else
     nopt = 0
  End If
! - ensemble forecast probabilities -
  If (ianalc == ia_gcm) Then
     If (iensc > 1) Then
        Allocate (efps(mya,nu,ng), Stat=ifail)
        If (ifail /= 0) GoTo 1
     End If
  End If
! - retroactive predictions -
  If (iretro == 1) Then
     Allocate (yret(mya,nur), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - transformed retroactive predictions -
     If (igauss > 0) Then
        Allocate (yrett(mya,nur), Stat=ifail)
        If (ifail /= 0) GoTo 1
     End If
! - retroactive prediction limits -
     Allocate (yrpls(mya,nur), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - retroactive forecast probabilities -
     Allocate (rfps(mya,nur,ng), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - retroactive odds -
     Allocate (rodds(mya,nur,ng,2), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
! - SVDs -
  Select Case (ianalc)
   Case (ia_cca)
     ifail = init_svds(ianalc)
     If (ifail /= 0) GoTo 1
   Case (ia_pcr, ia_mlr)
     ifail = init_svds(ianalc)
     If (ifail /= 0) GoTo 1
! - regression constants -
     Allocate (b0(mya,1), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - regression coefficients -
     Allocate (b(mxa,mya), Stat=ifail)
     If (ifail /= 0) GoTo 1
   Case (ia_gcm)
     If (ndom > 1) Then
        ig0 = 0
     Else
        ig0 = 1
     End If
! - interpolated X values -
     Allocate (xiny(mya,nu,ig0:ndom), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - regression constants -
     Allocate (b0(ig0:ndom,mya), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - regression coefficients -
     Allocate (b(ig0:ndom,mya), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - nearest latitudes and longitudes -
     ifail = init_space(nlt,nlg,ndom)
     If (ifail /= 0) GoTo 1
  End Select
! - copy of training period independent data -
  Select Case (ianalc)
   Case (ia_cca, ia_pcr, ia_mlr)
     Allocate (xc(mxa,nu,1), Stat=ifail)
   Case (ia_gcm)
     Allocate (xc(mya,nu,ig0:ndom), Stat=ifail)
  End Select
  If (ifail /= 0) GoTo 1
! - copy of training period dependent data -
  Allocate (yc(mya,nu), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - prediction error variance for retroactive forecasts -
  If (iretro == 1) Then
     Select Case (ianalc)
      Case (ia_cca, ia_pcr, ia_mlr)
        Allocate (xvp(1,nur), Stat=ifail)
      Case (ia_gcm)
        Allocate (xvp(mya,nur), Stat=ifail)
     End Select
     If (ifail /= 0) GoTo 1
  End If
! - SVD workspace -
  If (ianalc == ia_cca) ifail = init_svdwk()
  If (ifail /= 0) GoTo 1
!
! Calculate required workspace
  Select Case (ianalc)
   Case (ia_cca)
     mxyt = Max(mxa, mya, 2*nu - 3*lxyt)
     lrwk = lxyt*(6 + 4*lxyt) + mxyt
     liwk = 8*lxyt
   Case (ia_pcr, ia_mlr)
     mxt = Max(mxa, 2*nu - 3*lxt)
     lrwk = lxt*(6 + 4*lxt) + Max(mxt, mya)
     liwk = 8*lxt
   Case (ia_gcm)
     lxt = Min(ndom, nu)
     mxt = Max(mya, 2*nu - 3*lxt)
     lrwk = lxt*(6 + 4*lxt) + Max(mxt, mya)
     liwk = 8*lxt
     If (imc == 4) Then ! - add additional space for model combination -
        lnn = Min(ndom, nu)
        lrwk = Max(lrwk, 3*lnn + Max(Max(ndom, nu), 4*lnn*(lnn + 1)))
        liwk = Max(lrwk, 8*lnn)
     End If
  End Select
!
! Allocate workspace
  Allocate (iwk(liwk), Stat=ifail)
  If (ifail /= 0) GoTo 1
  Allocate (rwk(lrwk), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - hindcasts for optimization -
  If (lproj) Then
     mopt = 1
  Else
     mopt = Max (1, nopt)
  End If
  Select Case (igauss)
   Case (0)
     Allocate (yopt(mya,nu,mopt), Stat=ifail)
   Case Default
     Allocate (yopt(mya,nu,Max(2, mopt)), Stat=ifail)
  End Select
  If (ifail /= 0) GoTo 1
!
! Calculate half cross-validation window
  hcw = (lcw + 1)/2
!
! Identify initial grid point
  Call first_field ()
!
  Return
!
! Errors
1 ifail = 1
!
  Return
 End Subroutine init_analysis1
!
!
!
 Subroutine init_analysis2 (ifail, &
            nco)
!
! Initialises additional memory for analysis
!
! Modules
  Use arrays,   Only: yopt, yt
  Use iofiles,  Only: mya
  Use settings, Only: igauss, nu
  Use pcs,      Only: init_ccamaps
!
! Arguments
!
! Input scalars
  Integer, Intent(In), Optional :: nco ! - optimal number of CCA modes -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Allocated, Present
!
! Executable Statements
!
! Deallocate unneeded workspace
  ifail = 0
  If (Allocated(yopt)) Deallocate (yopt)
!
! Allocate additional memory
! - transformed response data -
  If (igauss > 0) Then
     Allocate (yt(mya,nu), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
!
! Allocate additional memory for CCA maps
  If (Present(nco)) Then
     ifail = init_ccamaps(nco)
     If (ifail /= 0) GoTo 1
  End If
!
  Return
!
! Errors
1 ifail = 1
!
  Return
 End Subroutine init_analysis2
!
!
!
 Subroutine init_results (ianalc, nlt, nlg, ifail)
!
! Initialises memory for results
!
! Modules
  Use data_cpt_constants, Only: ia_pfv
  Use arrays,             Only: awk, rwk, iwk, xc, xvp, yc, ypc
  Use iofiles,            Only: mya, nx, ny
  Use settings,           Only: iretro, nt, nu, nur
  Use pcs,                Only: ce
  Use categories,         Only: init_categories
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianalc ! - analysis identifier -
  Integer, Intent(In) :: nlt    ! - maximum number of latitudes -
  Integer, Intent(In) :: nlg    ! - maximum number of longitudes -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Allocated, Max
!
! Executable Statements
!
! Deallocate unneeded workspace
  If (Allocated(rwk)) Deallocate (rwk)
  If (Allocated(iwk)) Deallocate (iwk)
  If (Allocated(ce))  Deallocate (ce)
  If (Allocated(xvp)) Deallocate (xvp)
  If (Allocated(yc))  Deallocate (yc)
  If (Allocated(xc))  Deallocate (xc)
  If (Allocated(ypc)) Deallocate (ypc)
!
! Allocate additional memory
  Select Case (ianalc)
   Case Default
     ifail = init_categories(ianalc, iretro, mya, nu, nur)
   Case (ia_pfv)
     ifail = init_categories(ianalc, iretro, mya, nu, nt)
  End Select
  If (ifail /= 0) GoTo 1
! - workspace -
  Allocate (awk(Max(nx, ny, nlg, nlt)), Stat=ifail)
  If (ifail /= 0) GoTo 1
  Allocate (rwk(Max(2*nu, nx, ny, nlg, nlt)), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
  Return
!
! Errors
1 ifail = 1
!
  Return
 End Subroutine init_results
!
!
!
 Function perform_analysis() &
          Result (perform)
!
! Performs CCA, PCR, or GCM
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Memory allocation problem
!    ifail =  2 Calculation error
!    ifail = -1 Terminated
!
! Modules
  Use data_numbers,       Only: rp, one
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ng
  Use data_io_constants,  Only: ffmts, iout, lprd
  Use labels,             Only: cg_to_l, cg_tperiod
  Use time,               Only: pprd, &
                                Operator(+), &
                                get_cdate
  Use gui,                Only: print_text, print_warning
  Use progress_meter,     Only: lstop, &
                                set_progress_increment
  Use arrays,             Only: b, b0, iusex, iusey, x, xiny, xm, xvp, y, yhat, yhatt, ym, yopt, yret, yrett, yrpls, yt, &
                                convert_units
  Use iofiles,            Only: gifile, mxa, mya, xfile, yfile, &
                                file_reset
  Use settings,           Only: clf, iensc, iev, igauss, ig0, imc, intp, iretro, izero, lcw, lzero, ncv, nret, nretro, nt, nt1, &
                                nu, nur, nu1
  Use space,              Only: get_interpolated, set_nearest_grids
  Use fields,             Only: xfield, yfield
  Use season,             Only: lb, lprst
  Use pcs,                Only: bz, eofx, eofy, hx_map, hy_map, hx_ser, hy_ser, iec, ieofx, ieofy, iey, mcc, mxe, mye, mu, &
                                ncc, nco, nccu, npx, npy, nxe, nxo, nye, nyo, r, s, svx, svy, tsx, tsy
  Use missing,            Only: ktuse
  Use categories,         Only: pobs, rfps, rodds, tobst, &
                                calc_reg_probs
  Use climate,            Only: ave, sdev
  Use ensemble,           Only: close_combine, get_combined, init_combine
  Use cca,                Only: full_cca, cv_cca, cca_predict
  Use pcr,                Only: full_pcr, cv_pcr, pcr_predict
  Use gcm,                Only: full_gcm, cv_gcm, gcm_predict
!
! Function type
  Integer :: perform
!
! Locals
!
! Local scalars
  Integer :: iml   ! - model index -
  Integer :: ir    ! - index of current retroactive step -
  Integer :: irl   ! - index of last retroactive training-period date -
  Integer :: ir1   ! - index of first retroactive date -
  Integer :: irn   ! - index of last retroactive date -
  Integer :: iru   ! - index of last new date in retroactive training period -
  Integer :: mlt   ! - maximum number of latitudes -
  Integer :: mlg   ! - maximum number of longitudes -
  Integer :: nr    ! - total number of retroaoctive steps in current cycle -
  Integer :: ifail ! - error indicator -
!
  Character(Len=lprd) :: cprd1 ! - period -
  Character(Len=lprd) :: cprd2 ! - period -
  Character(Len=lprd) :: cprd3 ! - period -
!
  Type(pprd) :: prdc ! - current period -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count, MatMul, MaxVal, Min, Real, Trim
!
! Executable Statements
!
! Check for consistent units
  perform = -1
  lconvert = .false.
  Select Case (ianal)
   Case (ia_gcm)
     Do iml = 1, xfile%nms
        If (Trim(xfield(iml)%var) == Trim(yfield(1)%var)) Then
           If (Trim(xfield(iml)%unit) /= 'none' .and. Trim(yfield(1)%unit) /= 'none') Then
              If (Trim(xfield(iml)%unit) /= Trim(yfield(1)%unit)) Then
                 Call convert_units (xfield(iml)%unitc, yfield(1)%unit, x(1:mxa,:,1), ifail)
                 If (ifail == 0) Then
                    lconvert = .true.
                 Else
                    Call print_warning ('Unable to convert from '//Trim(xfield(iml)%unit)//' to '//Trim(yfield(1)%unit), &
                         lcursor=.true.)
                 End If
              End If
           End If
        End If
     End Do
     Select Case (yfile%idstr)
      Case (1)
        mlt = MaxVal(yfield(:)%region%nlts)
        mlg = MaxVal(yfield(:)%region%nlgs)
      Case (2)
        mlt = MaxVal(yfield(:)%nva)
        mlg = MaxVal(yfield(:)%nva)
     End Select
   Case Default
     mlt = 0
     mlg = 0
  End Select
!
! Initialise analysis
  Call init_analysis1 (ianal, .false., ng, nret, nopt, mlt, mlg, xfile%nms, ifail)
  If (ifail /= 0) Then
     perform = 1
     Return
  End If
!
! Initialise combination
  Select Case (ianal)
   Case (ia_gcm)
     If (imc == 4 .and. xfile%nms > 1) Then
        Call init_combine (xfile%nms, nt, ifail)
        If (ifail /= 0) Then
           perform = 1
           GoTo 1
        End If
     End If
  End Select
!
! Estimate number of calculations
  If (iev /= 3) Then
     Call set_progress_increment &
          (one/Real(num_calcs(ianal, iretro, igauss, izero, 0, mya, nu, nret, nretro, nu1, nur, nopt), Kind=rp), .false.)
  Else
     Call set_progress_increment &
          (one/Real(num_calcs(ianal, iretro, igauss, izero, 1, mya, nu, nret, nretro, nu1, nur, nopt), Kind=rp), .false.)
  End If
!
! Calculate optimization settings
  Select Case (ianal)
   Case (ia_cca) ! - CCA -
     If (nopt == 1) Then
        nopt = 0
        nxo = nxe
        nyo = nye
        nco = ncc
     End If
     iey = 1
     iec = 1
   Case (ia_pcr) ! - PCR -
     If (nopt == 1) Then
        nopt = 0
        nxo = nxe
     End If
   Case (ia_mlr) ! - MLR -
     nopt = 0
     nxe = Min(mxa, nu-1)
     mxe = nxe
   Case (ia_gcm) ! - GCM -
     nopt = 0
!
! Identify nearest grids
     Call set_nearest_grids (yfile%idstr, intp, xfile%nms, yfile%nfl, xfield(:)%nlt, xfield(:)%nlg, yfield(1:)%nva, &
          xfield(:)%region, yfield(1:)%region, ifail)
     Select Case (ifail)
      Case (0)
        Continue
      Case (-1)
        GoTo 1
      Case Default  
        perform = ifail + 2
        GoTo 1
     End Select
!
! Calculate interpolated values
     Call get_interpolated (ianal, xfile%idstr, yfile%idstr, intp, mya, iusex, iusey, xfile%nms, yfile%nfl, nu, xfile%ngs, &
          xfield(:)%nva, yfield(1:)%nva, x, xfield(:)%rmiss, xfield(:)%region, yfield(1:)%region, xiny(:,:,1:))
  End Select
!
! Open goodness index file
  If (gifile%cfile%lset) Then
     Open (Unit=iout, File=gifile%cfile%ffile, Access='sequential', Action='write', Form=ffmts(gifile%cfile%ffmt%ifmt), &
           Iostat=ifail, Status='unknown')
     If (ifail /= 0) Call file_reset (gifile%cfile%lset, gifile%cfile%fname)
  End If
!
! Perform retroactive calculations
  prdc = yfile%prd1 + (yfile%it1 - 1)
  cprd1 = get_cdate(prdc, 1, 2)
! - optimise model for successive training periods -
  ir = 0
  If (iretro == 1) Then
     Print *, 'Calculating retroactive forecasts ...'
     ir1 = 1
     irn = 0
     iru = nu1
     retro: Do irl = nt1, nt-1, nretro
        ir = ir + 1
        irn = irn + Count(ktuse(irl+1:Min(irl+nretro, nt)))
        nr = irn + 1 - ir1
        If (nr > 0) Then
           ncv = iru - lcw
           prdc = yfile%prd1 + (yfile%it1 + irl - 2)
           cprd2 = get_cdate(prdc, 1, 2)
           Print *, Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2)
           Print *, ' '
           If (gifile%cfile%lset) Then
              Call print_text (Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2), &
                   iout=iout)
              Call print_text (' ', &
                   iout=iout)
           End If
           If (nopt > 0) Then
              Select Case (igauss)
               Case (0)
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cv_cca (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
                         yhat=yhat, nxo=nxo, nyo=nyo, nco=nco)
                  Case (ia_pcr, ia_mlr)
                    Call cv_pcr (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ir, lopt, &
                         yhat=yhat, nxo=nxo)
                 End Select
               Case Default
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cv_cca (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
                         yhatt=yhat, nxo=nxo, nyo=nyo, nco=nco)
                  Case (ia_pcr, ia_mlr)
                    Call cv_pcr (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ir, lopt, &
                         yhatt=yhat, nxo=nxo)
                 End Select
              End Select
           Else
              Select Case (igauss)
               Case (0)
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cv_cca (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
                         yhat=yhat) 
                  Case (ia_pcr, ia_mlr)
                    Call cv_pcr (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ir, lopt, &
                         yhat=yhat)
                  Case (ia_gcm)
                    Call cv_gcm (imc, xfile%nms, iru, ncv, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:), &
                         yhat=yhat)
                 End Select
               Case Default
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cv_cca (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
                         yhatt=yhat)
                  Case (ia_pcr, ia_mlr)
                    Call cv_pcr (iru, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), &
                         ieofx, nxe, mxe, ir, lopt, &
                         yhatt=yhat)
                  Case (ia_gcm)
                    Call cv_gcm (imc, xfile%nms, iru, ncv, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:), &
                         yhatt=yhatt)
                 End Select
              End Select
              If (ianal == ia_mlr) nxo = nxe
           End If
           If (lstop) GoTo 1
! - construct full model -
           Select Case (igauss)
            Case (0)
              Select Case (ianal)
               Case (ia_cca)
                 Call full_cca (iru, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, ieofy, nyo, svy, eofy, tsy, &
                      nco, mu, r, s, ifail)
               Case (ia_pcr, ia_mlr)
                 Call full_pcr (iru, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, bz, b, npx, ifail)
               Case (ia_gcm)
                 Call full_gcm (imc, xfile%nms, iru, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:))
              End Select
            Case Default
              Select Case (ianal)
               Case (ia_cca)
                 Call full_cca (iru, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, ieofy, nyo, svy, eofy, tsy, &
                      nco, mu, r, s, ifail, &
                      yt=yopt(:,:,2))
               Case (ia_pcr, ia_mlr)
                 Call full_pcr (iru, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, bz, b, npx, ifail, &
                      yt=yopt(:,:,2))
               Case (ia_gcm)
                 Call full_gcm (imc, xfile%nms, iru, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:), &
                      yt=yopt(:,:,2))
              End Select
           End Select
           If (lstop) GoTo 1
           If (ifail /= 0) Then
              perform = 2
              GoTo 1
           End If
           If (npx < nxo) nxo = npx
           If (npy < nyo) nyo = npy
! - produce retroactive forecast -
           prdc = yfile%prd1 + (yfile%it1 + irl - 1)
           cprd3 = get_cdate(prdc, 1, 2)
           If (nr > 1) Then
              prdc = yfile%prd1 + (yfile%it1 + irl + nr - 2)
              cprd2 = get_cdate(prdc, 1, 2)
              Print *, 'Calculating retroactive forecasts for '//Trim(cprd3)//' - '//Trim(cprd2)
           Else
              Print *, 'Calculating retroactive forecasts for '//Trim(cprd3)
           End If
           If (xfile%nms > 1) Call get_combined (mya, xfile%nms, b0(1:,:), b(1:,:), xiny(1:mya,nu1+ir1:nu1+irn,1:), &
                                   xiny(1:mya,nu1+ir1:nu1+irn,0))
           Select Case (igauss)
            Case (0)
              If (.not.lprst) Then
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cca_predict (nr, mxa, mya, nxo, nyo, nco, ieofx, lzero, &
                         x(1:mxa,nu1+ir1:nu1+irn,1), ieofy, yret(:,ir1:irn),     &
                         nt=iru, xvp=xvp(1,ir1:irn))
                  Case (ia_pcr, ia_mlr)
                    Call pcr_predict (nr, mxa, mya, ieofx, lzero, x(1:mxa,nu1+ir1:nu1+irn,1), b, yret(:,ir1:irn), &
                         nt=iru, nxe=nxo, xvp=xvp(1,ir1:irn))
                  Case (ia_gcm)
                    Call gcm_predict (nr, mya, lzero, xiny(1:mya,nu1+ir1:nu1+irn,ig0), b0(ig0,:), b(ig0,:), yret(:,ir1:irn), &
                         nt=iru, xvp=xvp(:,ir1:irn))
                 End Select
              Else
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cca_predict (nr, mxa, mya, nxo, nyo, nco, ieofx, lzero, &
                         x(1:mxa,nu1+ir1:nu1+irn,1), ieofy, yret(:,ir1:irn),     &
                         nt=iru, xvp=xvp(1,ir1:irn), yprst=y(1:mya,nu1+ir1:nu1+irn,0))
                  Case (ia_pcr, ia_mlr)
                    Call pcr_predict (nr, mxa, mya, ieofx, lzero, x(1:mxa,nu1+ir1:nu1+irn,1), b, yret(:,ir1:irn), &
                         nt=iru, nxe=nxo, xvp=xvp(1,ir1:irn), yprst=y(1:mya,nu1+ir1:nu1+irn,0))
                  Case (ia_gcm)
                    Call gcm_predict (nr, mya, lzero, xiny(1:mya,nu1+ir1:nu1+irn,ig0), b0(ig0,:), b(ig0,:), yret(:,ir1:irn), &
                         nt=iru, xvp=xvp(:,ir1:irn), yprst=y(1:mya,nu1+ir1:nu1+irn,0))
                 End Select
              End If
            Case Default
              If (.not.lprst) Then
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cca_predict (nr, mxa, mya, nxo, nyo, nco, ieofx, lzero, &
                         x(1:mxa,nu1+ir1:nu1+irn,1), ieofy, yret(:,ir1:irn),     &
                         nt=iru, xvp=xvp(1,ir1:irn), nc=iru, tclim1=y(1:mya,1:iru,1), fcastt=yrett(:,ir1:irn))
                  Case (ia_pcr, ia_mlr)
                    Call pcr_predict (nr, mxa, mya, ieofx, lzero, x(1:mxa,nu1+ir1:nu1+irn,1), b, yret(:,ir1:irn), &
                         nt=iru, nxe=nxo, xvp=xvp(1,ir1:irn), nc=iru, tclim1=y(1:mya,1:iru,1), fcastt=yrett(:,ir1:irn))
                  Case (ia_gcm)
                    Call gcm_predict (nr,mya,lzero,xiny(1:mya,nu1+ir1:nu1+irn,ig0),b0(ig0,:),b(ig0,:),yret(:,ir1:irn), &
                         nt=iru, xvp=xvp(:,ir1:irn), nc=iru, tclim1=y(1:mya,1:iru,1), fcastt=yrett(:,ir1:irn))
                 End Select
              Else
                 Select Case (ianal)
                  Case (ia_cca)
                    Call cca_predict (nr, mxa, mya, nxo, nyo, nco, ieofx, lzero, &
                         x(1:mxa,nu1+ir1:nu1+irn,1), ieofy, yret(:,ir1:irn),     &
                         nt=iru, xvp=xvp(1,ir1:irn), nc=iru, tclim1=y(1:mya,1:iru,1), tclim0=y(1:mya,1:iru,0), &
                         fcastt=yrett(:,ir1:irn), yprst=y(1:mya,nu1+ir1:nu1+irn,0))
                  Case (ia_pcr, ia_mlr)
                    Call pcr_predict (nr, mxa, mya, ieofx, lzero, x(1:mxa,nu1+ir1:nu1+irn,1), b, yret(:,ir1:irn),      &
                         nt=iru, nxe=nxo, xvp=xvp(1,ir1:irn), nc=iru,tclim1=y(1:mya,1:iru,1), tclim0=y(1:mya,1:iru,0), &
                         fcastt=yrett(:,ir1:irn), yprst=y(1:mya,nu1+ir1:nu1+irn,0))
                  Case (ia_gcm)
                    Call gcm_predict (nr,mya,lzero,xiny(1:mya,nu1+ir1:nu1+irn,ig0),b0(ig0,:),b(ig0,:),yret(:,ir1:irn), &
                         nt=iru, xvp=xvp(:,ir1:irn), nc=iru, tclim1=y(1:mya,1:iru,1), tclim0=y(1:mya,1:iru,0),         &
                         fcastt=yrett(:,ir1:irn), yprst=y(1:mya,nu1+ir1:nu1+irn,0))
                 End Select
              End If
           End Select
           If (lstop) GoTo 1
           If (iev == 3) Then
              Select Case (ianal)
               Case (ia_cca)
                 Call cca_predict (iru, mxa, mya, nxo, nyo, nco, ieofx, lzero, x(1:mxa,1:iru,1), ieofy, yhat(:,:))
               Case (ia_pcr, ia_mlr)
                 Call pcr_predict (iru, mxa, mya, ieofx, lzero, x(1:mxa,1:iru,1), b, yhat(:,:))
               Case (ia_gcm)
                 Call gcm_predict (iru, mya, lzero, xiny(1:mya,1:iru,ig0), b0(ig0,:), b(ig0,:), yhat(:,:))
              End Select
              If (lstop) GoTo 1
           End If
! - calculate probabilities -
           If (iensc == 1) Then
              Select Case (igauss)
               Case (0)
                 Call calc_reg_probs (ianal, igauss, iru, nr, ng, nxo, mya, lb, ave(:,lb), sdev(:,lb), y(1:mya,:,lb), yhat, &
                      xvp(:,ir1:irn), yret(:,ir1:irn), tobst, pobs, clf, rfps(:,ir1:irn,:), rodds(:,ir1:irn,:,:), yrpls(:,ir1:irn))
               Case Default
! add persistence component to yopt?
                 Call calc_reg_probs (ianal, igauss, iru, nr, ng, nxo, mya, lb,                   &
                      ave(:,lb), sdev(:,lb), yopt(:,:,2), yhat, xvp(:,ir1:irn), yrett(:,ir1:irn), &
                      tobst, pobs, clf, rfps(:,ir1:irn,:),rodds(:,ir1:irn,:,:), yrpls(:,ir1:irn), &
                      nc=iru, clim=y(1:mya,1:iru,lb))
              End Select
           Else
!              Call calc_ens_probs (iensc,iru,nr,ng,mya,xfile%nem,x(1:mxa,:,1),yhat,pobs,clf,rfps(:,ir1:irn,:), &
!                   rodds(:,ir1:irn,:,:),yrpls(:,ir1:irn))
           End If
           If (lstop) GoTo 1
           Print *, ' '
           If (gifile%cfile%lset) Call print_text (' ', &
                                 iout=iout)
        End If
! - update indices -
        iru = nu1 + irn
        ir1 = irn + 1
     End Do retro
  End If
!
! Optimise using full training period
  If (ianal == ia_mlr) nxo = nxe
  prdc = yfile%prd1 + (yfile%it1 + nt - 2)
  cprd2 = get_cdate(prdc, 1, 2)
  ncv = nu - lcw
  ir = ir + 1
  If (nopt > 0) Then
     Print *, 'Optimizing cross-validated performance ...'
     Print *, Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2)
     Print *, ' '
     If (gifile%cfile%lset) Then
        Call print_text (Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2), &
             iout=iout)
        Call print_text (' ', &
             iout=iout)
     End If
     Select Case (igauss)
      Case (0)
        Select Case (ianal)
         Case (ia_cca)
           Call cv_cca (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
                yhat=yhat, nxo=nxo, nyo=nyo, nco=nco)
         Case (ia_pcr, ia_mlr)
           Call cv_pcr (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxe, mxe, ir, lopt, &
                yhat=yhat, nxo=nxo)
        End Select
      Case Default
        Select Case (ianal)
         Case (ia_cca)
           Call cv_cca (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
                yhat=yhat, yhatt=yhatt, nxo=nxo, nyo=nyo, nco=nco)
         Case (ia_pcr, ia_mlr)
           Call cv_pcr (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxe, mxe, ir, lopt, &
                yhat=yhat, yhatt=yhatt, nxo=nxo)
        End Select
     End Select
!
! Cross-validate optimal / chosen model
  Else
     Print *, 'Cross-validating model ...'
     Print *, Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2)
     Print *, ' '
     If (gifile%cfile%lset) Then
        Call print_text (Trim(cg_tperiod)//': '//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2), &
             iout=iout)
        Call print_text (' ', &
             iout=iout)
     End If
     Select Case (igauss)
      Case (0)
        Select Case (ianal)
         Case (ia_cca)
           Call cv_cca (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, nxo, ieofy, nyo, nyo, nco, nco, ir, lopt, &
                yhat=yhat)
         Case (ia_pcr, ia_mlr)
           Call cv_pcr (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, nxo, ir, lopt, &
                yhat=yhat)
         Case (ia_gcm)
           Call cv_gcm (imc, xfile%nms, nu, ncv, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:), &
                yhat=yhat)
        End Select
      Case Default
        Select Case (ianal)
         Case (ia_cca)
           Call cv_cca (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, nxo, ieofy, nyo, nyo, nco, nco, ir, lopt, &
                yhat=yhat, yhatt=yhatt)
         Case (ia_pcr, ia_mlr)
           Call cv_pcr (nu, ncv, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, nxo, ir, lopt, &
                yhat=yhat, yhatt=yhatt)
         Case (ia_gcm)
           Call cv_gcm (imc, xfile%nms, nu, ncv, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:), &
                yhat=yhat, yhatt=yhatt)
        End Select
     End Select
  End If
  If (lstop) GoTo 1
!
! Close goodness file
  If (gifile%cfile%lset) Then
     Close (Unit=iout)
     Call file_reset (gifile%cfile%lset, gifile%cfile%fname)
  End If
!
! Allocate additional memory
  Select Case (ianal)
   Case (ia_cca)
     nccu = Min(nxo, nyo, mcc)
     Call init_analysis2 (ifail, &
          nco=nccu)
   Case (ia_pcr, ia_mlr, ia_gcm)
     Call init_analysis2 (ifail)
  End Select
  If (ifail /= 0) Then
     perform = 1
     GoTo 1
  End If
!
! Fit model using all data
  Print *, 'Constructing model using full training period ('//Trim(cprd1)//' '//Trim(cg_to_l)//' '//Trim(cprd2)//') ...'
  If (ianal == ia_gcm .and. xfile%nms > 1) Call get_combined (mya, xfile%nms, b0(1:,:), b(1:,:), xiny(:,:,1:), xiny(:,:,0))
  Select Case (igauss)
   Case (0)
     Select Case (ianal)
      Case (ia_cca)
        Call full_cca (nu, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, &
             ieofy, nyo, svy, eofy, tsy, nccu, mu, r, s, ifail, &
             hx_map=hx_map, hy_map=hy_map, hx_ser=hx_ser, hy_ser=hy_ser)
      Case (ia_pcr, ia_mlr)
        Call full_pcr (nu, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, bz, b, npx, ifail)
      Case (ia_gcm)
        Call full_gcm (imc, xfile%nms, nu, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:))
     End Select
   Case Default
     Select Case (ianal)
      Case (ia_cca)
        Call full_cca (nu, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, &
             ieofy, nyo, svy, eofy, tsy, nccu, mu, r, s, ifail, &
             yt=yt, hx_map=hx_map, hy_map=hy_map, hx_ser=hx_ser, hy_ser=hy_ser)
      Case (ia_pcr, ia_mlr)
        Call full_pcr (nu, mxa, x(1:mxa,:,1), mya, y(1:mya,:,1), ieofx, nxo, svx, eofx, tsx, bz, b, npx, ifail, &
             yt=yt)
      Case (ia_gcm)
        Call full_gcm (imc, xfile%nms, nu, mya, ig0, xiny(:,:,ig0:), y(1:mya,:,1), b0(ig0:,:), b(ig0:,:), &
             yt=yt)
     End Select
  End Select
  If (lstop) GoTo 1
  If (ifail /= 0) Then
     perform = 2
     GoTo 1
  End If
!
! Calculate regression constants
  Select Case (ianal)
   Case (ia_pcr, ia_mlr)
     b0(1:mya,1) = ym(1:mya) - MatMul(xm(1:mxa), b(1:mxa,1:mya))
  End Select
  perform = 0
!
! Free memory
1 If (ianal == ia_gcm) Then
     If (imc == 4 .and. xfile%nms > 1) Call close_combine ()
  End If
!
  Return
 End Function perform_analysis
!
!
!
 Function get_nopt(nxe, mxe, &
          nye, mye, ncc, mcc) &
          Result (nopt)
!
! Calculates number of optimization steps
!
! Function type
  Integer :: nopt
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nxe ! - minimum number of X EOF modes -
  Integer, Intent(In) :: mxe ! - maximum number of X EOF modes -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: nye ! - minimum number of Y EOF modes -
  Integer, Intent(In), Optional :: mye ! - maximum number of Y EOF modes -
  Integer, Intent(In), Optional :: ncc ! - minimum number of CCA modes -
  Integer, Intent(In), Optional :: mcc ! - maximum number of CCA modes -
!
! Locals
!
! Local scalars
  Integer :: ic ! - CCA mode index -
  Integer :: ix ! - X EOF index -
  Integer :: iy ! - Y EOF index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Min, Present
!
! Executable Statements
!
! Calculate number of retroactive cross-validations
  If (Present(nye) .and. Present(mye) .and. Present(ncc) .and. Present(mcc)) Then
     nopt = 0
     Do ix = nxe, mxe
        Do iy = nye, mye
           Do ic = ncc, Min(mcc, ix, iy)
              nopt = nopt + 1
           End Do
        End Do
     End Do
  Else
     nopt = mxe + 1 - nxe
  End If
!
  Return
 End Function get_nopt
!
!
!
 Function num_calcs(ianal, iretro, igauss, izero, ifit, mya, nt, nret, nretro, nt1, ntr, nopt)
!
! Calculates total number of calculations to perform for progress meter
!
! Modules
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm, ia_pfv
!
! Function type
  Integer :: num_calcs
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal  ! - analysis method -
  Integer, Intent(In) :: iretro ! - retroactive forecast flag -
  Integer, Intent(In) :: igauss ! - transform to gaussian flag -
  Integer, Intent(In) :: izero  ! - zero-bound flag -
  Integer, Intent(In) :: ifit   ! - calculate fitted values flag -
  Integer, Intent(In) :: mya    ! - number of available gridpoints -
  Integer, Intent(In) :: nt     ! - number of cases in training period -
  Integer, Intent(In) :: nret   ! - number of retroactive iterations -
  Integer, Intent(In) :: nretro ! - model update interval -
  Integer, Intent(In) :: nt1    ! - initial number of training cases -
  Integer, Intent(In) :: ntr    ! - number of cases in retroactive period -
  Integer, Intent(In) :: nopt   ! - number of optimization steps -
!
! Locals
!
! Local scalars
  Integer :: ncvr ! - total number of retroactive cross-validations -
  Integer :: itg  ! - transformation flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Max, Nint, Real, Sqrt
!
! Executable Statements
!
! Calculate number of retroactive cross-validations
  ncvr = nret*nt1 + nret*(nret - 1)*nretro/2
!
! Calculate number of calculations to perform
  Select Case (igauss)
   Case (0)
     itg = 0
   Case Default
     itg = 1
  End Select
  num_calcs = 1
  Select Case (ianal)
   Case (ia_cca) ! - CCA -
     If (iretro == 1) Then
        num_calcs = num_calcs + ncvr*(1 + 2*nopt) + nret*(1 + (1 - itg)*izero) ! - cv_analysis -
        num_calcs = num_calcs + 6*nret                                         ! - full_analysis -
        num_calcs = num_calcs + 2*ntr                                          ! - analysis_prediction -
        num_calcs = num_calcs + ifit*2*ncvr                                    ! - fit_analysis -
        num_calcs = num_calcs + nret*(Max(mya, ntr) + itg + izero)             ! - calc_probs -
     End If
     num_calcs = num_calcs + nt*(1 + 2*nopt) + 1 + itg + izero                 ! - cv_analysis -
     num_calcs = num_calcs + 8                                                 ! - full_analysis -
   Case (ia_pcr, ia_mlr) ! - PCR, MLR -
     If (iretro == 1) Then
        num_calcs = num_calcs + ncvr*(1 + nopt) + nret*(1 + (1 - itg)*izero)   ! - cv_analysis -
        num_calcs = num_calcs + 2*nret                                         ! - full_analysis -
        num_calcs = num_calcs + ntr                                            ! - analysis_prediction -
        num_calcs = num_calcs + ifit*ncvr                                      ! - fit_analysis -
        num_calcs = num_calcs + nret*(Max(mya, ntr) + itg + izero)             ! - calc_probs -
     End If
     num_calcs = num_calcs + nt*(1 + nopt) + 1 + itg + izero                   ! - cv_analysis -
     num_calcs = num_calcs + 5                                                 ! - full_analysis -
   Case (ia_gcm) ! - GCM -
     num_calcs = num_calcs + mya                                               ! - get_nearest_grids -
     num_calcs = num_calcs + mya                                               ! - get_interpolated -
     If (iretro == 1) Then
        num_calcs = num_calcs + ncvr + nret*(1 + (1 - itg)*izero)              ! - cv_analysis -
        num_calcs = num_calcs + 2*nret                                         ! - full_analysis -
        num_calcs = num_calcs + ntr                                            ! - analysis_prediction -
        num_calcs = num_calcs + ifit*ncvr                                      ! - fit_analysis -
        num_calcs = num_calcs + nret*(Max(mya, ntr) + itg + izero)             ! - calc_probs -
     End If
     num_calcs = num_calcs + nt + 1 + itg + izero                              ! - cv_analysis -
     num_calcs = num_calcs + 2                                                 ! - full_analysis -
   Case (ia_pfv) ! - PFV -
     num_calcs = num_calcs + mya                                               ! - get_nearest_grids -
     num_calcs = num_calcs + mya                                               ! - get_interpolated -
  End Select
  If (ianal /= ia_pfv) num_calcs = num_calcs + Nint(Sqrt(Real(mya*nt)))        ! - set_cv_categories -
  If (iretro == 1) num_calcs = num_calcs + Nint(Sqrt(Real(mya*ntr)))           ! - set_ra_categories -
!
  Return
 End Function num_calcs
!
!
!
 Subroutine proj_read_analysis (iin, ianalp, pver, lres, ifail)
!
! Reads program analysis from project file
!
! Modules
  Use data_numbers, Only: sp, tolv
  Use scores,       Only: gm
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin ! - input file unit number -
!
  Real(Kind=sp), Intent(In) :: pver ! - project version -
!
  Logical, Intent(In) :: lres ! - read results? -
!
! Input/output scalars
  Integer, Intent(InOut) :: ianalp ! - project analysis -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Read analysis settings
  If (.not.lres) Then
     Read (Unit=iin, Err=1, End=2) ianalp, icalc, lconvert
!
! Read optimization settings
  Else
     If (pver > 15.05_sp-tolv .and. nopt > 1) Then
        Read (Unit=iin, Err=1, End=2) lopt
        Read (Unit=iin, Err=1, End=2) gm
     Else
        If (Allocated(lopt)) Deallocate (lopt)
        If (Allocated(gm))   Deallocate (gm)
     End If
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
 End Subroutine proj_read_analysis
!
!
!
 Subroutine proj_write_analysis (iout, lres, ifail)
!
! Writes analysis settings to project file
!
! Modules
  Use scores, Only: gm
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output file unit number -
!
  Logical, Intent(In) :: lres ! - read results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write analysis settings
  If (.not.lres) Then
     Write (Unit=iout, Err=1) ianal, icalc, lconvert
!
! Write optimization settings
  Else If (nopt > 0) Then
     Write (Unit=iout, Err=1) lopt
     Write (Unit=iout, Err=1) gm
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
 End Subroutine proj_write_analysis
!
!
!
 Function close_analysis()
!
! Frees memory allocated for analysis
!
! Modules
  Use arrays,     Only: awk, b, b0, rwk, iwk, x, xc, xiny, xm, xsd, xvp, y, yc, yhat, yhatt, ym, yopt, ypc, yret, yrett, yrpls, &
                        ysd, yt
  Use space,      Only: close_space
  Use fields,     Only: imfs
  Use scores,     Only: gm, &  
                        close_scores
  Use pcs,        Only: close_svds
  Use missing,    Only: close_missing
  Use categories, Only: efps, rfps, rodds, &
                        close_categories
  Use climate,    Only: close_climate
!
! Function type
  Integer :: close_analysis
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory allocated for scoring
  close_analysis = close_scores()
!
! Free other memory
  If (Allocated(rwk))   Deallocate (rwk)
  If (Allocated(awk))   Deallocate (awk)
  If (Allocated(iwk))   Deallocate (iwk)
  If (Allocated(yt))    Deallocate (yt)
  If (Allocated(yopt))  Deallocate (yopt)
  If (Allocated(xvp))   Deallocate (xvp)
  If (Allocated(yc))    Deallocate (yc)
  If (Allocated(xc))    Deallocate (xc)
  If (Allocated(ypc))   Deallocate (ypc)
  If (Allocated(b))     Deallocate (b)
  If (Allocated(b0))    Deallocate (b0)
  close_analysis = close_space()
  If (Allocated(b))     Deallocate (b)
  If (Allocated(b0))    Deallocate (b0)
  If (Allocated(xiny))  Deallocate (xiny)
  close_analysis = close_space()
  If (Allocated(rodds)) Deallocate (rodds)
  If (Allocated(rfps))  Deallocate (rfps)
  If (Allocated(imfs))  Deallocate (imfs)
  If (Allocated(yrpls)) Deallocate (yrpls)
  If (Allocated(yret))  Deallocate (yret)
  If (Allocated(yrett)) Deallocate (yrett)
  If (Allocated(efps))  Deallocate (efps)
  If (Allocated(gm))    Deallocate (gm)
  If (Allocated(lopt))  Deallocate (lopt)
  If (Allocated(yhatt)) Deallocate (yhatt)
  If (Allocated(yhat))  Deallocate (yhat)
  If (Allocated(ysd))   Deallocate (ysd)
  If (Allocated(xsd))   Deallocate (xsd)
  If (Allocated(ym))    Deallocate (ym)
  If (Allocated(xm))    Deallocate (xm)
!
! Free memory allocated for climatologies
  close_analysis = close_climate()
!
! Free memory allocated for categories
  close_analysis = close_categories()
!
! Free memory allocated for SVDs
  close_analysis = close_svds()
!
! Free memory allocated for missing values
  close_analysis = close_missing()
!
! Free other memory
  If (Allocated(y)) Deallocate (y)
  If (Allocated(x)) Deallocate (x)
!
  close_analysis = 0
!
  Return
 End Function close_analysis
End Module analysis
