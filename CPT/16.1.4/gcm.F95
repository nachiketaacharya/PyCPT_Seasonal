! Author: Simon Mason
Module gcm
!
! Modules
  Use data_numbers, Only: rp
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: cv_gcm, full_gcm, gcm_predict
!
Contains
!
!
 Subroutine full_gcm (imc, nms, nt, ny, ig0, x, y, b0, b1, &
            yt)
!
! Performs GCM using all data
!
! Modules
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: xc, yc
  Use settings,       Only: igauss, istd
  Use distribs,       Only: gaussian, set_gamma
  Use climate,        Only: ave, clim, nc, sdev
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imc ! - model combination identifier -
  Integer, Intent(In) :: nms ! - number of models -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: ny  ! - number of y spatial points -
  Integer, Intent(In) :: ig0 ! - lower index limit -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: y ! - response variables -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:,ig0:), Intent(InOut) :: x ! - explanatory variables -
!
! Output arrays
  Real(Kind=rp), Dimension(ig0:,:), Intent(Out) :: b0 ! - regression constants -
  Real(Kind=rp), Dimension(ig0:,:), Intent(Out) :: b1 ! - regression coefficients -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: yt ! - transformed response variables -
!
! Executable Statements
!
! Backup data
  xc(1:ny,1:nt,ig0:) = x(1:ny,1:nt,ig0:)
  yc(1:ny,1:nt) = y(1:ny,1:nt)
  If (igauss > 0) Then
     If (istd == 3) Call set_gamma (.false.)
     Call gaussian (igauss, 1, ny, nt, yc(1:ny,1:nt), &
          nr=nc, ref=clim(1:ny,1:nc,1), ave=ave(:,1), sdev=sdev(:,1))
  End If
!
! Copy transformed data
  If (Present(yt)) yt(:,1:nt) = yc(:,1:nt)
  If (update_progress_meter(.false.) /= 0) Return
!
! Calculate combination/calibration parameters
  Call do_gcm (imc, nms, nt, ny, ig0, xc(:,:,ig0:), yc, b0(ig0:,:), b1(ig0:,:))
  If (update_progress_meter(.false.) /= 0) Return
!
  Return
 End Subroutine full_gcm
!
!
!
 Subroutine cv_gcm (imc, nms, nt, ncv, nv, ig0, x, y, b0, b1, &
            yhat,yhatt)
!
! Performs cross-validated GCM
!
! Modules
  Use data_io_constants, Only: iout
  Use labels,            Only: l_goodness
  Use progress_meter,    Only: update_progress_meter
  Use arrays,            Only: xc, yc, yopt, &
                               get_cvsample
  Use iofiles,           Only: gifile
  Use settings,          Only: hcw, igauss, igood, istd, lzero
  Use iofiles,           Only: gifile
  Use distribs,          Only: gaussian, set_gamma, transform_cv
  Use scores,            Only: goodness
  Use climate,           Only: ave, sdev
  Use ensemble,          Only: get_combined
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imc ! - model combination identifier -
  Integer, Intent(In) :: nms ! - number of models -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: ncv ! - length of cross-validated training period -
  Integer, Intent(In) :: nv  ! - number of y spatial points -
  Integer, Intent(In) :: ig0 ! - lower index limit -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: y ! - response variables -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:,ig0:), Intent(InOut) :: x ! - explanatory variables -
!
! Output arrays
  Real(Kind=rp), Dimension(ig0:,:), Intent(Out) :: b0 ! - regression constants -
  Real(Kind=rp), Dimension(ig0:,:), Intent(Out) :: b1 ! - regression coefficients -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: yhat  ! - cross-validated hindcasts -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: yhatt ! - transformed cross-validated hindcasts -
!
! Locals
!
! Local scalars
  Integer :: it ! - cross-validation time-step -
!
  Real(Kind=rp) :: gm ! - goodness metric -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Trim
!
! Executable Statements
!
! Cross-validate
  time_step: Do it = 1, nt
     Call get_cvsample (hcw, it, y(:,1:nt), yc(:,1:nt), &
          v3=x(:,1:nt,ig0:), v3cv=xc(:,1:nt,ig0:))
!
! Transform
     If (igauss > 0) Then
        If (istd == 3) Call set_gamma (.false.)
        Call gaussian (igauss, 1, nv, ncv, yc(1:nv,1:ncv), &
             nr=ncv, ave=ave(:,1), sdev=sdev(:,1))
     End If
     If (update_progress_meter(.false.) /= 0) Return
!
! Calculate combination/calibration parameters
     Call do_gcm (imc, nms, ncv, nv, ig0, xc(:,:,ig0:), yc, b0(ig0:,:), b1(ig0:,:))
     If (nms>1) Call get_combined (nv, nms, b0(1:,:), b1(1:,:), x(:,it:it,1:), x(:,it:it,0))
     Call gcm_predict (1, nv, lzero, x(:,it:it,ig0), b0(ig0,:), b1(ig0,:), yopt(:,it:it,1))
     If (update_progress_meter(.false.) /= 0) Return
  End Do time_step
!
! Calculate and print goodness metric
  gm = goodness(igood, nt, nv, yopt(:,:,1), y(:,:))
  Write (Unit=*, Fmt='(1X,A,F18.3)') Trim(l_goodness%c)//': ', gm
  If (gifile%cfile%lset) Write (Unit=iout, Fmt='(1X,A,F18.3)') Trim(l_goodness%c)//': ', gm
  If (Present(yhatt)) yhatt(:,1:nt) = yopt(:,1:nt,1)
  If (Present(yhat)) Then
     yhat(:,1:nt) = yopt(:,1:nt,1)
     If (istd /= 3) Call transform_cv (1, nv, nt, ncv, ave(:,1), sdev(:,1), yhat) ! - transform -
  End If
!
  Return
 End Subroutine cv_gcm
!
!
!
 Subroutine do_gcm (imc, nms, nt, ny, ig0, x, yc, b0, b1)
!
! Performs GCM
!
! Modules
  Use data_numbers, Only: zero, eps, one
  Use arrays,       Only: xm, xsd, ym, ysd, &
                          moments
  Use settings,     Only: igcms, igood, lxt
  Use scores,       Only: get_regr, goodness
  Use pcs,          Only: bz, eofx, ieofx, npx, svx, tsx, xwk
  Use ensemble,     Only: get_combined
  Use pcr,          Only: do_pcr
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imc ! - model combination identifier -
  Integer, Intent(In) :: nms ! - number of models -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: ny  ! - number of y spatial points -
  Integer, Intent(In) :: ig0 ! - lower index limit -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: yc ! - response variables -
!
  Real(Kind=rp), Dimension(:,:,ig0:), Intent(InOut) :: x ! - explanatory variables -
!
! Output arrays
  Real(Kind=rp), Dimension(ig0:,:), Intent(Out) :: b0 ! - regression constants -
  Real(Kind=rp), Dimension(ig0:,:), Intent(Out) :: b1 ! - regression coefficients -
!
! Locals
!
! Local scalars
  Integer :: i     ! - location index -
  Integer :: ibest ! - best model -
  Integer :: imdl  ! - current model -
  Integer :: ifail ! - error indicator -
!
  Real(Kind=rp) :: gm    ! - goodness metric -
  Real(Kind=rp) :: gbest ! - best score -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MatMul, Real, Sum, Transpose
!
! Executable Statements
!
! Select model combination approach
  Call moments (ny, nt, yc, ym, ysd)
  If (nms > 1) Then
     Select Case (imc)
      Case (1) ! - uncalibrated average -
        b1(1:,:) = one/Real(nms, Kind=rp)
        b0(1:,:) = zero
      Case (2) ! - calibrated average -
        Select Case (igcms)
         Case (0, 1) ! - correct mean biases -
           b1(1:,:) = one/Real(nms, Kind=rp)
           Do imdl = 1, nms
              Call moments (ny, nt, x(:,:,imdl), xm, xsd)
              b0(imdl,:) = ym(:) - xm(:)
           End Do
           b0(1:,:) = zero
         Case Default ! - correct mean and variance biases -
           Do imdl = 1, nms
              Call moments (ny, nt, x(:,:,imdl), xm, xsd)
              b0(imdl,:) = ym(:) - xm(:)
              Do i = 1, ny
                 If (xsd(i) > eps) Then
                    b1(imdl,i) = ysd(i)/xsd(i)
                 Else
                    b1(imdl,i) = zero
                 End If
                 b0(imdl,i) = ym(i) - xm(i)*b1(imdl,i)
              End Do
           End Do
        End Select
      Case (3) ! - recalibrated average -
        Do imdl = 1, nms
           Do i = 1, ny
              Call get_regr (nt, x(i,:,imdl), yc(i,:), b0(imdl,i), b1(imdl,i))
           End Do
        End Do
      Case (4) ! - skill-weighted average -
        Do i = 1, ny
           xwk(:,1:nt) = Transpose(x(i,:,1:))
           Call do_pcr (nt, nms, xwk(:,1:nt), 1, yc(:,1:nt), ieofx, nms, svx, eofx, tsx, bz(:,i:i), .true., lxt, npx, ifail)
           b1(1:nms,i:i) = MatMul(eofx(1:nms,1:npx), bz(1:npx,1:1))
           b0(1:nms,i) = (ym(i) - Sum(xm(1:nms)*b1(1:nms,i)))/Real(nms, Kind=rp)
        End Do
      Case (5) ! - best model by location -
        Do i = 1, ny
           gbest = -one
           ibest = 0
           Do imdl = 1, nms
              gm=goodness(igood, nt, 1, x(imdl:imdl,:,1), yc(imdl:imdl,:))
              If (gm > gbest) Then
                 gbest = gm
                 ibest = imdl
              End If
           End Do
           b1(1:,i) = zero
           b1(ibest,i) = one
           b0(1:,i) = zero
        End Do
      Case (6) ! - best overall model -
        gbest = -one
        ibest = 0
        Do imdl = 1, nms
           gm = goodness(igood, nt, ny,x(:,:,1), yc(:,:))
           If (gm > gbest) Then
              gbest = gm
              ibest = imdl
           End If
        End Do
        b1(1:,:) = zero
        b1(ibest,:)=one
        b0(1:,:)=zero
     End Select
!
! Combine models
     Call get_combined (ny, nms, b0(1:,:), b1(1:,:), x(:,:,1:), x(:,:,0))
  End If
!
! Recalibrate combined forecasts if necessary
  Call moments (ny, nt, x(:,:,ig0), xm, xsd)
  Select Case (igcms)
   Case (0) ! - no correction -
     b0(ig0,:) = zero
     b1(ig0,:) = one
   Case (1) ! - correct mean biases -
     b1(ig0,:) = one
     b0(ig0,:) = ym(:) - xm(:)
   Case (2) ! - correct mean and variance biases -
     Do i = 1, ny
        If (xsd(i) > eps) Then
           b1(ig0,i) = ysd(i)/xsd(i)
        Else
           b1(ig0,i) = zero
        End If
        b0(ig0,i) = ym(i) - xm(i)*b1(ig0,i)
     End Do
   Case (3) ! - correct for skill -
     Do i = 1, ny
        Call get_regr (nt, x(i,:,ig0), yc(i,:), b0(ig0,i), b1(ig0,i))
     End Do
  End Select
!
  Return
 End Subroutine do_gcm
!
!
!
 Subroutine gcm_predict (nf, ny, lzero, x, b0, b1, fcast, &
            nt, xvp, nc, tclim1, tclim0, fcastt, yprst)
!
! Calculates predictions given new predictor values
!
! Modules
  Use data_numbers,   Only: zero, one
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: xm, xsd
  Use settings,       Only: igauss, istd, istdo
  Use distribs,       Only: gaussian, gaussian_inv, stdize, ustdize, zero_bound
  Use climate,        Only: ave, sdev
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nf ! - number of forecasts -
  Integer, Intent(In) :: ny ! - number of Y variables -
!
  Logical, Intent(In) :: lzero ! - apply zero-bound? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: nt ! - number of cases in training period -
  Integer, Intent(In), Optional :: nc ! - number of cases in climatology -
!
! Arrays,
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: b0 ! - regression constants -
  Real(Kind=rp), Dimension(:), Intent(In) :: b1 ! - regression coefficients -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: x ! - new predictor values -
!
! - optional input arrays -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: tclim1 ! - climatological data for transformation (not for anomalies) -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: tclim0 ! - climatological data for transformation (with persistence) -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: yprst  ! - persistence component -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: fcast ! - forecast values -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: xvp    ! - predictors over forecast period -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: fcastt ! - transformed forecast values -
!
! Locals
!
! Local scalars
  Integer :: i ! - variable index -
  Integer :: k ! - case index -
!
  Real(Kind=rp) :: df ! - number of cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present, Real
!
! Executable Statements
!
! Predict anomalies
  If (Present(xvp)) df=Real(nt, Kind=rp)
  Do k = 1, nf
     fcast(:,k) = b0(:) + b1(:)*x(:,k)
!
! Calculate prediction error variance
     If (Present(xvp)) Then
        Do i = 1, ny
           If (xsd(i) > zero) Then
              xvp(i,k) = (one + ((x(i,k) - xm(i))**2)/xsd(i))/df
           Else
              xvp(i,k) = zero
           End If
        End Do
     End If
!
! Update progress meter
     If (update_progress_meter(.false.) /= 0) Return
  End Do
!
! Transform
! - transform from gaussian -
  If (Present(tclim1) .and. Present(nc)) Then
     If (Present(fcastt) .and. .not.Present(yprst)) fcastt(:,:) = fcast(:,:)
     Call gaussian_inv (igauss,istd,1,ny,nf,fcast(:,:), &
          nr=nc, ref=tclim1(:,1:nc), lgam=.false., ave=ave(:,1), sdev=sdev(:,1))
     If (.not.Present(yprst) .and. istd == 3) Call stdize (istd, igauss, 1, ny, nf, ave(:,1), sdev(:,1), fcast(:,:))
  End If
! - add persistence component -
  If (Present(yprst)) Then
     If (Present(fcastt) .and. (istd /= 0 .and. istd /= 3)) Call ustdize (istd, igauss, 0, ny, nf, ave(:,1), sdev(:,1), fcast(:,:))
     fcast(:,:) = fcast(:,:) + yprst(:,:)
     If (istd /= 0) Call stdize (istd, igauss, 0, ny, nf, ave(:,0), sdev(:,0), fcast(:,:))
     If (Present(fcastt)) Then
        fcastt(:,:) = fcast(:,:)
        If (istd /= 3) Then
           Call gaussian (igauss, 0, ny, nf, fcastt(:,:), &
                nr=nc, ref=tclim0(:,1:nc), ave=ave(:,0), sdev=sdev(:,0))
        End If
     End If
! - apply zero-transform -
     If (lzero) Call zero_bound (istd, ny, nf, ave(:,0), sdev(:,0), fcast)
  Else
     If (lzero) Call zero_bound (istdo, ny, nf, ave(:,1), sdev(:,1), fcast(:,:))
  End If
!
  Return
 End Subroutine gcm_predict
End Module gcm
