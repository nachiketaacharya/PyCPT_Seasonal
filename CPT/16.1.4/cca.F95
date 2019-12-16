! Author: Simon Mason
Module cca
!
! Modules
  Use data_numbers, Only: rp
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: cca_predict, cv_cca, full_cca
!
Contains
!
!
 Subroutine full_cca (nt, nx, x, ny, y, ieofx, nxe, svx, eofx, tsx, ieofy, nye, svy, eofy, tsy, ncc, mu, r, s, ifail, &
            yt, hx_map, hy_map, hx_ser, hy_ser)
!
! Performs CCA using all data
!
! Modules
  Use data_numbers,   Only: zero, tol
#if DP == 1
  Use cpt_lapack,     Only: gesdd=>dgesdd
#else
  Use cpt_lapack,     Only: gesdd=>sgesdd
#endif
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: iusex, iusey, iwk, lrwk, rwk, xc, xm, xsd, yc, ym, ysd, &
                            moments
  Use iofiles,        Only: xfile, yfile
  Use settings,       Only: igauss, istd, lxt, lxyt, lyt
  Use space,          Only: rlatx, rlaty, &
                            latitude_weight
  Use fields,         Only: xfield, yfield
  Use distribs,       Only: gaussian, set_gamma,stdize
  Use pcs,            Only: ce, lc1, lr1, mxe, npx, npy
  Use climate,        Only: ave, clim, nc, sdev
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt    ! - number of cases -
  Integer, Intent(In) :: nx    ! - number of x spatial points -
  Integer, Intent(In) :: ny    ! - number of y spatial points -
  Integer, Intent(In) :: nxe   ! - number of X EOF modes -
  Integer, Intent(In) :: nye   ! - number of Y EOF modes -
  Integer, Intent(In) :: ncc   ! - number of CCA modes -
  Integer, Intent(In) :: ieofx ! - X EOF option -
  Integer, Intent(In) :: ieofy ! - Y EOF option -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: x ! - explanatory variables -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: y ! - response variables -
!
! Output arrays
  Real(Kind=rp), Dimension(:), Intent(Out) :: mu  ! - canonical correlations -
  Real(Kind=rp), Dimension(:), Intent(Out) :: svx ! - singular values of x -
  Real(Kind=rp), Dimension(:), Intent(Out) :: svy ! - singular values of y -
!
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: eofx ! - x EOF patterns -
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: eofy ! - y EOF patterns -
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: tsx  ! - time-series of x EOFs (transposed) -
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: tsy  ! - time-series of y EOFs (transposed) -
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: r    ! - canonical Y EOF weights -
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: s    ! - canonical X EOF weights (transposed) -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: yt     ! - transformed response variables -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: hx_map ! - X homogeneous covariance maps -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: hy_map ! - Y homogeneous covariance maps -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: hx_ser ! - X homogeneous covariance maps time series -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: hy_ser ! - Y homogeneous covariance maps time series -
!
! Locals
!
! Local scalars
  Integer :: k   ! - case index -
  Integer :: ie  ! - EOF mode index -
  Integer :: mnc ! - maximum number of CCA modes to be counted -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Abs, Count, MatMul, MaxVal, Min, MinVal, Present, Transpose
!
! Executable Statements
!
! Backup data
  ifail = 1
  xc(1:nx,1:nt,1) = x(1:nx,1:nt)
  yc(1:ny,1:nt) = y(1:ny,1:nt)
  If (igauss > 0) Then
     If (istd == 3) Call set_gamma (.false.)
     Call gaussian (igauss, 1, ny, nt, yc(1:ny,1:nt), &
          nr=nc, ref=clim(1:ny,1:nc,1), ave=ave(:,1), sdev=sdev(:,1))
  End If
!
! Copy transformed data
  If (Present(yt)) yt(:,1:nt) = yc(:,1:nt)
!
! Calculate means and standard deviations
  Call moments (nx, nt, xc(:,:,1), xm, xsd)
  Call moments (ny, nt, yc, ym, ysd)
! - update progress meter -
  If (update_progress_meter(.false.) /= 0) Return
!
! Calculate anomalies according to PCA type
  Call stdize (3-ieofx, 0, 1, nx, nt, xm, xsd, xc(:,1:nt,1))
  Call stdize (3-ieofy, 0, 1, ny, nt, ym, ysd, yc(:,1:nt))
!
! Scale by latitude
  If (xfile%idstr == 1) Call latitude_weight (xfile%nfl, xfield(1:)%nlt, xfield(1:)%region, rlatx(:,:), iusex, nt, xc(:,:,1))
  If (yfile%idstr == 1) Call latitude_weight (yfile%nfl, yfield(1:)%nlt, yfield(1:)%region, rlaty(:,:), iusey, nt, yc(:,:))
! - update progress meter -
  If (update_progress_meter(.false.) /= 0) Return
!
! Perform X EOF prefiltering
  Call gesdd ('S', nx, nt, xc(:,:,1), nx, svx, eofx, nx, tsx, lxt, rwk, lrwk, iwk, ifail)
  If (ifail /= 0) Return
! - determine number of non-zero eigenvalues -
  npx = Count(svx(1:Min(lxt, nt-1)) > zero)
  If (npx > 2) Then
     If (svx(npx)*svx(npx - 2)/svx(npx - 1)**2 < tol) npx = npx - 1 ! - check for probable rounding errors -
  End If
! - rescale loadings by latitude -
  If (xfile%idstr == 1) Call latitude_weight (xfile%nfl, xfield(:)%nlt, xfield(:)%region, rlatx(:,:), iusex, npx, eofx(:,:))
! - ensure that largest absolute loading is positive -
  Do ie = 1, npx
     If (Abs(MaxVal(eofx(1:nx,ie))) < Abs(MinVal(eofx(1:nx,ie)))) Then
        eofx(1:nx,ie) = -eofx(1:nx,ie)
        tsx(ie,1:nt) = -tsx(ie,1:nt)
     End If
  End Do
  If (npx < lxt) Then
     svx(npx+1:lxt) = zero
     tsx(npx+1:lxt,1:nt) = zero
     eofx(1:nx,npx+1:lxt) = zero
  End If
! - update progress meter -
  If (update_progress_meter(.false.) /= 0) Return
!
! Perform Y EOF prefiltering
  Call gesdd ('S', ny, nt, yc, ny, svy, eofy, ny, tsy, lyt, rwk, lrwk, iwk, ifail)
  If (ifail /= 0) Return
! - determine number of non-zero eigenvalues -
  npy = Count(svy(1:Min(lyt, nt-1)) > zero)
  If (npy > 2) Then
     If (svy(npy)*svy(npy - 2)/svy(npy - 1)**2 < tol) npy = npy - 1 ! - check for probable rounding errors -
  End If
! - rescale loadings by latitude -
  If (yfile%idstr == 1) Call latitude_weight (yfile%nfl, yfield(1:)%nlt, yfield(1:)%region, rlaty(:,:), iusey, npy, eofy(:,:))
! - ensure that largest absolute loading is positive -
  Do ie = 1, npy
     If (Abs(MaxVal(eofy(1:ny,ie))) < Abs(MinVal(eofy(1:ny,ie)))) Then
        eofy(1:ny,ie) = -eofy(1:ny,ie)
        tsy(ie,1:nt) = -tsy(ie,1:nt)
     End If
  End Do
  If (npy < lyt) Then
     svy(npy+1:lyt) = zero
     tsy(npy+1:lyt,1:nt) = zero
     eofy(1:ny,npy+1:lyt) = zero
  End If
! - update progress meter -
  If (update_progress_meter(.false.) /= 0) Return
!
! Compute CCA
  ce(1:nye,1:nxe) = MatMul(tsy(1:nye,1:nt), Transpose(tsx(1:nxe,1:nt)))
  Call gesdd ('S', nye, nxe, ce, lc1, mu, r, lr1, s, lxyt, rwk, lrwk, iwk, ifail)
  If (ifail /= 0) Return
  mnc = Min(nxe, nye)
  If (lxyt > mnc) Then
     mu(mnc+1:) = zero
     r(:,mnc+1:) = zero
     If (lr1 > nye) r(nye+1:,1:mnc) = zero
     s(mnc+1:,:) = zero
     If (mxe > nxe) s(1:mnc,nxe:) = zero
  End If
! - update progress meter -
  If (update_progress_meter(.false.) /= 0) Return
!
! Compute homogeneous covariance maps
! - X maps -
  If (Present(hx_map)) Then
     Call hcov_maps (nx, nt, nxe, ncc, svx, eofx, tsx, s, hx_map, hx_ser, rwk(1:nxe*ncc))
! - update progress meter -
     If (update_progress_meter(.false.) /= 0) Return
  End If
! - Y maps -
  If (Present(hy_map)) Then
     Call hcov_maps (ny, nt, nye, ncc, svy, eofy, tsy, Transpose(r(1:nye,1:ncc)), hy_map, hy_ser, rwk(1:nye*ncc))
! - update progress meter -
     If (update_progress_meter(.false.) /= 0) Return
  End If
!
! Rescale EOF time series
! - X EOFs -
  Do k = 1, nt
     Do ie = 1, npx
        tsx(ie,k) = tsx(ie,k)*svx(ie)
     End Do
! - Y EOFs -
     Do ie = 1, npy
        tsy(ie,k) = tsy(ie,k)*svy(ie)
     End Do
  End Do
! - update progress meter -
  If (update_progress_meter(.false.) /= 0) Return
  ifail = 0
!
  Return
 End Subroutine full_cca
!
!
!
 Subroutine hcov_maps (nv, nt, ne, nc, sv, eof, ts, rs, hmap, hser, rwk)
!
! Constructs CCA maps
!
! Arguments
!
! Inpur scalars
  Integer, Intent(In) :: nv ! - number of variables -
  Integer, Intent(In) :: nt ! - number of cases -
  Integer, Intent(In) :: ne ! - number of EOF modes -
  Integer, Intent(In) :: nc ! - number of CCA modes -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: sv ! - singular values -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: eof ! - spatial loadings -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: ts  ! - temporal scores -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: rs  ! - canonical EOF weights -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: hmap ! - homogeneous covariance maps -
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: hser ! - homogeneous covariance map time series -
!
! Workspace
!
! Workspace arrays
  Real(Kind=rp) :: rwk(:) ! - real workspace -
!
! Locals
!
! Local scalars
  Integer :: i1, i2 ! - indices -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MatMul,Reshape
!
! Executable Statements
!
! Calculate homogeneous covariance maps
  ForAll (i2 = 1:nc)
     ForAll (i1 = 1:ne)
        rwk((i2-1)*ne+i1) = sv(i1)*rs(i2,i1)
     End ForAll
  End ForAll
  hmap(1:nv,1:nc) = MatMul(eof(1:nv,1:ne), Reshape(rwk(1:ne*nc), (/ne,nc/)))
!
! Calculate corresponding time series
  hser(1:nc,1:nt) = MatMul(rs(1:nc,1:ne), ts(1:ne,1:nt))
!
  Return
 End Subroutine hcov_maps
!
!
!
 Subroutine cv_cca (nt, ncv, nx, x, ny, y, ieofx, nxe, mxe, ieofy, nye, mye, ncc, mcc, ir, lopt, &
            yhat, yhatt, nxo, nyo, nco)
!
! Performs cross-validated CCA
!
! Modules
  Use data_numbers,      Only: zero, one, tol
  Use data_io_constants, Only: iout
#if DP == 1
  Use cpt_lapack,        Only: gesdd=>dgesdd
#else
  Use cpt_lapack,        Only: gesdd=>sgesdd
#endif
  Use labels,            Only: l_goodness
  Use gui,               Only: print_text
  Use progress_meter,    Only: lstop, &
                               update_progress_meter
  Use arrays,            Only: iusex, iusey, iwk, lrwk, rwk, xc, xm, xsd, yc, ym, yopt, ysd, &
                               get_cvsample,moments
  Use iofiles,           Only: gifile, xfile, yfile
  Use settings,          Only: igauss, igood, istd, hcw, lxt, lxyt, lyt, lzero, sw_elev
  Use space,             Only: rlatx, rlaty, &
                               latitude_weight
  Use fields,            Only: xfield, yfield
  Use distribs,          Only: gaussian, set_gamma, stdize, transform_cv
  Use scores,            Only: gm, &
                               goodness
  Use pcs,               Only: ce, eofx, eofy, icco, lc1, lr1, npx, npy, mu, r, s, svx, svy, tsx, tsy
  Use climate,           Only: ave, sdev
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt    ! - number of cases -
  Integer, Intent(In) :: ir    ! - retroactive step -
  Integer, Intent(In) :: ncv   ! - length of cross-validated training period -
  Integer, Intent(In) :: nx    ! - number of x spatial points -
  Integer, Intent(In) :: ny    ! - number of y spatial points -
  Integer, Intent(In) :: nxe   ! - minimum number of X EOF modes -
  Integer, Intent(In) :: mxe   ! - maximum number of X EOF modes -
  Integer, Intent(In) :: nye   ! - minimum number of Y EOF modes -
  Integer, Intent(In) :: mye   ! - maximum number of Y EOF modes -
  Integer, Intent(In) :: ncc   ! - minimum number of CCA modes -
  Integer, Intent(In) :: mcc   ! - maximum number of CCA modes -
  Integer, Intent(In) :: ieofx ! - X EOF option -
  Integer, Intent(In) :: ieofy ! - Y EOF option -
!
! Output scalars
! - optional output scalars -
  Integer, Intent(Out), Optional :: nxo ! - optimal number of X EOF modes -
  Integer, Intent(Out), Optional :: nyo ! - optimal number of Y EOF modes -
  Integer, Intent(Out), Optional :: nco ! - optimal number of CCA modes -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: x ! - explanatory variables -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: y ! - response variables -
!
! Output arrays
  Logical, Dimension(:,:), Intent(Out) :: lopt ! - optimised? -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: yhat  ! - cross-validated hindcasts -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: yhatt ! - transformed cross-validated hindcasts -
!
! Locals
!
! Local scalars
  Integer :: it    ! - cross-validation time-step -
  Integer :: ixe   ! - X EOF mode index -
  Integer :: iye   ! - Y EOF mode index -
  Integer :: icc   ! - current number of CCA modes -
  Integer :: ixu   ! - used X EOF mode index -
  Integer :: iyu   ! - used Y EOF mode index -
  Integer :: icu   ! - current number of used CCA modes -
  Integer :: iopt  ! - optimization step index -
  Integer :: jopt  ! - index of optimized settings -
  Integer :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic AdjustR, Any, Count, MatMul, Min, Present, Transpose, Trim
!
! Executable Statements
!
! Set optimization checks
  lopt(:,ir) = .true.
!
! Cross-validate
  time_step: Do it = 1, nt
     Call get_cvsample (hcw, it, y(:,1:nt), yc(:,1:nt), &
          v2=x(:,1:nt), v2cv=xc(:,1:nt,1))
!
! Transform
     If (igauss > 0) Then
        If (istd == 3) Call set_gamma (.false.)
        Call gaussian (igauss, 1, ny, ncv, yc(1:ny,1:ncv), &
             nr=ncv, ave=ave(:,1), sdev=sdev(:,1))
     End If
!
! Calculate means and standard deviations
     Call moments (nx, ncv, xc(:,1:ncv,1), xm, xsd)
     Call moments (ny, ncv, yc(:,1:ncv), ym, ysd)
!
! Calculate anomalies according to PC type
     Call stdize (3-ieofx, 0, 1, nx, ncv, xm, xsd, xc(:,1:ncv,1))
     Call stdize (3-ieofy, 0, 1, ny, ncv, ym, ysd, yc(:,1:ncv))
!
! Scale by latitude
     If (xfile%idstr == 1) Call latitude_weight (xfile%nfl, xfield(1:)%nlt, xfield(1:)%region, rlatx(:,:), iusex, ncv, xc(:,:,1))
     If (yfile%idstr == 1) Call latitude_weight (yfile%nfl, yfield(1:)%nlt, yfield(1:)%region, rlaty(:,:), iusey, ncv, yc(:,:))
!
! Perform EOF prefiltering
! - X variables -
     Call gesdd ('S', nx, ncv, xc(:,:,1), nx, svx, eofx, nx, tsx, lxt, rwk, lrwk, iwk, ifail)
! - Y variables -
     If (ifail == 0) Call gesdd ('S', ny, ncv, yc, ny, svy, eofy, ny, tsy, lyt, rwk, lrwk, iwk, ifail)
     If (update_progress_meter(.false.) /= 0) Return
!
! Determine number of non-zero eigenvalues
     iopt = 0
1    If (ifail == 0) Then
        npx = Count(svx(1:Min(mxe,ncv-1)) > zero)
        npy = Count(svy(1:Min(mye,ncv-1)) > zero)
! - check for near-zero eigenvalues -
        If (npx == 0) Then
           ifail = 1
           GoTo 1
        Else If (npx > 2) Then
           If (svx(npx)*svx(npx - 2)/svx(npx - 1)**2 < tol) npx = npx - 1
        End If
        If (npy == 0) Then
           ifail = 1
           GoTo 1
        Else If (npy > 2) Then
           If (svy(npy)*svy(npy - 2)/svy(npy - 1)**2 < tol) npy = npy - 1
        End If
! - rescale time scores -
        If (npx < mxe) Then
           svx(npx+1:mxe) = zero
           tsx(npx+1:mxe,1:ncv) = zero
           eofx(1:nx,npx+1:mxe) = zero
        End If
        If (npy < mye) Then
           svy(npy+1:mye) = zero
           tsy(npy+1:mye,1:ncv) = zero
           eofy(1:ny,npy+1:mye) = zero
        End If
!
! Rescale by latitude
        If (xfile%idstr == 1) Call latitude_weight (xfile%nfl, xfield(1:)%nlt, xfield(1:)%region, rlatx(:,:), iusex, npx, eofx(:,:))
        If (yfile%idstr == 1) Call latitude_weight (yfile%nfl, yfield(1:)%nlt, yfield(1:)%region, rlaty(:,:), iusey, npy, eofy(:,:))
!
! Compute CCA
        Do ixe = nxe, mxe
           ixu = Min(ixe, npx)
           Do iye = nye, mye
              iyu = Min(iye, npy)
              ce(1:iyu,1:ixu) = MatMul(tsy(1:iyu,1:ncv), Transpose(tsx(1:ixu,1:ncv)))
              Call gesdd ('S', iyu, ixu, ce, lc1, mu, r, lr1, s, lxyt, rwk, lrwk, iwk, ifail)
!
! Check for over-fitting
              If (icco == 0 .and. Any(mu(1:Min(ixu, iyu)) > one - tol)) Then
                 Do icc = ncc, Min(mcc, ixe, iye)
                    iopt = iopt + 1
                    lopt(iopt,ir) = .false.
                    yopt(1:ny,it,iopt) = ym(1:ny)
                    If (update_progress_meter(.false.) /= 0) Return
                    If (update_progress_meter(.false.) /= 0) Return
                 End Do
                 Cycle
              Else
!
! Predict anomaly
                 Do icc = ncc,Min(mcc, ixe, iye)
                    icu = Min(icc, ixu, iyu)
                    iopt = iopt + 1
                    If (icco == 1 .or. lopt(iopt,ir)) Then
                       Call cca_predict (1, nx, ny, ixu, iyu, icu, ieofx, lzero, x(1:nx,it:it), ieofy, yopt(1:ny,it:it,iopt))
                    Else
                       yopt(1:ny,it,iopt) = ym(1:ny)
                       If (update_progress_meter(.false.) /= 0) Return
                       If (update_progress_meter(.false.) /= 0) Return
                    End If
                 End Do
              End If
           End Do
        End Do
!
! Supply mean if CCA failed
     Else
        Do ixe = nxe, mxe
           Do iye = nye, mye
              Do icc = ncc, Min(mcc, ixe, iye)
                 iopt = iopt + 1
                 yopt(1:ny,it,iopt) = ym(1:ny)
                 If (update_progress_meter(.false.) /= 0) Return
                 If (update_progress_meter(.false.) /= 0) Return
              End Do
           End Do
        End Do
     End If
     If (lstop) Return
  End Do time_step
!
! Standardise
  If (igauss > 0) Then
     yc(1:ny,1:nt) = y(1:ny,1:nt)
     If (istd /= 3) Then
        Call gaussian (igauss, 1, ny, nt, yc(:,:), & ! - NB strictly, this step should be cross-validated (as below) -
             nr=nt, ave=ave(:,1), sdev=sdev(:,1))
     End If
  End If
!
! Calculate and print goodness metric
  If ((Present(nxo) .or. Present(nyo) .or. Present(nco)) .and. Any(lopt(:,ir))) Then
     If (sw_elev%on) Then
        Print *, '              CURRENT                         OPTIMUM'
        Print *, ' '
        Print *, '     Number of Modes    Goodness      Number of Modes    Goodness'
        Print *, '        X    Y   CCA       Index         X    Y   CCA       Index'
        Print *, ' '
     End If
     If (gifile%cfile%lset) Then
        Call print_text ('               CURRENT                         OPTIMUM', &
             iout=iout)
        Call print_text (' ', &
             iout=iout)
        Call print_text ('      Number of Modes    Goodness      Number of Modes    Goodness', &
             iout=iout)
        Call print_text ('         X    Y   CCA       Index         X    Y   CCA       Index', &
             iout=iout)
        Call print_text (' ', &
             iout=iout)
     End If
     iopt = 0
     jopt = 1
     Do ixe = nxe, mxe
        Do iye = nye, mye
           Do icc = ncc, Min(mcc, ixe, iye)
              iopt = iopt + 1
              If (lopt(iopt,ir)) Then
                 Select Case (igauss)
                  Case (0)
                    gm(iopt,ir) = goodness(igood, nt, ny, yopt(:,:,iopt), y)
                  Case Default
                    gm(iopt,ir) = goodness(igood, nt, ny, yopt(:,:,iopt), yc)
                 End Select
                 If (iopt == 1 .or. gm(iopt,ir) > gm(jopt,ir)) Then
                    jopt = iopt
                    nxo = ixe
                    nyo = iye
                    nco = icc
                 End If
                 If (sw_elev%on) Write (Unit=*, Fmt='(2(5X,2I5,I6,F12.3))') &
                                   ixe, iye, icc, gm(iopt,ir), nxo, nyo, nco, gm(jopt,ir)
                 If (gifile%cfile%lset) Write (Unit=iout, Fmt='(2(5X,2I5,I6,F12.3))') &
                                              ixe, iye, icc, gm(iopt,ir), nxo, nyo, nco, gm(jopt,ir)
              Else
                 If (sw_elev%on) Write (Unit=*, Fmt='(5X,2I5,I6,A12,5X,2I5,I6,F12.3)') &
                                   ixe, iye, icc, AdjustR('N/A'), nxo, nyo, nco, gm(jopt,ir)
                 If (gifile%cfile%lset) Write (Unit=iout, Fmt='(5X,2I5,I6,A12,5X,2I5,I6,F12.3)') &
                                              ixe, iye, icc, AdjustR('N/A'), nxo, nyo, nco, gm(jopt,ir)
              End If
           End Do
        End Do
        If (lstop) Return
     End Do
     If (.not.sw_elev%on) Then
        Write (Unit=*, Fmt='(1X,A,F18.3,A,I0,A,I0,A,I0,A)') &
              Trim(l_goodness%c)//': ', gm(jopt,ir), ' (using ', nxo, ' X,', nyo, ' Y, and ', nco, ' CCA modes)'
     End If
  Else
     jopt = 1
     gm(jopt,ir) = goodness(igood, nt, ny, yopt(:,:,jopt), y)
     Write (Unit=*, Fmt='(1X,A,F18.3)') Trim(l_goodness%c)//': ',gm(jopt,ir)
     If (gifile%cfile%lset) Write (Unit=iout, Fmt='(1X,A,F18.3)') Trim(l_goodness%c)//': ', gm(jopt,ir)
! - set optimized values if there were no successful calculations -
     If (Present(nxo)) nxo = 1
     If (Present(nyo)) nyo = 1
     If (Present(nco)) nco = 1
  End If
  If (lstop) Return
  If (update_progress_meter(.false.) /= 0) Return
  If (Present(yhatt)) yhatt(:,1:nt) = yopt(:,1:nt,jopt)
  If (Present(yhat)) Then
     yhat(:,1:nt) = yopt(:,1:nt,jopt)
     If (istd /= 3) Call transform_cv (1, ny, nt, ncv, ave(:,1), sdev(:,1), yhat) ! - transform -
  End If
!
  Return
 End Subroutine cv_cca
!
!
!
 Subroutine cca_predict (nf, nx, ny, nxe, nye, ncc, ieofx, lzero, x, ieofy, fcast, &
            nt, xvp, fsx, fsxc, nc, tclim1, tclim0, fcastt, yprst)
!
! Calculates predictions given new predictor values
!
! Modules
  Use data_numbers,   Only: one
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: rwk, xc, xm, xsd, ym, ysd
  Use settings,       Only: igauss, istd, istdo
  Use distribs,       Only: gaussian, gaussian_inv, stdize, ustdize, zero_bound
  Use pcs,            Only: eofx, eofy, mxe, mu, prjc, r, s, svx, svy
  Use climate,        Only: ave, sdev
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nf    ! - number of forecasts -
  Integer, Intent(In) :: nx    ! - number of X variables -
  Integer, Intent(In) :: ny    ! - number of Y variables -
  Integer, Intent(In) :: nxe   ! - number of X EOF modes -
  Integer, Intent(In) :: nye   ! - number of Y EOF modes -
  Integer, Intent(In) :: ncc   ! - number of CCA modes -
  Integer, Intent(In) :: ieofx ! - X EOF option -
  Integer, Intent(In) :: ieofy ! - Y EOF option -
!
  Logical, Intent(In) :: lzero ! - apply zero-bound? -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: nt  ! - number of cases in training period -
  Integer, Intent(In), Optional :: nc  ! - number of cases in transformation climatology -
!
! Arrays,
!
! Input arrays
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
  Real(Kind=rp), Dimension(:), Intent(Out), Optional :: xvp ! - predictors over forecast period -
!
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: fsx    ! - predictor time scores -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: fsxc   ! - predictor CCA time scores -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: fcastt ! - transformed forecast values -
!
! Locals
!
! Local scalars
  Integer :: i ! - EOF mode index -
  Integer :: k ! - case index -
!
  Real(Kind=rp) :: df ! - number of cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic MatMul, Present, Real, Sum, Transpose
!
! Executable Statements
!
! Standardize predictors
  xc(:,1:nf,1) = x(:,:)
  Call stdize (3-ieofx, 0, 1, nx, nf, xm, xsd, xc(:,1:nf,1))
!
! Predict anomalies
  If (Present(xvp)) df = Real(nt, Kind=rp)
  Do k = 1, nf
!
! Project X anomalies onto the X EOFs and X CCA modes
     rwk(1:nxe) = MatMul(Transpose(eofx(:,1:nxe)), xc(:,k,1))
     If (Present(fsx)) Then
        fsx(1:nxe,k) = rwk(1:nxe)
        Do i = nxe+1, mxe
           fsx(i,k) = Sum(eofx(:,i)*xc(:,k,1))
        End Do
     End If
     If (Present(fsxc)) fsxc(1:ncc,k) = MatMul(s(1:ncc,1:nxe), fsx(1:nxe,k)/svx(1:nxe))
! - standardise to unit variance by scaling by inverse of singular values -
     rwk(1:nxe) = rwk(1:nxe)/svx(1:nxe)
!
! Project X EOFs onto the CCA modes
     prjc(1:ncc) = MatMul(s(1:ncc,1:nxe), rwk(1:nxe))
     If (Present(xvp)) xvp(k) = one/df + Sum(prjc(1:ncc))**2
! - scale by canonical correlations -
     prjc(1:ncc) = prjc(1:ncc)*mu(1:ncc)
     If (update_progress_meter(.false.) /= 0) Return
!
! Project CCA modes back onto the Y EOF modes
     rwk(1:nye) = MatMul(r(1:nye,1:ncc), prjc(1:ncc))
! - scale by singular values -
     rwk(1:nye) = rwk(1:nye)*svy(1:nye)
!
! Project Y EOF modes back onto the Y anomalies
     fcast(:,k) = MatMul(eofy(:,1:nye), rwk(1:nye))
!
! Rescale
     Select Case (ieofy)
      Case (1)
        fcast(:,k) = fcast(:,k)*ysd(:) + ym(:)
      Case (2)
        fcast(:,k) = fcast(:,k) + ym(:)
     End Select
     If (update_progress_meter(.false.) /= 0) Return
  End Do
!
! Transform
! - transform from gaussian -
  If (Present(tclim1) .and. Present(nc)) Then
     If (Present(fcastt) .and. .not.Present(yprst)) fcastt(:,:) = fcast(:,:)
     Call gaussian_inv (igauss, istd, 1, ny, nf, fcast(:,:), &
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
 End Subroutine cca_predict
End Module cca
