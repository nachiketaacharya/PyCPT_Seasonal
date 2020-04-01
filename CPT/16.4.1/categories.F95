! Author: Simon Mason
Module categories
!
! Modules
  Use data_numbers,       Only: rp
  Use data_cpt_constants, Only: ng, nts
  Use time,               Only: pdate
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: calc_odds, calc_reg_probs, calc_thresholds, change_cat, check_refdates, close_categories, init_refdates,     &
            init_categories, init_icat, init_thresholds, proj_read_categories, proj_write_categories, set_cv_categories, &
            set_percentiles, set_ra_categories, set_thresholds
!
! Scalars
!
! Integer scalars
  Integer, Public :: ithr    ! - threshold type -
  Integer, Public :: ithr_bk ! - backup threshold type -
  Integer, Public :: icatc   ! - current threshold -
  Integer, Public :: icatp   ! - previous threshold -
!
! Arrays
!
! Integer arrays
  Integer, Dimension(nts), Public :: irefs ! - reference-date indices -
!
  Integer, Dimension(:,:), Allocatable, Public :: ifor  ! - forecast categories -
  Integer, Dimension(:,:), Allocatable, Public :: iobs  ! - observed categories -
  Integer, Dimension(:,:), Allocatable, Public :: irfor ! - retroactive forecast categories -
  Integer, Dimension(:,:), Allocatable, Public :: irobs ! - retroactive observed categories -
!
! Real arrays
  Real(Kind=rp), Dimension(nts), Public :: pthr ! - percentile thresholds -
  Real(Kind=rp), Dimension(nts), Public :: thr  ! - absolute thresholds -
!
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: hps   ! - highest forecast probabilities -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: pobs  ! - climatological percentile thresholds -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: tobs  ! - climatological absolute thresholds -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: tobst ! - transformed climatological absolute thresholds -
  Real(Kind=rp), Dimension(:,:), Allocatable, Public :: tobsx ! - model absolute thresholds -
!
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: efps  ! - ensemble forecast probabilities -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: fps   ! - forecast probabilities -
  Real(Kind=rp), Dimension(:,:,:), Allocatable, Public :: rfps  ! - retroactive forecast probabilities -
!
  Real(Kind=rp), Dimension(:,:,:,:), Allocatable, Public :: odds  ! - odds -
  Real(Kind=rp), Dimension(:,:,:,:), Allocatable, Public :: rodds ! - retroactive odds -
!
! Derived-Type arrays
  Type(pdate), Dimension(nts), Public :: refs ! - reference dates -
!
Contains
!
!
 Subroutine init_icat ()
!
! Initialises category indicator
!
! Executable Statements
!
! Initialize category indicator
  icatc = 1
  icatp = 1
!
  Return
 End Subroutine init_icat
!
!
!
 Subroutine init_refdates ()
!
! Initialises reference dates
!
! Modules
  Use time, Only: Assignment(=)
!
! Executable Statements
!
! Initialize reference dates
  refs(:) = 0
  irefs(:) = 0
!
  Return
 End Subroutine init_refdates
!
!
!
 Function check_refdates(ydate1, nt) &
          Result (check)
!
! Checks reference dates
!
! Modules
  Use time,     Only: pdate, &
                      Operator(==), Operator(<), Operator(>), Operator(+), &
                      iseq, &
                      date_diff
  Use missing,  Only: ktuse
!
! Function Type
  Integer :: check
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt ! - number of cases -
!
  Type(pdate), Intent(In) :: ydate1 ! - first date of training period -
!
! Locals
!
! Local scalars
  Integer :: j  ! - threshold index -
  Integer :: iy ! - year index -
!
  Type(pdate) :: ydate2 ! - last date of training period -
!
! Functions and Subroutines
!
! Intrinsic Functions
  Intrinsic Count
!
! Executable Statements
!
! Initialize reference dates
  check = 0
  ydate2 = ydate1 + (nt - 1)
  If (.not.refs(1) == 0) Then
     Do j = 1, nts
! - year is too early -
        If (refs(j) < ydate1) Then
           iy = 1
           irefs(j) = 0
           Call find_iref (1)
           refs(j) = ydate1 + (iy - 1)
           Select Case (check)
            Case (0)
              check = j
            Case Default
              check = 3
           End Select
        Else If (refs(j) > ydate2) Then ! - year is too late -
           iy = nt
           irefs(j) = Count(ktuse(1:nt)) + 1
           Call find_iref (-1)
           refs(j) = ydate1 + (iy - 1)
           Select Case (check)
            Case (0)
              check = j
            Case Default
              check = 3
           End Select
        Else ! - year is within limits -
           irefs(j) = date_diff(ydate1, refs(j), iseq) + 1
           If (ktuse(irefs(j))) Then ! - year is available -
              irefs(j) = Count(ktuse(1:irefs(j)))
           Else ! - identify nearest available year -
              Select Case (j)
               Case (1)
                 iy = 1
                 irefs(j) = 0
                 Call find_iref (1)
               Case (2)
                 iy = nt
                 irefs(j) = Count(ktuse(1:nt)) + 1
                 Call find_iref (-1)
              End Select
              refs(j) = ydate1 + (iy - 1)
              If (j == 1) Then
                 check = 4
              Else
                 Select Case (check)
                  Case (0)
                    check = 5
                  Case (4)
                    check = 6
                  Case Default
                    check = 3
                 End Select
              End If
           End If
        End If
     End Do
  Else
     iy = 0
     Do j = 1, nts
        iy = iy + 1
        Do
           If (ktuse(iy)) Exit
           iy = iy + 1
        End Do
        irefs(j) = j
        refs(j) = ydate1 + (iy - 1)
     End Do
  End If
!
  Return
!
 Contains
!
!
  Subroutine find_iref (iinc)
!
! Identifies next available reference
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iinc ! - increment -
!
! Executble Statements
!
! Identify next available reference
  Do
     Do
        If (ktuse(iy)) Exit
        iy = iy + iinc
     End Do
     irefs(j) = irefs(j) + iinc
     If (irefs(1) /= irefs(2)) Exit
  End Do
!
  Return
  End Subroutine find_iref
 End Function check_refdates
!
!
!
 Function init_categories(ianalc, iretro, m, n, nr) &
          Result (init)
!
! Initialises memory for categories
!
! Modules
  Use data_cpt_constants, Only: ia_pfv
!
! Function Type
  Integer :: init
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianalc ! - analysis identifier -
  Integer, Intent(In) :: iretro ! - retroactive verification flag -
  Integer, Intent(In) :: m      ! - total number of variables -
  Integer, Intent(In) :: n      ! - number of Cases -
  Integer, Intent(In) :: nr     ! - number of retroactive Cases -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Allocate memory
! - observed categories -
  If (ianalc /= ia_pfv) Then
     Allocate (iobs(m,n), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - hindcast categories -
     Allocate (ifor(m,n), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
! - retroactive observed categories -
  If (iretro == 1) Then
     Allocate (irobs(m,nr), Stat=ifail)
     If (ifail /= 0) GoTo 1
! - retroactive forecast categories -
     Allocate (irfor(m,nr), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
!
  init = 0
  Return
! 
! Errors
1 init = 1
  Return
 End Function init_categories
!
!
!
 Function init_thresholds(nsn) &
          Result (init)
!
! Initialises memory for thresholds
!
! Modules
  Use data_cpt_constants, Only: ng, nts
  Use iofiles,            Only: mya
  Use settings,           Only: iensc
!
! Function Type
  Integer :: init
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nsn ! - number of seasons -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Allocate memory
! - percentiles / climatological probabilities -
  Allocate (pobs(nsn*mya,ng), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - absolute thresholds -
  Allocate (tobs(nsn*mya,nts), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - transformed absolute thresholds -
  Allocate (tobst(nsn*mya,nts), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - model absolute thresholds -
  If (iensc > 1) Then
     Allocate (tobsx(nsn*mya,nts), Stat=ifail)
     If (ifail /= 0) GoTo 1
  End If
!
  init = 0
  Return
!
! Errors
1 init = 1
!
  Return
 End Function init_thresholds
!
!
!
 Subroutine calc_thresholds (igauss, mya, nsn, nc, clim, tobs, tobst, ifail)
!
! Sets climatological thresholds and percentiles
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: igauss ! - transform Y data flag -
  Integer, Intent(In) :: mya    ! - number of available variables -
  Integer, Intent(In) :: nsn    ! - number of seasons -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nc ! - length of climatological period -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: clim ! - climatological data -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: tobs  ! - climatological absolute thresholds -
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: tobst ! - transformed climatological absolute thresholds -
!
! Executable Statements
!
! Sets climatological thresholds and percentiles
  Select Case (ithr)
   Case (1, 2)
     ifail = set_percentiles(nsn, nc(:), clim(:,:,1)) ! - clim is not used for PFV -
     If (ifail /= 0) Return
     Call set_thresholds(igauss, mya, nsn, nc(:), clim(:,:,:), tobs, tobst, ifail)
     If (ifail /= 0) Return
   Case (3)
     Call set_thresholds(igauss, mya, nsn, nc(:), clim(:,:,:), tobs, tobst, ifail)
     If (ifail /= 0) Return
     ifail = set_percentiles(nsn, nc(:), clim(:,:,1))
     If (ifail /= 0) Return
  End Select
!
  Return
 End Subroutine calc_thresholds
!
!
!
 Subroutine set_thresholds (igauss, mya, nsn, nc, clim, tobs, tobst, ifail)
!
! Sets category thresholds
!
! Modules
  Use data_numbers,   Only: zero, one
  Use progress_meter, Only: set_progress_increment, update_progress_meter
  Use arrays,         Only: y, &
                            insertion_sort
  Use distribs,       Only: q_empirical, q_normal
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: igauss ! - transform Y data flag -
  Integer, Intent(In) :: mya    ! - number of available variables -
  Integer, Intent(In) :: nsn    ! - number of seasons -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nc ! - length of climatological period -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: clim ! - climatological data -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: tobs  ! - climatological absolute thresholds -
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: tobst ! - transformed climatological absolute thresholds -
!
! Locals
!
! Local scalars
  Integer :: i   ! - variable index -
  Integer :: j   ! - category index -
  Integer :: isn ! - season index -
!
  Real(Kind=rp) :: sp ! - cumulative probability -
!
! Executable Statements
!
! Set thresholds by percentile
  ifail=-1
  Select Case (ithr)
   Case (1)
     Do isn = 1, nsn
        Do i = 1, mya
           Do j = 1, nts
              tobs((isn-1)*mya+i,j) = q_empirical(clim(i,:,isn), nc(isn), pthr(j))
           End Do
           If (update_progress_meter(.false.) /= 0) Return
        End Do
     End Do
!
! Set thresholds by absolute value
   Case (2)
     Do i = 1, mya*nsn
        tobs(i,:) = thr(:)
        If (update_progress_meter(.false.) /= 0) Return
     End Do
!
! Set thresholds by reference dates
   Case (3)
     If (nsn > 0) Then
        Do isn = 1, nsn
           Do i = 1, mya
              Do j = 1, nts
                 tobs((isn-1)*mya+i,j) = y(i,irefs(j),isn)
              End Do
           End Do
        End Do
     Else
        Do i = 1, mya
           Do j = 1, nts
              tobs(i,j) = y(i,irefs(j),0)
           End Do
        End Do
     End If
     Call insertion_sort (mya*nsn, nts, tobs(:,:), 'a')
     Call set_progress_increment (Real(mya, Kind=rp), .true.)
     If (update_progress_meter(.false.) /= 0) Return
     Call set_progress_increment (one/Real(mya, Kind=rp), .true.)
  End Select
!
! Transform thresholds
  Select Case (igauss)
   Case (0)
     tobst(:,:) = tobs(:,:)
   Case Default
     Do i = 1, mya*nsn
        sp = zero
        Do j = 1, nts
           sp = sp + pobs(i,j)
           tobst(i,j) = q_normal(sp)
        End Do
     End Do
     If (update_progress_meter(.false.) /= 0) Return
  End Select
  ifail = 0
!
  Return
 End Subroutine set_thresholds
!
!
!
 Function set_percentiles(nsn, nc, clim)
!
! Sets category thresholds as percentiles
!
! Modules
  Use data_numbers,       Only: zero, one
  Use data_cpt_constants, Only: ng, nts
  Use progress_meter,     Only: update_progress_meter
  Use iofiles,            Only: mya
  Use settings,           Only: igauss
  Use distribs,           Only: cdf_empirical, q_normal
!
! Function Type
  Integer :: set_percentiles
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nsn ! - number of seasons -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nc ! - length of climatological period -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: clim ! - climatological data -
!
! Locals
!
! Local scalars
  Integer :: i   ! - variable index -
  Integer :: ii  ! - variable index -
  Integer :: j   ! - category index -
  Integer :: isn ! - season index -
!
  Real(Kind=rp) :: sp ! - cumulative probability -
!
! Executable Statements
!
! Set percentiles
  set_percentiles = -1
  Select Case (ithr)
   Case (1)
     ii = 0
     Do isn = 1, nsn
        Do i = 1, mya
           ii = ii + 1
           pobs(ii,1:nts) = pthr(:)
           pobs(ii,ng) = one
           Do j = ng, 2, -1
              pobs(ii,j) = pobs(ii,j) - pobs(ii,j-1)
           End Do
           If (update_progress_meter(.false.) /= 0) Return
        End Do
     End Do
!
! Set percentiles from absolute value
   Case (2)
     ii = 0
     Do isn = 1, nsn
        Do i = 1, mya
           ii = ii + 1
           Do j = 1, nts
              pobs(ii,j) = cdf_empirical(clim(ii,:), nc(isn), thr(j))
           End Do
           pobs(ii,ng) = one
           Do j = ng, 2, -1
              pobs(ii,j) = pobs(ii,j) - pobs(ii,j-1)
           End Do
           If (update_progress_meter(.false.) /= 0) Return
        End Do
     End Do
!
! Set percentiles from reference dates
   Case (3)
     Do i = 1, mya
        Do j = 1, nts
           pobs(i,j) = cdf_empirical(clim(i,:), nc(1), tobs(i,j))
        End Do
        pobs(i,ng) = one
        Do j = ng, 2, -1
           pobs(i,j) = pobs(i,j) - pobs(i,j-1)
        End Do
        If (update_progress_meter(.false.) /= 0) Return
     End Do
     If (igauss > 0) Then
        Do i = 1, mya
           sp = zero
           Do j = 1, nts
              sp = sp + pobs(i,j)
              tobst(i,j) = q_normal(sp)
           End Do
           pobs(i,ng) = one
           Do j = ng, 2, -1
              pobs(i,j) = pobs(i,j) - pobs(i,j-1)
           End Do
        End Do
     End If
  End Select
  set_percentiles = 0
!
  Return
 End Function set_percentiles
!
!
!
 Function set_cv_categories(lb)
!
! Calculate cross-validated categories
!
! Modules
  Use data_numbers,   Only: one
  Use progress_meter, Only: set_progress_increment, update_progress_meter
  Use arrays,         Only: rwk, y, yhat, &
                            insertion_sort
  Use iofiles,        Only: mya
  Use settings,       Only: hcw, ncv, nu
!
! Function Type
  Integer :: set_cv_categories
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: lb ! - lower bound -
!
! Locals
!
! Local scalars
  Integer :: i  ! - Case index -
  Integer :: j  ! - threshold index -
  Integer :: k  ! - time index -
  Integer :: i1 ! - start of cross-validation sample -
!
! Functions and Subroutines
!
! Intrinsic Functions
  Intrinsic Cshift, Real
!
! Executable Statements
!
! Set categories by percentile
  set_cv_categories = -1
  Select Case (ithr)
   Case (1)
     Call set_progress_increment (one/Real(nu, Kind=rp), .true.)
! - cross-validate -
     Do k = 1, nu
        i1 = hcw + k - 1
        Do i = 1, mya
           rwk(1:nu) = Cshift(y(i,1:nu,lb), i1)
           Call insertion_sort (ncv, rwk, 'a')
! - categorise observations -
           Call categorise_per (1, nts, y(i,k:k,lb), pthr(:), ncv, rwk(:), iobs(i,k:k))
! - categorise forecasts -
           Call categorise_per (1, nts, yhat(i,k:k), pthr(:), ncv, rwk(:), ifor(i,k:k))
! - update progress meter -
           If (update_progress_meter(.false.) /= 0) Return
        End Do
     End Do
     Call set_progress_increment (Real(nu, Kind=rp), .true.)
!
! Set categories by absolute value
   Case (2)
     Do i = 1, mya
        Call categorise_abs (nu, nts, y(i,:,lb), thr(:), iobs(i,:))  ! - categorise observations -
        Call categorise_abs (nu, nts, yhat(i,:), thr(:), ifor(i,:)) ! - categorise forecasts -
        If (update_progress_meter(.false.) /= 0) Return
     End Do
!
! Set categories by reference dates
   Case (3)
     Do i = 1, mya
        Do j = 1, nts
           rwk(j) = y(i,irefs(j),lb)
        End Do
        Call insertion_sort (nts, rwk(:), 'a')
        Call categorise_abs (nu, nts, y(i,:,lb), rwk(:), iobs(i,:)) ! - categorise observations -
        Call categorise_abs (nu, nts, yhat(i,:), rwk(:), ifor(i,:)) ! - categorise forecasts -
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Select
  set_cv_categories = 0
!
  Return
 End Function set_cv_categories
!
!
!
 Function set_ra_categories(ianal, lb, nsn, is0, nc, clim)
!
! Calculate retroactive categories
!
! Modules
  Use data_numbers,       Only: zero, one
  Use data_cpt_constants, Only: ia_cca, ia_pcr, ia_mlr, ia_gcm
  Use progress_meter,     Only: set_progress_increment, update_progress_meter
  Use arrays,             Only: rwk, y, yret, &
                                insertion_sort
  Use iofiles,            Only: mya
  Use settings,           Only: nret, nretro, nt, ntr, nt1, nu, nur, nu1, nyr
  Use missing,            Only: ktuse
!
! Function Type
  Integer :: set_ra_categories
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis identifier -
  Integer, Intent(In) :: lb    ! - lower bound -
  Integer, Intent(In) :: nsn   ! - number of seasons -
  Integer, Intent(In) :: is0   ! - zeroth season -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nc ! - length of climatological period -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: clim ! - climatological data -
!
! Locals
!
! Local scalars
  Integer :: i   ! - grid / station index -
  Integer :: j   ! - category / threshold index -
  Integer :: k   ! - Case index -
  Integer :: irl ! - index of last retroactive training-period date -
  Integer :: ir1 ! - index of first retroactive date -
  Integer :: irn ! - index of last retroactive date -
  Integer :: ir0 ! - index of first new date in retroactive training period -
  Integer :: iru ! - index of last new date in retroactive training period -
  Integer :: isn ! - number of seasons -
  Integer :: isi ! - last season index -
  Integer :: is1 ! - season index -
  Integer :: nr  ! - total number of retroaoctive steps in current cycle -
!
  Real(Kind=rp) :: fmax ! - maximum forecast probability -
!
! Functions and Subroutines
!
! Intrinsic Functions
  Intrinsic Count, Min, Real, Sum
!
! Executable Statements
!
! Set categories by percentile for all but PFV
  set_ra_categories = -1
  Select Case (ithr)
   Case (1)
     Select Case (ianal)
      Case (ia_cca, ia_pcr, ia_mlr, ia_gcm)
        Call set_progress_increment (one/Real(nret, Kind=rp), .true.)
        Do i = 1, mya
           ir1 = 1
           irn = 0
           ir0 = 1
           iru = nu1
           retro: Do irl = nt1, nt-1, nretro
              irn = irn + Count(ktuse(irl+1:Min(irl+nretro, nt)))
              nr = irn + 1 - ir1
              If (nr > 0) Then
                 rwk(ir0:iru) = y(i,ir0:iru,lb)
                 Call insertion_sort (iru, rwk, 'a')
! - categorise observations -
                 Call categorise_per (nr, nts, y(i,nu1+ir1:nu1+irn,lb), pthr(:), iru, rwk, irobs(i,ir1:irn))
! - categorise forecasts -
                 Call categorise_per (nr, nts, yret(i,ir1:irn), pthr(:), iru, rwk, irfor(i,ir1:irn))
              End If
! - update indices -
              ir1 = irn + 1
              ir0 = iru + 1
              iru = nu1 + irn
! - update progress meter -
              If (update_progress_meter(.false.) /= 0) Return
           End Do retro
        End Do
        Call set_progress_increment (Real(nret, Kind=rp), .true.)
!
! Set categories for PFV
     Case (5)
        Call set_progress_increment (one/Real(nur, Kind=rp), .true.)
        Do i = 1, mya
! - categorise observations -
           Do is1 = 1, nsn
              If (is1 > ntr) Exit
              isn = is0 + is1
              If (isn > nsn) isn = isn - nsn
              isi = nsn*(nyr - 1) + is1
              If (isi <= ntr) Then
                 Call categorise_per (nyr, nts, y(i,1:nyr,isn), pthr(:), nc(isn), clim(i,1:nc(isn),isn), irobs(i,is1:isi:nsn))
              Else
                 isi = isi - nsn
                 Call categorise_per (nyr-1, nts, y(i,1:nyr-1,isn), pthr(:), nc(isn), clim(i,1:nc(isn),isn), irobs(i,is1:isi:nsn))
              End If
           End Do
! - identify category with highest probability -
           Do k = 1, nur
              If (Sum(rfps(i,k,:)) > zero) Then
                 irfor(i,k) = 1
                 fmax = rfps(i,k,1)
                 Do j = 2, ng
                    If (rfps(i,k,j) > fmax) Then
                       irfor(i,k) = j
                       fmax = rfps(i,k,j)
                    End If
                 End Do
              Else
                 irfor(i,k) = 0
              End If
              If (update_progress_meter(.false.) /= 0) Return
           End Do
        End Do
        Call set_progress_increment (Real(nur, Kind=rp), .true.)
     End Select
!
! Set categories by absolute value
   Case (2)
     irl = nu1 + 1
     Do i = 1, mya
        Call categorise_abs (nur, nts, y(i,irl:nu,lb), thr, irobs(i,:)) ! - categorise observations -
        Call categorise_abs (nur, nts, yret(i,:), thr, irfor(i,:))      ! - categorise forecasts -
        If (update_progress_meter(.false.) /= 0) Return
     End Do
!
! Set categories by reference dates
   Case (3)
     irl = nu1 + 1
     Do i = 1, mya
        Do j = 1, nts
           rwk(j) = y(i,irefs(j),1)
        End Do
        Call insertion_sort (nts, rwk(:), 'a')
        Call categorise_abs (nur, nts, y(i,irl:nu,lb), rwk(:), irobs(i,:)) ! - categorise observations -
        Call categorise_abs (nur, nts, yret(i,:), rwk(:), irfor(i,:))      ! - categorise forecasts -
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Select
  set_ra_categories = 0
!
  Return
 End Function set_ra_categories
!
!
!
 Subroutine categorise_abs (n, nts, x, t, icat)
!
! Categorises data given threshold values
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n   ! - number of Cases -
  Integer, Intent(In) :: nts ! - number of thresholds -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: x ! - data to be categorised -
  Real(Kind=rp), Dimension(:), Intent(In) :: t ! - absolute thresholds -
!
! Output arrays
  Integer, Dimension(:), Intent(Out) :: icat ! - categories -
!
! Locals
!
! Local scalars
  Integer :: k ! - Case index -
  Integer :: j ! - category index -
!
! Executable Statements
!
! Determine category
  Do k = 1, n
     icat(k)=1
     Do j = 1, nts
        If (t(j) > x(k)) Exit
        icat(k) = icat(k) + 1
     End Do
  End Do
!
  Return
 End Subroutine categorise_abs
!
!
!
 Subroutine categorise_per (n, nts, x, p, nc, c, icat)
!
! Categorises data given percentile thresholds
!
! Modules
  Use distribs, Only: q_empirical
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n   ! - number of Cases -
  Integer, Intent(In) :: nts ! - number of thresholds -
  Integer, Intent(In) :: nc  ! - length of reference climatology -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: x ! - data to be categorised -
  Real(Kind=rp), Dimension(:), Intent(In) :: p ! - percentiles -
  Real(Kind=rp), Dimension(:), Intent(In) :: c ! - reference climatology (sorted in ascending order) -
!
! Output arrays
  Integer, Dimension(:), Intent(Out) :: icat ! - categories -
!
! Locals
!
! Local scalars
  Integer :: j ! - category index -
!
! Local arrays
  Real(Kind=rp), Dimension(nts) :: t ! - thresholds -
!
! Executable Statements
!
! Set thresholds
  Do j = 1, nts
     t(j) = q_empirical(c, nc, p(j))
  End Do
!
! Determine category
  Call categorise_abs (n, nts, x, t, icat)
!
  Return
 End Subroutine categorise_per
!
!
!
 Subroutine calc_reg_probs (ianal, igauss, nt, nf, ng, nx, ny, lb, ave, sdev, y, yhat, xvp, fcast, tobst, &
            pobs, clf, fps, odds, fpls, &
            nenf, fens, nc, clim, pev)
!
! Calculates forecast probabilities and prediction intervals
!
! Modules
  Use data_numbers,       Only: zero, eps, one, two, oneh
  Use data_cpt_constants, Only: ia_gcm
  Use maths,              Only: interval
  Use progress_meter,     Only: set_progress_increment, update_progress_meter
  Use settings,           Only: dofr, istd, izero
  Use distribs,           Only: gaussian_inv, s_student_t, q_student_t_2, zero_bound
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal  ! - analysis identifier -
  Integer, Intent(In) :: igauss ! - transform to gaussian flag -
  Integer, Intent(In) :: nt     ! - number of Cases in training period -
  Integer, Intent(In) :: nf     ! - number of forecasts -
  Integer, Intent(In) :: ng     ! - number of categories -
  Integer, Intent(In) :: nx     ! - number of predictors -
  Integer, Intent(In) :: ny     ! - number of predictands -
  Integer, Intent(In) :: lb     ! - lower-bound -
!
  Real(Kind=rp), Intent(In) :: clf ! - forecast confidence level -
!
! - optional input scalars -
  Integer, Intent(In), Optional :: nc   ! - number of years in climatology -
  Integer, Intent(In), Optional :: nenf ! - number of ensemble members -
!
! Input arrays
  Real(Kind=rp), Dimension(:), Intent(In) :: ave  ! - climatological averages -
  Real(Kind=rp), Dimension(:), Intent(In) :: sdev ! - climatological standard deviations -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: y     ! - training period data -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: yhat  ! - estimated values -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: xvp   ! - predictors over forecast period -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: fcast ! - forecast values -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: tobst ! - thresholds -
  Real(Kind=rp), Dimension(:,:), Intent(In) :: pobs  ! - climatological percentile thresholds -
!
! - optional input arrays -
  Real(Kind=rp), Dimension(:,:), Intent(In), Optional :: clim ! - climatology -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:,:), Intent(Out) :: fps ! - forecast probabilities -
!
  Real(Kind=rp), Dimension(:,:,:,:), Intent(Out) :: odds ! - odds -
!
  Type(interval), Dimension(:,:), Intent(Out) :: fpls ! - prediction limits -
!
! - optional output arrays -
  Real(Kind=rp), Dimension(:,:), Intent(Out), Optional :: pev ! - prediction error standard deviations -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(Out), Optional :: fens ! - ensemble forecasts -
!
! Locals
!
! Local scalars
  Integer :: i  ! - variable index -
  Integer :: j  ! - category index -
  Integer :: k  ! - forecast index -
  Integer :: l  ! - category index -
  Integer :: ii ! - variable index -
  Integer :: ie ! - ensemble member index -
!
  Real(Kind=rp) :: t     ! - Student's t-deviate -
  Real(Kind=rp) :: t0    ! - lower tail Student's t-deviate -
  Real(Kind=rp) :: p     ! - lower tail Student's t-deviate -
  Real(Kind=rp) :: s2    ! - estimated mean squared error -
  Real(Kind=rp) :: psd   ! - prediction error standard deviation -
!
! Functions and Subroutines
!
! Intrinsic Functions
  Intrinsic Min, Present, Real, Sqrt, Sum
!
! Executable Statements
!
! Calculate critical t statistic
  dofr = Real(nt-nx-1, Kind=rp)
  If (dofr > zero) Then
     t0 = q_student_t_2(one-clf/oneh, dofr)
  Else
     t0 = zero
  End If
!
! Calculate goodness of fit statistics
  Call set_progress_increment (one/Real(Min(nf, ny), Kind=rp), .true.)
  ii = 1
  Do i = 1, ny
     If (dofr > zero) s2 = Sum((y(i,1:nt) - yhat(i,1:nt))**2)/dofr
!
! Calculate forecast probabilities
     Do k = 1, nf
        If (dofr > zero) Then
           psd = Sqrt(s2*(one + xvp(ii,k)))
        Else
           psd = zero
        End If
        If (psd>eps) Then
           Do j = 1, ng-1
              t = (tobst(i,j) - fcast(i,k))/psd
              fps(i,k,j+1) = s_student_t(t, dofr)*oneh
           End Do
           fps(i,k,1) = oneh
           Do j = 1, ng-1
              fps(i,k,j) = fps(i,k,j) - fps(i,k,j+1)
           End Do
        Else
           l = 1
           Do j = 1, ng-1
              If (fcast(i,k) > tobst(i,j)) Then
                 l = l + 1
              Else
                 Exit
              End If
           End Do
           fps(i,k,l) = oneh
        End If
!
! Calculate odds and relative odds
        Call calc_odds (1, ng, 1, pobs(i:i,:), fps(i:i,k:k,:), odds(i:i,k:k,:,:))
!
! Calculate prediction intervals
        t = t0*psd
        fpls(i,k)%lower = fcast(i,k) - t
        fpls(i,k)%upper = fcast(i,k) + t
        If (Present(pev)) pev(i,k) = psd
        If (Present(nenf)) Then
           Do ie=nenf/2+1,nenf
              p = one - Real(ie, Kind=rp)/Real(nenf+1, Kind=rp)
              If (dofr > zero) Then
                 t = q_student_t_2(two*p, dofr)*psd
              Else
                 t = zero
              End If
              fens(i,k,ie) = fcast(i,k) + t
              fens(i,k,nenf+1-ie) = fcast(i,k) - t
           End Do
        End If
!
! Update progress meter
        If (update_progress_meter(.false.) /= 0) Return
     End Do
     If (ianal == ia_gcm) ii = ii + 1
  End Do
  Call set_progress_increment (Real(Min(nf, ny), Kind=rp), .true.)
!
! Transform If necessary
  If (igauss > 0 .and. istd /= 3) Then
     Call gaussian_inv (igauss, istd, lb, ny, nf, fpls(:,:)%lower, &
          nr=nc, ref=clim(:,1:nc), lgam=.true., ave=ave(:), sdev=sdev(:))
     Call gaussian_inv (igauss, istd, lb, ny, nf, fpls(:,:)%upper, &
          nr=nc, ref=clim(:,1:nc), lgam=.true., ave=ave(:), sdev=sdev(:))
     If (Present(fens)) Then
        Do ie = 1, nenf
           Call gaussian_inv (igauss, istd, lb, ny, nf, fens(:,:,ie), &
                nr=nc, ref=clim(:,1:nc), lgam=.true., ave=ave(:), sdev=sdev(:))
        End Do
     End If
     If (update_progress_meter(.false.) /= 0) Return
  End If
!
! Apply zero-bound
  If (izero == 1) Then
     Call zero_bound (istd, ny, nf, ave, sdev, fpls(:,:)%lower)
     Call zero_bound (istd, ny, nf, ave, sdev, fpls(:,:)%upper)
     If (Present(fens)) Then
        Do ie = 1, nenf
           Call zero_bound (istd, ny, nf, ave, sdev, fens(:,:,ie))
        End Do
     End If
     If (update_progress_meter(.false.) /= 0) Return
  End If
!
  Return
 End Subroutine calc_reg_probs
!
!
!
 Subroutine calc_odds (nf, ng, nv, pobs, fps, odds)
!
! Calculates odds and relative odds
!
! Modules
  Use data_numbers, Only: zero, one, oneh, eps
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nf ! - number of forecasts -
  Integer, Intent(In) :: ng ! - number of categories -
  Integer, Intent(In) :: nv ! - number of variables -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(In) :: pobs ! - climatological probabilities -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: fps ! - forecast probabilities -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:,:,:), Intent(Out) :: odds ! - odds -
!
! Locals
!
! Local scalars
  Integer :: i ! - variable index -
  Integer :: j ! - category index -
  Integer :: k ! - forecast index -
!
! Local arrays
  Real(Kind=rp), Dimension(ng) :: codds ! - climatological odds -
!
! Executable Statements
!
! Calculate odds
  Do i = 1, nv
     Where (pobs(i,:) > zero)
        codds(:) = (one - pobs(i,:))/pobs(i,:) ! - inverted so that relative odds are calculated by multiplication rather than division -
     Else Where
        codds(:) = zero
     End Where
     Do k = 1, nf
        Do j = 1, ng
           If (fps(i,k,j) > eps .and. fps(i,k,j) < oneh - eps) Then
              odds(i,k,j,1) = fps(i,k,j)/(oneh - fps(i,k,j))
           Else
              odds(i,k,j,1) = -one
           End If
!
! Calculate odds relative to climatology
           If (pobs(i,j) > eps .and. pobs(i,j) < one - eps) Then
              odds(i,k,j,2) = odds(i,k,j,1)*codds(j)
           Else
              odds(i,k,j,2) = -one
           End If
        End Do
     End Do
  End Do
!
  Return
 End Subroutine calc_odds
!
!
!
 Function change_cat() RESULT (change)
!
! Changes category
!
! Modules
  Use settings, Only: change_setting
!
! Function Type
  Integer :: change
!
! Executable Statements
!
! Change category
  change=change_setting(icatc, icatp, 1, ng)
!
  Return
 End Function change_cat
!
!
!
 Subroutine proj_read_categories (iin, ihind, iretro, nsn, lres, ifail)
!
! Reads category settings from project file
!
! Modules
  Use time, Only: Operator(+)
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin    ! - input file unit number -
  Integer, Intent(In) :: ihind  ! - hindcasts flag -
  Integer, Intent(In) :: iretro ! - retroactive forecast flag -
  Integer, Intent(In) :: nsn    ! - number of seasons -
!
  Logical, Intent(In) :: lres ! - read results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Read category settings
  ithr_bk = ithr
  If (.not.lres) Then
     Read (Unit=iin, Err=1, End=2) refs, irefs, ithr, pthr, thr
!
! Read results
  Else
     If (init_thresholds(nsn) /= 0) GoTo 4
     Read (Unit=iin, Err=1, End=2) pobs
     Read (Unit=iin, Err=1, End=2) tobs
     Read (Unit=iin, Err=1, End=2) tobst
     If (ihind == 1) Then
        Read (Unit=iin, Err=1, End=2) ifor
        Read (Unit=iin, Err=1, End=2) iobs
     End If
     If (iretro == 1) Then                   
        Read (Unit=iin, Err=1, End=2) irobs
        Read (Unit=iin, Err=1, End=2) irfor
        Read (Unit=iin, Err=1, End=2) rfps
        Read (Unit=iin, Err=1, End=2) rodds(:,:,:,2)
        Read (Unit=iin, Err=1, End=2) rodds(:,:,:,1)
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
! - End of file -
2 ifail = 2
  Return
!
! - problem allocating memory -
4 ifail = 4
  Return
!
 End Subroutine proj_read_categories
!
!
!
 Subroutine proj_write_categories (iout, ihind, iretro, lres, ifail)
!
! Writes category settings to project file
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout   ! - output file unit number -
  Integer, Intent(In) :: ihind  ! - hindcasts flag -
  Integer, Intent(In) :: iretro ! - retroactive forecast flag -
!
  Logical, Intent(In) :: lres ! - read results? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write category settings
  If (.not.lres) Then
     Write (Unit=iout, Err=1) refs, irefs, ithr, pthr, thr
!
! Read results
  Else
     Write (Unit=iout, Err=1) pobs
     Write (Unit=iout, Err=1) tobs
     Write (Unit=iout, Err=1) tobst
     If (ihind == 1) Then
        Write (Unit=iout, Err=1) ifor
        Write (Unit=iout, Err=1) iobs
     End If
     If (iretro == 1) Then
        Write (Unit=iout, Err=1) irobs
        Write (Unit=iout, Err=1) irfor
        Write (Unit=iout, Err=1) rfps
        Write (Unit=iout, Err=1) rodds(:,:,:,2)
        Write (Unit=iout, Err=1) rodds(:,:,:,1)
     End If
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
 End Subroutine proj_write_categories
!
!
!
 Function close_categories()
!
! Frees memory allocated for categories
!
! Function Type
  Integer :: close_categories
!
! Functions and Subroutines
!
! Intrinsic Functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory allocated for categories
  If (Allocated(irfor)) Deallocate (irfor)
  If (Allocated(irobs)) Deallocate (irobs)
  If (Allocated(ifor))  Deallocate (ifor)
  If (Allocated(iobs))  Deallocate (iobs)
  If (Allocated(tobsx)) Deallocate (tobsx)
  If (Allocated(tobst)) Deallocate (tobst)
  If (Allocated(tobs))  Deallocate (tobs)
  If (Allocated(pobs))  Deallocate (pobs)
  close_categories = 0
!
  Return
 End Function close_categories
End Module categories
