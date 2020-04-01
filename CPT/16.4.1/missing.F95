! Author: Simon Mason
Module missing
!
! Modules
  Use data_numbers, Only: rp
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: close_missing, compress, init_kav, init_missing, proj_read_missing1, proj_read_missing2, &
            proj_write_missing1, proj_write_missing2, replace_missing_xy, replace_missing_z, replace_missing_zyp
!
! Arrays
!
! Integer arrays
  Integer, Dimension(:), Allocatable, Private :: iwk ! - integer workspace -
!
! Real arrays
  Real(Kind=rp), Dimension(:), Allocatable, Private :: xave ! - means of explanatory variables -
  Real(Kind=rp), Dimension(:), Allocatable, Private :: rwk  ! - double precision workspace -
!
! Logical arrays
  Logical, Dimension(:), Allocatable, Public :: kfuse ! - used forecasts? -
  Logical, Dimension(:), Allocatable, Public :: ktuse ! - used training cases? -
  Logical, Dimension(:), Allocatable, Public :: kcuse ! - used climatological and training cases? -
!
  Logical, Dimension(:,:), Pointer, Public :: kavx => Null() ! - X cases available? -
  Logical, Dimension(:,:), Pointer, Public :: kavy => Null() ! - Y cases available? -
  Logical, Dimension(:,:), Pointer, Public :: kavz => Null() ! - Z cases available? -
!
! Scalars
!
! Integer scalars
  Integer, Public :: immx ! - method for replacing missing X data flag -
  Integer, Public :: immy ! - method for replacing missing Y data flag -
  Integer, Public :: ipmx ! - maximum % of missing X data -
  Integer, Public :: ipmy ! - maximum % of missing X data -
  Integer, Public :: ipvx ! - maximum % of missing X variables -
  Integer, Public :: ipvy ! - maximum % of missing X variables -
  Integer, Public :: nnsx ! - number of X near-neighbours -
  Integer, Public :: nnsy ! - number of Y near-neighbours -
!
  Integer, Private :: lnn  ! - minimum of nns and nt -
  Integer, Private :: liwk ! - integer workspace dimensions -
  Integer, Private :: lrwk ! - real workspace dimensions -
!
! Real scalars
  Real(Kind=rp), Public :: xmiss ! - X missing values -
  Real(Kind=rp), Public :: ymiss ! - Y missing values -
!
Contains
!
!
 Subroutine init_kav (nt, nfl, nss, kav, ifail)
!
! Initialises memory for cases available
!
! On exit, returns:
!    0 if no errors
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl ! - number of fields and lagged fields -
  Integer, Intent(In) :: nss ! - number of seasons -
  Integer, Intent(In) :: nt  ! - number of cases -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Logical, Dimension(:,:), Pointer :: kav ! - cases available? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Associated, Max
!
! Executable Statements
!
! Alloculate memory for cases available
  If (Associated(kav)) Then
     Deallocate (kav)
     Nullify (kav)
  End If
  Allocate (kav(nt,nfl/nss), Stat=ifail)
  If (ifail /= 0) Then
     ifail = 1
     Return
  End If
  kav(1:Max(nss,1),:) = .true.
  If (nt > nss) kav(nss+1:,:) = .false.
  ifail = 0
!
  Return
 End Subroutine init_kav
!
!
!
 Function init_missing(lktuse, nt, nx, ny) &
          Result (init)
!
! Initialises memory and settings for reading data
!
! On exit, returns:
!    0 if no errors
!    1 if problem allocating memory
!
! Modules
  Use arrays, Only: iusex, iusey
!
! Function type
  Integer :: init
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nt ! - number of cases -
  Integer, Intent(In) :: nx ! - number of X-variables -
  Integer, Intent(In) :: ny ! - number of Y-variables -
!
  Logical, Intent(In) :: lktuse ! - initialise ktuse? -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Allocate data space
  init = 1
! - X variables flags -
  Allocate (iusex(nx), Stat=ifail)
  If (ifail /= 0) Return
! - Y variables flags -
  Allocate (iusey(ny), Stat=ifail)
  If (ifail /= 0) Return
! - available cases flags -
  Allocate (kcuse(nt), Stat=ifail)
  If (ifail /= 0) Return
! - available training cases flags -
  If (lktuse) Then
     Allocate (ktuse(nt))
     If (ifail /= 0) Return
  End If
!
  init = 0
!
  Return
 End Function init_missing
!
!
!
 Function replace_missing_xy(ianal, icex, nt, ntc, ntcy, k1, nsn, lensn, lb, lssr, linit, lcsv, lprst, nu) &
          Result (replace)
!
! Replaces missing X and Y data and compresses arrays
!
! Modules
  Use data_numbers,        Only: one
  Use data_cpt_constants,  Only: ia_pfv
  Use data_time_constants, Only: isq_sn, nmn
  Use time,                Only: pprd, &
                                 Operator(+), Operator(<), &
                                 nss, &
                                 get_cdate
  Use screen,              Only: window_update
  Use progress_meter,      Only: set_progress_increment
  Use errors,              Only: cpt_error
  Use arrays,              Only: iusex, iusey, x, y
  Use iofiles,             Only: mfile, mxa, mya, nx, ny, xfile
  Use settings,            Only: izero
  Use fields,              Only: mfield, xfield, yfield
!
! Function type
  Integer :: replace
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis indicator -
  Integer, Intent(In) :: icex  ! - climatological period lies outside the range of the training period flag -
  Integer, Intent(In) :: nsn   ! - number of seasons -
  Integer, Intent(In) :: lb    ! - lower-bound -
  Integer, Intent(In) :: lensn ! - total length of season -
  Integer, Intent(In) :: k1    ! - index of first case in training period -
  Integer, Intent(In) :: nt    ! - number of cases -
  Integer, Intent(In) :: ntc   ! - number of cases within climatological and training periods -
  Integer, Intent(In) :: ntcy  ! - number of years within climatological and training periods -
!
  Logical, Intent(In) :: lcsv  ! - calculate seasonal values? -
  Logical, Intent(In) :: linit ! - initialise? -
  Logical, Intent(In) :: lprst ! - include persistence component? -
!
! Output scalars
  Integer, Intent(Out) :: nu ! - number of used cases -
!
! Input arrays
  Logical, Dimension(lb:), Intent(In) :: lssr ! - season has shortened record? -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
  Integer :: l     ! - month index -
  Integer :: l0    ! - initial index -
  Integer :: l1    ! - initial index -
  Integer :: l2    ! - initial index -
  Integer :: na    ! - adjustment to n -
  Integer :: nact  ! - number of actions -
  Integer :: nrep  ! - number of repeats -
  Integer :: ns    ! - number of seasons to replace at the same time -
  Integer :: nsd   ! - number of seasons to replace individually -
  Integer :: nsi   ! - season increment -
!
  Logical :: lcomp ! - compress? -
!
  Type(pprd) :: error_prd ! - period of erroneous data -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count, Min, Real
!
! Executable Statements
!
! Initialise missing value replacement
  replace = -1
  If (linit) Then
     If (init_missing(.false., ntc, nx, ny) /= 0) Return
     l0 = 1
  Else
     l0 = lb
  End If
!
! Initialise progress increment for replacing missing values
  nact = 2*(ny + nx + ntc) + 3
  Call set_progress_increment (one/Real(nact, Kind=rp),.false.)
!
! Replace missing Y values and check for missing variables
  kcuse(:) = .true.
  If (lcsv) Then ! - monthly data -
     ns = 1
     nrep = Min(lensn, nmn)
     nsi = 1
     nsd = 1
     lcomp = .false.
  Else If (lprst) Then ! - persistence component -
     ns = 2
     nrep = 0
     nsi = 1
     nsd = 1
     lcomp = .true.
  Else
     ns = 1
     nrep = nsn
     nsi = nrep
     If (mfile%iseq /= isq_sn) Then
        nsd = 1
     Else
        nsd = nss
     End If
     lcomp = .true.
  End If
  l1 = 1
  l2 = 1
  Do l = l0, nrep
     If (.not.lssr(l)) Then
        na = 0
     Else
        na = 1
     End If
     If (nsi > 1) l1 = l ! - seasons are collated in kcuse - 
     If (nsd > 1) l2 = l
     Call replace_missing (lcomp, mfile%nfl/nsd, yfield(1:)%nv, ntcy-na, ns, y(:,:,l:l+ns-1), izero, immy, ipmy, ipvy, &
          nnsy, yfield(l2:)%rmiss, yfield(1:)%nva, mya, iusey(:), kcuse(l1::nsi), ifail)
     Select Case (ifail)
      Case (0)
        Call window_update (mya)
      Case (1, 2)
        Call cpt_error ('replace_missing_xy', ifail, .false., &
             c_arg1='Y', c_arg2=Trim(mfile%cfile%ffile))
        Return
      Case (3, 4, 5)
        Call cpt_error ('replace_missing_xy', ifail, .true.)
     End Select
     If (lcsv) kcuse(:) = .true.
  End Do
  If (lprst) yfield(0)%nva = yfield(1)%nva
  mfield(1)%nva = yfield(1)%nva
!
! Replace missing X values and check for missing variables
  Select Case (ianal)
   Case Default
     Call replace_missing (lcomp, xfile%nfl, xfield(:)%nv, nt, 1, x(:,:,:), 0, immx, ipmx, ipvx, nnsx, &
          xfield(:)%rmiss, xfield(:)%nva, mxa, iusex(:), kcuse(k1:), ifail)
   Case (ia_pfv)
     If (.not.lcsv) Then ! - postpone until after calculation of seasonal values because of data compression -
        Call identify_missing_fps (xfile%nfl/nss, xfield(:)%nv, nt, xfile%ngs, x(:,:,:), ipvx, xfield(:)%rmiss, &
             xfield(:)%nva, mxa, iusex(:), kcuse(k1:), ifail)
     End If
  End Select
  Select Case (ifail)
   Case (0)
     Call window_update (mxa)
   Case (1, 2)
     Call cpt_error ('replace_missing_xy', ifail, .false., &
          c_arg1='X', c_arg2=Trim(xfile%cfile%ffile))
     Return
   Case (3, 4, 5)
     Call cpt_error ('replace_missing_xy', ifail, .true.)
   Case (:-1)
     error_prd = xfile%prd1 + ((xfile%it1 - ifail) - 2)
     ifail = 6
     Call cpt_error ('replace_missing_xy', ifail, .false., &
          c_arg1=Trim(xfile%cfile%ffile), c_arg2=get_cdate(error_prd, 1, 2))
     replace = 1
     Return
  End Select
!
! Compress to eliminate unused cases
  If ((.not.lcsv) .and. nsn == 1) Then
     Select Case (icex)
      Case (0)
        Call compress (ntc, kcuse(:), nu, &
             aa=x(1:mxa,:,:), ba=y(1:mya,:,:))
      Case Default
        Call compress (ntc, kcuse(:), nu, &
             ba=y(1:mya,:,:))
        Call compress (nt, kcuse(k1:), nu, &
             aa=x(1:mxa,:,:))
     End Select
     Select Case (icex)
      Case (-1, 1, 3)
        nu = Count(kcuse(k1:k1+nt-1))
      Case (-2, 2)
        nu = Count(kcuse(1:nt))
     End Select
  Else ! - postpone compression for monthly and seasonal data -
     nu = ntc
  End If
!
  replace = 0
!
  Return
!
 Contains
!
!
  Subroutine replace_missing (lcomp, nfl, nv, n, ns, v, izero, imm, ipm, ipv, nns, vmiss, nva, mva, iuse, kuse, ifail)
!
! Replaces missing values and removes missing variables
!
! On exit:
!    ifail = 0 No errors
!    ifail = 1 all fields are constant
!    ifail = 2 too many missing values in all fields
!    ifail = 3 at least one of the fields has too many missing variables (missing values are replaced in the other fields) -
!    ifail = 4 insufficient variables for near-neighbours
!    ifail = 5 insufficient memory for near-neighbours
!
! Modules
  Use data_numbers,   Only: zero, half, one, oneh, onem
  Use maths,          Only: norm_rand
  Use progress_meter, Only: set_progress_increment, update_progress_meter
  Use arrays,         Only: rwk, vcopy, &
                            insertion_sort
  Use distribs,       Only: q_empirical
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl   ! - number of fields and lagged fields -
  Integer, Intent(In) :: n     ! - number of cases -
  Integer, Intent(In) :: ns    ! - number of simultaneities (for persistence component) -
  Integer, Intent(In) :: izero ! - zero-bound flag -
  Integer, Intent(In) :: ipm   ! - maximum % of missing values -
  Integer, Intent(In) :: ipv   ! - maximum % of missing variables -
  Integer, Intent(In) :: nns   ! - number of near-neighbours -
!
  Logical, Intent(In) :: lcomp ! - compress? -
!
! Input/output scalars
  Integer, Intent(InOut) :: imm ! - missing value replacement method -
!
! Output scalars
  Integer, Intent(Out) :: mva   ! - total number of non-missing variables -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nv ! - number of variables by field -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: vmiss ! - missing value indicators -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:,:), Intent(InOut) :: v ! - variables (missing variables deleted on output if lcomp) -
!
  Logical, Dimension(:), Intent(InOut) :: kuse ! - used cases? -
!
! Output arrays
  Integer, Dimension(:), Intent(Out) :: iuse ! - indices of used variables -
  Integer, Dimension(:), Intent(Out) :: nva  ! - number of non-missing variables by field -
!
! Locals
!
! Local arrays
  Integer, Dimension(ns) :: nmiss ! - actual number of missing values -
!
! Local scalars
  Integer :: i     ! - variable index -
  Integer :: ii    ! - non-missing variable index -
  Integer :: j     ! - variable index by field -
  Integer :: k     ! - time index -
  Integer :: l     ! - field index -
  Integer :: isi   ! - simultaneity index -
  Integer :: imiss ! - maximum number of missing cases permitted -
  Integer :: kmiss ! - maximum number of missing variables permitted -
  Integer :: mv    ! - total number of variables -
  Integer :: nnm   ! - number of non-missing values -
!
  Real(Kind=rp) :: vbar ! - mean -
  Real(Kind=rp) :: vvar ! - variance -
  Real(Kind=rp) :: vmed ! - median -
  Real(Kind=rp) :: vm   ! - current missing value -
  Real(Kind=rp) :: tol  ! - missing value rounding tolerance -
!
  Logical :: lconstant ! - constant fields? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Abs, All, Any, Count, Nint, Real, Sum
!
! Executable Statements
!
! Identify which variables to use
  mva = 0
  nva(:) = 0
  iuse(:) = 0
  lconstant = .true.
  i = 0
  imiss = Nint(Real(n*ipm, Kind=rp)/oneh) ! - maximum number of missing values per series -
  lconstant = .true.
  Do l = 1, nfl
     tol = Abs(vmiss(l))/onem
     vm = vmiss(l)
     Do j = 1, nv(l)
        i = i + 1
! - avoid rounding errors on missing values -
        Where (v(i,1:n,1:ns) > vm-tol .and. v(i,1:n,1:ns) < vm+tol) v(i,1:n,1:ns) = vm
! - count number of missing values -
        Do isi = 1, ns
           If (Any(v(i,2:n,isi) /= v(i,1,isi))) Then
              lconstant = .false.
              nmiss(isi) = Count(v(i,1:n,isi) == vm)
           Else
              nmiss(isi) = n
           End If
        End Do
! - retain variables with sufficient observations -
        If (Any(nmiss(:) <= imiss)) Then
           nva(l) = nva(l) + 1
           mva = mva + 1
           iuse(mva) = i
        End If
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Do
  mv = i
!
! Check number of available variables
  If (All(nva(1:nfl) > 0)) Then
     ifail = 0
  Else
     If (lconstant) Then ! - check for constant fields -
        ifail = 1
        Return
     Else If (All(nva(1:nfl) == 0)) Then ! - too many missing values in all fields -
        ifail = 2
        Return
     Else If (Any(nva(1:nfl) == 0)) Then ! - at least one field with too many missing values -
        ifail = 3
     End If
  End If
  If (update_progress_meter(.false.) /= 0) Return
!
! Identify which cases to use
  kmiss = Nint(Real(mva*ipv, Kind=rp)/oneh) ! - maximum number of missing series per case -
  Do k = 1, n
     i = 0
     ii = 1
     nmiss(:) = 0
     Do l = 1, nfl
        Do j = 1, nva(l)
           i = i + 1
           If (iuse(ii) /= i) Cycle
           ii = ii + 1
           Do isi = 1, ns
              If (v(i,k,isi) == vmiss(l)) nmiss(isi) = nmiss(isi) + 1
           End Do
        End Do
     End Do
     If (Any(nmiss(:) > kmiss)) kuse(k) = .false.
     If (update_progress_meter(.false.) /= 0) Return
  End Do
!
! Replace missing values
  If (imiss > 0) Then
     Call set_progress_increment (Real(mv, Kind=rp)/Real(mva, Kind=rp), .true.)
! - initialise for near-neighbours -
     If (imm == 4) Then
        If (mva > 1) Then
           Call init_neighbours (nns, mva, n, ns, ifail)
           If (ifail /= 0) Then
              imm = 1
              ifail = 5
           End If
        Else
           imm = 1
           ifail = 4
        End If
        If (ifail /= 0) Return
     End If
! - replace missing values -
     i = 0
     Do l = 1, nfl
        Do j = 1, nva(l)
           i = i + 1
           Do isi = 1, ns
              nmiss(isi) = Count(v(iuse(i),1:n,isi) == vmiss(l))
              If (nmiss(isi) > 0) Then
                 Select Case (imm)
                  Case (1) ! - long-term mean -
                    vbar = Sum(v(iuse(i),1:n,isi), Mask=(v(iuse(i),1:n,isi) /= vmiss(l)))/Real(n-nmiss(isi), Kind=rp)
                    Do k = 1, n
                       If (.not.kuse(k)) Cycle
                       If (v(iuse(i),k,isi) == vmiss(l)) v(iuse(i),k,isi) = vbar
                    End Do
                  Case (2) ! - long-term median -
                    nnm = 0
                    Do k = 1, n
                       If (.not.kuse(k)) Cycle
                       If (v(iuse(i),k,isi) /= vmiss(l)) Then
                          nnm = nnm + 1
                          rwk(nnm) = v(iuse(i),k,isi)
                       End If
                    End Do
                    Call insertion_sort (nnm, rwk(:), 'a')
                    vmed = q_empirical(rwk(:), nnm, half)
                    Where (v(iuse(i),1:n,isi) == vmiss(l)) v(iuse(i),1:n,isi) = vmed
                  Case (3) ! - random numbers -
                    vbar = Sum(v(iuse(i),1:n,isi), Mask=(v(iuse(i),1:n,isi) /= vmiss(l)))/Real(n-nmiss(isi), Kind=rp)
                    vvar = Sum((v(iuse(i),1:n,isi)-vbar)**2, Mask=(v(iuse(i),1:n,isi) /= vmiss(l)))/Real(n-nmiss(isi)-1, Kind=rp)
                    Do k = 1, n
                       If (.not.kuse(k)) Cycle
                       If (v(iuse(i),k,isi) /= vmiss(l)) Then
                          v(iuse(i),k,isi) = norm_rand(vbar, vvar)
                          If (izero == 1) Then
                             If (v(iuse(i),k,isi) < zero) v(iuse(i),k,isi) = zero
                          End If
                       End If
                    End Do
                  Case (4) ! - near-neighbours -
                    Do k = 1, n
                       If (kuse(k)) Then
                          If (v(iuse(i),k,isi) == vmiss(l)) Then
                             vcopy(i,k,isi) = &
                                  neighbour(nns, nfl, nva, n, v(iuse(1:mva),1:n,isi), i, l, v(iuse(1:mva),k,isi), vmiss(:))
                             If (izero == 1) Then
                                If (vcopy(i,k,isi) < zero) vcopy(i,k,isi) = zero
                             End If
                          Else
                             vcopy(i,k,isi) = v(iuse(i),k,isi)
                          End If
                       Else
                          vcopy(i,k,isi) = v(iuse(i),k,isi)
                       End If
                    End Do
                 End Select
              Else
                 If (imm == 4) vcopy(i,1:n,isi) = v(iuse(i),1:n,isi)
              End If
           End Do
           If (update_progress_meter(.false.) /= 0) Return
        End Do
     End Do
! - uninitialise for near-neighbours -
     If (imm == 4) Then
        v(iuse(1:mva),1:n,1:ns) = vcopy(1:mva,1:n,1:ns)
     End If
     Call uninit_neighbours ()
     Call set_progress_increment (Real(mva, Kind=rp)/Real(mv, Kind=rp), .true.)
  Else
     Call set_progress_increment (Real(mv, Kind=rp), .true.)
     If (update_progress_meter(.false.) /= 0) Return
     Call set_progress_increment (one/Real(mv, Kind=rp), .true.)
  End If
  ifail = 0
!
! Remove unused variables
  If (lcomp) Then
     Do i = 1, mva
        If (iuse(i) > i) v(i,1:n,1:ns) = v(iuse(i),1:n,1:ns)
     End Do
     If (mva < mv) v(mva+1:,:,:) = vmiss(1)
  End If
!
  Return
  End Subroutine replace_missing
!
!
!
  Subroutine identify_missing_fps (nfl, nv, n, ngs, fps, ipv, vmiss, nva, mva, iuse, kuse, ifail)
!
! Checks that probabilities sum to one, and removes variables with too many missing forecast probabilities
!
! On exit:
!    ifail = 0 no errors
!    ifail = 1 mva=0
!    ifail = 2 at least one of the fields have no missing variables
!    ifail =-1 probabilities do not add to 1 or 100
!
! Modules
  Use data_numbers,   Only: zero, one, oneh, tolp
  Use maths,          Only: approx
  Use progress_meter, Only: update_progress_meter
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl ! - number of fields and lagged fields -
  Integer, Intent(In) :: n   ! - number of cases -
  Integer, Intent(In) :: ngs ! - number of categories -
  Integer, Intent(In) :: ipv ! - maximum % of missing variables -
!
! Output scalars
  Integer, Intent(Out) :: mva   ! - total number of non-missing variables -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nv ! - number of variables by field -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: vmiss ! - missing value indicators -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:,:), Intent(InOut) :: fps ! - forecast probabilities (missing variables deleted on output) -
!
  Logical, Dimension(:), Intent(InOut) :: kuse ! - used cases? -
!
! Output arrays
  Integer, Dimension(:), Intent(Out) :: nva  ! - number of non-missing variables by field -
  Integer, Dimension(:), Intent(Out) :: iuse ! - indices of used variables -
!
! Locals
!
! Local scalars
  Integer :: i     ! - variable index -
  Integer :: j     ! - variable index by field -
  Integer :: ij    ! - available variable index -
  Integer :: k     ! - time index -
  Integer :: l     ! - field index -
  Integer :: m     ! - number of variables -
  Integer :: kmiss ! - maximum number of missing variables permitted -
  Integer :: nmiss ! - actual number of missing values -
!
  Real(Kind=rp) :: sp ! - sum of probabilities -
!
  Logical :: luse ! - used variable? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic All, Any, Count, Real, Sum
!
! Executable Statements
!
! Identify which variables to use
  mva = 0
  nva(:) = 0
  iuse(:) = 0
  i = 0
  Do l = 1, nfl
     Do j = 1, nv(l)
        i = i + 1
        Do k = 1, n
           sp = Sum(fps(i,k,1:ngs))
           If (approx(sp, one, utol=tolp) .or. approx(sp/oneh, one, utol=tolp)) Then
              Continue
           Else If (All(fps(i,k,1:ngs) > vmiss(l)-tolp) .and. All(fps(i,k,1:ngs) < vmiss(l)+tolp)) Then
              fps(i,k,1:ngs) = -one
           Else If (All(fps(i,k,1:ngs) < zero)) Then
              fps(i,k,1:ngs) = -one
           Else If (approx(sp, zero, utol=tolp)) Then
              fps(i,k,1:ngs) = -one
           Else ! - trap invalid probability sums -
              ifail = -k
              Return
           End If
        End Do
        If (Any(fps(i,1:n,1:ngs) > -one)) Then
           nva(l) = nva(l) + 1
           mva = mva + 1
           iuse(mva) = i
        End If
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Do
!
! Check number of available variables
  If (All(nva(1:nfl) > 0)) Then
     ifail = 0
  Else
     If (All(nva(1:nfl) == 0)) Then
        ifail = 1
     Else
        ifail = 2
     End If
     Return
  End If
!
! Compress data array
  Do i = 1, mva
     If (iuse(i) > i) fps(i,1:n,1:ngs) = fps(iuse(i),1:n,1:ngs)
  End Do
  If (update_progress_meter(.false.) /= 0) Return
!
! Determine maximum number of missing series per case
  kmiss = Nint(Real(mva*ipv, Kind=rp)/oneh)
!
! Identify which cases to use
  nmiss = 0
  Do k = 1, n
     If (kuse(k)) Then
        nmiss = Count(fps(1:mva,k,1) < zero)
        If (nmiss > kmiss) kuse(k) = .false.
     End If
     If (update_progress_meter(.false.) /= 0) Return
  End Do
!
! Recheck which variables to use
  m = mva
  mva = 0
  nva(1:nfl) = 0
  i = 0
  ij = 1
  Do l = 1, nfl
     Do j = 1, nv(l)
        i = i + 1
        If (iuse(ij) == i) Then
           luse = .false.
           Do k = 1, n
              If (.not.kuse(k)) Then
                 fps(ij,k,1:ngs) = -one
              Else If (fps(ij,k,1) > -one) Then
                 luse = .true.
              End If
           End Do
           ij = ij + 1
           If (luse) Then
              nva(l) = nva(l) + 1
              mva = mva + 1
              iuse(mva) = i
           Else If (ij <= m) Then
              iuse(ij-1:m-1) = iuse(ij:m)
              iuse(m) = 0
              fps(ij-1:m-1,1:n,1:ngs) = fps(ij:m,1:n,1:ngs)
              fps(m,1:n,1:ngs) = -one
              m = m - 1
           End If
        End If
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Do
!
  Return
  End Subroutine identify_missing_fps
 End Function replace_missing_xy
!
!
!
 Function replace_missing_z(ianal,nsn) &
          Result (replace)
!
! Replaces missing values and removes missing variables in forecast data
!
! Modules
  Use data_cpt_constants, Only: ng
  Use screen,             Only: window_update
  Use arrays,             Only: iusex, z
  Use iofiles,            Only: mxa, mza, zfile
  Use settings,           Only: nf
  Use fields,             Only: xfield, zfield
!
! Function type
  Integer :: replace
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ianal ! - analysis indicator -
  Integer, Intent(In) :: nsn   ! - number of seasons -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Replace missing values
  mza = mxa
  zfield(:)%nva = xfield(:)%nva
  kfuse(:) = .true.
  Select Case (ianal)
   Case (1:4)
     Call replace_missing (zfile%nfl, zfield(:)%nva, iusex(:), nf, nsn, immx, zfield(:)%rmiss, z(:,:,1), kfuse)
     ifail = 0
   Case (5)
     Call identify_missing_fps_z (zfile%nfl, zfield(:)%nva, iusex(:), nf, ng, zfield(:)%rmiss, z(:,:,:), ifail)
  End Select
  Call window_update (mza)
  replace = ifail
!
  Return
!
 Contains
!
!
  Subroutine replace_missing (nfl, nva, iuse, nf, nsn, imm, zmiss, z, kuse)
!
! Replaces missing values and removes missing variables in forecast data
!
! Modules
  Use data_numbers,   Only: tol, oneh
  Use maths,          Only: norm_rand
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: x
  Use iofiles,        Only: mxa
  Use settings,       Only: nu
  Use scores,         Only: get_mean, get_median, get_var
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl ! - number of fields and lagged fields -
  Integer, Intent(In) :: nf  ! - number of forecasts -
  Integer, Intent(In) :: nsn ! - number of seasons -
  Integer, Intent(In) :: imm ! - missing value replacement method -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nva  ! - number of used variables per field -
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: zmiss ! - missing value indicators -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: z ! - variables (missing variables deleted on output) -
!
  Logical, Dimension(:), Intent(InOut) :: kuse ! - used cases? -
!
! Locals
!
! Local scalars
  Integer :: i      ! - variable index -
  Integer :: j      ! - variable index by field -
  Integer :: k      ! - time index -
  Integer :: l      ! - field index -
  Integer :: i1, i2 ! - indices -
  Integer :: immu   ! - used missing value replacement option -
  Integer :: kmiss  ! - maximum number of missing variables permitted -
  Integer :: nmiss  ! - actual number of missing values -
  Integer :: ifail  ! - error indicator -
!
  Real(Kind=rp) :: vbar ! - mean -
  Real(Kind=rp) :: vvar ! - variance -
  Real(Kind=rp) :: vmed ! - median -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count
!
! Executable Statements
!
! Initalise memory
  immu = imm
  If (immu == 4) Then
     Call init_neighbours (nnsx, mxa, nu, nsn, ifail)
     If (ifail /= 0) immu = 1
  End If
!
! Remove unused variables and replace missing values
  i = 0
  Do l = 1, nfl
     Do j = 1, nva(l)
        i = i + 1
        If (i < iuse(i)) z(i,1:nf) = z(iuse(i),1:nf)
        Where (z(i,1:nf) > zmiss(l)-tol .and. z(i,1:nf) < zmiss(l)+tol) z(i,1:nf) = zmiss(l)
        nmiss = Count(z(i,1:nf) == zmiss(l))
        If (nmiss /= 0) Then
           Select Case (immu)
            Case (1) ! - long-term mean -
              vbar = get_mean(nu, x(i,1:nu,1), zmiss(l))
              Where (z(i,1:nf) == zmiss(l)) z(i,1:nf) = vbar
            Case (2) ! - long-term median -
              vmed = get_median(nu, x(i,1:nu,1), zmiss(l))
              Where (z(i,1:nf)==zmiss(l)) z(i,1:nf) = vmed
            Case (3) ! - random numbers -
              vbar = get_mean(nu, x(i,1:nu,1), zmiss(l))
              vvar = get_var(nu, x(i,1:nu,1), vbar, zmiss(l))
              Where (z(i,1:nf) == zmiss(l)) z(i,1:nf) = norm_rand(vbar, vvar)
            Case (4) ! - near-neighbours -
              Do k = 1, nf
                 If (z(i,k) == zmiss(l)) z(i,k) = neighbour(nnsx, nfl, nva, nu, x(:,:,1), i, l, z(:,k), zmiss(:))
              End Do
           End Select
        End If
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Do
  Call uninit_neighbours ()
!
! Identify which cases to use
  kmiss = Nint(Real(mxa*ipvx, Kind=rp)/oneh)
  Do k = 1, nf
     nmiss = 0
     i1 = 1
     i2 = 0
     Do l = 1, nfl
        i2 = i2 + nva(l)
        nmiss = nmiss + Count(z(i1:i2,k) == zmiss(l))
        i1 = i1 + nva(l)
     End Do
     If (nmiss > kmiss) kuse(k) = .false.
  End Do
!
  Return
  End Subroutine replace_missing
!
!
!
  Subroutine identify_missing_fps_z (nfl, nva, iuse, n, ngs, zmiss, fps, ifail)
!
! Checks that probabilities sum to one, and removes variables in forecast data
!
! Modules
  Use data_numbers,   Only: one, oneh, tolp
  Use maths,          Only: approx
  Use progress_meter, Only: update_progress_meter
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl ! - number of fields and lagged fields -
  Integer, Intent(In) :: n   ! - number of cases -
  Integer, Intent(In) :: ngs ! - number of categories -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nva  ! - number of used variables per field -
  Integer, Dimension(:), Intent(In) :: iuse ! - indices of used variables -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: zmiss ! - missing value indicators -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:,:), Intent(InOut) :: fps ! - forecast probabilities (missing variables deleted on output) -
!
! Locals
!
! Local scalars
  Integer :: i ! - variable index -
  Integer :: j ! - variable index by field -
  Integer :: k ! - time index -
  Integer :: l ! - field index -
!
  Real(Kind=rp) :: sp ! - sum of probabilities -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Sum
!
! Executable Statements
!
! Remove unused variables
  i = 0
  Do l = 1, nfl
     Do j = 1, nva(l)
        i = i + 1
        If (i < iuse(i)) fps(i,1:n,1:ngs) = fps(iuse(i),1:n,1:ngs)
!
! Remove unused cases
        Do k = 1, n
           sp = Sum(fps(i,k,1:ngs))
           If (approx(sp, one, utol=tolp) .or. approx(sp/oneh, one, utol=tolp)) Then
              Continue
           Else If (All(fps(i,k,1:ngs) > zmiss(l)-tolp) .and. All(fps(i,k,1:ngs) < zmiss(l)+tolp)) Then
              fps(i,k,1:ngs) = zmiss(l)
           Else ! - trap invalid probability sums -
              ifail = -k
              Return
           End If
        End Do
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Do
  ifail = 0
!
  Return
  End Subroutine identify_missing_fps_z
 End Function replace_missing_z
!
!
!
 Function replace_missing_zyp(ave, sdev) &
          Result (replace)
!
! Replaces missing values and removes missing variables in forecast persistence data
!
! Modules
  Use arrays,   Only: iusey, zyp
  Use iofiles,  Only: mfile
  Use settings, Only: nf
  Use fields,   Only: yfield
!
! Function type
  Integer :: replace
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: ave  ! - climatological averages -
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: sdev ! - climatological standard deviations -
!
! Executable Statements
!
! Replace missing values
  Call replace_missing (mfile%nfl, yfield(1:)%nva, nf, zyp(:,:,0), ave, sdev, immy, yfield(1:)%rmiss, iusey)
  replace = 0
!
  Return
!
 Contains
!
!
  Subroutine replace_missing (nfl, nva, nf, p, ave, sdev, imm, pmiss, iuse)
!
! Modules
  Use data_numbers,   Only: tol
  Use maths,          Only: norm_rand
  Use progress_meter, Only: update_progress_meter
  Use arrays,         Only: y
  Use iofiles,        Only: mya
  Use settings,       Only: igauss, istd, istdo, nu
  Use distribs,       Only: stdize, ustdize
  Use scores,         Only: get_mean, get_median, get_var
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfl ! - number of fields and lagged fields -
  Integer, Intent(In) :: nf  ! - number of forecasts -
  Integer, Intent(In) :: imm ! - missing value replacement method -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nva  ! - number of used variables per field -
  Integer, Dimension(:), Intent(In) :: iuse ! - used variables indices -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: pmiss ! - missing value indicator -
!
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: ave  ! - climatological averages -
  Real(Kind=rp), Dimension(:,0:), Intent(In) :: sdev ! - climatological standard deviations -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: p ! - variables (missing variables deleted on output) -
!
! Locals
!
! Local scalars
  Integer :: i     ! - variable index -
  Integer :: j     ! - variable index by field -
  Integer :: k     ! - time index -
  Integer :: l     ! - field index -
  Integer :: immu  ! - used missing value replacement option -
  Integer :: isn   ! - season index -
  Integer :: nmiss ! - actual number of missing values -
  Integer :: ifail ! - error indicator -
!
  Real(Kind=rp) :: vbar ! - mean -
  Real(Kind=rp) :: vvar ! - variance -
  Real(Kind=rp) :: vmed ! - median -
!
! Local arrays
  Real(Kind=rp), Dimension(:,:), Allocatable :: pt ! - temporary copy of data -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Count
!
! Executable Statements
!
! Initalise memory
  immu = imm
  If (immu == 4) Then
     Call init_neighbours (nnsy, mya, nu, 1, ifail)
     If (ifail /= 0) immu = 1
     Allocate (pt(mya,nf), Stat=ifail)
  End If
! - temporarily isolate persistence component -
  If ((istd /= istdo .and. istd /= 0) .or. istdo /= 0) Then
     Do isn = 0, 1
        Call ustdize (istd, igauss, isn, mya, nu, ave(:,isn), sdev(:,isn), y(1:mya,1:nu,isn))
     End Do
  End If
  y(1:mya,1:nu,0) = y(1:mya,1:nu,0) - y(1:mya,1:nu,1)
!
! Remove unused variables and replace missing values
  i = 0
  Do l = 1, nfl
     Do j = 1, nva(l)
        i = i + 1
        If (i < iuse(i)) p(i,1:nf) = p(iuse(i),1:nf)
        If (immu == 4) pt(i,:) = p(i,1:nf)
        Where (p(i,1:nf) > pmiss(l)-tol .and. p(i,1:nf) < pmiss(l)+tol) p(i,1:nf) = pmiss(l)
        nmiss = Count(p(i,1:nf) == pmiss(l))
        If (nmiss /= 0) Then
           Select Case (immu)
            Case (1) ! - long-term mean -
              vbar = get_mean(nu, y(i,1:nu,0), pmiss(l))
              Where (p(i,1:nf) == pmiss(l)) p(i,1:nf) = vbar
            Case (2) ! - long-term median -
              vmed = get_median(nu, y(i,1:nu,0), pmiss(l))
              Where (p(i,1:nf) == pmiss(l)) p(i,1:nf) = vmed
            Case (3) ! - random numbers -
              vbar = get_mean(nu, y(i,1:nu,0), pmiss(l))
              vvar = get_var(nu, y(i,1:nu,0), vbar, pmiss(l))
              Where (p(i,1:nf) == pmiss(l)) p(i,1:nf) = norm_rand(vbar, vvar)
            Case (4) ! - near-neighbours -
              Do k = 1, nf
                 If (p(i,k) == pmiss(l)) pt(i,k) = neighbour(nnsy, nfl, nva, nu, y(:,:,0), i, l, p(:,k), pmiss(:))
              End Do
           End Select
        End If
        If (update_progress_meter(.false.) /= 0) Return
     End Do
  End Do
  Call uninit_neighbours ()
!
! Replace forecast component
  If (immu == 4) p(1:mya,1:nf) = pt(:,:)
  y(1:mya,1:nu,0) = y(1:mya,1:nu,0) + y(1:mya,1:nu,1)
  If ((istd /= istdo .and. istd /= 0) .or. istdo /= 0) Then
     Do isn = 0, 1
        Call stdize (istd, igauss, isn, mya, nu, ave(:,isn), sdev(:,isn), y(1:mya,1:nu,isn))
     End Do
  End If
!
  Return
  End Subroutine replace_missing
 End Function replace_missing_zyp
!
!
!
 Subroutine compress (n, kuse, nu, &
            ia, aa, ba)
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - total number of cases -
!
! Output scalars
  Integer, Intent(Out) :: nu ! - number of used cases -
!
! Input arrays
  Logical, Dimension(:), Intent(In) :: kuse ! - used cases -
!
! Input/output arrays
  Integer, Dimension(:,:), Intent(InOut), Optional :: ia ! - data to be compressed -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(InOut), Optional :: aa ! - data to be compressed -
  Real(Kind=rp), Dimension(:,:,:), Intent(InOut), Optional :: ba ! - data to be compressed -
!
! Locals
!
! Local scalars
  Integer :: k ! - case index -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Compress arrays
  nu = 0
  Do k = 1, n
     If (.not.kuse(k)) Cycle
     nu = nu + 1
     If (nu < k) Then
        If (Present(ia)) ia(:,nu) = ia(:,k)
        If (Present(aa)) aa(:,nu,:) = aa(:,k,:)
        If (Present(ba)) ba(:,nu,:) = ba(:,k,:)
     End If
  End Do
!
  Return
 End Subroutine compress
!
!
!
 Function neighbour(nns, nfl, nva, n, v, i, l, w, vmiss)
!
! Estimates a missing value from near-neighbours
!
! Modules
  Use data_numbers, Only: zero, eps
!
! Function type
  Real(Kind=rp) :: neighbour
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nns ! - number of near-neighbours -
  Integer, Intent(In) :: nfl ! - number of fields and lagged fields -
  Integer, Intent(In) :: n   ! - number of cases -
  Integer, Intent(In) :: i   ! - current series -
  Integer, Intent(In) :: l   ! - current field / lagged-field -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nva ! - number of series per field -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: w     ! - data for missing case -
  Real(Kind=rp), Dimension(:), Intent(In) :: vmiss ! - missing value flags -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v  ! - data -
!
! Locals
!
! Local arrays
  Integer, Dimension(nns) :: ibest ! - near-neighbours -
!
  Real(Kind=rp), Dimension(nns) :: rbest ! - best correlation -
  Real(Kind=rp), Dimension(nns) :: vmu   ! - used missing value flags -
!
! Local scalars
  Integer :: j   ! - variable by field index -  
  Integer :: ij  ! - variable index -
  Integer :: ik  ! - case index -
  Integer :: il  ! - field / lagged-field index -
  Integer :: ir  ! - near-neighbour index -
  Integer :: nnm ! - number of non-missing cases -
  Integer :: nnu ! - number of used near-neighbours -
!
  Real(Kind=rp) :: xbar ! - mean -
  Real(Kind=rp) :: ybar ! - mean -
  Real(Kind=rp) :: df   ! - degrees of freedom -
  Real(Kind=rp) :: sxx  ! - sum of squares of independent variable -
  Real(Kind=rp) :: syy  ! - sum of squares of dependent variable -
  Real(Kind=rp) :: sxy  ! - sum of cross products -
  Real(Kind=rp) :: r    ! - correlation -
  Real(Kind=rp) :: b0   ! - regression constant -
  Real(Kind=rp) :: b1   ! - regression coefficient -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Any, Count, Min, Real, Sqrt
!
! Executable Statements
!
! Fit regression line
  ibest(:) = 0
  rbest(:) = zero
  vmu(:) = zero
  ij = 0
  Do il = 1, nfl
     Do j = 1, nva(il)
        ij = ij + 1
        If (ij == i) Cycle
        If (w(ij) == vmiss(il)) Cycle
! - calculate means -
        nnm = 0
        xbar = zero
        ybar = zero
        Do ik = 1, n
           If (v(i,ik) /= vmiss(l) .and. v(ij,ik) /= vmiss(il)) Then
              nnm = nnm + 1
              xbar = xbar + v(ij,ik)
              ybar = ybar + v(i,ik)
           End If
        End Do
        If (nnm <= 2) Cycle
        df = Real(nnm, Kind=rp)
        xbar = xbar/df
        ybar = ybar/df
! - calculate sum of squares and cross-products -
        sxx = zero
        syy = zero
        sxy = zero
        Do ik = 1, n
           If (v(i,ik) /= vmiss(l) .and. v(ij,ik) /= vmiss(il)) Then
              sxx = sxx + (v(ij,ik) - xbar)**2
              syy = syy + (v(i,ik) - ybar)**2
              sxy = sxy + (v(ij,ik) - xbar)*(v(i,ik) - ybar)
           End If
        End Do
! - calculate correlation -
        If (sxx > eps .and. syy > eps) Then
           r = sxy/Sqrt(sxx*syy)
           Do ir = 1, Min(nns, ij)
              If (r > rbest(ir)) Then
! - calculate coefficients -
                 If (ir < nns) Then
                    ibest(ir+1:nns) = ibest(ir:nns-1)
                    rbest(ir+1:nns) = rbest(ir:nns-1)
                    vmu(ir+1:nns) = vmu(ir:nns-1)
                 End If
                 ibest(ir) = ij
                 rbest(ir) = r
                 vmu(ir) = vmiss(il)
                 If (nns == 1 .or. Count(ibest(:) /= 0) == 0) Then
                    b1 = sxy/sxx
                    b0 = ybar - b1*xbar
                 End If
                 Exit
              End If
           End Do
        End If
     End Do
  End Do
!
! Estimate missing value
  If (Any(ibest(:) /= 0)) Then
     nnu = Count(ibest(:) > 0)
     If (nnu == 1) Then
        neighbour = b0 + b1*w(ibest(1))
     Else
        neighbour = get_neighbour(nns, nnu, n, v, w, vmu, i, l, ibest)
     End If
!
! Calculate mean if no close neighbours
  Else
     nnm = 0
     ybar = zero
     Do ik = 1, n
        If (v(i,ik) /= vmiss(l)) Then
           nnm = nnm + 1
           ybar = ybar + v(i,ik)
        End If
     End Do
     df = Real(nnm, Kind=rp)
     ybar = ybar/df
     neighbour = ybar
  End If
!
  Return
 End Function neighbour
!
!
!
 Subroutine init_neighbours (nns, nv, n, nsn, ifail)
!
! Initialises memory and settings to perform MLR
!
! Modules
  Use arrays, Only: bs=>b2, vcopy, xk, yc
  Use pcs,    Only: bz=>bz2, eof=>eofx2, sv=>svx2, ts=>tsx2, &
                    init_pc
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nns ! - number of used near-neighbours -
  Integer, Intent(In) :: nv  ! - total number of variables -
  Integer, Intent(In) :: n   ! - number of cases -
  Integer, Intent(In) :: nsn ! - number of seasons -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - memory allocation status -
!
! Locals
!
! Local scalars
  Integer :: mxt ! - maximum of mxa, and nu -
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Max
!
! Executable Statements
!
! Allocated additional memory
  Allocate (vcopy(nv,n,nsn), Stat=ifail)
  If (ifail /= 0) Return
  If (nns == 1) Return
!
! Allocate additional memory
! - means -
  Allocate (xave(nns), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - independent variables principal components -
  If (init_pc(nns, n, lnn, sv, eof, ts) /= 0) GoTo 1
! - principal component regression coefficients -
  Allocate (bz(lnn,1), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - regression coefficients -
  Allocate (bs(nv,1), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - copy of training period independent data -
  Allocate (xk(nns,n,1), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - copy of training period dependent data -
  Allocate (yc(n,1), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
! Allocate workspace
  mxt = Max(nns, n)
  lrwk = 3*lnn + Max(mxt, 4*lnn*(lnn+1))
  liwk = 8*lnn
  Allocate (iwk(liwk), Stat=ifail)
  If (ifail /= 0) GoTo 1
  Allocate (rwk(lrwk), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
  ifail = 0
  Return
!
! Errors
1 ifail = 1
!
  Return
 End Subroutine init_neighbours
!
!
!
 Function get_neighbour(nns, nnu, n, v, w, vmiss, i, l, ibest)
!
! Modules
  Use data_numbers, Only: zero, tol
#if DP==1
  Use cpt_lapack,   Only: gesdd=>dgesdd
#else
  Use cpt_lapack,   Only: gesdd=>sgesdd
#endif
  Use arrays,       Only: bs=>b2, xk, yc
  Use pcs,          Only: bz=>bz2, eof=>eofx2, sv=>svx2, ts=>tsx2
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nns ! - number of near-neighbours -
  Integer, Intent(In) :: n   ! - number of cases -
  Integer, Intent(In) :: i   ! - current series -
  Integer, Intent(In) :: l   ! - current field / lagged-field -
!
! Input/output scalars
  Integer, Intent(InOut) :: nnu ! - number of used near-neighbours -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: ibest ! - near-neighbours -
!
  Real(Kind=rp), Dimension(:), Intent(In) :: w     ! - data for missing case -
  Real(Kind=rp), Dimension(:), Intent(In) :: vmiss ! - missing value flags -
!
  Real(Kind=rp), Dimension(:,:), Intent(In) :: v ! - data -
!
! Function type
  Real(Kind=rp) :: get_neighbour
!
! Locals
!
! Local scalars
  Integer :: ie    ! - EOF index -
  Integer :: ik    ! - case index -
  Integer :: nnm   ! - number of non-missing cases -
  Integer :: npe   ! - number of positive singular values -
  Integer :: ifail ! - error indicator -
!
  Real(Kind=rp) :: ybar ! - average -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic All, Count, MatMul, Min, Real
!
! Executable Statements
!
! Standardise
1 nnm = 0
  xave(:) = zero
  ybar = zero
  Do ik = 1, n
     If (v(i,ik) /= vmiss(l) .and. All(v(ibest(1:nnu),ik) /= vmiss(1:nnu))) Then
        nnm = nnm + 1
        xk(1:nnu,nnm,1) = v(ibest(1:nnu),ik)
        xave(1:nnu) = xave(1:nnu) + v(ibest(1:nnu),ik)
        yc(nnm,1) = v(i,ik)
        ybar = ybar + v(i,ik)
     End If
  End Do
  If (nnm == 0) Then
     If (nnu > 1) Then
        nnu = nnu - 1
        GoTo 1
     Else
        nnm = 0
        ybar = zero
        Do ik = 1, n
           If (v(i,ik) /= vmiss(l)) Then
              nnm = nnm + 1
              ybar = ybar + v(i,ik)
           End If
        End Do
        get_neighbour = ybar
        Return
     End If
  End If
  xave(:) = xave(:)/Real(nnm, Kind=rp)
  ybar = ybar/Real(nnm, Kind=rp)
  Do ie = 1, nnu
     xk(ie,1:nnm,1) = xk(ie,1:nnm,1) - xave(ie)
  End Do
  Call gesdd ('S', nnu, nnm, xk(:,:,1), nns, sv, eof, nns, ts, lnn, rwk, lrwk, iwk, ifail)
  If (ifail /= 0) Then
     get_neighbour = ybar
     Return
  End If
  npe = Count(sv(1:Min(nnu, nnm)) > tol)
  Do ik = 1, nnm
     Do ie = 1, npe
        ts(ie,ik) = ts(ie,ik)/sv(ie)
     End Do
  End Do
  yc(1:nnm,1) = yc(1:nnm,1) - ybar
  bz(1:nnu,1) = MatMul(ts(1:nnu,1:nnm), yc(1:nnm,1))
  bs(1:nnu,1) = MatMul(eof(1:nnu,1:nnu), bz(1:nnu,1))
  get_neighbour = ybar
  Do ie = 1, nnu
     get_neighbour = get_neighbour + (w(ibest(ie)) - xave(ie))*bs(ie,1)
  End Do
!
  Return
 End Function get_neighbour
!
!
!
 Subroutine proj_read_missing1 (iin, nss, lsetx, ntx, nflx, lsety, nty, nfly, lsetz, ntz, nflz, ifail)
!
! Reads missing value settings from project file
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin  ! - input file unit number -
  Integer, Intent(In) :: nss  ! - number of seasons -
  Integer, Intent(In) :: nflx ! - number of X fields -
  Integer, Intent(In) :: nfly ! - number of Y fields -
  Integer, Intent(In) :: nflz ! - number of forecast fields -
  Integer, Intent(In) :: ntx  ! - number of X cases -
  Integer, Intent(In) :: nty  ! - number of Y cases -
  Integer, Intent(In) :: ntz  ! - number of forecasts -
!
  Logical, Intent(In) :: lsetx ! - X file set? -
  Logical, Intent(In) :: lsety ! - Y file set? -
  Logical, Intent(In) :: lsetz ! - forecast file set? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Read missing value settings
  Read (Unit=iin, Err=1, End=2) immx, immy, ipmx, ipmy, ipvx, ipvy, nnsx, nnsy
  ifail = 0
  If (lsetx) Call proj_read_amissing (iin, ntx, nflx, nss, kavx, ifail)
  If (ifail /= 0) Return
  If (lsety) Call proj_read_amissing (iin, nty, nfly, 1, kavy, ifail)
  If (ifail /= 0) Return
  If (lsetz) Call proj_read_amissing (iin, ntz, nflz, nss, kavz, ifail)
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
 Contains
!
!
 Subroutine proj_read_amissing (iin, nt, nfl, nss, kav, ifail)
!
! Reads missing value settings from project file
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin ! - input file unit number -
  Integer, Intent(In) :: nt  ! - number of cases -
  Integer, Intent(In) :: nfl ! - number of fields -
  Integer, Intent(In) :: nss ! - number of seasons -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
  Logical, Dimension(:,:), Pointer :: kav ! - cases available? -
!
! Executable Statements
!
! Read available cases flags
  Call init_kav (nt, nfl, nss, kav, ifail)
  If (ifail /= 0) GoTo 3
  Read (Unit=iin, Err=1, End=2) kav
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
! - problem initializing memory -
3 ifail = 3
  Return
!
  End Subroutine proj_read_amissing
 End Subroutine proj_read_missing1
!
!
!
 Subroutine proj_write_missing1 (iout, lsetx, lsety, lsetz, ifail)
!
! Writes missing value settings to project file
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output file unit number -
!
  Logical, Intent(In) :: lsetx ! - X file set? -
  Logical, Intent(In) :: lsety ! - Y file set? -
  Logical, Intent(In) :: lsetz ! - forecast file set? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write domain settings
  Write (Unit=iout, Err=1) immx, immy, ipmx, ipmy, ipvx, ipvy, nnsx, nnsy
  If (lsetx) Write (Unit=iout, Err=1) kavx
  If (lsety) Write (Unit=iout, Err=1) kavy
  If (lsetz) Write (Unit=iout, Err=1) kavz
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
 End Subroutine proj_write_missing1
!
!
!
 Subroutine proj_read_missing2 (iin, ifail)
!
! Reads missing value indicators from project file
!
! Modules
  Use arrays,   Only: iusex, iusey
  Use iofiles,  Only: nx, ny
  Use settings, Only: nt
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin  ! - input file unit number -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Read missing value indicators
  If (init_missing(.true., nt, nx, ny) /= 0) GoTo 4
  Read (Unit=iin, Err=1, End=2) iusex
  Read (Unit=iin, Err=1, End=2) iusey
  Read (Unit=iin, Err=1, End=2) kcuse
  Read (Unit=iin, Err=1, End=2) ktuse
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
! - end of file -
4 ifail = 4
  Return
 End Subroutine proj_read_missing2
!
!
!
 Subroutine proj_write_missing2 (iout, ifail)
!
! Writes missing value indicators to project file
!
! Modules
  Use arrays, Only: iusex, iusey
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iout ! - output file unit number -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Executable Statements
!
! Write missing value indicators
  Write (Unit=iout, Err=1) iusex
  Write (Unit=iout, Err=1) iusey
  Write (Unit=iout, Err=1) kcuse
  Write (Unit=iout, Err=1) ktuse
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
 End Subroutine proj_write_missing2
!
!
!
 Subroutine uninit_neighbours ()
!
! Frees memory
!
! Modules
  Use arrays, Only: bs=>b2, vcopy, xk, yc
  Use pcs,    Only: bz=>bz2, eof=>eofx2, sv=>svx2, ts=>tsx2
!
! Functions and Subroutines
!
! Intrinsic routines
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory
  If (Allocated(rwk))   Deallocate (rwk)
  If (Allocated(iwk))   Deallocate (iwk)
  If (Allocated(yc))    Deallocate (yc)
  If (Allocated(xk))    Deallocate (xk)
  If (Allocated(bs))    Deallocate (bs)
  If (Allocated(bz))    Deallocate (bz)
  If (Allocated(ts))    Deallocate (ts)
  If (Allocated(eof))   Deallocate (eof)
  If (Allocated(sv))    Deallocate (sv)
  If (Allocated(xave))  Deallocate (xave)
  If (Allocated(vcopy)) Deallocate (vcopy)
!
  Return
 End Subroutine uninit_neighbours
!
!
!
 Function close_missing()
!
! Frees memory allocated for missing values
!
! Modules
  Use arrays, Only: iusex, iusey
!
! Function type
  Integer :: close_missing
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory
  If (Allocated(ktuse)) Deallocate (ktuse)
  If (Allocated(kcuse)) Deallocate (kcuse)
  If (Allocated(iusey)) Deallocate (iusey)
  If (Allocated(iusex)) Deallocate (iusex)
!
  close_missing = 0
!
  Return
 End Function close_missing
End Module missing
