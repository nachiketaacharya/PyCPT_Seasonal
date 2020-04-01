Module ensemble
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: close_combine, ensemble_ave, get_combined, init_combine
!
Contains
!
!
 Subroutine init_combine (m, n, ifail)
!
! Initialises memory and settings to combine model predictions
!
! Modules
  Use iofiles, Only: mya
  Use pcs,     Only: bz, eofx, svx, tsx, xwk, &
                     init_pc
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: m ! - number of models -
  Integer, Intent(In) :: n ! - number of cases -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: lmn ! - minimum of the two input arguments -
!
! Executable Statements
!
! Allocate memory
! - independent variables principal components -
  If (init_pc(m, n, lmn, svx, eofx, tsx) /= 0) GoTo 1
! - principal component regression coefficients -
  Allocate (bz(lmn,mya), Stat=ifail)
  If (ifail /= 0) GoTo 1
! - workspace -
  Allocate (xwk(m,n), Stat=ifail)
  If (ifail /= 0) GoTo 1
!
  Return
!
! Errors
1 Call close_combine ()
!
  Return
 End Subroutine init_combine
!
!
!
 Subroutine get_combined (m, nms, b0, b1, x, c)
!
! Combines model predictions
!
! Modules
  Use data_numbers, Only: rp, zero
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: m   ! - number of locations -
  Integer, Intent(In) :: nms ! - number of models -
!
! Input arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: b0 ! - regression constants -
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: b1 ! - regression coefficients -
!
  Real(Kind=rp), Dimension(:,:,:), Intent(In) :: x ! - model predictions -
!
! Output arrays
  Real(Kind=rp), Dimension(:,:), Intent(Out) :: c ! - combined predictions -
!
! Locals
!
! Local scalars
  Integer :: i ! - location index -
  Integer :: j ! - model index -
!
! Executable Statements
!
! Combine models
  Do i = 1, m
     c(i,:) = zero
     Do j = 1, nms
        c(i,:) = c(i,:) + b0(j,i) + b1(j,i)*x(i,:,j)
     End Do
  End Do
!
  Return
 End Subroutine get_combined
!
!
!
 Subroutine close_combine ()
!
! Frees memory for combining model predictions
!
! Modules
  Use pcs, Only: bz, eofx, svx, tsx, xwk
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Allocate memory
  If (Allocated(xwk))  Deallocate (xwk)
  If (Allocated(bz))   Deallocate (bz)
  If (Allocated(tsx))  Deallocate (tsx)
  If (Allocated(eofx)) Deallocate (eofx)
  If (Allocated(svx))  Deallocate (svx)
!
  Return
 End Subroutine close_combine
!
!
!
 Subroutine ensemble_ave (nfs, nls, n, nv, nem, v)
!
! Averages ensemble members
!
! Modules
  Use data_numbers, Only: rp, zero
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfs ! - number of fields -
  Integer, Intent(In) :: nls ! - number of lagged-fields -
  Integer, Intent(In) :: n   ! - number of cases -
  Integer, Intent(In) :: nem ! - number of ensemble members -
!
! Input arrays
  Integer, Dimension(:), Intent(In) :: nv ! - number of variables -
!
! Input/output arrays
  Real(Kind=rp), Dimension(:,:), Intent(InOut) :: v ! - data -
!
! Locals
!
! Local scalars
  Integer :: i   ! - variable index -
  Integer :: ie  ! - ensemble member index -
  Integer :: ifd ! - field index -
  Integer :: ilf ! - lagged field index -
  Integer :: iv  ! - current variable -
  Integer :: iv1 ! - variable offset -
  Integer :: k   ! - time index -
  Integer :: l   ! - field/lagged field index -
!
! Local arrays
  Real(Kind=rp), Dimension(n) :: eave ! - ensemble average -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Real
!
! Executable Statements
!
! Average ensemble members
  iv1 = 0
  Do ifd = 1, nfs
     Do ilf = 1, nls
        l = (ifd - 1)*nem*nls + ilf
        Do i = 1, nv(l)
           eave(:) = zero
           Do ie = 1, nem
              iv = iv1 + ((ie - 1)*nls + (ilf - 1))*nv(l) + i
              Do k = 1, n
                 eave(k) = eave(k) + v(iv,k)
              End Do
           End Do
           eave(:) = eave(:)/Real(nem, Kind=rp)
           Do ie = 1, nem
              iv = iv1 + ((ie - 1)*nls + (ilf - 1))*nv(l) + i
              Do k = 1, n
                 v(iv,k) = eave(k)
              End Do
           End Do
        End Do
     End Do
     l = ifd*nem*nls
     iv1 = iv1 + nem*nls*nv(l)
  End Do
!
  Return
 End Subroutine ensemble_ave
End Module ensemble
