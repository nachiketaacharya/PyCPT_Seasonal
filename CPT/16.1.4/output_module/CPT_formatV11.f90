!
! CPT_formatV11 module of subroutines is used to save data in CPT V11 formats. 
! The output data are ASCII values. The user can call the following subroutines 
! to create appropriate data file.
! 
! 1. write_cpt_grid_v11  -->  gridded data
! 2. write_cpt_stns_v11  -->  station data
! 3. write_cpt_unrf_v11  -->  unreferenced data
!
! The data can be provided in single or double precision. Different versions of
! each routine are provided depending upon whether the data contain multiple fields
! and or lagged fields.
!
! Use of these subroutines:
! Include the following statement at the beginning of subroutines or modules.
!
!   USE CPT_formatV11
!
! Each routine requests information about the date for the first data in the file.
! In case of seasonal averages the date format involves a start and end date. For
! example, if the data are DJF averages from 1971/72 to 2000/01, then the first 
! data are for DJF 1971/72. The start date for the first data is therefore December 
! 1971 and the end date is Feburary 1972. The following values would then be set:
! 
!  period1_sdate_iyr=1971
!  period1_sdate_imn=12
!  period1_sdate_idy=0
!  period1_edate_iyr=1972
!  period1_edate_imn=2
!  period1_edate_idy=0
!
! Because it is implicit that the whole of December and the whole of February
! are included the days of the month are not needed and so can be set to 0.
!
! See the CPT help pages for more details on the dates.
!
! If there are dates with only missing values, then the output files can be saved
! with these dates omitted to save space. In this case, the optional argument, kuse,
! should be passed, and the argument nt set to the total number of cases including
! the ones to be omitted. The input data, v, should be compressed (i.e., it should 
! not contain these missing cases), but kuse should be set for each case (.true. for
! the cases that are present, and .false. for those that are missing and should be
! omitted). Note that for daily data, leap years are counted according to the 
! British implementation of the Gregorian calendar (see function ndays).
!
! These subroutines and functions are written by Simon Mason and modified by Lulin
! Song.
!
! This file was first created by Lulin Song on April 1, 2010.
! Modified by Simon Mason on April 20, 2010.
! Modified by Simon Mason on May 13, 2010.
! Modified by Simon Mason on May 16, 2011.
! Modified by Simon Mason on October 03, 2011.
! Modified by Simon Mason on May 11, 2011.
!
! $Id$
MODULE numbers
!
! Implicit declarations
  IMPLICIT NONE
!
! Parameters
!
! - private parameters -
  INTEGER, PARAMETER, PUBLIC :: sp=KIND(1.0e0)          ! - single precision kind parameter -
  INTEGER, PARAMETER, PUBLIC :: dp=KIND(1.0d0)          ! - double precision kind parameter -
!
! Real parameters
  REAL(KIND=sp), PARAMETER, PUBLIC :: zero_sp=  0.0_sp  ! - zero -
  REAL(KIND=sp), PARAMETER, PUBLIC ::  one_sp=  1.0_sp  ! - one -
  REAL(KIND=sp), PARAMETER, PUBLIC ::  ten_sp= 10.0_sp  ! - ten -
  REAL(KIND=sp), PARAMETER, PUBLIC :: r360_sp=360.0_sp  ! - 360 degrees -
!
  REAL(KIND=dp), PARAMETER, PUBLIC :: zero_dp=  0.0_dp  ! - zero -
  REAL(KIND=dp), PARAMETER, PUBLIC ::  one_dp=  1.0_dp  ! - one -
  REAL(KIND=dp), PARAMETER, PUBLIC ::  ten_dp= 10.0_dp  ! - ten -
  REAL(KIND=dp), PARAMETER, PUBLIC :: r360_dp=360.0_dp  ! - 360 degrees -
!
END MODULE numbers
!
!
!
MODULE CPT_formatV11
!
! Modules
  USE numbers
!
! Implicit declarations
  IMPLICIT NONE
!
! Parameters
!
! Integer parameters
! - public parameters -
  INTEGER, PARAMETER, PUBLIC :: lvar=32 ! - maximum length of variable name -
  INTEGER, PARAMETER, PUBLIC :: lstn=16 ! - maximum length of station name -
!
! - private parameters -
  INTEGER, PARAMETER, PRIVATE :: ioutd=21  ! - default output unit number -
  INTEGER, PARAMETER, PRIVATE ::   mwu=29  ! - maximum width of output field -
  INTEGER, PARAMETER, PRIVATE ::   nmn=12  ! - number of months -
  INTEGER, PARAMETER, PRIVATE ::  lfil=800 ! - maximum length of file with full path -
  INTEGER, PARAMETER, PRIVATE ::  ldat=12  ! - maximum length of date -
  INTEGER, PARAMETER, PRIVATE ::  lprd=25  ! - maximum length of period -
!
! Character parameters
  CHARACTER(LEN=32), PARAMETER, PRIVATE :: cxmlns_cpt= & ! - CPT XML namespace
     'http://iri.columbia.edu/CPT/v10/'
!
! Scalars
!
! Integer scalars
  INTEGER, PRIVATE :: iseq ! - time sequence identifier -
  INTEGER, PRIVATE :: iout ! - output unit number -
!
! Character scalars
  CHARACTER(LEN=mwu), PRIVATE :: cout ! - output field -
!
! Logical scalars
  LOGICAL, PRIVATE :: lopen ! - file opened flag -
!
! Derived type definitions
!
! - date -
  TYPE date
     INTEGER :: iyr ! - year -
     INTEGER :: imn ! - month -
     INTEGER :: idy ! - day -
  END TYPE date
!
! - period -
  TYPE period
     TYPE(date) :: sdate ! - start date -
     TYPE(date) :: edate ! - end date -
  END TYPE period
!
! - level -
  TYPE level
     REAL(KIND=dp) :: hght    ! - height -
!
     CHARACTER(LEN=5) :: unit ! - units -
  END TYPE level
!
! Interfaces
!
! Interface operators
  INTERFACE OPERATOR(==)
     MODULE PROCEDURE same_date
     MODULE PROCEDURE equal_date
  END INTERFACE
!
  INTERFACE OPERATOR(<)
     MODULE PROCEDURE lt_date
  END INTERFACE
!
  INTERFACE OPERATOR(>)
     MODULE PROCEDURE gt_date
  END INTERFACE
!
  INTERFACE OPERATOR(+)
     MODULE PROCEDURE add_date
  END INTERFACE
!
  INTERFACE OPERATOR(+)
     MODULE PROCEDURE add_period
  END INTERFACE
!
! Generic interfaces
  INTERFACE write_cpt_grid_v11
    MODULE PROCEDURE write_cpt_grid_v11_sp
    MODULE PROCEDURE write_cpt_grid_v11_dp
    MODULE PROCEDURE write_cpt_grid_lags_v11_sp
    MODULE PROCEDURE write_cpt_grid_lags_v11_dp
    MODULE PROCEDURE write_cpt_grid_fields_v11_sp
    MODULE PROCEDURE write_cpt_grid_fields_v11_dp
    MODULE PROCEDURE write_cpt_grid_fields_lags_v11_sp
    MODULE PROCEDURE write_cpt_grid_fields_lags_v11_dp
  END INTERFACE
!
  INTERFACE write_cpt_stns_v11
    MODULE PROCEDURE write_cpt_stns_v11_sp
    MODULE PROCEDURE write_cpt_stns_v11_dp
    MODULE PROCEDURE write_cpt_stns_lags_v11_sp
    MODULE PROCEDURE write_cpt_stns_lags_v11_dp
    MODULE PROCEDURE write_cpt_stns_fields_v11_sp
    MODULE PROCEDURE write_cpt_stns_fields_v11_dp
    MODULE PROCEDURE write_cpt_stns_fields_lags_v11_sp
    MODULE PROCEDURE write_cpt_stns_fields_lags_v11_dp
  END INTERFACE
!
  INTERFACE write_cpt_unrf_v11
    MODULE PROCEDURE write_cpt_unrf_v11_sp
    MODULE PROCEDURE write_cpt_unrf_v11_dp
    MODULE PROCEDURE write_cpt_unrf_lags_v11_sp
    MODULE PROCEDURE write_cpt_unrf_lags_v11_dp
  END INTERFACE
!
  INTERFACE magnitude
    MODULE PROCEDURE magnitude_int
    MODULE PROCEDURE magnitude_sp
    MODULE PROCEDURE magnitude_dp
  END INTERFACE magnitude
!
  INTERFACE iprec
    MODULE PROCEDURE iprec_sp
    MODULE PROCEDURE iprec_dp
  END INTERFACE iprec
!
  INTERFACE get_cdate
    MODULE PROCEDURE get_cdate_date
    MODULE PROCEDURE get_cdate_period
  END INTERFACE get_cdate
!
CONTAINS
!
!
 SUBROUTINE write_cpt_grid_v11_sp ( outfile,nt,nlt,nlg,tseq,v,miss,rlat,rlng,var,unit,&
                                    period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                    period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                    ifail,&
                                    mdate_iyr,mdate_imn,mdate_idy,&
                                    z_hght,z_unit,&
                                    kuse&
                                  )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp,r360=>r360_sp
!
! Outputs gridded data
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr ! - start year of first period -
  INTEGER, INTENT(IN) :: period1_sdate_imn ! - start month of first period -
  INTEGER, INTENT(IN) :: period1_sdate_idy ! - start day of first period -
  INTEGER, INTENT(IN) :: period1_edate_iyr ! - end year of first period -
  INTEGER, INTENT(IN) :: period1_edate_imn ! - end month of first period -
  INTEGER, INTENT(IN) :: period1_edate_idy ! - end day of first period -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! - optional input scalars -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr ! - year made ('start date' for model forecasts) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn ! - month made ('start date' for model forecasts) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy ! - day made ('start date' for model forecasts) -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght    ! - atmoshperic level -- height -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit ! - atmoshperic level -- units -
!
! Input arrays
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:) ! - data (minimum dimensions: nlg, nlt, nn) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)  ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)  ! - longitudes (minimum dimensions: nlg) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  INTEGER :: i  ! - latitude index -
  INTEGER :: j  ! - longitude index -
  INTEGER :: k  ! - time index -
  INTEGER :: nn ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
  TYPE(date) :: mdate1 ! - first forecast date -
!
  TYPE(level) ::  z ! - level -
!
  TYPE(period) :: period1 ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  period1%sdate%iyr=period1_sdate_iyr
  period1%sdate%imn=period1_sdate_imn
  period1%sdate%idy=period1_sdate_idy
  period1%edate%iyr=period1_edate_iyr
  period1%edate%imn=period1_edate_imn
  period1%edate%idy=period1_edate_idy
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     mdate1%iyr=mdate_iyr
     mdate1%imn=mdate_imn
     mdate1%idy=mdate_idy
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     z%hght=z_hght
     z%unit=z_unit
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  nn=0
  DO k=1,nt
     IF (PRESENT(kuse)) THEN
        IF (.NOT.kuse(k)) CYCLE
     END IF
     nn=nn+1
     IF (lmdate) THEN
        IF (lz) THEN
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_z=z,cpt_s=mdate1+(k-1),cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
        ELSE
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_s=mdate1+(k-1),cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
        END IF
     ELSE
        IF (lz) THEN
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_z=z,cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
        ELSE
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
        END IF
     END IF
     IF (ifail/=0) GOTO 1
     DO j=1,nlg
        IF (rlng(j)>r360) THEN
           WRITE (cout,FMT=*) rlng(j)-r360
        ELSE
           WRITE (cout,FMT=*) rlng(j)
        END IF
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
     DO i=1,nlt
        WRITE (cout,FMT=*) rlat(i)
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
        DO j=1,nlg
           WRITE (cout,FMT=*) v(j,i,nn)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_v11_sp
!
!
!
 SUBROUTINE write_cpt_grid_v11_dp ( outfile,nt,nlt,nlg,tseq,v,miss,rlat,rlng,var,unit,&
                                    period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                    period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                    ifail,&
                                    mdate_iyr,mdate_imn,mdate_idy,&
                                    z_hght,z_unit,&
                                    kuse&
                                  )
!
! Modules
  USE numbers, ONLY: rp=>dp,r360=>r360_dp
!
! Outputs gridded data
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr ! - start year of first period -
  INTEGER, INTENT(IN) :: period1_sdate_imn ! - start month of first period -
  INTEGER, INTENT(IN) :: period1_sdate_idy ! - start day of first period -
  INTEGER, INTENT(IN) :: period1_edate_iyr ! - end year of first period -
  INTEGER, INTENT(IN) :: period1_edate_imn ! - end month of first period -
  INTEGER, INTENT(IN) :: period1_edate_idy ! - end day of first period -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! - optional input scalars -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr ! - year made ('start date' for model forecasts) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn ! - month made ('start date' for model forecasts) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy ! - day made ('start date' for model forecasts) -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght    ! - atmoshperic level -- height -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit ! - atmoshperic level -- units -
!
! Input arrays
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:) ! - data (minimum dimensions: nlg, nlt, nn) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)  ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)  ! - longitudes (minimum dimensions: nlg) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  INTEGER :: i  ! - latitude index -
  INTEGER :: j  ! - longitude index -
  INTEGER :: k  ! - time index -
  INTEGER :: nn ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
  TYPE(date) :: mdate1 ! - first forecast date -
!
  TYPE(level) ::  z ! - level -
!
  TYPE(period) :: period1 ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  period1%sdate%iyr=period1_sdate_iyr
  period1%sdate%imn=period1_sdate_imn
  period1%sdate%idy=period1_sdate_idy
  period1%edate%iyr=period1_edate_iyr
  period1%edate%imn=period1_edate_imn
  period1%edate%idy=period1_edate_idy
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     mdate1%iyr=mdate_iyr
     mdate1%imn=mdate_imn
     mdate1%idy=mdate_idy
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     z%hght=z_hght
     z%unit=z_unit
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  nn=0
  DO k=1,nt
     IF (PRESENT(kuse)) THEN
        IF (.NOT.kuse(k)) CYCLE
     END IF
     nn=nn+1
     IF (lmdate) THEN
        IF (lz) THEN
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_z=z,cpt_s=mdate1+(k-1),cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=miss)
        ELSE
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_s=mdate1+(k-1),cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=miss)
        END IF
     ELSE
        IF (lz) THEN
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_z=z,cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=miss)
        ELSE
           CALL write_tag (iout,ifail, &
                           cpt_field=TRIM(var),cpt_t=period1+(k-1), &
                           cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                           cpt_units=TRIM(unit),cpt_missing=miss)
        END IF
     END IF
     IF (ifail/=0) GOTO 1
     DO j=1,nlg
        IF (rlng(j)>r360) THEN
           WRITE (cout,FMT=*) rlng(j)-r360
        ELSE
           WRITE (cout,FMT=*) rlng(j)
        END IF
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
     DO i=1,nlt
        WRITE (cout,FMT=*) rlat(i)
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
        DO j=1,nlg
           WRITE (cout,FMT=*) v(j,i,nn)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_v11_dp
!
!
!
 SUBROUTINE write_cpt_grid_lags_v11_sp ( outfile,nt,nlt,nlg,nls,tseq,v,miss,rlat,rlng,var,unit,&
                                         period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                         period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                         ifail,&
                                         mdate_iyr,mdate_imn,mdate_idy,&
                                         z_hght,z_unit,&
                                         kuse&
                                       )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp,r360=>r360_sp
!
! Outputs gridded data with additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nls  ! - number of lags -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! - optional input scalars -
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght    ! - atmoshperic level -- height -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit ! - atmoshperic level -- units -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each lag (minimum dimensions: nls) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:) ! - data (minimum dimensions: nlg, nlt, nn, nls) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)    ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)    ! - longitudes (minimum dimensions: nlg) -
!
! - optional input arrays -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr(:) ! - year made ('start date' for model forecasts) (minimum dimensions: nls) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn(:) ! - month made ('start date' for model forecasts) (minimum dimensions: nls) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy(:) ! - day made ('start date' for model forecasts) (minimum dimensions: nls) -
!
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(date) :: mdate1(nls) ! - first forecast date -
!
  TYPE(period) :: period1(nls) ! - first period -
!
! Local scalars
  INTEGER :: i  ! - latitude index -
  INTEGER :: j  ! - longitude index -
  INTEGER :: k  ! - time index -
  INTEGER :: l  ! - lag index -
  INTEGER :: nn ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
  TYPE(level) ::  z ! - level -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO l=1,nls
     period1(l)%sdate%iyr=period1_sdate_iyr(l)
     period1(l)%sdate%imn=period1_sdate_imn(l)
     period1(l)%sdate%idy=period1_sdate_idy(l)
     period1(l)%edate%iyr=period1_edate_iyr(l)
     period1(l)%edate%imn=period1_edate_imn(l)
     period1(l)%edate%idy=period1_edate_idy(l)
  END DO
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     DO l=1,nls
        mdate1(l)%iyr=mdate_iyr(l)
        mdate1(l)%imn=mdate_imn(l)
        mdate1(l)%idy=mdate_idy(l)
     END DO
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     z%hght=z_hght
     z%unit=z_unit
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  nn=0
  DO k=1,nt
     DO l=1,nls
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        IF (lmdate) THEN
           IF (lz) THEN
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var),cpt_z=z, &
                              cpt_s=mdate1(l)+(k-1),cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
           ELSE
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var), &
                              cpt_s=mdate1(l)+(k-1),cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
           END IF
        ELSE
           IF (lz) THEN
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var),cpt_z=z, &
                              cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
           ELSE
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var), &
                              cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
           END IF
        END IF
        IF (ifail/=0) GOTO 1
        DO j=1,nlg
           IF (rlng(j)>r360) THEN
              WRITE (cout,FMT=*) rlng(j)-r360
           ELSE
              WRITE (cout,FMT=*) rlng(j)
           END IF
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        DO i=1,nlt
           WRITE (cout,FMT=*) rlat(i)
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
           DO j=1,nlg
              WRITE (cout,FMT=*) v(j,i,nn,l)
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_lags_v11_sp
!
!
!
 SUBROUTINE write_cpt_grid_lags_v11_dp ( outfile,nt,nlt,nlg,nls,tseq,v,miss,rlat,rlng,var,unit,&
                                         period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                         period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                         ifail,&
                                         mdate_iyr,mdate_imn,mdate_idy,&
                                         z_hght,z_unit,&
                                         kuse&
                                       )
!
! Modules
  USE numbers, ONLY: rp=>dp,r360=>r360_dp
!
! Outputs gridded data with additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nls  ! - number of lags -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! - optional input scalars -
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght    ! - atmoshperic level -- height -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit ! - atmoshperic level -- units -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each lag (minimum dimensions: nls) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:) ! - data (minimum dimensions: nlg, nlt, nn, nls) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)    ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)    ! - longitudes (minimum dimensions: nlg) -
!
! - optional input arrays -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr(:) ! - year made ('start date' for model forecasts) (minimum dimensions: nls) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn(:) ! - month made ('start date' for model forecasts) (minimum dimensions: nls) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy(:) ! - day made ('start date' for model forecasts) (minimum dimensions: nls) -
!
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(date) :: mdate1(nls) ! - first forecast date -
!
  TYPE(period) :: period1(nls) ! - first period -
!
! Local scalars
  INTEGER :: i  ! - latitude index -
  INTEGER :: j  ! - longitude index -
  INTEGER :: k  ! - time index -
  INTEGER :: l  ! - lag index -
  INTEGER :: nn ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
  TYPE(level) ::  z ! - level -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO l=1,nls
     period1(l)%sdate%iyr=period1_sdate_iyr(l)
     period1(l)%sdate%imn=period1_sdate_imn(l)
     period1(l)%sdate%idy=period1_sdate_idy(l)
     period1(l)%edate%iyr=period1_edate_iyr(l)
     period1(l)%edate%imn=period1_edate_imn(l)
     period1(l)%edate%idy=period1_edate_idy(l)
  END DO
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     DO l=1,nls
        mdate1(l)%iyr=mdate_iyr(l)
        mdate1(l)%imn=mdate_imn(l)
        mdate1(l)%idy=mdate_idy(l)
     END DO
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     z%hght=z_hght
     z%unit=z_unit
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  nn=0
  DO k=1,nt
     DO l=1,nls
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        IF (lmdate) THEN
           IF (lz) THEN
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var),cpt_z=z, &
                              cpt_s=mdate1(l)+(k-1),cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=miss)
           ELSE
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var), &
                              cpt_s=mdate1(l)+(k-1),cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=miss)
           END IF
        ELSE
           IF (lz) THEN
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var), &
                              cpt_z=z,cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=miss)
           ELSE
              CALL write_tag (iout,ifail, &
                              cpt_field=TRIM(var), &
                              cpt_t=period1(l)+(k-1), &
                              cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                              cpt_units=TRIM(unit),cpt_missing=miss)
           END IF
        END IF
        IF (ifail/=0) GOTO 1
        DO j=1,nlg
           IF (rlng(j)>r360) THEN
              WRITE (cout,FMT=*) rlng(j)-r360
           ELSE
              WRITE (cout,FMT=*) rlng(j)
           END IF
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        DO i=1,nlt
           WRITE (cout,FMT=*) rlat(i)
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
           DO j=1,nlg
              WRITE (cout,FMT=*) v(j,i,nn,l)
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_lags_v11_dp
!
!
!
 SUBROUTINE write_cpt_grid_fields_v11_sp ( outfile,nt,nlt,nlg,nfs,tseq,v,miss,rlat,rlng,var,unit,&
                                           period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                           period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                           lensemble,ifail,&
                                           mdate_iyr,mdate_imn,mdate_idy,&
                                           z_hght,z_unit,&
                                           kuse&
                                         )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp,r360=>r360_sp
!
! Outputs gridded data with multiple fields
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
  LOGICAL, INTENT(IN) :: lensemble ! - set to .true. if fields represent ensemble members -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each field (minimum dimensions: nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:) ! - data (minimum dimensions: nlg, nlt, nn, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)    ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)    ! - longitudes (minimum dimensions: nlg) -
!
  CHARACTER(LEN=*), INTENT(IN) :: var(:)  ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:) ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr(:) ! - year made ('start date' for model forecasts) (minimum dimensions: nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn(:) ! - month made ('start date' for model forecasts) (minimum dimensions: nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy(:) ! - day made ('start date' for model forecasts) (minimum dimensions: nfs) -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght(:) ! - atmoshperic level -- height (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit(:) ! - atmoshperic level -- units (minimum dimensions: nfs) -
!
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(date) :: mdate1(nfs) ! - first forecast date -
!
  TYPE(period) :: period1(nfs) ! - first period -
!
  TYPE(level) ::  z(nfs) ! - level -
!
! Local scalars
  INTEGER :: i   ! - latitude index -
  INTEGER :: j   ! - longitude index -
  INTEGER :: k   ! - time index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     period1(ifd)%sdate%iyr=period1_sdate_iyr(ifd)
     period1(ifd)%sdate%imn=period1_sdate_imn(ifd)
     period1(ifd)%sdate%idy=period1_sdate_idy(ifd)
     period1(ifd)%edate%iyr=period1_edate_iyr(ifd)
     period1(ifd)%edate%imn=period1_edate_imn(ifd)
     period1(ifd)%edate%idy=period1_edate_idy(ifd)
  END DO
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     DO ifd=1,nfs
        mdate1(ifd)%iyr=mdate_iyr(ifd)
        mdate1(ifd)%imn=mdate_imn(ifd)
        mdate1(ifd)%idy=mdate_idy(ifd)
     END DO
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     DO ifd=1,nfs
        z(ifd)%hght=z_hght(ifd)
        z(ifd)%unit=z_unit(ifd)
     END DO
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  DO ifd=1,nfs
     nn=0
     DO k=1,nt
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        IF (lmdate) THEN
           IF (lz) THEN
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           ELSE
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)), &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           END IF
        ELSE
           IF (lz) THEN
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           ELSE
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)), &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           END IF
        END IF
        IF (ifail/=0) GOTO 1
        DO j=1,nlg
           IF (rlng(j)>r360) THEN
              WRITE (cout,FMT=*) rlng(j)-r360
           ELSE
              WRITE (cout,FMT=*) rlng(j)
           END IF
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        DO i=1,nlt
           WRITE (cout,FMT=*) rlat(i)
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
           DO j=1,nlg
              WRITE (cout,FMT=*) v(j,i,nn,ifd)
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_fields_v11_sp
!
!
!
 SUBROUTINE write_cpt_grid_fields_v11_dp ( outfile,nt,nlt,nlg,nfs,tseq,v,miss,rlat,rlng,var,unit,&
                                           period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                           period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                           lensemble,ifail,&
                                           mdate_iyr,mdate_imn,mdate_idy,&
                                           z_hght,z_unit,&
                                           kuse&
                                         )
!
! Modules
  USE numbers, ONLY: rp=>dp,r360=>r360_dp
!
! Outputs gridded data with multiple fields
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
  LOGICAL, INTENT(IN) :: lensemble ! - set to .true. if fields represent ensemble members -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each field (minimum dimensions: nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:) ! - data (minimum dimensions: nlg, nlt, nn, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)    ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)    ! - longitudes (minimum dimensions: nlg) -
!
  CHARACTER(LEN=*), INTENT(IN) :: var(:)  ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:) ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr(:) ! - year made ('start date' for model forecasts) (minimum dimensions: nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn(:) ! - month made ('start date' for model forecasts) (minimum dimensions: nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy(:) ! - day made ('start date' for model forecasts) (minimum dimensions: nfs) -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght(:) ! - atmoshperic level -- height (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit(:) ! - atmoshperic level -- units (minimum dimensions: nfs) -
!
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(date) :: mdate1(nfs) ! - first forecast date -
!
  TYPE(period) :: period1(nfs) ! - first period -
!
  TYPE(level) ::  z(nfs) ! - level -
!
! Local scalars
  INTEGER :: i   ! - latitude index -
  INTEGER :: j   ! - longitude index -
  INTEGER :: k   ! - time index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     period1(ifd)%sdate%iyr=period1_sdate_iyr(ifd)
     period1(ifd)%sdate%imn=period1_sdate_imn(ifd)
     period1(ifd)%sdate%idy=period1_sdate_idy(ifd)
     period1(ifd)%edate%iyr=period1_edate_iyr(ifd)
     period1(ifd)%edate%imn=period1_edate_imn(ifd)
     period1(ifd)%edate%idy=period1_edate_idy(ifd)
  END DO
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     DO ifd=1,nfs
        mdate1(ifd)%iyr=mdate_iyr(ifd)
        mdate1(ifd)%imn=mdate_imn(ifd)
        mdate1(ifd)%idy=mdate_idy(ifd)
     END DO
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     DO ifd=1,nfs
        z(ifd)%hght=z_hght(ifd)
        z(ifd)%unit=z_unit(ifd)
     END DO
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  DO ifd=1,nfs
     nn=0
     DO k=1,nt
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        IF (lmdate) THEN
           IF (lz) THEN
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           ELSE
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)), &
                                 cpt_s=mdate1(ifd)+(k-1),cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           END IF
        ELSE
           IF (lz) THEN
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           ELSE
              IF (lensemble) THEN
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              ELSE
                 CALL write_tag (iout,ifail, &
                                 cpt_field=TRIM(var(ifd)), &
                                 cpt_t=period1(ifd)+(k-1), &
                                 cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                 cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
              END IF
           END IF
        END IF
        IF (ifail/=0) GOTO 1
        DO j=1,nlg
           IF (rlng(j)>r360) THEN
              WRITE (cout,FMT=*) rlng(j)-r360
           ELSE
              WRITE (cout,FMT=*) rlng(j)
           END IF
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        DO i=1,nlt
           WRITE (cout,FMT=*) rlat(i)
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
           DO j=1,nlg
              WRITE (cout,FMT=*) v(j,i,nn,ifd)
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_fields_v11_dp
!
!
!
 SUBROUTINE write_cpt_grid_fields_lags_v11_sp ( outfile,nt,nlt,nlg,nls,nfs,tseq,v,miss,rlat,rlng,var,unit,&
                                                period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                                period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                                lensemble,ifail,&
                                                mdate_iyr,mdate_imn,mdate_idy,&
                                                z_hght,z_unit,&
                                                kuse&
                                              )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp,r360=>r360_sp
!
! Outputs gridded data with multiple fields and additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nls  ! - number of lags -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
!
  LOGICAL, INTENT(IN) :: lensemble ! - set to .true. if fields represent ensemble members -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:,:) ! - start year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:,:) ! - start month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:,:) ! - start day of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:,:) ! - end year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:,:) ! - end month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:,:) ! - end day of first period for each field (minimum dimensions: nls, nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:,:) ! - data (minimum dimensions: nlg, nlt, nn, nls, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)      ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)      ! - longitudes (minimum dimensions: nlg) -
!
  CHARACTER(LEN=*), INTENT(IN) :: var(:)  ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:) ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr(:,:) ! - year made ('start date' for model forecasts) (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn(:,:) ! - month made ('start date' for model forecasts) (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy(:,:) ! - day made ('start date' for model forecasts) (minimum dimensions: nls, nfs) -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght(:) ! - atmoshperic level -- height (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit(:) ! - atmoshperic level -- units (minimum dimensions: nfs) -
!
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(date) :: mdate1(nls,nfs) ! - first forecast date -
!
  TYPE(period) :: period1(nls,nfs) ! - first period -
!
  TYPE(level) ::  z(nfs) ! - level -
!
! Local scalars
  INTEGER :: i   ! - latitude index -
  INTEGER :: j   ! - longitude index -
  INTEGER :: k   ! - time index -
  INTEGER :: l   ! - lag index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     DO l=1,nls
        period1(l,ifd)%sdate%iyr=period1_sdate_iyr(l,ifd)
        period1(l,ifd)%sdate%imn=period1_sdate_imn(l,ifd)
        period1(l,ifd)%sdate%idy=period1_sdate_idy(l,ifd)
        period1(l,ifd)%edate%iyr=period1_edate_iyr(l,ifd)
        period1(l,ifd)%edate%imn=period1_edate_imn(l,ifd)
        period1(l,ifd)%edate%idy=period1_edate_idy(l,ifd)
     END DO
  END DO
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     DO ifd=1,nfs
        DO l=1,nls
           mdate1(l,ifd)%iyr=mdate_iyr(l,ifd)
           mdate1(l,ifd)%imn=mdate_imn(l,ifd)
           mdate1(l,ifd)%idy=mdate_idy(l,ifd)
        END DO
     END DO
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     DO ifd=1,nfs
        z(ifd)%hght=z_hght(ifd)
        z(ifd)%unit=z_unit(ifd)
     END DO
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  DO ifd=1,nfs
     nn=0
     DO k=1,nt
        DO l=1,nls
           IF (PRESENT(kuse)) THEN
              IF (.NOT.kuse(k)) CYCLE
           END IF
           nn=nn+1
           IF (lmdate) THEN
              IF (lz) THEN
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              ELSE
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)), &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              END IF
           ELSE
              IF (lz) THEN
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              ELSE
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)), &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              END IF
           END IF
           IF (ifail/=0) GOTO 1
           DO j=1,nlg
              IF (rlng(j)>r360) THEN
                 WRITE (cout,FMT=*) rlng(j)-r360
              ELSE
                 WRITE (cout,FMT=*) rlng(j)
              END IF
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
           DO i=1,nlt
              WRITE (cout,FMT=*) rlat(i)
              WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
              DO j=1,nlg
                 WRITE (cout,FMT=*) v(j,i,nn,l,ifd)
                 WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
              END DO
              WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
           END DO
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_fields_lags_v11_sp
!
!
!
 SUBROUTINE write_cpt_grid_fields_lags_v11_dp ( outfile,nt,nlt,nlg,nls,nfs,tseq,v,miss,rlat,rlng,var,unit,&
                                                period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                                period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                                lensemble,ifail,&
                                                mdate_iyr,mdate_imn,mdate_idy,&
                                                z_hght,z_unit,&
                                                kuse&
                                              )
!
! Modules
  USE numbers, ONLY: rp=>dp,r360=>r360_dp
!
! Outputs gridded data with multiple fields and additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nlt  ! - total number of latitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nlg  ! - total number of longitudes of each field and lead-time -
  INTEGER, INTENT(IN) :: nls  ! - number of lags -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
  LOGICAL, INTENT(IN) :: lensemble ! - set to .true. if fields represent ensemble members -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:,:) ! - start year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:,:) ! - start month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:,:) ! - start day of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:,:) ! - end year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:,:) ! - end month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:,:) ! - end day of first period for each field (minimum dimensions: nls, nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:,:) ! - data (minimum dimensions: nlg, nlt, nn, nls, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)      ! - latitudes (minimum dimensions: nlt) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)      ! - longitudes (minimum dimensions: nlg) -
!
  CHARACTER(LEN=*), INTENT(IN) :: var(:)  ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:) ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_iyr(:,:) ! - year made ('start date' for model forecasts) (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_imn(:,:) ! - month made ('start date' for model forecasts) (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN), OPTIONAL :: mdate_idy(:,:) ! - day made ('start date' for model forecasts) (minimum dimensions: nls, nfs) -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: z_hght(:) ! - atmoshperic level -- height (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: z_unit(:) ! - atmoshperic level -- units (minimum dimensions: nfs) -
!
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(date) :: mdate1(nls,nfs) ! - first forecast date -
!
  TYPE(period) :: period1(nls,nfs) ! - first period -
!
  TYPE(level) ::  z(nfs) ! - level -
!
! Local scalars
  INTEGER :: i   ! - latitude index -
  INTEGER :: j   ! - longitude index -
  INTEGER :: k   ! - time index -
  INTEGER :: l   ! - lag index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
  LOGICAL :: lmdate ! - forecast date flag -
  LOGICAL :: lz     ! - atmospheric level flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     DO l=1,nls
        period1(l,ifd)%sdate%iyr=period1_sdate_iyr(l,ifd)
        period1(l,ifd)%sdate%imn=period1_sdate_imn(l,ifd)
        period1(l,ifd)%sdate%idy=period1_sdate_idy(l,ifd)
        period1(l,ifd)%edate%iyr=period1_edate_iyr(l,ifd)
        period1(l,ifd)%edate%imn=period1_edate_imn(l,ifd)
        period1(l,ifd)%edate%idy=period1_edate_idy(l,ifd)
     END DO
  END DO
!
! Initialise optional field variables
  IF (PRESENT(mdate_iyr).AND.PRESENT(mdate_imn).AND.PRESENT(mdate_idy)) THEN
     DO ifd=1,nfs
        DO l=1,nls
           mdate1(l,ifd)%iyr=mdate_iyr(l,ifd)
           mdate1(l,ifd)%imn=mdate_imn(l,ifd)
           mdate1(l,ifd)%idy=mdate_idy(l,ifd)
        END DO
     END DO
     lmdate=.true.
  ELSE
     lmdate=.false.
  END IF
  IF(PRESENT(z_hght) .AND. PRESENT(z_unit) ) THEN
     DO ifd=1,nfs
        z(ifd)%hght=z_hght(ifd)
        z(ifd)%unit=z_unit(ifd)
     END DO
     lz=.true.
  ELSE
     lz=.false.
  END IF
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',outfile, 'failed. Exit.'
     CLOSE(iout)
     RETURN
  END IF
       
! Print gridded data
  DO ifd=1,nfs
     nn=0
     DO k=1,nt
        DO l=1,nls
           IF (PRESENT(kuse)) THEN
              IF (.NOT.kuse(k)) CYCLE
           END IF
           nn=nn+1
           IF (lmdate) THEN
              IF (lz) THEN
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              ELSE
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)), &
                                    cpt_s=mdate1(l,ifd)+(k-1),cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              END IF
           ELSE
              IF (lz) THEN
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd),cpt_m=ifd, &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_z=z(ifd), &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              ELSE
                 IF (lensemble) THEN
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)),cpt_m=ifd, &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 ELSE
                    CALL write_tag (iout,ifail, &
                                    cpt_field=TRIM(var(ifd)), &
                                    cpt_t=period1(l,ifd)+(k-1), &
                                    cpt_nrow=nlt,cpt_ncol=nlg,cpt_row='Y',cpt_col='X', &
                                    cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss,KIND=dp))
                 END IF
              END IF
           END IF
           IF (ifail/=0) GOTO 1
           DO j=1,nlg
              IF (rlng(j)>r360) THEN
                 WRITE (cout,FMT=*) rlng(j)-r360
              ELSE
                 WRITE (cout,FMT=*) rlng(j)
              END IF
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
           DO i=1,nlt
              WRITE (cout,FMT=*) rlat(i)
              WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) TRIM(ADJUSTL(cout))
              DO j=1,nlg
                 WRITE (cout,FMT=*) v(j,i,nn,l,ifd)
                 WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
              END DO
              WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1)
           END DO
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE(iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_grid_fields_lags_v11_dp
!
!
!
 SUBROUTINE write_cpt_stns_v11_sp ( outfile,nv,nt,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                    period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                    period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                    ifail,&
                                    kuse&
                                  )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp
!
! Outputs station data with no additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nv   ! - number of stations -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr ! - start year of first period -
  INTEGER, INTENT(IN) :: period1_sdate_imn ! - start month of first period -
  INTEGER, INTENT(IN) :: period1_sdate_idy ! - start day of first period -
  INTEGER, INTENT(IN) :: period1_edate_iyr ! - end year of first period -
  INTEGER, INTENT(IN) :: period1_edate_imn ! - end month of first period -
  INTEGER, INTENT(IN) :: period1_edate_idy ! - end day of first period -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! Input arrays
  REAL(KIND=rp), INTENT(IN) :: v(:,:)  ! - data (minimum dimensions: nv, nn) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:) ! - latitudes (minimum dimensions: nv) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:) ! - longitudes (minimum dimensions: nv) -
!
  CHARACTER(LEN=lstn), INTENT(IN) :: cstn(:) ! - station names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  INTEGER :: i  ! - station index -
  INTEGER :: k  ! - time index -
  INTEGER :: nn ! - number of non-missing cases -
!
  TYPE(period) :: period1 ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  period1%sdate%iyr=period1_sdate_iyr
  period1%sdate%imn=period1_sdate_imn
  period1%sdate%idy=period1_sdate_idy
  period1%edate%iyr=period1_edate_iyr
  period1%edate%imn=period1_edate_imn
  period1%edate%idy=period1_edate_idy
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                                &
                  cpt_field=TRIM(var),cpt_nrow=nn,cpt_ncol=nv,cpt_row='T',cpt_col='station', &
                  cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
  IF (ifail/=0) GOTO 1
!
! Print station information
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
  DO i=1,nv
     WRITE (cout,FMT=*) rlat(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
  DO i=1,nv
     WRITE (cout,FMT=*) rlng(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
  nn=0
  DO k=1,nt
     IF (PRESENT(kuse)) THEN
        IF (.NOT.kuse(k)) CYCLE
     END IF
     nn=nn+1
     cout=get_cdate(period1+(k-1))
     WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
     DO i=1,nv
        WRITE (cout,FMT=*) v(i,nn)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_v11_sp
!
!
!
 SUBROUTINE write_cpt_stns_v11_dp ( outfile,nv,nt,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                    period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                    period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                    ifail,&
                                    kuse&
                                  )
!
! Modules
  USE numbers, ONLY: rp=>dp
!
! Outputs station data with no additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nv   ! - number of stations -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr ! - start year of first period -
  INTEGER, INTENT(IN) :: period1_sdate_imn ! - start month of first period -
  INTEGER, INTENT(IN) :: period1_sdate_idy ! - start day of first period -
  INTEGER, INTENT(IN) :: period1_edate_iyr ! - end year of first period -
  INTEGER, INTENT(IN) :: period1_edate_imn ! - end month of first period -
  INTEGER, INTENT(IN) :: period1_edate_idy ! - end day of first period -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! Input arrays
  REAL(KIND=rp), INTENT(IN) :: v(:,:)  ! - data (minimum dimensions: nv, nn) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:) ! - latitudes (minimum dimensions: nv) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:) ! - longitudes (minimum dimensions: nv) -
!
  CHARACTER(LEN=lstn), INTENT(IN) :: cstn(:) ! - station names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  INTEGER :: i  ! - station index -
  INTEGER :: k  ! - time index -
  INTEGER :: nn ! - number of non-missing cases -
!
  TYPE(period) :: period1 ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  period1%sdate%iyr=period1_sdate_iyr
  period1%sdate%imn=period1_sdate_imn
  period1%sdate%idy=period1_sdate_idy
  period1%edate%iyr=period1_edate_iyr
  period1%edate%imn=period1_edate_imn
  period1%edate%idy=period1_edate_idy
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                                &
                  cpt_field=TRIM(var),cpt_nrow=nn,cpt_ncol=nv,cpt_row='T',cpt_col='station', &
                  cpt_units=TRIM(unit),cpt_missing=miss)
  IF (ifail/=0) GOTO 1
!
! Print station information
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
  DO i=1,nv
     WRITE (cout,FMT=*) rlat(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
  DO i=1,nv
     WRITE (cout,FMT=*) rlng(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
  nn=0
  DO k=1,nt
     IF (PRESENT(kuse)) THEN
        IF (.NOT.kuse(k)) CYCLE
     END IF
     nn=nn+1
     cout=get_cdate(period1+(k-1))
     WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
     DO i=1,nv
        WRITE (cout,FMT=*) v(i,nn)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_v11_dp
!
!
!
 SUBROUTINE write_cpt_stns_lags_v11_sp ( outfile,nv,nt,nls,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                         period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                         period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                         ifail,&
                                         kuse&
                                       )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp
!
! Outputs station data with additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nv   ! - number of stations -
  INTEGER, INTENT(IN) :: nls  ! - number of lead-times -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each lag (minimum dimensions: nls) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:) ! - data (minimum dimensions: nv, nn, nls) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)  ! - latitudes (minimum dimensions: nv) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)  ! - longitudes (minimum dimensions: nv) -
!
  CHARACTER(LEN=lstn), INTENT(IN) :: cstn(:) ! - station names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nls) ! - first period -
!
! Local scalars
  INTEGER :: i  ! - station index -
  INTEGER :: k  ! - time index -
  INTEGER :: l  ! - lag index -
  INTEGER :: nn ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO l=1,nls
     period1(l)%sdate%iyr=period1_sdate_iyr(l)
     period1(l)%sdate%imn=period1_sdate_imn(l)
     period1(l)%sdate%idy=period1_sdate_idy(l)
     period1(l)%edate%iyr=period1_edate_iyr(l)
     period1(l)%edate%imn=period1_edate_imn(l)
     period1(l)%edate%idy=period1_edate_idy(l)
  END DO
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                                    &
                  cpt_field=TRIM(var),cpt_nrow=nn*nls,cpt_ncol=nv,cpt_row='T',cpt_col='station', &
                  cpt_units=TRIM(unit),cpt_missing=REAL(miss,KIND=dp))
  IF (ifail/=0) GOTO 1
!
! Print station information
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
  DO i=1,nv
     WRITE (cout,FMT=*) rlat(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
  DO i=1,nv
     WRITE (cout,FMT=*) rlng(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
  nn=0
  DO k=1,nt
     DO l=1,nls
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        cout=get_cdate(period1(l)+(k-1))
        WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
        DO i=1,nv
           WRITE (cout,FMT=*) v(i,nn,l)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_lags_v11_sp
!
!
!
 SUBROUTINE write_cpt_stns_lags_v11_dp ( outfile,nv,nt,nls,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                         period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                         period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                         ifail,&
                                         kuse&
                                       )
!
! Modules
  USE numbers, ONLY: rp=>dp
!
! Outputs station data with additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nv   ! - number of stations -
  INTEGER, INTENT(IN) :: nls  ! - number of lead-times -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
  CHARACTER(LEN=*), INTENT(IN) :: unit    ! - field units -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each lag (minimum dimensions: nls) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:) ! - data (minimum dimensions: nv, nn, nls) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:)  ! - latitudes (minimum dimensions: nv) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:)  ! - longitudes (minimum dimensions: nv) -
!
  CHARACTER(LEN=lstn), INTENT(IN) :: cstn(:) ! - station names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nls) ! - first period -
!
! Local scalars
  INTEGER :: i  ! - station index -
  INTEGER :: k  ! - time index -
  INTEGER :: l  ! - lag index -
  INTEGER :: nn ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO l=1,nls
     period1(l)%sdate%iyr=period1_sdate_iyr(l)
     period1(l)%sdate%imn=period1_sdate_imn(l)
     period1(l)%sdate%idy=period1_sdate_idy(l)
     period1(l)%edate%iyr=period1_edate_iyr(l)
     period1(l)%edate%imn=period1_edate_imn(l)
     period1(l)%edate%idy=period1_edate_idy(l)
  END DO
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                                    &
                  cpt_field=TRIM(var),cpt_nrow=nn*nls,cpt_ncol=nv,cpt_row='T',cpt_col='station', &
                  cpt_units=TRIM(unit),cpt_missing=miss)
  IF (ifail/=0) GOTO 1
!
! Print station information
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
  DO i=1,nv
     WRITE (cout,FMT=*) rlat(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
  DO i=1,nv
     WRITE (cout,FMT=*) rlng(i)
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
  nn=0
  DO k=1,nt
     DO l=1,nls
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        cout=get_cdate(period1(l)+(k-1))
        WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
        DO i=1,nv
           WRITE (cout,FMT=*) v(i,nn,l)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_lags_v11_dp
!
!
!
 SUBROUTINE write_cpt_stns_fields_v11_sp ( outfile,nv,nt,nfs,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                           period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                           period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                           ifail,&
                                           kuse&
                                         )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp
!
! Outputs station data with multiple fields, but no additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
! Input arrays
  INTEGER, INTENT(IN) :: nv(:) ! - number of stations per field (minimum dimensions: nfs) -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each field (minimum dimensions: nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:)  ! - data (minimum dimensions: nv, nn, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:,:) ! - latitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:,:) ! - longitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: miss(:)   ! - missing value flags (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN) :: cstn(:,:) ! - station names (minimum dimensions: nv, nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: var(:)    ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:)   ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nfs) ! - first period -
!
! Local scalars
  INTEGER :: i   ! - station index -
  INTEGER :: k   ! - time index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     period1(ifd)%sdate%iyr=period1_sdate_iyr(ifd)
     period1(ifd)%sdate%imn=period1_sdate_imn(ifd)
     period1(ifd)%sdate%idy=period1_sdate_idy(ifd)
     period1(ifd)%edate%iyr=period1_edate_iyr(ifd)
     period1(ifd)%edate%imn=period1_edate_imn(ifd)
     period1(ifd)%edate%idy=period1_edate_idy(ifd)
  END DO

! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  DO ifd=1,nfs
     CALL write_tag (iout,ifail,                                                                          &
                     cpt_field=TRIM(var(ifd)),cpt_nrow=nn,cpt_ncol=nv(ifd),cpt_row='T',cpt_col='station', &
                     cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss(ifd),KIND=dp))
     IF (ifail/=0) GOTO 1
!
! Print station information
     DO i=1,nv(ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i,ifd)))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlat(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlng(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
     nn=0
     DO k=1,nt
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        cout=get_cdate(period1(ifd)+(k-1))
        WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
        DO i=1,nv(ifd)
           WRITE (cout,FMT=*) v(i,nn,ifd)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_fields_v11_sp
!
!
!
 SUBROUTINE write_cpt_stns_fields_v11_dp ( outfile,nv,nt,nfs,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                           period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                           period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                           ifail,&
                                           kuse&
                                         )
!
! Modules
  USE numbers, ONLY: rp=>dp
!
! Outputs station data with multiple fields, but no additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
! Input arrays
  INTEGER, INTENT(IN) :: nv(:) ! - number of stations per field (minimum dimensions: nfs) -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period for each field (minimum dimensions: nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period for each field (minimum dimensions: nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:)  ! - data (minimum dimensions: nv, nn, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:,:) ! - latitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:,:) ! - longitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: miss(:)   ! - missing value flags (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN) :: cstn(:,:) ! - station names (minimum dimensions: nv, nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: var(:)    ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:)   ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nfs) ! - first period -
!
! Local scalars
  INTEGER :: i   ! - station index -
  INTEGER :: k   ! - time index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period

!
! Initialise sequencing and field information
  iseq=tseq
  DO ifd=1,nfs
     period1(ifd)%sdate%iyr=period1_sdate_iyr(ifd)
     period1(ifd)%sdate%imn=period1_sdate_imn(ifd)
     period1(ifd)%sdate%idy=period1_sdate_idy(ifd)
     period1(ifd)%edate%iyr=period1_edate_iyr(ifd)
     period1(ifd)%edate%imn=period1_edate_imn(ifd)
     period1(ifd)%edate%idy=period1_edate_idy(ifd)
  END DO

! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  DO ifd=1,nfs
     CALL write_tag (iout,ifail,                                                                          &
                     cpt_field=TRIM(var(ifd)),cpt_nrow=nn,cpt_ncol=nv(ifd),cpt_row='T',cpt_col='station', &
                     cpt_units=TRIM(unit(ifd)),cpt_missing=miss(ifd))
     IF (ifail/=0) GOTO 1
!
! Print station information
     DO i=1,nv(ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i,ifd)))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlat(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlng(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
     nn=0
     DO k=1,nt
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        cout=get_cdate(period1(ifd)+(k-1))
        WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
        DO i=1,nv(ifd)
           WRITE (cout,FMT=*) v(i,nn,ifd)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_fields_v11_dp
!
!
!
 SUBROUTINE write_cpt_stns_fields_lags_v11_sp ( outfile,nv,nt,nls,nfs,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                                period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                                period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                                ifail,&
                                                kuse&
                                              )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp
!
! Outputs station data with multiple fields, but no additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nls  ! - number of lags -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
! Input arrays
  INTEGER, INTENT(IN) :: nv(:) ! - number of stations per field (minimum dimensions: nfs) -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:,:) ! - start year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:,:) ! - start month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:,:) ! - start day of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:,:) ! - end year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:,:) ! - end month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:,:) ! - end day of first period for each field (minimum dimensions: nls, nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:) ! - data (minimum dimensions: nv, nn, nls, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:,:)  ! - latitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:,:)  ! - longitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: miss(:)    ! - missing value flags (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN) :: cstn(:,:) ! - station names (minimum dimensions: nv, nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: var(:)    ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:)   ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nls,nfs) ! - first period -
!
! Local scalars
  INTEGER :: i   ! - station index -
  INTEGER :: k   ! - time index -
  INTEGER :: l   ! - lag index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     DO l=1,nls
        period1(l,ifd)%sdate%iyr=period1_sdate_iyr(l,ifd)
        period1(l,ifd)%sdate%imn=period1_sdate_imn(l,ifd)
        period1(l,ifd)%sdate%idy=period1_sdate_idy(l,ifd)
        period1(l,ifd)%edate%iyr=period1_edate_iyr(l,ifd)
        period1(l,ifd)%edate%imn=period1_edate_imn(l,ifd)
        period1(l,ifd)%edate%idy=period1_edate_idy(l,ifd)
     END DO
  END DO

! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  DO ifd=1,nfs
     CALL write_tag (iout,ifail,                                                                              &
                     cpt_field=TRIM(var(ifd)),cpt_nrow=nn*nls,cpt_ncol=nv(ifd),cpt_row='T',cpt_col='station', &
                     cpt_units=TRIM(unit(ifd)),cpt_missing=REAL(miss(ifd),KIND=dp))
     IF (ifail/=0) GOTO 1
!
! Print station information
     DO i=1,nv(ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i,ifd)))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlat(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlng(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
     nn=0
     DO k=1,nt
        DO l=1,nls
           IF (PRESENT(kuse)) THEN
              IF (.NOT.kuse(k)) CYCLE
           END IF
           nn=nn+1
           cout=get_cdate(period1(l,ifd)+(k-1))
           WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
           DO i=1,nv(ifd)
              WRITE (cout,FMT=*) v(i,nn,l,ifd)
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_fields_lags_v11_sp
!
!
!
 SUBROUTINE write_cpt_stns_fields_lags_v11_dp ( outfile,nv,nt,nls,nfs,tseq,v,miss,rlat,rlng,cstn,var,unit,&
                                                period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                                period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                                ifail,&
                                                kuse&
                                              )
!
! Modules
  USE numbers, ONLY: rp=>dp
!
! Outputs station data with multiple fields, but no additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nls  ! - number of lags -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  CHARACTER(LEN=   *), INTENT(IN) :: outfile ! - file name with full path -
!
! Input arrays
  INTEGER, INTENT(IN) :: nv(:) ! - number of stations per field (minimum dimensions: nfs) -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:,:) ! - start year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:,:) ! - start month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:,:) ! - start day of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:,:) ! - end year of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:,:) ! - end month of first period for each field (minimum dimensions: nls, nfs) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:,:) ! - end day of first period for each field (minimum dimensions: nls, nfs) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:,:) ! - data (minimum dimensions: nv, nn, nls, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlat(:,:)  ! - latitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: rlng(:,:)  ! - longitudes (minimum dimensions: nv, nfs) -
  REAL(KIND=rp), INTENT(IN) :: miss(:)    ! - missing value flags (minimum dimensions: nfs) -
!
  CHARACTER(LEN=*), INTENT(IN) :: cstn(:,:) ! - station names (minimum dimensions: nv, nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: var(:)    ! - field variable (minimum dimensions: nfs) -
  CHARACTER(LEN=*), INTENT(IN) :: unit(:)   ! - field units (minimum dimensions: nfs) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nls,nfs) ! - first period -
!
! Local scalars
  INTEGER :: i   ! - station index -
  INTEGER :: k   ! - time index -
  INTEGER :: l   ! - lag index -
  INTEGER :: ifd ! - field index -
  INTEGER :: nn  ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO ifd=1,nfs
     DO l=1,nls
        period1(l,ifd)%sdate%iyr=period1_sdate_iyr(l,ifd)
        period1(l,ifd)%sdate%imn=period1_sdate_imn(l,ifd)
        period1(l,ifd)%sdate%idy=period1_sdate_idy(l,ifd)
        period1(l,ifd)%edate%iyr=period1_edate_iyr(l,ifd)
        period1(l,ifd)%edate%imn=period1_edate_imn(l,ifd)
        period1(l,ifd)%edate%idy=period1_edate_idy(l,ifd)
     END DO
  END DO

! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),nfs,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  DO ifd=1,nfs
     CALL write_tag (iout,ifail,                                                                              &
                     cpt_field=TRIM(var(ifd)),cpt_nrow=nn*nls,cpt_ncol=nv(ifd),cpt_row='T',cpt_col='station', &
                     cpt_units=TRIM(unit(ifd)),cpt_missing=miss(ifd))
     IF (ifail/=0) GOTO 1
!
! Print station information
     DO i=1,nv(ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cstn(i,ifd)))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:Y'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlat(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='no',ERR=1) 'cpt:X'
     DO i=1,nv(ifd)
        WRITE (cout,FMT=*) rlng(i,ifd)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print station data
     nn=0
     DO k=1,nt
        DO l=1,nls
           IF (PRESENT(kuse)) THEN
              IF (.NOT.kuse(k)) CYCLE
           END IF
           nn=nn+1
           cout=get_cdate(period1(l,ifd)+(k-1))
           WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
           DO i=1,nv(ifd)
              WRITE (cout,FMT=*) v(i,nn,l,ifd)
              WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
           END DO
           WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
        END DO
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_stns_fields_lags_v11_dp
!
!
!
 SUBROUTINE write_cpt_unrf_v11_sp ( outfile,nv,nt,tseq,v,miss,var,cvar,&
                                    period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                    period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                    ifail,&
                                    kuse&
                                  )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp
!
! Outputs unreferenced data with no additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nv   ! - number of series -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr ! - start year of first period -
  INTEGER, INTENT(IN) :: period1_sdate_imn ! - start month of first period -
  INTEGER, INTENT(IN) :: period1_sdate_idy ! - start day of first period -
  INTEGER, INTENT(IN) :: period1_edate_iyr ! - end year of first period -
  INTEGER, INTENT(IN) :: period1_edate_imn ! - end month of first period -
  INTEGER, INTENT(IN) :: period1_edate_idy ! - end day of first period -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Input arrays
  REAL(KIND=rp), INTENT(IN) :: v(:,:) ! - data (minimum dimensions: nv, nn)
!
  CHARACTER(LEN=*), INTENT(IN) :: cvar(:) ! - variable names (minimum dimensions: nv)
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Locals
!
! Local scalars
  INTEGER :: i  ! - series index -
  INTEGER :: k  ! - time index -
  INTEGER :: nn ! - number of non-missing cases -
!
  TYPE(period) :: period1 ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  period1%sdate%iyr=period1_sdate_iyr
  period1%sdate%imn=period1_sdate_imn
  period1%sdate%idy=period1_sdate_idy
  period1%edate%iyr=period1_edate_iyr
  period1%edate%imn=period1_edate_imn
  period1%edate%idy=period1_edate_idy
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                              &
                  cpt_field=TRIM(var),cpt_nrow=nn,cpt_ncol=nv,cpt_row='T',cpt_col='index', &
                  cpt_missing=REAL(miss,KIND=dp))
  IF (ifail/=0) GOTO 1
!
! Print index names
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cvar(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print unreferenced data
  nn=0
  DO k=1,nt
     IF (PRESENT(kuse)) THEN
        IF (.NOT.kuse(k)) CYCLE
     END IF
     nn=nn+1
     cout=get_cdate(period1+(k-1))
     WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
     DO i=1,nv
        WRITE (cout,FMT=*) v(i,nn)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_unrf_v11_sp
!
!
!
 SUBROUTINE write_cpt_unrf_v11_dp ( outfile,nv,nt,tseq,v,miss,var,cvar,&
                                    period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                    period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                    ifail,&
                                    kuse&
                                  )
!
! Modules
  USE numbers, ONLY: rp=>dp
!
! Outputs unreferenced data with additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nv   ! - number of series -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  INTEGER, INTENT(IN) :: period1_sdate_iyr ! - start year of first period -
  INTEGER, INTENT(IN) :: period1_sdate_imn ! - start month of first period -
  INTEGER, INTENT(IN) :: period1_sdate_idy ! - start day of first period -
  INTEGER, INTENT(IN) :: period1_edate_iyr ! - end year of first period -
  INTEGER, INTENT(IN) :: period1_edate_imn ! - end month of first period -
  INTEGER, INTENT(IN) :: period1_edate_idy ! - end day of first period -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Input arrays
  REAL(KIND=rp), INTENT(IN) :: v(:,:) ! - data (minimum dimensions: nv, nn) -
!
  CHARACTER(LEN=*), INTENT(IN) :: cvar(:) ! - index names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Locals
!
! Local scalars
  INTEGER :: i  ! - series index -
  INTEGER :: k  ! - time index -
  INTEGER :: nn ! - number of non-missing cases -
!
  TYPE(period) :: period1 ! - first period -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  period1%sdate%iyr=period1_sdate_iyr
  period1%sdate%imn=period1_sdate_imn
  period1%sdate%idy=period1_sdate_idy
  period1%edate%iyr=period1_edate_iyr
  period1%edate%imn=period1_edate_imn
  period1%edate%idy=period1_edate_idy
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                              &
                  cpt_field=TRIM(var),cpt_nrow=nn,cpt_ncol=nv,cpt_row='T',cpt_col='index', &
                  cpt_missing=miss)
  IF (ifail/=0) GOTO 1
!
! Print index names
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cvar(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print unreferenced data
  nn=0
  DO k=1,nt
     IF (PRESENT(kuse)) THEN
        IF (.NOT.kuse(k)) CYCLE
     END IF
     nn=nn+1
     cout=get_cdate(period1+(k-1))
     WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
     DO i=1,nv
        WRITE (cout,FMT=*) v(i,nn)
        WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
     END DO
     WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_unrf_v11_dp
!
!
!
 SUBROUTINE write_cpt_unrf_lags_v11_sp ( outfile,nv,nt,nls,tseq,v,miss,var,cvar,&
                                         period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                         period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                         ifail,&
                                         kuse&
                                       )
!
! Modules
  USE numbers, ONLY: rp=>sp,dp
!
! Outputs unreferenced data with additional lags
! Single precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nv   ! - number of series -
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nls  ! - number of lead-times -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period of each lag (minimum dimensions: nls) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:) ! - data (minimum dimensions: nv, nn, nls) -
!
  CHARACTER(LEN=*),INTENT(IN) :: cvar(:) ! - variable names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nls) ! - first period -
!
! Local scalars
  INTEGER :: i  ! - series index -
  INTEGER :: k  ! - time index -
  INTEGER :: l  ! - lag index -
  INTEGER :: nn ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC REAL
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO l=1,nls
     period1(l)%sdate%iyr=period1_sdate_iyr(l)
     period1(l)%sdate%imn=period1_sdate_imn(l)
     period1(l)%sdate%idy=period1_sdate_idy(l)
     period1(l)%edate%iyr=period1_edate_iyr(l)
     period1(l)%edate%imn=period1_edate_imn(l)
     period1(l)%edate%idy=period1_edate_idy(l)
  END DO
!
! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                                   &
                  cpt_field=TRIM(var),cpt_nrow=nn*nls,cpt_ncol=nv,cpt_row='T',cpt_col='index', &
                  cpt_missing=REAL(miss,KIND=dp))
  IF (ifail/=0) GOTO 1
!
! Print index names
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cvar(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print unreferenced data
  nn=0
  DO k=1,nt
     DO l=1,nls
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        cout=get_cdate(period1(l)+(k-1))
        WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
        DO i=1,nv
           WRITE (cout,FMT=*) v(i,nn,l)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_unrf_lags_v11_sp
!
!
!
 SUBROUTINE write_cpt_unrf_lags_v11_dp ( outfile,nv,nt,nls,tseq,v,miss,var,cvar,&
                                         period1_sdate_iyr,period1_sdate_imn,period1_sdate_idy,&
                                         period1_edate_iyr,period1_edate_imn,period1_edate_idy,&
                                         ifail,&
                                         kuse&
                                       )
!
! Modules
  USE numbers, ONLY: rp=>dp
!
! Outputs unreferenced data with additional lags
! Double precision version
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: nv   ! - number of series -
  INTEGER, INTENT(IN) :: nt   ! - number of cases -
  INTEGER, INTENT(IN) :: nls  ! - number of lead-times -
  INTEGER, INTENT(IN) :: tseq ! - time sequence identifier -
                              ! - 1 -> time sequence increase by year -
                              ! - 2 -> time sequence increase by month -
                              ! - 3 -> time sequence increase daily -
!
  REAL(KIND=rp), INTENT(IN) :: miss ! - missing value flag -
!
  CHARACTER(LEN=*), INTENT(IN) :: outfile ! - file name with full path -
  CHARACTER(LEN=*), INTENT(IN) :: var     ! - field variable -
!
! Input arrays
  INTEGER, INTENT(IN) :: period1_sdate_iyr(:) ! - start year of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_imn(:) ! - start month of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_sdate_idy(:) ! - start day of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_iyr(:) ! - end year of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_imn(:) ! - end month of first period of each lag (minimum dimensions: nls) -
  INTEGER, INTENT(IN) :: period1_edate_idy(:) ! - end day of first period of each lag (minimum dimensions: nls) -
!
  REAL(KIND=rp), INTENT(IN) :: v(:,:,:) ! - data (minimum dimensions: nv, nn, nls) -
!
  CHARACTER(LEN=*),INTENT(IN) :: cvar(:) ! - variable names (minimum dimensions: nv) -
!
! - optional input arrays -
  LOGICAL, INTENT(IN), OPTIONAL :: kuse(:) ! - missing cases indicator (minimum dimensions: nt)
                                           ! - The second dimension of v can be nn where nn==COUNT(kuse(:))
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local arrays
  TYPE(period) :: period1(nls) ! - first period -
!
! Local scalars
  INTEGER :: i  ! - series index -
  INTEGER :: k  ! - time index -
  INTEGER :: l  ! - lag index -
  INTEGER :: nn ! - number of non-missing cases -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ACHAR
  INTRINSIC ADJUSTL
  INTRINSIC COUNT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Initialise time sequencing
  iseq=tseq
!
! Initialise starting period
  DO l=1,nls
     period1(l)%sdate%iyr=period1_sdate_iyr(l)
     period1(l)%sdate%imn=period1_sdate_imn(l)
     period1(l)%sdate%idy=period1_sdate_idy(l)
     period1(l)%edate%iyr=period1_edate_iyr(l)
     period1(l)%edate%imn=period1_edate_imn(l)
     period1(l)%edate%idy=period1_edate_idy(l)
  END DO

! Open file and write out first two lines
  iout=ioutd
  DO
     INQUIRE (UNIT=iout,OPENED=lopen)
     IF (.NOT.lopen) EXIT
     iout=iout+1
  END DO
  CALL open_output (iout,TRIM(outfile),1,ifail)
  IF (ifail/=0) THEN
     PRINT *,'OPEN file ',TRIM(outfile),' failed. Exit.'
     CLOSE (UNIT=iout)
     RETURN
  END IF
!
! Print tag line
  IF (PRESENT(kuse)) THEN
     nn=COUNT(kuse(1:nt))
  ELSE
     nn=nt
  END IF
  CALL write_tag (iout,ifail,                                                                  &
                  cpt_field=TRIM(var),cpt_nrow=nn*nls,cpt_ncol=nv,cpt_row='T',cpt_col='index', &
                  cpt_missing=miss)
  IF (ifail/=0) GOTO 1
!
! Print index names
  DO i=1,nv
     WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cvar(i)))
  END DO
  WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
!
! Print unreferenced data
  nn=0
  DO k=1,nt
     DO l=1,nls
        IF (PRESENT(kuse)) THEN
           IF (.NOT.kuse(k)) CYCLE
        END IF
        nn=nn+1
        cout=get_cdate(period1(l)+(k-1))
        WRITE (UNIT=iout,ADVANCE='no',FMT='(A)',ERR=1) TRIM(ADJUSTL(cout))
        DO i=1,nv
           WRITE (cout,FMT=*) v(i,nn,l)
           WRITE (UNIT=iout,FMT='(2A)',ADVANCE='no',ERR=1) ACHAR(9),TRIM(ADJUSTL(cout))
        END DO
        WRITE (UNIT=iout,FMT='(A)',ADVANCE='yes',ERR=1) ' '
     END DO
  END DO
  ifail=0
!
  CLOSE (UNIT=iout)
  RETURN
!
! Error
1 ifail=3
!
  CLOSE(iout)
  RETURN
 END SUBROUTINE write_cpt_unrf_lags_v11_dp
!
!
!
 FUNCTION same_date(d1,d2)
!
! Function type
  LOGICAL :: same_date
!
! Arguments
!
! Input scalars
  TYPE(date), INTENT(IN) :: d1 ! - first date -
  TYPE(date), INTENT(IN) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  same_date=.false.
  IF (d1%iyr/=d2%iyr) RETURN
  IF (d1%imn/=d2%imn) RETURN
  IF (d1%idy/=d2%idy) RETURN
  same_date=.true.
!
  RETURN
 END FUNCTION same_date
!
!
!
 FUNCTION equal_date(d1,i1)
!
! Function type
  LOGICAL :: equal_date
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: i1 ! - constant -
!
  TYPE(date), INTENT(IN) :: d1 ! - first date -
!
! Executable Statements
!
! Compare dates
  equal_date=.false.
  IF (d1%iyr/=i1) RETURN
  IF (d1%imn/=i1) RETURN
  IF (d1%idy/=i1) RETURN
  equal_date=.true.
!
  RETURN
 END FUNCTION equal_date
!
!
!
 FUNCTION lt_date(d1,d2)
!
! Function type
  LOGICAL :: lt_date
!
! Arguments
!
! Input scalars
  TYPE(date), INTENT(IN) :: d1 ! - first date -
  TYPE(date), INTENT(IN) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  lt_date=.true.
  IF (d1%iyr<d2%iyr) THEN
     RETURN
  ELSE IF (d1%iyr==d2%iyr) THEN
     IF (d1%imn<d2%imn) THEN
        RETURN
     ELSE IF (d1%imn==d2%imn) THEN
        IF (d1%idy<d2%idy) RETURN
     END IF
  END IF
  lt_date=.false.
!
  RETURN
 END FUNCTION lt_date
!
!
!
 FUNCTION gt_date(d1,d2)
!
! Function type
  LOGICAL :: gt_date
!
! Arguments
!
! Input scalars
  TYPE(date), INTENT(IN) :: d1 ! - first date -
  TYPE(date), INTENT(IN) :: d2 ! - second date -
!
! Executable Statements
!
! Compare dates
  gt_date=.true.
  IF (d1%iyr>d2%iyr) THEN
     RETURN
  ELSE IF (d1%iyr==d2%iyr) THEN
     IF (d1%imn>d2%imn) THEN
        RETURN
     ELSE IF (d1%imn==d2%imn) THEN
        IF (d1%idy>d2%idy) RETURN
     END IF
  END IF
  gt_date=.false.
!
  RETURN
 END FUNCTION gt_date
!
!
!
 FUNCTION add_date(d,i)
!
! Function type
  TYPE(date) :: add_date
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: i ! - increment -
!
  TYPE(date), INTENT(IN) :: d ! - date -
!
! Executable Statements
!
! Increment date
  add_date=d
  SELECT CASE (iseq)
   CASE (1) ! - year increment -
     add_date%iyr=d%iyr+i
   CASE (3) ! - day increment -
     add_date%idy=d%idy+i
     IF (add_date%idy<1) THEN
        DO
           add_date%idy=add_date%idy+ndays(add_date%iyr,add_date%imn)
           add_date%imn=add_date%imn-1
           IF (add_date%imn<1) THEN
              add_date%iyr=add_date%iyr-1
              add_date%imn=add_date%imn+nmn
           END IF
           IF (add_date%idy>=1) EXIT
        END DO
     ELSE IF (add_date%idy>ndays(add_date%iyr,add_date%imn)) THEN
        DO
           add_date%idy=add_date%idy-ndays(add_date%iyr,add_date%imn)
           add_date%imn=add_date%imn+1
           IF (add_date%imn>nmn) THEN
              add_date%iyr=add_date%iyr+1
              add_date%imn=add_date%imn-nmn
           END IF
           IF (add_date%idy<=ndays(add_date%iyr,add_date%imn)) EXIT
        END DO
     END IF
  END SELECT
!
  RETURN
 END FUNCTION add_date
!
!
!
 FUNCTION add_period(p,i)
!
! Function type
  TYPE(period) :: add_period
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: i ! - increment -
!
  TYPE(period), INTENT(IN) :: p ! - period -
!
! Executable Statements
!
! Increment period
  add_period%sdate=p%sdate+i
  add_period%edate=p%edate+i
!
  RETURN
 END FUNCTION add_period
!
!
!
 FUNCTION ndays(iyr,imn)
!
! Calculates number of days in the month
! NB - assumes the Gregorian calendard as implemented by Britain and the British Empire
!
! Function type
  INTEGER :: ndays
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: iyr ! - year -
  INTEGER, INTENT(IN) :: imn ! - month -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC MOD
!
! Executable Statements
!
! Define number of days
  SELECT CASE (imn)
   CASE (1)  ! January
     ndays=31
   CASE (3)  ! March
     ndays=31
   CASE (4)  ! April
     ndays=30
   CASE (5)  ! May
     ndays=31
   CASE (6)  ! June
     ndays=30
   CASE (7)  ! July
     ndays=31
   CASE (8)  ! August
     ndays=31
   CASE (9)  ! September
     IF (iyr/=1752) THEN
        ndays=30
     ELSE
        ndays=19
     END IF
   CASE (10) ! October
     ndays=31
   CASE (11) ! November
     ndays=30
   CASE (12) ! December
     ndays=31
!
! Check for leap years
   CASE (2)  ! February
     IF (MOD(iyr,4)==0) THEN
        IF (MOD(iyr,100)==0) THEN
           IF ((MOD(iyr,400)==0).AND.(iyr>1752)) THEN
              ndays=29
           ELSE
              ndays=28
           END IF
        ELSE
           ndays=29
        END IF
     ELSE
        ndays=28
     END IF
  END SELECT
!
  RETURN
 END FUNCTION ndays
!
!
!
 FUNCTION date_diff(d1,d2,iseq)
!
! Function type
  INTEGER :: date_diff
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: iseq ! - date sequence -
!
  TYPE(date), INTENT(IN) :: d1 ! - first date -
  TYPE(date), INTENT(IN) :: d2 ! - second date -
!
! Locals
!
! Local scalars
  INTEGER :: iya ! - current year -
  INTEGER :: ima ! - current month -
  INTEGER :: id  ! - direction -
!
! Executable Statements
!
! Compare dates
  SELECT CASE (iseq)
   CASE (1) ! - yearly -
     date_diff=d2%iyr-d1%iyr
     IF ((d1%imn>0).AND.(d2%imn>0)) THEN
        IF (d2%imn-d1%imn<0) THEN
           date_diff=date_diff-1
        ELSE IF ((d1%idy>0).AND.(d2%idy>0)) THEN
           IF (d2%imn-d1%imn<0) THEN
              IF (d2%idy-d1%idy<0) date_diff=date_diff-1
           END IF
        END IF
     END IF
   CASE (3) ! - daily -
     date_diff=d2%idy-d1%idy
     IF (d2>d1) THEN
        id=1
     ELSE IF (d2<d1) THEN
        id=-1
     ELSE
        RETURN
     END IF
     iya=d1%iyr
     ima=d1%imn
     DO
        date_diff=date_diff+ndays(iya,ima)
        ima=ima+id
        IF ((ima>nmn).OR.(ima<1)) THEN
           iya=iya+id
           ima=ima-nmn*id
        END IF
        IF ((iya==d2%iyr).AND.(ima==d2%imn)) EXIT
     END DO
   CASE DEFAULT
     date_diff=0
  END SELECT
!
  RETURN
 END FUNCTION date_diff
!
!
!
 FUNCTION get_cdate_date(adate) RESULT(get_cdate)
!
! Creates date as a character string
!
! ISO format
!
! Function type
  CHARACTER(LEN=ldat) :: get_cdate
!
! Arguments
!
  TYPE(date), INTENT(IN) :: adate ! - date -
!
! Executable Statements
!
! Create date
  IF (adate%imn>0) THEN
     IF (adate%idy>0) THEN
           WRITE (UNIT=get_cdate,FMT='(I4,A,I2.2,A,I2.2)') adate%iyr,'-',adate%imn,'-',adate%idy
     ELSE
           WRITE (UNIT=get_cdate,FMT='(I4,A,I2.2)') adate%iyr,'-',adate%imn
     END IF
  ELSE
     WRITE (UNIT=get_cdate,FMT='(I4)') adate%iyr
  END IF
!
  RETURN
 END FUNCTION get_cdate_date
!
!
!
 FUNCTION get_cdate_period(aperiod) RESULT (get_cdate)
!
! Creates period as a character string
!
! ISO format
!
! Function type
  CHARACTER(LEN=lprd) :: get_cdate
!
! Arguments
!
  TYPE(period), INTENT(IN) :: aperiod ! - date -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC TRIM
!
! Executable Statements
!
! Create date
  IF (aperiod%edate==aperiod%sdate) THEN
     get_cdate=get_cdate_date(aperiod%sdate)
  ELSE
     IF (aperiod%edate%iyr>aperiod%sdate%iyr) THEN
        get_cdate=TRIM(get_cdate_date(aperiod%sdate))//'/'//get_cdate_date(aperiod%edate)
     ELSE IF (aperiod%edate%imn/=aperiod%sdate%imn) THEN
        IF (aperiod%sdate%idy>0) THEN
           WRITE (UNIT=get_cdate,FMT='(I4,4(A,I2.2))') &
              aperiod%sdate%iyr,'-',aperiod%sdate%imn,'-',aperiod%sdate%idy,'/',aperiod%edate%imn,'-',aperiod%edate%idy
        ELSE
           WRITE (UNIT=get_cdate,FMT='(I4,2(A,I2.2))') aperiod%sdate%iyr,'-',aperiod%sdate%imn,'/',aperiod%edate%imn
        END IF
     ELSE
        WRITE (UNIT=get_cdate,FMT='(I4,3(A,I2.2))') &
           aperiod%sdate%iyr,'-',aperiod%sdate%imn,'-',aperiod%sdate%idy,'/',aperiod%edate%idy
     END IF
  END IF
!
  RETURN
 END FUNCTION get_cdate_period
!
!
!
 SUBROUTINE open_output (iout,afile,nfs,ifail)
!
! Opens CPT output file and prints XMLNS header
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: iout ! - output unit number -
  INTEGER, INTENT(IN) :: nfs  ! - number of fields -
!
  CHARACTER(LEN=*), INTENT(IN) :: afile ! - output file name with full path -
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC TRIM
!
! Open output file as 'sequential'
  OPEN (UNIT=iout,FILE=TRIM(afile),ACCESS='sequential',ACTION='write',FORM='formatted', &
        IOSTAT=ifail,STATUS='unknown')
!
! Error
  IF (ifail/=0) THEN
     ifail=1
     RETURN
  END IF
!
! Print XMLNS header 'formatted'
  WRITE (UNIT=iout,FMT='(A)',IOSTAT=ifail) 'xmlns:cpt='//cxmlns_cpt
  IF (ifail/=0) GOTO 1
!
! Print number of fields
  CALL write_tag (iout,ifail, &
                  cpt_nfields=nfs)
!
! Error
1 IF (ifail/=0) ifail=3
!
  RETURN
 END SUBROUTINE open_output
!
!
!
 SUBROUTINE write_tag (iout,ifail,                                                                             &
                       cpt_nfields,cpt_ncats,cpt_name,cpt_field,cpt_c,cpt_prob,cpt_cmode,cpt_mode,cpt_t,cpt_s, &
                       cpt_z,cpt_m,cpt_clev,cpt_limit,cpt_nrow,cpt_ncol,cpt_row,cpt_col,cpt_units,cpt_missing)
!
! Modules
  USE numbers, ONLY: rp=>dp,one=>one_dp
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: iout ! - output unit number -
!
! - optional input scalars -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_nfields ! - number of fields -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_ncats   ! - number of categories -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_c       ! - current category -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_mode    ! - current mode -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_m       ! - ensemble member -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_nrow    ! - number of rows -
  INTEGER, INTENT(IN), OPTIONAL :: cpt_ncol    ! - number of columns -
!
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: cpt_prob    ! - climatological probability -
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: cpt_clev    ! - confidence level -
  REAL(KIND=rp), INTENT(IN), OPTIONAL :: cpt_missing ! - missing values flag -
!
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_field ! - field -
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_name  ! - name -
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_cmode ! - current mode -
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_row   ! - rows -
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_limit ! - confidence limit -
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_col   ! - columns -
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cpt_units ! - units -
!
  TYPE(level), INTENT(IN), OPTIONAL :: cpt_z ! - level -
!
  TYPE(period), INTENT(IN), OPTIONAL :: cpt_t ! - date -
!
  TYPE(date), INTENT(IN), OPTIONAL :: cpt_s ! - start date -
!
! Output scalars
  INTEGER, INTENT(OUT) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  CHARACTER(LEN=    15) :: cfmt ! - format statement -
  CHARACTER(LEN=lvar+1) :: cout ! - output field -
!
  LOGICAL :: lfirst ! - first output field -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ADJUSTL
  INTRINSIC NINT
  INTRINSIC PRESENT
  INTRINSIC TRIM
!
! Executable Statements
!
! Print nfields
  IF (PRESENT(cpt_nfields)) THEN
     WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_nfields),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1) 'cpt:nfields=',cpt_nfields
     RETURN
  END IF
!
! Print ncats
  IF (PRESENT(cpt_ncats)) THEN
     WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_ncats),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1) 'cpt:ncats=',cpt_ncats
     RETURN
  END IF
!
! Print name
  IF (PRESENT(cpt_name)) THEN
     WRITE (UNIT=iout,FMT='(A)',ERR=1) 'cpt:Name='//TRIM(cpt_name)
     RETURN
  END IF
!
! Print tags line
  lfirst=.true.
  IF (PRESENT(cpt_field)) THEN
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:field='//TRIM(cpt_field)
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_c)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_c),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') 'cpt:C=',cpt_c
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_prob)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=cout,FMT=*) cpt_prob
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:clim_prob='//TRIM(ADJUSTL(cout))
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_cmode)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:mode='//TRIM(cpt_cmode)
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_mode)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_mode),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') 'cpt:mode=',cpt_mode
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_m)) THEN
     IF (cpt_m>0) THEN
        IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
        WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_m),')'
        WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') 'cpt:M=',cpt_m
        lfirst=.false.
     END IF
  END IF
  IF ((PRESENT(cpt_clev)).AND.(PRESENT(cpt_limit))) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=iout,FMT='(3A)',ERR=1,ADVANCE='no') 'cpt:climit=',TRIM(cpt_limit),' ('
     WRITE (UNIT=cfmt,FMT='(A,2(I1,A))') '(F',iprec(cpt_clev,3)+3,'.',iprec(cpt_clev,3),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') cpt_clev
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') '%)'
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_z)) THEN
     IF (TRIM(cpt_z%unit)/='none') THEN
        IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
        IF (cpt_z%hght>one) THEN
           WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(NINT(cpt_z%hght)),')'
           WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') 'cpt:Z=',cpt_z
        ELSE
           WRITE (UNIT=cout,FMT=*) cpt_z%hght
           WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:Z='//ADJUSTL(cout)//' '//TRIM(cpt_z%unit)
        END IF
        lfirst=.false.
     END IF
  END IF
  IF (PRESENT(cpt_t)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     cout=get_cdate(cpt_t)
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:T='//TRIM(cout)
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_s)) THEN
     IF (.NOT.(cpt_s==0)) THEN
        IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
        cout=get_cdate(cpt_s)
        WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:S='//TRIM(cout)
        lfirst=.false.
     END IF
  END IF
  IF (PRESENT(cpt_nrow)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_nrow),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') 'cpt:nrow=',cpt_nrow
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_ncol)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=cfmt,FMT='(A,I1,A)') '(A,I',magnitude(cpt_ncol),')'
     WRITE (UNIT=iout,FMT=cfmt,ERR=1,ADVANCE='no') 'cpt:ncol=',cpt_ncol
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_row)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:row='//cpt_row
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_col)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:col='//cpt_col
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_units)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:units='//cpt_units
     lfirst=.false.
  END IF
  IF (PRESENT(cpt_missing)) THEN
     IF (.NOT.lfirst) WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') ', '
     WRITE (UNIT=cout,FMT=*) cpt_missing
     WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='no') 'cpt:missing='//ADJUSTL(cout)
     lfirst=.false.
  END IF
  WRITE (UNIT=iout,FMT='(A)',ERR=1,ADVANCE='yes') ' '
!
  ifail=0
  RETURN
!
1 ifail=1
  RETURN
 END SUBROUTINE write_tag
!
!
!
 FUNCTION magnitude_int(ix) RESULT (magnitude)
!
! Calculates order of magnitude of an integer value
!
! Function type
  INTEGER :: magnitude
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: ix
!
! Locals
!
! Local scalars
  INTEGER :: iax ! - absolute value of argument -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ABS
!
! Executable Statements
!
! Identify order of magnitude
  iax=ABS(ix)
  IF (iax>0) THEN
     magnitude=1
     DO
        IF (iax<10**magnitude) EXIT
        magnitude=magnitude+1
     END DO
  ELSE 
     magnitude=0
  END IF
!
  RETURN
 END FUNCTION magnitude_int
!
!
!
 FUNCTION magnitude_sp(x) RESULT (magnitude)
!
! Calculates order of magnitude of a single precision value
! Single precision version
!
! Modules
  USE numbers, ONLY: rp=>sp,zero=>zero_sp,one=>one_sp,ten=>ten_sp
!
! Function type
  INTEGER :: magnitude
!
! Arguments
!
! Input scalars
  REAL(KIND=rp), INTENT(IN) :: x
!
! Locals
!
! Local scalars
  REAL(KIND=rp) :: ax ! - absolute value of argument -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ABS
!
! Executable Statements
!
! Identify order of magnitude
  ax=ABS(x)
  IF (ax<one) THEN
     IF (ax>zero) THEN
        magnitude=-1
        DO
           IF (ax>ten**magnitude) EXIT
           magnitude=magnitude-1
        END DO
     ELSE
        magnitude=0
     END IF
  ELSE
     magnitude=1
     DO
        IF (ax<ten**magnitude) EXIT
        magnitude=magnitude+1
     END DO
  END IF
!
  RETURN
 END FUNCTION magnitude_sp
!
!
!
 FUNCTION magnitude_dp(x) RESULT (magnitude)
!
! Calculates order of magnitude of a double precision value
! Double precision version
!
! Modules
  USE numbers, ONLY: rp=>dp,zero=>zero_dp,one=>one_dp,ten=>ten_dp
!
! Function type
  INTEGER :: magnitude
!
! Arguments
!
! Input scalars
  REAL(KIND=rp), INTENT(IN) :: x
!
! Locals
!
! Local scalars
  REAL(KIND=rp) :: ax ! - absolute value of argument -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC ABS
!
! Executable Statements
!
! Identify order of magnitude
  ax=ABS(x)
  IF (ax<one) THEN
     IF (ax>zero) THEN
        magnitude=-1
        DO
           IF (ax>ten**magnitude) EXIT
           magnitude=magnitude-1
        END DO
     ELSE
        magnitude=0
     END IF
  ELSE
     magnitude=1
     DO
        IF (ax<ten**magnitude) EXIT
        magnitude=magnitude+1
     END DO
  END IF
!
  RETURN
 END FUNCTION magnitude_dp
!
!
!
 FUNCTION iprec_sp (r,mprec) RESULT (iprec)
!
! Returns number of decimal places
! Single precision version
!
! Modules
  USE numbers, ONLY: rp=>sp,ten=>ten_sp
!
! Function type
  INTEGER :: iprec
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: mprec ! - maximum precision required -
!
  REAL(KIND=rp), INTENT(IN) :: r ! - value -
!
! Locals
!
! Local scalars
  INTEGER :: ip ! - current precision -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC MOD
  INTRINSIC NINT
!
! Executable Statements
!
! Identify precision
  iprec=0
  DO ip=mprec,1,-1
     IF (MOD(NINT(r*ten**mprec),10**ip)==0) EXIT
     iprec=iprec+1
  END DO
!
  RETURN
 END FUNCTION iprec_sp
!
!
!
 FUNCTION iprec_dp (r,mprec) RESULT (iprec)
!
! Returns number of decimal places
! Double precision version
!
! Modules
  USE numbers, ONLY: rp=>dp,ten=>ten_dp
!
! Function type
  INTEGER :: iprec
!
! Arguments
!
! Input scalars
  INTEGER, INTENT(IN) :: mprec ! - maximum precision required -
!
  REAL(KIND=rp), INTENT(IN) :: r ! - value -
!
! Locals
!
! Local scalars
  INTEGER :: ip ! - current precision -
!
! Functions and Subroutines
!
! Intrinsic functions
  INTRINSIC MOD
  INTRINSIC NINT
!
! Executable Statements
!
! Identify precision
  iprec=0
  DO ip=mprec,1,-1
     IF (MOD(NINT(r*ten**mprec),10**ip)==0) EXIT
     iprec=iprec+1
  END DO
!
  RETURN
 END FUNCTION iprec_dp
END MODULE CPT_formatV11
