! Author: Simon Mason
Module iofiles
!
! Modules
#if GUI == 1
  Use data_numbers,      Only: rp
#endif
  Use data_io_constants, Only: ldir, ldir_v15, ldsc, lext, lfil, lfil_v15, lnam, lnam_v15, lprd, lprd_v15, lstr
  Use time,              Only: pdate, pprd
!
! Declarations
!
! Implicit declarations
  Implicit None
!
! Accessibility
  Private
  Public :: file_mtime, file_reset, file_version, files_reset, get_new_file, get_old_file, geto_gen, init_dirs, init_ifile, &
            init_iofiles, init_ofiles, max_width, open_iofile, parse_file, proj_read_iofiles, proj_write_iofiles,           &
            prompt_results_file, reset_iofiles, set_fmt, set_grey_files, set_mfile
#if GUI == 1
  Public :: check_new_file, get_fmt, geto_file
#else
  Public :: output_format, output_precision
#endif
!
! Derived type definitions
!
! - file format -
  Public :: fformat
  Type fformat
     Sequence
     Integer :: iver ! - CPT version number -
     Integer :: ifmt ! - file format -
     Integer :: iacc ! - file access -
     Integer :: iprc ! - data precision -
     Integer :: lrec ! - record length -
  End Type fformat
!
! - generic file -
  Public cptfile
  Type cptfile
     Sequence
     Character(Len=lfil) :: ffile ! - file -
     Character(Len=ldir) :: fdir  ! - file directory -
     Character(Len=lnam) :: fname ! - file name -
     Character(Len=lext) :: fext  ! - file extension -
     Character(Len=ldsc) :: desc  ! - file description -
!
     Type(fformat) :: ffmt        ! - file format -
!
     Logical :: lset              ! - file set? -
  End Type cptfile
!
! - version 15 generic file -
  Public cptfile_v15
  Type cptfile_v15
     Sequence
     Character(Len=lfil_v15) :: ffile ! - file -
     Character(Len=ldir_v15) :: fdir  ! - file directory -
     Character(Len=lnam_v15) :: fname ! - file name -
     Character(Len=    lext) :: fext  ! - file extension -
     Character(Len=    ldsc) :: desc  ! - file description -
!
     Type(fformat) :: ffmt        ! - file format -
!
     Logical :: lset              ! - file set? -
  End Type cptfile_v15
!
! - input files -
  Public :: ifile
  Type ifile
     Sequence
     Type(cptfile) :: cfile
!
     Character(Len=lstr) :: cgss  ! - field structure -
     Character(Len=lprd) :: cprd1 ! - first period in file -
     Character(Len=lprd) :: cprdn ! - last period in file -
     Character(Len=lprd) :: cprdb ! - beginning date of interest -
     Character(Len=  21) :: cssn  ! - season -
!
     Integer :: idstr             ! - data structure flag -
     Integer :: igeog             ! - geographical reference flag -
     Integer :: iseq              ! - time sequence flag -
     Integer :: nms               ! - number of models -
     Integer :: nem               ! - number of ensemble members -
     Integer :: nfs               ! - number of fields (including ensemble members) -
     Integer :: nls               ! - number of lagged fields -
     Integer :: nfl               ! - total number of fields and lagged fields -
     Integer :: ngs               ! - number of categories -
     Integer :: nse               ! - number of sub-seasons -
     Integer :: nt                ! - number of time steps -
     Integer :: ntm               ! - number of time steps (monthly) -
     Integer :: nat               ! - number of available time steps -
     Integer :: it1               ! - index of first date of interest -
     Integer :: ntag              ! - number of additional tag lines -
     Integer :: mtime             ! - date and time modified -
     Integer :: lmax              ! - maximum length of line -
!
     Type(pprd) :: prd1           ! - period of first data -
     Type(pprd) :: prdn           ! - period of last data -
     Type(pdate) :: fdate         ! - first date of interest -
     Type(pdate) :: bdate         ! - beginning date of interest (including persistence component) -
!
     Logical :: lstack            ! - stacked fields? -
     Logical :: lensemble         ! - ensemble fields? -
  End Type ifile
!
! - version 15 input files -
  Public :: ifile_v15
  Type ifile_v15
     Sequence
     Type(cptfile_v15) :: cfile
!
     Character(Len=    lstr) :: cgss  ! - field structure -
     Character(Len=lprd_v15) :: cprd1 ! - first period in file -
     Character(Len=lprd_v15) :: cprdn ! - last period in file -
     Character(Len=      21) :: cssn  ! - season -
!
     Integer :: idstr                 ! - data structure flag -
     Integer :: igeog                 ! - geographical reference flag -
     Integer :: iseq                  ! - time sequence flag -
     Integer :: nms                   ! - number of models -
     Integer :: nem                   ! - number of ensemble members -
     Integer :: nfs                   ! - number of fields (including ensemble members) -
     Integer :: nls                   ! - number of lagged fields -
     Integer :: nfl                   ! - total number of fields and lagged fields -
     Integer :: ngs                   ! - number of categories -
     Integer :: nt                    ! - number of time steps -
     Integer :: nat                   ! - number of available time steps -
     Integer :: it1                   ! - index of first date of interest -
     Integer :: ntag                  ! - number of additional tag lines -
     Integer :: mtime                 ! - date and time modified -
     Integer :: lmax                  ! - maximum length of line -
!
     Type(pprd) :: prd1               ! - period of first data -
     Type(pprd) :: prdn               ! - period of last data -
     Type(pdate) :: fdate             ! - first date of interest -
     Type(pdate) :: bdate             ! - beginning date of interest (including persistence component) -
!
     Logical :: lstack                ! - stacked fields? -
     Logical :: lensemble             ! - ensemble fields? -
  End Type ifile_v15
!
! - version 14 input files -
  Public :: ifile_v14
  Type ifile_v14
     Sequence
     Character(Len=    lfil) :: ffile ! - file -
     Character(Len=ldir_v15) :: fdir  ! - file directory -
     Character(Len=lnam_v15) :: fname ! - file name -
     Character(Len=    lext) :: fext  ! - file extension -
     Character(Len=    lstr) :: cgss  ! - field structure -
     Character(Len=lprd_v15) :: cprd1 ! - first period in file -
     Character(Len=lprd_v15) :: cprdn ! - last period in file -
     Character(Len=      21) :: cssn  ! - season -
!
     Integer :: idstr                 ! - data structure flag -
     Integer :: igeog                 ! - geographical reference flag -
     Integer :: iseq                  ! - time sequence flag -
     Integer :: nms                   ! - number of models -
     Integer :: nem                   ! - number of ensemble members -
     Integer :: nfs                   ! - number of fields (including ensemble members) -
     Integer :: nls                   ! - number of lagged fields -
     Integer :: nfl                   ! - total number of fields and lagged fields -
     Integer :: ngs                   ! - number of categories -
     Integer :: nt                    ! - number of time steps -
     Integer :: nat                   ! - number of available time steps -
     Integer :: it1                   ! - index of first date of interest -
     Integer :: ntag                  ! - number of additional tag lines -
     Integer :: mtime                 ! - date and time modified -
!
     Type(fformat) :: ffmt            ! - file format -
!
     Type(pprd) :: prd1               ! - period of first data -
     Type(pprd) :: prdn               ! - period of last data -
     Type(pdate) :: fdate             ! - first date of interest -
     Type(pdate) :: bdate             ! - beginning date of interest (including persistence component) -
!
     Logical :: lset                  ! - file set? -
     Logical :: lstack                ! - stacked fields? -
     Logical :: lensemble             ! - ensemble fields? -
  End Type ifile_v14
!
! - output files -
  Public ofile
  Type ofile
     Sequence
     Type(cptfile) :: cfile ! - file format -
!
     Integer :: nfile       ! - file number -
  End Type ofile
!
#if GUI == 1
! - graphics files -
  Public gfile
  Type gfile
     Sequence
     Type(cptfile) :: cfile ! - file format -
!
     Real(Kind=rp) :: qual  ! - JPEG quality -
  End Type gfile
!
#endif
! Scalars
!
! Integer scalars
  Integer, Target, Public :: nx ! - total number of X variables -
  Integer, Target, Public :: ny ! - total number of Y variables -
  Integer, Target, Public :: nz ! - total number of Z variables -
!
  Integer, Pointer, Public :: nv ! - total number of variables -
!
  Integer, Public :: im_fopen ! - file open grey menu flag -
  Integer, Public :: im_sat   ! - start-at grey menu flag -
  Integer, Public :: im_sax   ! - X start-at menu grey flag -
  Integer, Public :: im_say   ! - Y start-at menu grey flag -
  Integer, Public :: im_saz   ! - Z start-at menu grey flag -
  Integer, Public :: ipic     ! - graphics type identifier -
  Integer, Public :: itags    ! - include CPT tags in output file -
  Integer, Public :: mxa      ! - total number of available X variables -
  Integer, Public :: mya      ! - total number of available Y variables -
  Integer, Public :: mza      ! - total number of available Z variables -
  Integer, Public :: nfile    ! - file number -
!
  Integer, Private :: icd       ! - delimiter indicator -
  Integer, Private :: ifmtd = 2 ! - file format identifier -
  Integer, Private :: iprcd = 2 ! - data precision identifier -
!
#if GUI == 1
  Integer, Private :: ipg ! - gray data precision flag -
!
#endif
! Character scalars
  Character(Len=   1), Public :: cxy       ! - current input file -
  Character(Len=   1), Public :: cdelim    ! - delimiter -
  Character(Len= 150), Public :: idir      ! - default input file directory -
  Character(Len= 150), Public :: idir_old  ! - old data file directory -
  Character(Len= 150), Public :: odir      ! - default output file directory -
  Character(Len= 150), Public :: odir_old  ! - old output file directory -
  Character(Len=  36), Public :: fname     ! - file name -
  Character(Len=  36), Public :: fname_old ! - old file name -
  Character(Len= 186), Public :: ffile     ! - file -
  Character(Len= 186), Public :: ffile_old ! - old file -
  Character(Len=   4), Public :: fext      ! - file extension -
!
  Character(Len= 150), Private :: lidir ! - latest input directory -
  Character(Len= 150), Private :: lodir ! - latest output directory -
!
!
! Derived type scalars
  Type(fformat), Public :: cformat ! - current format -
!
  Type(ifile), Public :: bkfile ! - backup input file -
  Type(ifile), Public :: mfile  ! - monthly Y data input file -
  Type(ifile), Public :: pfile  ! - persistence data input file -
!
  Type(ifile), Public :: xfile  ! - X data input file -
  Type(ifile), Public :: yfile  ! - Y data input file -
  Type(ifile), Public :: zfile  ! - X forecast data input file -
!
  Type(ofile), Public :: atfile ! - attributes diagram output file -
  Type(ofile), Public :: avfile ! - averages output file -
  Type(ofile), Public :: ccfile ! - canonical correlations output file -
  Type(ofile), Public :: cofile ! - correlations output file -
  Type(ofile), Public :: exfile ! - exceedence probabilities output file -
  Type(ofile), Public :: fcfile ! - X CCA temporal scores for forecasts output file -
  Type(ofile), Public :: fefile ! - prediction error variance output file -
  Type(ofile), Public :: flfile ! - prediction limits output file -
  Type(ofile), Public :: fofile ! - forecast odds output file -
  Type(ofile), Public :: fpfile ! - forecast probabilities output file -
  Type(ofile), Public :: fsfile ! - forecast ensembles output file -
  Type(ofile), Public :: fvfile ! - forecasts output file -
  Type(ofile), Public :: fxfile ! - X temporal scores for forecasts output file -
  Type(ofile), Public :: gifile ! - goodness index output file -
  Type(ofile), Public :: pbfile ! - PC regression coefficients output file -
  Type(ofile), Public :: psfile ! - probabilistic scores output file -
  Type(ofile), Public :: pvfile ! - p-values output file -
  Type(ofile), Public :: rbfile ! - regression coefficients output file -
  Type(ofile), Public :: rcfile ! - retroactive categories output file -
  Type(ofile), Public :: rhfile ! - ranked-hits output file -
  Type(ofile), Public :: rlfile ! - retroactive prediction limits output file -
  Type(ofile), Public :: rofile ! - ROC output file (individual point) -
  Type(ofile), Public :: rpfile ! - retroactive forecast probabilities output file -
  Type(ofile), Public :: rrfile ! - ROC output file (all points) -
  Type(ofile), Public :: skfile ! - skill output file -
  Type(ofile), Public :: thfile ! - thresholds output file -
  Type(ofile), Public :: wrfile ! - weather roulette output file -
  Type(ofile), Public :: xefile ! - X eigenvalues output file -
  Type(ofile), Public :: xifile ! - interpolated X data output file -
  Type(ofile), Public :: xlfile ! - X spatial loadings output file -
  Type(ofile), Public :: xmfile ! - X homogeneous covariance maps output file -
  Type(ofile), Public :: xofile ! - X input data output file -
  Type(ofile), Public :: xsfile ! - X temporal scores output file -
  Type(ofile), Public :: xtfile ! - X homogeneous covariance maps time series output file -
  Type(ofile), Public :: x2file ! - merged input data output file -
  Type(ofile), Public :: ycfile ! - Y observed categories output file -
  Type(ofile), Public :: yefile ! - Y eigenvalues output file -
  Type(ofile), Public :: yhfile ! - cross-validated predictions output file -
  Type(ofile), Public :: ylfile ! - Y spatial loadings output file -
  Type(ofile), Public :: ymfile ! - Y homogeneous covariance maps output file -
  Type(ofile), Public :: yofile ! - Y input data output file -
  Type(ofile), Public :: yrfile ! - retroactive predictions output file -
  Type(ofile), Public :: ysfile ! - Y temporal scores output file -
  Type(ofile), Public :: ytfile ! - Y homogeneous covariance maps time series output file -
  Type(ofile), Public :: zifile ! - interpolated Z data output file -
!
#if GUI == 1
  Type(gfile), Public, Save :: pic     ! - graphics output file -
  Type(gfile), Public, Save :: pic_old ! - backup greaphics output file -
!
#else
! Arrays
!
! Integer arrays
  Integer, Dimension(42), Public :: new_ids = & ! - new file identifiers -
    (/ 101, 111, 102, 121, 201, 122, 202, 203, 204, 301, &
       302, 303, 311, 312, 313, 401, 411, 412, 421, 422, &
       431, 432, 501, 502, 511, 512, 514, 513, 531, 521, &
       711, 712, 713, 723, 724, 721, 725, 726, 701, 601, &
       602, 603/)
!
#endif
! Interfaces
!
! Generic interfaces
  Public convert_old_ifile
  Interface convert_old_ifile
   Module Procedure ifile_v14_to_v16
   Module Procedure ifile_v15_to_v16
  End Interface convert_old_ifile
!
Contains
!
!
 Subroutine init_dirs (ifail)
!
! Sets working directories
!
! Modules
  Use data_io_constants, Only: cdir, ddir, default_ini, rdir
!
! Arguments
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: l ! - length of string -
!
  Logical :: le ! - file exists? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim, Trim
  Intrinsic get_environment_variable
!
! Executable Statements
!
! Set root directory
  Call get_environment_variable ('CPT_BIN_DIR', rdir)
  l = Len_Trim(rdir)
  If (Len_Trim(rdir) > 0) Then
     If (rdir(l:l) /= cdir) rdir = Trim(rdir)//cdir
     Inquire (File=Trim(rdir)//default_ini, Exist=le)
     If (.not.le) Then
        ifail = 1
        Return
     End If
  Else
     Inquire (File=default_ini, Exist=le)
     If (le) Then
        rdir = '.'//cdir
     Else
        Inquire (File='..'//cdir//default_ini, Exist=le)
        If (le) Then
           rdir = '..'//cdir
        Else
           Inquire (File='..'//cdir//'..'//cdir//default_ini, Exist=le)
           If (le) Then
              rdir = '..'//cdir//'..'//cdir
           Else
              ifail = 2
              Return
           End If
        End If
     End If
  End If
!
! Set data directory
  ddir = Trim(rdir)//'data'//cdir
  ifail = 0
!
  Return
 End Subroutine init_dirs
!
!
!
 Subroutine reset_iofiles ()
!
! Resets input and output files
!
! Modules
  Use data_io_constants, Only: ddir
#if GUI == 1
  Use data_gui_constants, Only: djqual
!
#endif
! Executable Statements
!
! X input file settings
  Call init_ifile (xfile, &
       dir=idir)
!
! Y input file settings
  Call init_ifile (yfile, &
       dir=idir)
  Call init_ifile (mfile, &
       dir=idir)
  Call init_ifile (pfile, &
       dir=idir)
!
! Backup file settings
  Call init_ifile (bkfile, &
       dir=ddir)
!
! Forecasts data file settings
  Call init_ifile (zfile, &
       dir=idir)
!
! Reset output-file settings
  Call init_ofiles ()
!
! Reset delimiter
  icd = 1
  If (cf_delim()==0) Return
!
! Reset inclusion of CPT output tags
  itags = 1
!
! Reset latest directory
  lidir = ' '
  lidir = ' '
!
#if GUI == 1
! Reset graphics file quality
  pic%qual = djqual
!
#endif
  Return
 End Subroutine reset_iofiles
!
!
!
 Function cf_delim()
!
! Sets delimiting character
!
! Function type
  Integer :: cf_delim
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Achar
!
! Executable Statements
!
! Reset delimiter
  Select Case (icd)
   Case (1) ! - tab -
     cdelim = Achar(9)
   Case (2) ! - space -
     cdelim = ' '
   Case (3) ! - comma -
     cdelim = ','
  End Select
  cf_delim = 1
!
  Return
 End Function cf_delim
!
!
!
 Subroutine init_iofiles ()
!
! Initialises data structures and output file numbers and descriptions
!
! Modules
  Use labels, Only: l_averages, l_correls, l_fcasts, l_goodness, l_scoresp, l_threshs, l_wrlt
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Initialise default file format
  Call init_dformat ()
!
! Initialise input data output files
! - observations -
  xofile%nfile = 101
  xofile%cfile%desc = 'X Input Data'
  yofile%nfile = 102
  yofile%cfile%desc = 'Y Input Data'
! - interpolations -
  xifile%nfile = 111
  xifile%cfile%desc = 'Interpolated X Data'
! - categories -
  ycfile%nfile = 121
  ycfile%cfile%desc = 'Y Categories Data'
  rcfile%nfile = 122
  rcfile%cfile%desc = 'Retroactive Categories'
! - merged data -
  x2file%nfile = 131
  x2file%cfile%desc = 'Merged Input Data'
!
! Initialise hindcasts output files
  yhfile%nfile = 201
  yhfile%cfile%desc = 'Cross-Validated Predictions'
  yrfile%nfile = 202
  yrfile%cfile%desc = 'Retroactive Predictions'
  rpfile%nfile = 203
  rpfile%cfile%desc = 'Retroactive Forecast Probabilities'
  rlfile%nfile = 204
  rlfile%cfile%desc = 'Retroactive Prediction Limits'
!
! Initialise modes output files
! - EOFs -
  xefile%nfile = 301
  xefile%cfile%desc = 'X Eigenvalues'
  xlfile%nfile = 302
  xlfile%cfile%desc = 'X Spatial Loadings'
  xsfile%nfile = 303
  xsfile%cfile%desc = 'X Temporal Scores'
  yefile%nfile = 311
  yefile%cfile%desc = 'Y Eigenvalues'
  ylfile%nfile = 312
  ylfile%cfile%desc = 'Y Spatial Loadings'
  ysfile%nfile = 313
  ysfile%cfile%desc = 'Y Temporal Scores'
! - CCA -
  ccfile%nfile = 401
  ccfile%cfile%desc = 'Canonical Correlations'
  xmfile%nfile = 411
  xmfile%cfile%desc = 'X CCA Map Loadings'
  xtfile%nfile = 412
  xtfile%cfile%desc = 'X CCA Map Series'
  ymfile%nfile = 421
  ymfile%cfile%desc = 'Y CCA Map Loadings'
  ytfile%nfile = 422
  ytfile%cfile%desc = 'Y CCA Map Series'
! - MLR and PCR -
  rbfile%nfile = 431
  rbfile%cfile%desc = 'Regression Coefficients'
  pbfile%nfile = 432
  pbfile%cfile%desc = 'PC Regression Coefficients'
!
! Initialise forecasts output files
! - probabilistic forecasts -
  fpfile%nfile = 501
  fpfile%cfile%desc = 'Forecast Probabilities'
  fofile%nfile = 502
  fofile%cfile%desc = 'Forecast Odds'
! - deterministic forecasts -
  fvfile%nfile = 511
  fvfile%cfile%desc = Trim(l_fcasts%c)
  fsfile%nfile = 512
  fsfile%cfile%desc = 'Forecast Ensembles'
  flfile%nfile = 513
  flfile%cfile%desc = 'Prediction Limits'
  fefile%nfile = 514
  fefile%cfile%desc = 'Prediction Error Variances'
! - flexible forecasts -
  exfile%nfile = 521
  exfile%cfile%desc = 'Exceedence Probabilities'
! - X and CCA mode predictors -
  fxfile%nfile = 531
  fxfile%cfile%desc = 'Predictor Time Scores'
  fcfile%nfile = 532
  fcfile%cfile%desc = 'Predictor CCA Time Scores'
! - interpolations -
  zifile%nfile = 541
  zifile%cfile%desc = 'Interpolated Z Data'
!
! Initialise climatologies output files
  thfile%nfile = 601
  thfile%cfile%desc = Trim(l_threshs%c)
  avfile%nfile = 602
  avfile%cfile%desc = Trim(l_averages%c)
  cofile%nfile = 603
  cofile%cfile%desc = Trim(l_correls%c)
!
! Initialise skill scores output files
! - goodness index -
  gifile%nfile = 701
  gifile%cfile%desc = Trim(l_goodness%c)
! - deterministic scores -
  skfile%nfile = 711
  skfile%cfile%desc = 'Skill Scores'
  pvfile%nfile = 712
  pvfile%cfile%desc = 'P-values'
  rofile%nfile = 713
  rofile%cfile%desc = 'ROC Results'
! - probabilistic scores -
  psfile%nfile = 721
  psfile%cfile%desc = Trim(l_scoresp%c)
  rrfile%nfile = 723
  rrfile%cfile%desc = 'ROC Results'
  atfile%nfile = 724
  atfile%cfile%desc = 'Reliability Results'
  rhfile%nfile = 725
  rhfile%cfile%desc = 'Ranked-hits Results'
  wrfile%nfile = 726
  wrfile%cfile%desc = Trim(l_wrlt%c)
!
  Return
 End Subroutine init_iofiles
!
!
!
 Subroutine init_dformat ()
!
! Initialises default file format
!
! Executable Statements
!
! Initialise default file format
  cformat%iver = 0
  cformat%ifmt = ifmtd
  cformat%iacc = 1
  cformat%iprc = iprcd
  cformat%lrec = 0
!
  Return
 End Subroutine init_dformat
!
!
!
 Subroutine init_ifile (afile, &
            dir)
!
! Initialises an input file
!
! Modules
  Use labels, Only: ca_dstruct_l
  Use time,   Only: Assignment(=)
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: dir ! - directory -
!
! Input/output scalars
  Type(ifile), Intent(InOut) :: afile ! - input file -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Present
!
! Executable Statements
!
! Set input file settings
  afile%cfile%ffile = ' '
  If (Present(dir)) afile%cfile%fdir=dir
  afile%cfile%fname = ' '
  afile%cfile%fext = ' '
  afile%cfile%desc = ' '
  afile%cfile%ffmt%iver = 0
  afile%cfile%ffmt%ifmt = 2
  afile%cfile%ffmt%iacc = 1
  afile%cfile%ffmt%iprc = 2
  afile%cfile%ffmt%lrec = 0
  afile%cfile%lset = .false.
  afile%cgss = ca_dstruct_l(0)(1:Len(afile%cgss))
  afile%cprd1 = 'N/A'
  afile%cprdn = 'N/A'
  afile%cprdb = ' '
  afile%cssn = ' '
  afile%idstr = 0
  afile%igeog = 0
  afile%iseq = 0
  afile%nms = 0
  afile%nem = 0
  afile%nfs = 0
  afile%nls = 0
  afile%nfl = 0
  afile%nse = 0
  afile%ngs = 0
  afile%nt = 0
  afile%ntm = 0
  afile%nat = 0
  afile%it1 = 0
  afile%ntag = 0
  afile%mtime = 0
  afile%lmax = 0
  afile%prd1 = 0
  afile%prdn = 0
  afile%fdate = 0
  afile%bdate = 0
  afile%lstack = .false.
  afile%lensemble = .false.
!
  Return
 End Subroutine init_ifile
!
!
!
 Subroutine init_ofiles ()
!
! Initialises output files
!
! Executable Statements
!
! Reset default format
  Call init_dformat ()
!
! Input files
! - X input data -
  Call init_file (xofile%cfile)
! - interpolated X data -
  Call init_file (xifile%cfile)
! - Y input data -
  Call init_file (yofile%cfile)
! - Y categories -
  Call init_file (ycfile%cfile)
! - merged input data -
  Call init_file (x2file%cfile)
!
! Prediction file settings
! - cross-validated predictions -
  Call init_file (yhfile%cfile)
! - retroactive categories -
  Call init_file (rcfile%cfile)
! - retroactive predictions -
  Call init_file (yrfile%cfile)
! - retroactive forecast probabilities -
  Call init_file (rpfile%cfile)
! - retroactive prediction limits -
  Call init_file (rlfile%cfile)
!
! X output file settings
! - eigenvalues -
  Call init_file (xefile%cfile)
! - spatial loadings -
  Call init_file (xlfile%cfile)
! - temporal scores -
  Call init_file (xsfile%cfile)
! - homogeneous maps -
  Call init_file (xmfile%cfile)
! - homogeneous time series -
  Call init_file (xtfile%cfile)
!
! Y output file settings
! - eigenvalues -
  Call init_file (yefile%cfile)
! - spatial loadings -
  Call init_file (ylfile%cfile)
! - temporal scores -
  Call init_file (ysfile%cfile)
! - homogeneous maps -
  Call init_file (ymfile%cfile)
! - homogeneous time series -
  Call init_file (ytfile%cfile)
!
! Canonical correlation output file settings
  Call init_file (ccfile%cfile)
!
! Regression output file settings
  Call init_file (rbfile%cfile)
!
! PC regression output file settings
  Call init_file (pbfile%cfile)
!
! Forecast files
! - forecast probabilities -
  Call init_file (fpfile%cfile)
! - forecast odds -
  Call init_file (fofile%cfile)
! - forecasts -
  Call init_file (fvfile%cfile)
! - forecast ensembles -
  Call init_file (fsfile%cfile)
! - prediction error variances -
  Call init_file (fefile%cfile)
! - prediction limits -
  Call init_file (flfile%cfile)
! - exceedence probabilities -
  Call init_file (exfile%cfile)
! - interpolated X data -
  Call init_file (zifile%cfile)
! - predictor time scores -
  Call init_file (fxfile%cfile)
! - predictor CCA time scores -
  Call init_file (fcfile%cfile)
!
! Skill scores and p-values files
  Call init_file (skfile%cfile)
  Call init_file (pvfile%cfile)
!
! Verification files
! - ROC files -
  Call init_file (rofile%cfile)
  Call init_file (rrfile%cfile)
! - attributes diagram file -
  Call init_file (atfile%cfile)
! - probabilistic scores file -
  Call init_file (psfile%cfile)
! - attributes diagram file -
  Call init_file (rhfile%cfile)
! - weather roulette file -
  Call init_file (wrfile%cfile)
! - goodness index file -
  Call init_file (gifile%cfile)
!
! Climatologies
! - thresholds file -
  Call init_file (thfile%cfile)
! - averages file -
  Call init_file (avfile%cfile)
! - correlations file -
  Call init_file (cofile%cfile)
!
  Return
 End Subroutine init_ofiles
!
!
!
 Subroutine init_file (f)
!
! Initialises a file
!
! Arguments
!
! Input/output scalars
  Type(cptfile), Intent(InOut) :: f ! - output file -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len, Repeat
!
! Executable Statements
  f%lset = .false.
  f%fname = Repeat(' ', Len(f%fname))
  f%ffile = Repeat(' ', Len(f%ffile))
  f%fdir = odir
  f%fext = Repeat(' ', Len(f%fext))
  f%ffmt%ifmt = cformat%ifmt
  f%ffmt%iacc = 1
  f%ffmt%iprc = cformat%iprc
  f%ffmt%lrec = 0
!
  Return
 End Subroutine init_file
!
!
!
 Subroutine file_reset (lfile, fname)
!
! Resets a file's name and its setting flag
!
! Arguments
!
! Output scalars
  Character(Len=*), Intent(Out) :: fname ! - file name -
!
  Logical, Intent(Out) :: lfile ! - file set? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len, Repeat
!
! Executable Statements
!
! Reset file flag
  lfile = .false.
!
! Reset file name
  fname = Repeat(' ', Len(fname))
!
  Return
 End Subroutine file_reset
!
!
!
 Subroutine files_reset ()
!
! Resets output files
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim
!
! Executable Statements
!
! Input data with missing values estimated
! - X input data -
  If (xofile%cfile%lset .or. Len_Trim(xofile%cfile%fname) > 0) Call file_reset (xofile%cfile%lset, xofile%cfile%fname)
! - interpolated X data -
  If (xifile%cfile%lset .or. Len_Trim(xifile%cfile%fname) > 0) Call file_reset (xifile%cfile%lset, xifile%cfile%fname)
! - Y input data -
  If (yofile%cfile%lset .or. Len_Trim(yofile%cfile%fname) > 0) Call file_reset (yofile%cfile%lset, yofile%cfile%fname)
! - Y categories -
  If (ycfile%cfile%lset .or. Len_Trim(ycfile%cfile%fname) > 0) Call file_reset (ycfile%cfile%lset, ycfile%cfile%fname)
! - merged input data -
  If (x2file%cfile%lset .or. Len_Trim(x2file%cfile%fname) > 0) Call file_reset (x2file%cfile%lset, x2file%cfile%fname)
!
! Predictions
! - cross-validated predictions -
  If (yhfile%cfile%lset .or. Len_Trim(yhfile%cfile%fname) > 0) Call file_reset (yhfile%cfile%lset, yhfile%cfile%fname)
! - retroactive categories -
  If (rcfile%cfile%lset .or. Len_Trim(rcfile%cfile%fname) > 0) Call file_reset (rcfile%cfile%lset, rcfile%cfile%fname)
! - retroactive predictions -
  If (yrfile%cfile%lset .or. Len_Trim(yrfile%cfile%fname) > 0) Call file_reset (yrfile%cfile%lset, yrfile%cfile%fname)
! - retroactive forecast probabilities -
  If (rpfile%cfile%lset .or. Len_Trim(rpfile%cfile%fname) > 0) Call file_reset (rpfile%cfile%lset, rpfile%cfile%fname)
! - retroactive prediction limits -
  If (rlfile%cfile%lset .or. Len_Trim(rlfile%cfile%fname) > 0) Call file_reset (rlfile%cfile%lset, rlfile%cfile%fname)
!
! Eigenvalues
! - X eigenvalues -
  If (xefile%cfile%lset .or. Len_Trim(xefile%cfile%fname) > 0) Call file_reset (xefile%cfile%lset, xefile%cfile%fname)
! - Y eigenvalues -
  If (yefile%cfile%lset .or. Len_Trim(yefile%cfile%fname) > 0) Call file_reset (yefile%cfile%lset, yefile%cfile%fname)
!
! Spatial loadings
! - X spatial loadings -
  If (xlfile%cfile%lset .or. Len_Trim(xlfile%cfile%fname) > 0) Call file_reset (xlfile%cfile%lset, xlfile%cfile%fname)
! - Y spatial loadings -
  If (ylfile%cfile%lset .or. Len_Trim(ylfile%cfile%fname) > 0) Call file_reset (ylfile%cfile%lset, ylfile%cfile%fname)
!
! Temporal scores
! - X scores -
  If (xsfile%cfile%lset .or. Len_Trim(xsfile%cfile%fname) > 0) Call file_reset (xsfile%cfile%lset, xsfile%cfile%fname)
! - Y scores -
  If (ysfile%cfile%lset .or. Len_Trim(ysfile%cfile%fname) > 0) Call file_reset (ysfile%cfile%lset, ysfile%cfile%fname)
!
! CCA results
! - canonical correlations -
  If (ccfile%cfile%lset .or. Len_Trim(ccfile%cfile%fname) > 0) Call file_reset (ccfile%cfile%lset, ccfile%cfile%fname)
! - X homogeneous maps -
  If (xmfile%cfile%lset .or. Len_Trim(xmfile%cfile%fname) > 0) Call file_reset (xmfile%cfile%lset, xmfile%cfile%fname)
! - Y homogeneous maps -
  If (ymfile%cfile%lset .or. Len_Trim(ymfile%cfile%fname) > 0) Call file_reset (ymfile%cfile%lset, ymfile%cfile%fname)
! - X homogeneous map series -
  If (xtfile%cfile%lset .or. Len_Trim(xtfile%cfile%fname) > 0) Call file_reset (xtfile%cfile%lset, xtfile%cfile%fname)
! - Y homogeneous map series -
  If (ytfile%cfile%lset .or. Len_Trim(ytfile%cfile%fname) > 0) Call file_reset (ytfile%cfile%lset, ytfile%cfile%fname)
!
! Regression coefficients
! - regression coefficients -
  If (rbfile%cfile%lset .or. Len_Trim(rbfile%cfile%fname) > 0) Call file_reset (rbfile%cfile%lset, rbfile%cfile%fname)
! - PC regression coefficients -
  If (pbfile%cfile%lset .or. Len_Trim(pbfile%cfile%fname) > 0) Call file_reset (pbfile%cfile%lset, pbfile%cfile%fname)
!
! Forecasts
! - forecast probabilities -
  If (fpfile%cfile%lset .or. Len_Trim(fpfile%cfile%fname) > 0) Call file_reset (fpfile%cfile%lset, fpfile%cfile%fname)
! - forecast odds -
  If (fofile%cfile%lset .or. Len_Trim(fofile%cfile%fname) > 0) Call file_reset (fofile%cfile%lset, fofile%cfile%fname)
! - forecast values -
  If (fvfile%cfile%lset .or. Len_Trim(fvfile%cfile%fname) > 0) Call file_reset (fvfile%cfile%lset, fvfile%cfile%fname)
! - forecast ensembles -
  If (fsfile%cfile%lset .or. Len_Trim(fsfile%cfile%fname) > 0) Call file_reset (fsfile%cfile%lset, fsfile%cfile%fname)
! - prediction error variances -
  If (fefile%cfile%lset .or. Len_Trim(fefile%cfile%fname) > 0) Call file_reset (fefile%cfile%lset, fefile%cfile%fname)
! - prediction limits -
  If (flfile%cfile%lset .or. Len_Trim(flfile%cfile%fname) > 0) Call file_reset (flfile%cfile%lset, flfile%cfile%fname)
! - exceedence probabilities -
  If (exfile%cfile%lset .or. Len_Trim(exfile%cfile%fname) > 0) Call file_reset (exfile%cfile%lset, exfile%cfile%fname)
! - predictor forecast scores -
  If (fxfile%cfile%lset .or. Len_Trim(fxfile%cfile%fname) > 0) Call file_reset (fxfile%cfile%lset, fxfile%cfile%fname)
! - predictor CCA forecast scores -
  If (fcfile%cfile%lset .or. Len_Trim(fcfile%cfile%fname) > 0) Call file_reset (fcfile%cfile%lset, fcfile%cfile%fname)
! - interpolated Z data -
  If (zifile%cfile%lset .or. Len_Trim(zifile%cfile%fname) > 0) Call file_reset (zifile%cfile%lset, zifile%cfile%fname)
!
! Skill scores and p-values
  If (skfile%cfile%lset .or. Len_Trim(skfile%cfile%fname) > 0) Call file_reset (skfile%cfile%lset, skfile%cfile%fname)
  If (pvfile%cfile%lset .or. Len_Trim(pvfile%cfile%fname) > 0) Call file_reset (pvfile%cfile%lset, pvfile%cfile%fname)
!
! Verification
! - ROC results -
  If (rofile%cfile%lset .or. Len_Trim(rofile%cfile%fname) > 0) Call file_reset (rofile%cfile%lset, rofile%cfile%fname)
  If (rrfile%cfile%lset .or. Len_Trim(rrfile%cfile%fname) > 0) Call file_reset (rrfile%cfile%lset, rrfile%cfile%fname)
! - reliability results -
  If (atfile%cfile%lset .or. Len_Trim(atfile%cfile%fname) > 0) Call file_reset (atfile%cfile%lset, atfile%cfile%fname)
! - probabilistic scores -
  If (psfile%cfile%lset .or. Len_Trim(psfile%cfile%fname) > 0) Call file_reset (psfile%cfile%lset, psfile%cfile%fname)
! - ranked hits -
  If (rhfile%cfile%lset .or. Len_Trim(rhfile%cfile%fname) > 0) Call file_reset (rhfile%cfile%lset, rhfile%cfile%fname)
! - weather roulette -
  If (wrfile%cfile%lset .or. Len_Trim(wrfile%cfile%fname) > 0) Call file_reset (wrfile%cfile%lset, wrfile%cfile%fname)
! - goodness index -
  If (gifile%cfile%lset .or. Len_Trim(gifile%cfile%fname) > 0) Call file_reset (gifile%cfile%lset, gifile%cfile%fname)
!
! Climatologies
! - thresholds -
  If (thfile%cfile%lset .or. Len_Trim(thfile%cfile%fname) > 0) Call file_reset (thfile%cfile%lset, thfile%cfile%fname)
! - averages -
  If (avfile%cfile%lset .or. Len_Trim(avfile%cfile%fname) > 0) Call file_reset (avfile%cfile%lset, avfile%cfile%fname)
! - correlations -
  If (cofile%cfile%lset .or. Len_Trim(cofile%cfile%fname) > 0) Call file_reset (cofile%cfile%lset, cofile%cfile%fname)
!
  Return
 End Subroutine files_reset
!
!
!
 Subroutine open_iofile (iin, ffile, lread, lfmt, lrecl, ifail)
!
! Opens an input file
!
! Modules
  Use data_io_constants, Only: lfli
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin ! - input file unit number -
!
  Character(Len=*), Intent(In) :: ffile ! - input file -
!
  Logical, Intent(In) :: lread ! - read access? -
  Logical, Intent(In) :: lfmt  ! - formatted file? -
  Logical, Intent(In) :: lrecl ! - record length? -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: ios ! - IO status -
!
  Character(Len=11) :: cfmt ! - formatting -
!
! Executabel Statements
!
! Set formatting
  If (lfmt) Then
     cfmt = 'formatted'
  Else
     cfmt = 'unformatted'
  End If
!
! Open file
  If (lread) Then
     If (lrecl) Then
        Open (Unit=iin, File=ffile, Access='sequential', Action='read', Form=cfmt, IOstat=ios, Recl=lfli, Status='old')
     Else
        Open (Unit=iin, File=ffile, Access='sequential', Action='read', Form=cfmt, IOstat=ios, Status='old')
     End If
  Else
     If (lrecl) Then
        Open (Unit=iin, File=ffile, Access='sequential', Action='write', Form=cfmt, IOstat=ios, Recl=lfli, Status='unknown')
     Else
        Open (Unit=iin, File=ffile, Access='sequential', Action='write', Form=cfmt, IOstat=ios, Status='unknown')
     End If
  End If
  Select Case (ios)
   Case (0)
     ifail = 0
   Case (128)
     ifail = 1
   Case (134)
     ifail = 2
   Case Default
     ifail = 3
  End Select
!
  Return
 End Subroutine open_iofile
!
!
!
 Function file_mtime(afile)
!
! Determines modification date and time of file
!
! Function type
  Integer :: file_mtime
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: afile ! - file -
!
! Locals
!
! Local scalars
#if FTN95 == 1
  Integer(Kind=2) :: mode, dev, rdev, nlink, ierr ! - dummy variables -
!
  Integer(Kind=3) :: isize ! - file size -
  Integer(Kind=3) :: atime ! - file access time -
  Integer(Kind=3) :: mtime ! - file modification time -
  Integer(Kind=3) :: ctime ! - file creation time -
#elif GFORTRAN == 1
  Integer :: ifail ! - error indicator -
!
! Local arrays
  Integer, Dimension(13) :: ifstat ! - file status values -
#endif
!
! Functions and Subroutines
!
#if FTN95 == 1
! Silverfrost library routines
  External :: fileinfo@
#elif GFORTRAN == 1
! Intrinsic routies
  Intrinsic stat
#endif
!
! Executable Statements
!
! Get file info
#if FTN95 == 1
  Call fileinfo@ (afile, mode, dev, rdev, nlink, isize, atime, mtime, ctime, ierr)
#elif GFORTRAN == 1
  Call stat (afile, ifstat, ifail)
#endif
!
! Convert to date format
#if FTN95 == 1
  file_mtime = mtime
#elif GFORTRAN == 1
  If (ifail == 0) Then
     file_mtime = ifstat(10)
  Else
     file_mtime = 0
  End If
#else
  file_mtime = 0
#endif
!
  Return
 End Function file_mtime
!
!
!
 Function file_version(afile)
!
! Determines CPT file format version
!
! Modules
  Use data_io_constants, Only: cxmlns, cxmlns_cpt, iin, ltag
!
! Function type
  Integer :: file_version
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: afile ! - filename -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
  Character(Len=ltag) :: ctag ! - CPT file tag -
!
  Logical :: lrecl ! - record length flag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Len
!
! Executable Statements
!
! Open file
  lrecl = .false.
1 Call open_iofile (iin, afile, .true., .true., lrecl, ifail)
  If (ifail /= 0) Then
     file_version = -ifail
     Return
  End If
!
! Identify file version
  Read (Unit=iin, Fmt='(A)', IOstat=ifail) ctag
  Close (Unit=iin)
  If (ifail /= 0) Then
     If (.not.lrecl) Then
        lrecl = .true.
        GoTo 1
     End If
     file_version=-4
     Return
  End If
!
! Check for CPT XML namespace flag
  If (Index(ctag,cxmlns//':cpt='//cxmlns_cpt(1:Len(cxmlns_cpt)-1)) > 0) Then
     file_version = 10
  Else
     file_version = 9
  End If
!
  Return
 End Function file_version
!
!
!
 Subroutine set_mfile (icsv, lnewy)
!
! Sets a monthly input file
!
! Modules
  Use time,   Only: get_cdate
  Use screen, Only: window_update
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: icsv ! - calculate seasonal values flag -
!
  Logical, Intent(In) :: lnewy ! - new Y-file? -
!
! Locals
!
! Local scalars
  Integer :: mseq ! - time-sequencing -
!
  Type(pprd) :: mprd1 ! - period of first data -
  Type(pprd) :: mprdn ! - period of last data -
!
! Executable Statements
!
! Set monthly input file settings
  If (icsv == 0) Then
     mfile = yfile
     mfile%cfile%desc = 'M Input File'
  Else
! - backup current settings -
     If (lnewy) Then
        mfile = yfile
     Else
        mseq = mfile%iseq
        mprd1 = mfile%prd1
        mprdn = mfile%prdn
        mfile = yfile
! - restore previous settings -
        mfile%iseq = mseq
        mfile%prd1 = mprd1
        mfile%prdn = mprdn
        mfile%cprd1 = get_cdate(mfile%prd1, 1, 2)
        mfile%cprdn = get_cdate(mfile%prdn, 1, 2)
     End If
  End If
!
! Update settings
  Call window_update (mfile%cprd1)
  Call window_update (mfile%cprdn)
!
  Return
 End Subroutine set_mfile
!
!
!
 Subroutine max_width (iin, lmax, ifail, &
            afile, idat, cline)
!
! Determines maximum width of file
!
! On exit:
!    ifail =  0 Successful
!    ifail =  1 Problem reading file
!    ifail =  2 File is too wide
!
! Modules
  Use data_io_constants, Only: lfli
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin  ! - input unit number -
!
! -optional input scalars -
  Integer, Intent(In), Optional :: idat ! - line mumber to return -
!
  Character(Len=*), Intent(In), Optional :: afile ! - file name -
!
! Output scalars
  Integer, Intent(Out) :: lmax  ! - maximum length of line -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! -optional output scalars -
  Character(Len=*), Intent(Out), Optional :: cline ! - last line of input -
!
! Locals
!
! Local scalars
  Integer :: i    ! - line counter -
  Integer :: l    ! - line number -
  Integer :: emax ! - estimated maximum length of line -
  Integer :: llen ! - current length of line -
!
  Character(Len=lfli) :: cl ! - line of input -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Len_Trim, Min, Present
!
! Executable Statements
!
! Open file
  If (Present(afile)) Then
     Call open_iofile (iin, afile, .true., .true., .true., ifail)
     If (ifail /= 0) GoTo 2
  End If
!
! Determine maximum length of line of data
  lmax = 0
  emax = 100
  i = 1
  If (Present(idat)) Then
     l = idat
  Else
     l = 2
  End If
  cl=' '
  Do
     Read (Unit=iin, Fmt='(A)', Err=2, End=1) cl(1:emax)
     llen = Len_Trim(cl(1:emax))
     If (llen <= lmax) Then ! - current line is shorter - 
        If (i == l) Then
           emax = Min(lmax+100, lfli)
           If (Present(cline)) cline = cl ! - return first line -
        End If
        i = i + 1
     Else ! - extend if line is too short -
        lmax = llen
        If (lmax >= lfli - 100) Then ! - return if line is too long -
           ifail = 4
           Return
        End If
        emax = Min(lmax+100, lfli)
        If (emax < 2**13) Then
           Backspace (Unit=iin)
        Else
           Rewind (Unit=iin)
        End If
     End If
  End Do
!
! Close file
1 If (Present(afile)) Close (Unit=iin)
  ifail = 0
  Return
!
! Errors
2 ifail = 1
  Return
 End Subroutine max_width
!
!
!
 Subroutine ifile_v15_to_v16 (afile_v15, afile_v16)
!
! Converts version 14 input files to version 16
!
! Modules
  Use time, Only: get_cdate
!
! Arguments
!
! Input scalars
  Type(ifile_v15), Intent(In) :: afile_v15 ! - version 14 file -
!
! Output scalars
  Type(ifile), Intent(Out) :: afile_v16 ! - version 16 file -
!
! Executable Statements
!
! Convert file
  afile_v16%cfile%ffile = afile_v15%cfile%ffile   ! - file -
  afile_v16%cfile%fdir = afile_v15%cfile%fdir     ! - file directory -
  afile_v16%cfile%fname = afile_v15%cfile%fname   ! - file name -
  afile_v16%cfile%fext = afile_v15%cfile%fext     ! - file extension -
  afile_v16%cfile%desc = afile_v15%cfile%desc     ! - file extension -
  afile_v16%cfile%ffmt = afile_v15%cfile%ffmt     ! - file format -
  afile_v16%cfile%lset = afile_v15%cfile%lset     ! - file flag -
  afile_v16%cgss = afile_v15%cgss                 ! - file structure -
  afile_v16%cprd1 = afile_v15%cprd1               ! - first date in file -
  afile_v16%cprdn = afile_v15%cprdn               ! - last date in file -
  afile_v16%cprdb = get_cdate(afile_v15%bdate, 2) ! - start-at-date -
  afile_v16%cssn = afile_v15%cssn                 ! - season -
  afile_v16%idstr = afile_v15%idstr               ! - data structure flag -
  afile_v16%igeog = afile_v15%igeog               ! - geographical reference flag -
  afile_v16%iseq = afile_v15%iseq                 ! - time sequence flag -
  afile_v16%nms = afile_v15%nms                   ! - number of models -
  afile_v16%nem = afile_v15%nem                   ! - number of ensemble members -
  afile_v16%nfs = afile_v15%nfs                   ! - number of fields -
  afile_v16%nls = afile_v15%nls                   ! - number of lagged fields -
  afile_v16%nfl = afile_v15%nfl                   ! - total number of fields and lagged fields -
  afile_v16%nse = 1                               ! - number of sub-seasons -
  afile_v16%ngs = afile_v15%ngs                   ! - number of categories -
  afile_v16%nt = afile_v15%nt                     ! - number of time steps -
  afile_v16%ntm = afile_v15%nt                    ! - number of time steps (monthly) -
  afile_v16%nat = afile_v15%nat                   ! - number of available time steps -
  afile_v16%it1 = afile_v15%it1                   ! - index of first date of interest -
  afile_v16%ntag = afile_v15%ntag                 ! - number of additional tag lines -
  afile_v16%mtime = afile_v15%mtime               ! - date last modified -
  afile_v16%lmax = afile_V15%lmax                 ! - maximum length of input line -
  afile_v16%prd1 = afile_v15%prd1                 ! - period of first data -
  afile_v16%prdn = afile_v15%prdn                 ! - period of last data -
  afile_v16%fdate = afile_v15%fdate               ! - parsed first date of interest -
  afile_v16%bdate = afile_v15%bdate               ! - parsed beginning date -
  afile_v16%lstack = afile_v15%lstack             ! - stacked fields flag -
  afile_v16%lensemble = afile_v15%lensemble       ! - ensemble fields flag -
!
  Return
 End Subroutine ifile_v15_to_v16
!
!
!
 Subroutine ifile_v14_to_v16 (afile_v14, afile_v16)
!
! Converts version 14 input files to version 16
!
! Modules
  Use time, Only: get_cdate
!
! Arguments
!
! Input scalars
  Type(ifile_v14), Intent(In) :: afile_v14 ! - version 14 file -
!
! Output scalars
  Type(ifile), Intent(Out) :: afile_v16 ! - version 16 file -
!
! Executable Statements
!
! Convert file
  afile_v16%cfile%ffile = afile_v14%ffile         ! - file -
  afile_v16%cfile%fdir = afile_v14%fdir           ! - file directory -
  afile_v16%cfile%fname = afile_v14%fname         ! - file name -
  afile_v16%cfile%fext = afile_v14%fext           ! - file extension -
  afile_v16%cfile%desc = ' '                      ! - file extension -
  afile_v16%cfile%ffmt = afile_v14%ffmt           ! - file format -
  afile_v16%cfile%lset = afile_v14%lset           ! - file flag -
  afile_v16%cgss = afile_v14%cgss                 ! - file structure -
  afile_v16%cprd1 = afile_v14%cprd1               ! - first date in file -
  afile_v16%cprdn = afile_v14%cprdn               ! - last date in file -
  afile_v16%cprdb = get_cdate(afile_v14%bdate, 2) ! - start-at-date -
  afile_v16%cssn = afile_v14%cssn                 ! - season -
  afile_v16%idstr = afile_v14%idstr               ! - data structure flag -
  afile_v16%igeog = afile_v14%igeog               ! - geographical reference flag -
  afile_v16%iseq = afile_v14%iseq                 ! - time sequence flag -
  afile_v16%nms = afile_v14%nms                   ! - number of models -
  afile_v16%nem = afile_v14%nem                   ! - number of ensemble members -
  afile_v16%nfs = afile_v14%nfs                   ! - number of fields -
  afile_v16%nls = afile_v14%nls                   ! - number of lagged fields -
  afile_v16%nfl = afile_v14%nfl                   ! - total number of fields and lagged fields -
  afile_v16%nse = 1                               ! - number of sub-seasons -
  afile_v16%ngs = afile_v14%ngs                   ! - number of categories -
  afile_v16%nt = afile_v14%nt                     ! - number of time steps -
  afile_v16%ntm = afile_v14%nt                    ! - number of time steps (monthly) -
  afile_v16%nat = afile_v14%nat                   ! - number of available time steps -
  afile_v16%it1 = afile_v14%it1                   ! - index of first date of interest -
  afile_v16%ntag = afile_v14%ntag                 ! - number of additional tag lines -
  afile_v16%mtime = afile_v14%mtime               ! - date last modified -
  afile_v16%lmax = 0                              ! - maximum length of input line -
  afile_v16%prd1 = afile_v14%prd1                 ! - period of first data -
  afile_v16%prdn = afile_v14%prdn                 ! - period of last data -
  afile_v16%fdate = afile_v14%fdate               ! - parsed first date of interest -
  afile_v16%bdate = afile_v14%bdate               ! - parsed beginning date -
  afile_v16%lstack = afile_v14%lstack             ! - stacked fields flag -
  afile_v16%lensemble = afile_v14%lensemble       ! - ensemble fields flag -
!
  Return
 End Subroutine ifile_v14_to_v16
!
!
!
 Subroutine get_old_file (ftype, ffile, fdir, ftypes, ffilts, nfilts, fname, ifail, &
            fext)
!
! Gets an existing file
!
! Modules
#if GUI == 1
  Use errors, Only: cpt_error
#endif
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfilts ! - number of file filters -
!
  Character(Len=*), Intent(In) :: ftype ! - file type -
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: fdir ! - file directory -
!
! - optional input/output scalars -
  Character(Len=*), Intent(InOut), Optional :: fext ! - selected file extension -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
  Character(Len=*), Intent(Out) :: ffile ! - selected file -
  Character(Len=*), Intent(Out) :: fname ! - selected file name -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: ftypes ! - file types -
  Character(Len=*), Dimension(:), Intent(In) :: ffilts ! - file filters -
!
#if GUI == 1
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Len_Trim
!
! ClearWin+ routines
  External :: get_filtered_file$
!
#endif
! Executable Statements
!
! Default to latest directory if directory is unset
!
! Create window
  ffile = ' '
#if GUI == 1
  If (Len_Trim(fdir) == 0 .and. Len_Trim(lidir) > 0) Then
     Call get_filtered_file$ ('Open '//ftype, ffile, lidir, ftypes, ffilts, nfilts, 1)
  Else
     Call get_filtered_file$ ('Open '//ftype, ffile, fdir, ftypes, ffilts, nfilts, 1)
  End If
#else
  Call get_filtered_file (ftype, ffile, 1)
#endif
!
#if GUI == 1
! Check for % sign in file name
  If (Index(ffile, '%') /= 0) Then
     ifail = 1
     Call cpt_error ('get_old_file', ifail, .false.)
     ifail = 1
     Return
  End If
!
#endif
! Isolate file name
  Call parse_file (ffile, fdir, fname, ifail, &
       fext_new=fext)
  If (ifail == 0) lidir = fdir
!
  Return
 End Subroutine get_old_file
!
!
!
 Subroutine get_new_file (ftype, lnewname, ffile, fdir, fext, ftypes, ffilts, nfilts, fname, ifail)
!
! Prompts for a new file name
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nfilts ! - number of file filters -
!
  Character(Len=*), Intent(In) :: ftype ! - file type -
  Character(Len=*), Intent(In) :: fext  ! - file extension -
!
  Logical, Intent(In) :: lnewname ! - new filename? -
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: fdir ! - file directory -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
  Character(Len=*), Intent(Out) :: ffile ! - selected file -
  Character(Len=*), Intent(Out) :: fname ! - selected file name -
!
! Input arrays
  Character(Len=*), Dimension(:), Intent(In) :: ftypes ! - file types -
  Character(Len=*), Dimension(:), Intent(In) :: ffilts ! - file filters -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Len, Len_Trim, Repeat
!
#if GUI == 1
! ClearWin+ routines
  External :: get_filtered_file$
!
#endif
! Executable Statements
!
! Select file window
  Do
     If (lnewname) ffile = Repeat(' ',Len(ffile))
#if GUI == 1
     If (Len_Trim(fdir)==0 .and. Len_Trim(lodir) > 0) Then
        Call get_filtered_file$ ('Define '//ftype, ffile, lodir, ftypes, ffilts, nfilts, 0)
     Else
        Call get_filtered_file$ ('Define '//ftype, ffile, fdir, ftypes, ffilts, nfilts, 0)
     End If
#else
     Call get_filtered_file (ftype, ffile, 0)
#endif
!
! Isolate file name
     If (Index(ffile, ffilts(1)(2:Len_Trim(ffilts(1))), Back=.true.) > 0) Then
        Call parse_file (ffile, fdir, fname, ifail, &
             fext_old=ffilts(1)(2:Len_Trim(ffilts(1))))
     Else
        Call parse_file (ffile, fdir, fname, ifail, &
             fext_old=fext)
     End If
     If (ifail == 1) Exit
     lodir = fdir
     ifail = check_new_file(ffile)
     If (ifail == 0) Exit
  End Do
!
  Return
 End Subroutine get_new_file
!
!
!
 Function check_new_file(ffile) &
          Result (check)
!
! Returns 0 if filename is OK to use.
!
! Modules
#if GUI == 1
  Use errors, Only: query
#else
  Use gui,    Only: print_warning
#endif
  Use errors, Only: cpt_error
!
! Function type
  Integer :: check
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: ffile ! - file -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
  Logical :: le ! - file exists? -
!
! Executable Statements
!
! Check file name
! - check for existence of file -
  Inquire (File=ffile, Exist=le, IOstat=ifail)
  If (ifail == 0) Then
     If (le) Then
#if GUI == 1
        check = 1 - query('check_new_file', 1, .false.)
     Else
        check = 0
     End If
#else
        Call print_warning ('File already exists. Overwriting old file', nopause=.true.)
     End If
     check = 0
#endif
!
! Invalid file name
  Else
     ifail = 1
     Call cpt_error ('check_new_file', ifail, .false.)
     check = 1
  End If
!
  Return
 End Function check_new_file
!
!
!
 Subroutine parse_file (ffile, fdir, fname, ifail, &
            fext_old, fext_new)
!
! Parses a file into its directory, name, and extension
!
! Modules
  Use data_io_constants, Only: cdir
!
! Arguments
!
! Input scalars
! - optional input scalars -
  Character(Len=*), Intent(In), Optional :: fext_old ! - file extension -
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: ffile ! - selected file -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
  Character(Len=*), Intent(Out) :: fdir  ! - file directory -
  Character(Len=*), Intent(Out) :: fname ! - selected file name -
!
! - optional output scalars -
  Character(Len=*), Intent(Out), Optional :: fext_new ! - file extension -
!
! Locals
!
! Local scalars
  Integer :: i1, i2, i3, i4 ! - position of filename -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Len_Trim, Present, Scan
!
! Executable Statements
!
! Isolate file name
  i3 = Len_Trim(ffile)
  If (i3 > 0) Then
! - identify directory -
     i1 = Scan(ffile, '/\', Back=.true.)
     If (i1 > 0) Then
        fdir = ffile(1:i1)
     Else
        fdir = '.'//cdir
     End If
! - ientify any file extension -
     If (.not.present(fext_old)) Then
        i2 = Index(ffile, '.', Back=.true.)
        If (i2 > i1) Then
           If (Present(fext_new)) Then
              fname = ffile(i1+1:i2-1)
              fext_new = ffile(i2+1:i3)
           Else
              fname = ffile(i1+1:i3)
           End If
        Else
           fname = ffile(i1+1:i3)
           If (Present(fext_new)) fext = ' '
        End If
     Else
! - remove file extension if present -
        i4 = Max(Index(fext_old, '.', Back=.true.), 1)
        i2 = Index(ffile(i1+1:i3), fext_old(i4:), Back=.true.)
        If (i2 > 0) Then
           i2 = Index(ffile(1:i3), fext_old(i4:), Back=.true.)
           i3 = i2 - 1
        Else
           i2 = Index(fext_old, '.', Back=.true.)
           If (i2 > 0) Then
              ffile = Trim(ffile)//fext_old(i2:)
           Else
              ffile = Trim(ffile)
           End If
        End If
        fname = ffile(i1+1:i3)
     Endif
     ifail = 0
  Else
     ifail = 1
  End If
!
  Return
 End Subroutine parse_file
!
!
!
 Subroutine prompt_results_file (rfile, lfirst)
!
! Prompts for a results file
!
! Modules
  Use gui, Only: gui_creturn, prompt_file
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lfirst ! - first file? -
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: rfile ! - results file -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for results file
  If (.not.lfirst) Call gui_creturn (.false.)
  Call prompt_file (Trim(rfile%cfile%desc), rfile%cfile%fname, rfile%nfile, nfile, geto_file, &
      itab=35)
!
  Return
 End Subroutine prompt_results_file
!
!
!
 Function geto_file()
!
! Prompts for output file name corresponding to output file type number nfile
!
! Function type
  Integer :: geto_file
!
! Executable Statements
!
! Define output file
#if GUI == 1
  If (nfile == xofile%nfile) Then ! - X input data output file -
     Call geto_gen (xofile)
  Else If (nfile == xifile%nfile) Then ! - interpolated X data output file -
     Call geto_gen (xifile)
  Else If (nfile == yofile%nfile) Then ! - Y input data output file -
     Call geto_gen (yofile)
  Else If (nfile == ycfile%nfile) Then ! - Y categories output file -
     Call geto_gen (ycfile)
  Else If (nfile == yhfile%nfile) Then ! - cross-validated predictions output file -
     Call geto_gen (yhfile)
  Else If (nfile == rcfile%nfile) Then ! - retroactive categories output file -
     Call geto_gen (rcfile)
  Else If (nfile == yrfile%nfile) Then ! - retroactive predictions output file -
     Call geto_gen (yrfile)
  Else If (nfile == rpfile%nfile) Then ! - retroactive forecast probabilities output file -
     Call geto_gen (rpfile)
  Else If (nfile == rlfile%nfile) Then ! - retroactive prediction limits output file -
     Call geto_gen (rlfile)
  Else If (nfile == xefile%nfile) Then ! - X eigenvalues output file -
     Call geto_gen (xefile)
  Else If (nfile == xlfile%nfile) Then ! - X spatial loadings output file -
     Call geto_gen (xlfile)
  Else If (nfile == xsfile%nfile) Then ! - X temporal scores output file -
     Call geto_gen (xsfile)
  Else If (nfile == yefile%nfile) Then ! - Y eigenvalues output file -
     Call geto_gen (yefile)
  Else If (nfile == ylfile%nfile) Then ! - Y spatial loadings output file -
     Call geto_gen (ylfile)
  Else If (nfile == ysfile%nfile) Then ! - Y temporal scores output file -
     Call geto_gen (ysfile)
  Else If (nfile == ccfile%nfile) Then ! - canonical correlations output file -
     Call geto_gen (ccfile)
  Else If (nfile == xmfile%nfile) Then ! - X homogeneous covariance maps output file -
     Call geto_gen (xmfile)
  Else If (nfile == xtfile%nfile) Then ! - X homogeneous covariance maps time series output file -
     Call geto_gen (xtfile)
  Else If (nfile == ymfile%nfile) Then ! - Y homogeneous covariance maps output file -
     Call geto_gen (ymfile)
  Else If (nfile == ytfile%nfile) Then ! - Y homogeneous covariance maps time series output file -
     Call geto_gen (ytfile)
  Else If (nfile == rbfile%nfile) Then ! - regression coefficients output file -
     Call geto_gen (rbfile)
  Else If (nfile == pbfile%nfile) Then ! - PC regression coefficients output file -
     Call geto_gen (pbfile)
  Else If (nfile == fpfile%nfile) Then ! - forecast probabilities output file -
     Call geto_gen (fpfile)
  Else If (nfile == fofile%nfile) Then ! - forecast odds output file -
     Call geto_gen (fofile)
  Else If (nfile == fvfile%nfile) Then ! - forecasts output file -
     Call geto_gen (fvfile)
  Else If (nfile == fsfile%nfile) Then ! - forecast ensembles output file -
     Call geto_gen (fsfile)
  Else If (nfile == fefile%nfile) Then ! - prediction error variances output file -
     Call geto_gen (fefile)
  Else If (nfile == flfile%nfile) Then ! - prediction limits output file -
     Call geto_gen (flfile)
  Else If (nfile == exfile%nfile) Then ! - exceedence probabilities output file -
     Call geto_gen (exfile)
  Else If (nfile == zifile%nfile) Then ! - interpolated Z data output file -
     Call geto_gen (zifile)
  Else If (nfile == fxfile%nfile) Then ! - predictor time scores output file -
     Call geto_gen (fxfile)
  Else If (nfile == fcfile%nfile) Then ! - predictor CCA time scores output file -
     Call geto_gen (fcfile)
  Else If (nfile == skfile%nfile) Then ! - skill output file -
     Call geto_gen (skfile)
  Else If (nfile == pvfile%nfile) Then ! - p-values output file -
     Call geto_gen (pvfile)
  Else If (nfile == rofile%nfile) Then ! - ROC output file -
     Call geto_gen (rofile)
  Else If (nfile == rrfile%nfile) Then ! - ROC output file -
     Call geto_gen (rrfile)
  Else If (nfile == atfile%nfile) Then ! - attributes diagram output file -
     Call geto_gen (atfile)
  Else If (nfile == psfile%nfile) Then ! - probabilistic scores output file -
     Call geto_gen (psfile)
  Else If (nfile == rhfile%nfile) Then ! - ranked hits output file -
     Call geto_gen (rhfile)
  Else If (nfile == wrfile%nfile) Then ! - weather roulette output file -
     Call geto_gen (wrfile)
  Else If (nfile == gifile%nfile) Then ! - goodness index output file -
     Call geto_gen (gifile)
  Else If (nfile == thfile%nfile) Then ! - thresholds output file -
     Call geto_gen (thfile)
  Else If (nfile == avfile%nfile) Then ! - averages output file -
     Call geto_gen (avfile)
  Else If (nfile == cofile%nfile) Then ! - correlations output file -
     Call geto_gen (cofile)
  Else If (nfile == x2file%nfile) Then ! - merged input data output file -
     Call geto_gen (x2file)
  End If
#endif
  geto_file=2
!
  Return
 End Function geto_file
!
!
!
#if GUI == 0
 Subroutine get_filtered_file (title, ffile, oldnew)
!
! Prompts for a file name, and checks whether it exists (corresponds to FTN95 file extension)
!
! Modules
  Use gui, Only: iprompt
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: oldnew ! - old file / new file flag -
! 
  Character(Len=*), Intent(In) :: title ! - file title -
! 
! Output scalars
  Character(Len=*), Intent(Out) :: ffile ! - file -
!
! Locals
!
! Local scalars
  Logical :: le ! - file exists? -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for file
1 Write (Unit=*, Fmt='(1X,2A)', Advance='no') title,': '
  Read (Unit=*, Fmt='(A)', Err=2, End=2) ffile
!
! Check for existence
  If (oldnew == 1) Then
     Inquire (File=Trim(ffile), Exist=le)
     If (.not.le) GoTo 2
  End If
!
  Return
!
2 Select Case (iprompt)
   Case (0)
     GoTo 1
   Case (1)
     Continue
   Case (2)
     Stop
  End Select
!
  Return
 End Subroutine get_filtered_file
!
!
!
#endif
 Subroutine geto_gen (afile)
!
! Prompts for an output file
!
! Modules
  Use data_io_constants, Only: ffmts
  Use screen,            Only: window_update
!
! Arguments
!
! Input/output scalars
  Type(ofile), Intent(InOut) :: afile ! - file -
!
! Locals
!
! Local scalars
  Integer :: ifail ! - error indicator -
!
  Character(Len=24) :: cfiles ! - error indicator -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Define output file
  Call set_fmt (0, afile%cfile%ffmt)
  Select Case (Trim(ffmts(afile%cfile%ffmt%ifmt)))
   Case ('formatted')
     afile%cfile%fext = '.txt'
     cfiles='Formatted output files'
   Case ('unformatted')
     afile%cfile%fext = '.dat'
     cfiles='Unformatted output files'
   Case ('GrADS')
     afile%cfile%fext = '.dat'
     cfiles='GrADS output files'
  End Select
  Call get_new_file (Trim(afile%cfile%desc)//' Output file', .false., afile%cfile%ffile, afile%cfile%fdir, afile%cfile%fext, &
       (/cfiles(:)/), (/'*'//afile%cfile%fext/), 1, afile%cfile%fname, ifail)
!
! Confirm selection
  If (ifail == 0) Then
     afile%cfile%lset = .true.
     Call window_update (afile%cfile%fname)
  Else
     Call init_file (afile%cfile)
  End If
!
  Return
 End Subroutine geto_gen
!
!
!
#if GUI == 1
 Function get_fmt()
!
! Sets file format options
!
! Modules
  Use data_io_constants, Only: cdel, cprcs, faccs, ffmts, nacc, ndel, nfmt, nprc
  Use gui,               Only: gui_creturn, init_win, print_text, prompt_list, win_prompt
!
! Function type
  Integer :: get_fmt
!
! Locals
!
! Local scalars
  Type(fformat) :: tformat ! - file format -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Backup settings
  tformat = cformat
!
! Set gray flags
  If (Trim(ffmts(cformat%ifmt)) == 'unformatted') Then
     ipg = 1
  Else If (Trim(ffmts(cformat%ifmt)) == 'formatted') Then
     ipg = 0
  Else If (Trim(ffmts(cformat%ifmt)) == 'GrADS') Then
     ipg = 0
     cformat%iprc = 1
  End If
!
! Create window
  Call init_win (title='File format')
!
! Indicate and prompt for file options
  Call gui_creturn (.true.)
  Call print_text ('The following format options will be used to write the data:')
  Call gui_creturn (.true.)
  Call prompt_list ('FORMAT', .true., ffmts, nfmt, cformat%ifmt, &
       iwdth=13, itab=15, cbf1=cf_ifmt)
  Call gui_creturn (.true.)
  Call prompt_list ('ACCESS', .true., faccs, nacc, cformat%iacc, &
       iwdth=13, itab=0)
  Call gui_creturn (.true.)
  Call prompt_list ('DATA PRECISION', .true., cprcs, nprc, cformat%iprc, &
       iwdth=13, itab=0, igrey=ipg)
  Call gui_creturn (.true.)
  Call print_text ('The following format character will be used to delimit items in formatted output files:')
  Call gui_creturn (.true.)
  Call prompt_list ('DELIMITER', .true., cdel, ndel, icd, &
       iwdth=13, itab=0, cbf1=cf_delim)
  get_fmt = win_prompt (.true.)
!
! Update settings
  If (get_fmt /= 1) cformat = tformat
  get_fmt = 2
!
  Return
 End Function get_fmt
!
!
!
 Function cf_ifmt()
!
! Sets file format identifier and corresponding precision, if appropriate
!
! Modules
  Use data_io_constants, Only: ffmts
  Use screen,            Only: window_update
!
! Function type
  Integer :: cf_ifmt
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Set gray data precision flag
  If (Trim(ffmts(cformat%ifmt)) == 'unformatted') Then
     ipg = 1
  Else If (Trim(ffmts(cformat%ifmt)) == 'formatted') Then
     cformat%iprc = 2
     ipg = 0
  Else If (Trim(ffmts(cformat%ifmt)) == 'GrADS') Then
     cformat%iprc = 1
     ipg = 0
  End If
  Call window_update (cformat%iprc)
  Call window_update (cformat%ifmt)
  cf_ifmt = 2
!
  Return
 End Function cf_ifmt
!
!
!
#else
 Function output_format()
!
! Prompts for output format
!
! Modules
  Use data_io_constants, Only: ffmts, nfmt
  Use gui,               Only: ganged_rbs, print_text
!
! Function type
  INTEGER :: output_format
!
! Executable Statements
!
! Get file format
  Call print_text ('Output format:')
  Call ganged_rbs (nfmt, cformat%ifmt, ffmts(:))
!
! Force single precision for GrADS output
  If (cformat%ifmt == 3) cformat%iprc = 1
  output_format = 1
!
  Return
 End Function output_format
!
!
!
 Function output_precision()
!
! Prompts for output precision
!
! Modules
  Use data_io_constants, Only: cprcs, nprc
  Use gui,               Only: ganged_rbs, print_text
!
! Function type
  Integer :: output_precision
!
! Executable Statements
!
! Get file precision
  Call print_text ('Output data precision:')
  Call ganged_rbs (nprc, cformat%iprc, cprcs(:))
  output_precision = 1
!
  Return
 End Function output_precision
!
!
!
#endif
 Subroutine set_fmt (nrec, ffmt)
!
! Sets file format
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nrec ! - number of records in one line -
!
! Input/output scalars
  Type(fformat), Intent(InOut) :: ffmt ! - file format -
!
! Executable Statements
!
! Apply file format settings
  ffmt%ifmt = cformat%ifmt
  ffmt%iacc = cformat%iacc
  ffmt%iprc = cformat%iprc
  ffmt%lrec = nrec
!
  Return
 End Subroutine set_fmt
!
!
!
 Subroutine set_grey_files (ixdstr, iydstr, iaction, igeog)
!
! Sets which menu items to display or grey
!
! Modules
  Use menus, Only: im_edit, inew
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: ixdstr ! - X file structure -
  Integer, Intent(In) :: iydstr ! - Y file structure -
!
! Output scalars
  Integer, Intent(Out) :: iaction ! - ready for calculation flag -
  Integer, Intent(Out) :: igeog   ! - spatial domain flag -
!
! Executable Statements
!
! Set file open and edit menu flag
  im_fopen = ixdstr + iydstr
  If (im_fopen > 0) Then
     inew = 0
     im_fopen = 1
     im_edit = 1
  Else
     inew = 1
     im_fopen = 0
     im_edit = 0
  End If
  If ((ixdstr == 1 .or. ixdstr == 2) .or. (iydstr == 1 .or. iydstr == 2)) Then
     igeog = 1
  Else
     igeog = 0
  End If
  If (ixdstr > 0 .and. iydstr > 0) Then
     iaction = 1
  Else
     iaction = 0
  End If
  im_sat = im_fopen
  im_sax = ixdstr
  im_say = iydstr
!
  Return
 End Subroutine set_grey_files
!
!
!
 Subroutine proj_read_iofiles (iin, pver, iold, ifail)
!
! Reads input file settings from project file
!
! Modules
  Use data_numbers,        Only: sp, tolv
  Use data_time_constants, Only: isq_mn
  Use time,                Only: iseq
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iin ! - input file unit number -
!
  Real(Kind=sp), Intent(In) :: pver ! - project file version -
!
! Output scalars
  Integer, Intent(Out) :: iold  ! - old monthly input file? -
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Locals
!
! Local scalars
  Type(ifile_v14) :: afile_v14 ! - version 14 file -
  Type(ifile_v15) :: afile_v15 ! - version 15 file -
!
! Executable Statements
!
! Read input file settings
  iold = 0
  If (pver > 16.00_sp) Then
     Read (Unit=iin, Err=1, End=2) xfile
     Read (Unit=iin, Err=1, End=2) yfile
     Read (Unit=iin, Err=1, End=2) zfile
     Read (Unit=iin, Err=1, End=2) mfile
     Read (Unit=iin, Err=1, End=2) icd, itags, mxa, mya, mza, nx, ny, nz
  Else If (pver > 15.03_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) afile_v15
     Call ifile_v15_to_v16 (afile_v15, xfile)
     Read (Unit=iin, Err=1, End=2) afile_v15
     Call ifile_v15_to_v16 (afile_v15, yfile)
     Read (Unit=iin, Err=1, End=2) afile_v15
     Call ifile_v15_to_v16 (afile_v15, zfile)
     Read (Unit=iin, Err=1, End=2) afile_v15
     Call ifile_v15_to_v16 (afile_v15, mfile)
     Read (Unit=iin, Err=1, End=2) icd, itags
  Else If (pver > 15.0_sp) Then
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, xfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, yfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, zfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, mfile)
     Read (Unit=iin, Err=1, End=2) icd, itags
  Else If (pver > 14.02_sp-tolv) Then
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, xfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, yfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, zfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, mfile)
     Read (Unit=iin, Err=1, End=2) icd, itags
  Else
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, xfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, yfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, zfile)
     Read (Unit=iin, Err=1, End=2) afile_v14
     Call ifile_v14_to_v16 (afile_v14, mfile)
     Read (Unit=iin, Err=1, End=2) icd
     itags=1
  End If
  If (xfile%nem == 1) xfile%nem = 0
  If (yfile%nem == 1) yfile%nem = 0
  If (zfile%nem == 1) zfile%nem = 0
  If (yfile%cfile%lset) iseq = yfile%iseq
!
! Reopen monthly and subseasonal input files if project version is pre-16.
  If (pver < 16.0_sp) Then
     If (mfile%iseq == isq_mn .or. mfile%iseq > 0) iold = mfile%iseq
  End If
!
! Reset maximum file widths
  If (xfile%cfile%lset .and. xfile%idstr > 1) Then
     Call max_width (iin+1, xfile%lmax, ifail, &
          afile=xfile%cfile%ffile)
     Select Case (ifail)
      Case (0)
        Continue
      Case (1)
        GoTo 3
      Case Default
        GoTo 5
     End Select
  Else
     xfile%lmax = 0
  End If
  If (yfile%cfile%lset .and. yfile%idstr > 1) Then
     Call max_width (iin+1, yfile%lmax, ifail, &
          afile=yfile%cfile%ffile)
     Select Case (ifail)
      Case (0)
        Continue
      Case (1)
        GoTo 3
      Case Default
        GoTo 5
     End Select
  Else
     yfile%lmax = 0
  End If
  If (zfile%cfile%lset .and. zfile%idstr > 1) Then
     Call max_width (iin+1, zfile%lmax, ifail, &
          afile=zfile%cfile%ffile)
     Select Case (ifail)
      Case (0)
        Continue
      Case (1)
        GoTo 3
      Case Default
        GoTo 5
     End Select
  Else
     zfile%lmax = 0
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
! - problem determining file width -
3 ifail = 3
  Return
!
! - problem opening file -
5 ifail = 5
  Return
!
 End Subroutine proj_read_iofiles
!
!
!
 Subroutine proj_write_iofiles (iout, ifail)
!
! Writes input file settings to project file
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
! Write input file settings
  Write (Unit=iout, Err=1) xfile
  Write (Unit=iout, Err=1) yfile
  Write (Unit=iout, Err=1) zfile
  Write (Unit=iout, Err=1) mfile
  Write (Unit=iout, Err=1) icd, itags, mxa, mya, mza, nx, ny, nz
!
  ifail = 0
  Return
!
! Errors in project file
! - problem writing file -
1 ifail = 1
  Return
!
 End Subroutine proj_write_iofiles
End Module iofiles
