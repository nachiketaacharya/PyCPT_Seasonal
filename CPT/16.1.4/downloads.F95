! Author: Simon Mason
Module downloads
!
! Implicit declarations
  Implicit None
!
! Declarations
!
! Accessibility
  Private
  Public :: download_data, init_download
!
! Scalars
!
! Integer scalars
  Integer, Public :: iddt = 0    ! - download data-type index -
  Integer, Public :: igrib = 0   ! - download variable GRIB index -
  Integer, Public :: isource = 0 ! - download dataset identifier -
!
  Integer, Private :: imd0  ! - month of forecast initialization -
  Integer, Private :: imd1  ! - first month of data to download -
  Integer, Private :: imd2  ! - last month of data to download -
  Integer, Private :: inino ! - NINO index -
  Integer, Private :: iyd1  ! - first year of data to download -
  Integer, Private :: iyd2  ! - last year of data to download -
  Integer, Private :: lsd   ! - length of season to average -
  Integer, Private :: nld   ! - number of lags -
  Integer, Private :: nvar  ! - number of variables -
!
  Logical, Private :: ldset   ! - domain set? -
  Logical, Private :: lnmme   ! - NMME dataset? -
  Logical, Private :: lyrange ! - Y-range?
!
! Character scalars
  Character(Len=  16), Private :: cg_data   ! - dataset -
  Character(Len=  16), Private :: csprev    ! - previous source -
  Character(Len=  64), Private :: url_dl    ! - IRI DL / portal home page -
  Character(Len=1024), Private :: url_dfile ! - download file URL -
!
Contains
!
!
 Subroutine init_download ()
!
! Initialises data download parameters
!
! Modules
  Use time, Only: dnow
!
! Executable Statements
!
! Initialise default dates
  iyd1 = 1
  iyd2 = dnow%iyr
  imd0 = 1
  imd1 = 1
  lsd = 1
  nld = 1
!
! Initialise indices
  inino = 3 ! - variable identifiers -
!
! Initialise source
  csprev = ' '
!
  Return
 End Subroutine init_download
!
!
!
 Function download_data()
!
! Provides access to the IRI Data Library
!
! Modules
  Use data_iridl, Only: id_gcm, id_obs, id_rea, &
                        ig_enso, ig_gphs, ig_mslp, ig_ssts, ig_tave, ig_tmax, ig_tmin, ig_prcp, ig_uvwd, ig_uwnd, ig_vwnd, ig_wetd
  Use labels,     Only: cg_nino_t
  Use time,       Only: pdate
  Use gui,        Only: use_url
#if GUI == 0
  Use iofiles,    Only: ofile, &
                        get_new_file
#endif
  Use fields,     Only: dfield=>tfield, &
                        init_field
  Use missing,    Only: xmiss
!
! Function type
  Integer :: download_data
!
! Locals
!
! Local scalars
  Integer :: ivar    ! - variable indicator -
  Integer :: iymn    ! - first year available -
  Integer :: iymx    ! - last year available -
#if GUI == 0
  Integer :: ierr    ! - error indicator -
#endif
  Integer :: ifail   ! - error indicator -
  Integer :: ndats   ! - number of datasets -
  Integer :: nem     ! - number of ensemble members -
  Integer :: nplv    ! - number of pressure levels -
  Integer :: nsource ! - dataset index -
!
  Character(Len= 1) :: cpgrd   ! - pressure level grid code -
  Character(Len=64) :: url_src ! - source URL -
!
  Logical :: lsinv          ! - date is in real-time variable URL? -
#if GUI == 1
  Logical :: ltsv = .false. ! - point to .tsv file?
#else
  Logical :: ltsv = .true.  ! - point to .tsv file?
#endif
!
  Type(pdate) :: ds2 ! - first date of second dataset -
!
#if GUI == 0
  Type(ofile) :: rdfile ! - renamed download file -
!
#endif
! Local arrays
  Integer, Dimension(:), Allocatable :: ivars ! - variables -
  Integer, Dimension(:), Allocatable :: iplvs ! - pressure level selections -
!
  Character(Len=  64), Dimension(:), Allocatable :: ccvt  ! - unit conversions -
  Character(Len=  64), Dimension(:), Allocatable :: cplvs ! - pressure levels -
  Character(Len=  16), Dimension(:), Allocatable :: cvars ! - variables -
!
  Character(Len= 256), Dimension(:,:), Allocatable :: url_vars ! - variables URLs -
!
  Logical, Dimension(:), Allocatable :: llsum ! - sum rather than average? -
!
! Functions and Subroutines
!
! Instrinsic functions
  Intrinsic Trim
#if GUI == 0
#if FTN95 == 1
  Integer, External :: start_process@
#else
  Intrinsic execute_command_line
#endif
#endif
!
! Executable Statements
!
! Initialise field
  download_data = 2
  Call init_field (dfield, 1, (/xmiss/), ifail)
!
! Read dataset configuration
  Call get_dataset (iddt, igrib, isource, ifail)
  If (ifail /= 0) GoTo 1
!
! Identify download file
  ivar = igrib
  Select Case (iddt)
! - observed data -
   Case (id_obs)
     Select Case (igrib)
      Case (ig_enso) ! - ENSO indices -
        Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, inino, &
             .false., .false., .false., ltsv, .false., .false., ifail,              &
             cbox=Trim(cg_nino_t), nvar=nvar, cvars=cvars)
      Case (ig_ssts, ig_tave, ig_tmax, ig_tmin, ig_prcp, ig_wetd) ! - SSTs, T mean, T max, T min, precipitation, wet-days -
        Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
             .false., .false., lsinv, ltsv, lyrange, .true., ifail,                &
             llsum=llsum, dfield=dfield)
     End Select
! - reanalysis data -
   Case (id_rea)
     Select Case (igrib)
      Case (ig_mslp) ! - MSLP -
        Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
             .false., .false., lsinv, ltsv, lyrange, .true., ifail,                &
             dfield=dfield)
      Case (ig_uvwd) ! - U and V by pressure level -
        Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
             .false., .true., lsinv, ltsv, lyrange, .true., ifail,                 &
             nplv=nplv, iplvs=iplvs, cplvs=cplvs, cpgrd=cpgrd, dfield=dfield)
      Case (ig_gphs, ig_uwnd, ig_vwnd) ! - other variables by pressure level -
        Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
             .false., .false., lsinv, ltsv, lyrange, .true., ifail,                &
             nplv=nplv, iplvs=iplvs, cplvs=cplvs, cpgrd=cpgrd, dfield=dfield)
     End Select
! - model output -
   Case (id_gcm)
     Select Case (ndats)
      Case (1) ! - models with a single dataset -
        Select Case (igrib)
         Case (ig_tave, ig_tmax, ig_tmin, ig_prcp, ig_ssts) ! - surface and near-surface variables -
           Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
                .true., .true., lsinv, ltsv, lyrange, .true., ifail,                  &
                ccvt=ccvt, nem=nem, dfield=dfield)
         Case (ig_gphs, ig_uwnd, ig_vwnd) ! - other variables by pressure level -
           Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
                .true., .true., lsinv, ltsv, lyrange, .true., ifail,                  &
                ccvt=ccvt, nem=nem, nplv=nplv, iplvs=iplvs, cplvs=cplvs, cpgrd=cpgrd, dfield=dfield)
        End Select
      Case (2) ! - models with separate hindcast and realtime datasets -
        Select Case (igrib)
         Case (ig_tave, ig_tmax, ig_tmin, ig_prcp, ig_ssts) ! - surface and near-surface variables -
           Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
                .true., .true., lsinv, ltsv, lyrange, .true., ifail,                  &
                ccvt=ccvt, nem=nem, ds2=ds2, dfield=dfield)
         Case (ig_gphs, ig_uwnd, ig_vwnd) ! - other variables by pressure level -
           Call get_download_url (Trim(cg_data), iymn, iymx, url_src, url_vars, ivar, &
                .true., .true., lsinv, ltsv, lyrange, .true., ifail,                  &
                ccvt=ccvt, nem=nem, ds2=ds2, nplv=nplv, iplvs=iplvs, cplvs=cplvs, cpgrd=cpgrd, dfield=dfield)
        End Select
     End Select
  End Select
!
! Open browser
  If (ifail == 0) Call use_url (Trim(url_dfile))
!
#if GUI == 0
! Download in batch versions
  If (ltsv) Then
     rdfile%cfile%fext = '.tsv'
! - prompt for download file -
     Call get_new_file ('Download file', .true., rdfile%cfile%ffile, rdfile%cfile%fdir, rdfile%cfile%fext, &
          (/'Download Files'/), (/rdfile%cfile%fext/), 1, rdfile%cfile%fname, ifail)
     If (ifail == 0) Then
! - replace special characters -
        Call encode_url (url_dfile)
#if FTN95 == 1
! - replace special characters -
        Call encode_url (rdfile%cfile%ffile)
        ifail = start_process@('bitsadmin', &
             '/transfer CPT_DL_Download /download /priority normal '//Trim(url_dfile)//' '//Trim(rdfile%cfile%ffile))
        ierr = ifail
#else
        Call execute_command_line ('wget '''//Trim(url_dfile)//''' -O '''//Trim(rdfile%cfile%ffile)//'''', &
             cmdstat=ifail, exitstat=ierr)
#endif
        If (ifail /= 0 .or. ierr /= 0) Print*, 'Try copying the URL into a browser and downloading manually.'
     End If
  End If
#endif
!
! Reset identifiers
1 iddt = 0
  igrib = 0
  isource = 0
!
! Free memory
  Call download_mem_dealloc ()
!
  Return
!
 Contains
!
!
  Subroutine get_dataset (iddt, ivar, isource, ifail)
!
! Reads dataset settings
!
! Modules
  Use data_io_constants, Only: ddir, iin
  Use data_iridl,        Only: c_ddt, dfile
  Use time,              Only: dnow
  Use gui,               Only: iw, &
                               init_win, prompt_list, win_prompt
  Use errors,            Only: cpt_error
  Use iofiles,           Only: open_iofile
  Use fields,            Only: dfield=>tfield
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iddt ! - download data-type index -
!
! Input/output scalars
  Integer, Intent(InOut) :: ivar ! - required variable -
!
! Output scalars
  Integer, Intent(Out) :: isource ! - download dataset identifier -
  Integer, Intent(Out) :: ifail   ! - error indicator -
!
! Locals
!
! Local scalars
  Integer :: i    ! - index -
  Integer :: j    ! - dataset index -
  Integer :: k    ! - variable index -
  Integer :: isrc ! - source index -
!
  Character(Len=  5) :: cds     ! - dataset structure -
  Character(Len=256) :: cline   ! - input line -
  Character(Len= 32) :: cpt_tag ! - CPT tag -
!
  Logical :: lpwd   ! - password protected? -
  Logical :: lpwdr  ! - password required? -
  Logical :: lnmmem ! - NMME model? -
!
! Local arrays
  Integer, Dimension(:,:), Allocatable :: isrcs ! - source identifiers -
!
  Character(Len=32), Dimension(:), Allocatable :: csrcs ! - sources -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Len, Trim
!
! Executable Statements
!
! Open dataset settings file
  Call open_iofile (iin, Trim(ddir)//dfile, .true., .true., .false., ifail)
  If (ifail /= 0) GoTo 3
  isource = 0
!
! Search for required data type
4 ifail = 6
  Do
     Read (Unit=iin, Fmt='(A)', Err=1, End=3) cline
     i = Index(cline, 'cpt:ddt='//Trim(c_ddt(iddt)))
     If (i > 0) Exit
  End Do
!
! Search for required dataset
! - construct tag -
  Write (Unit=cpt_tag, Fmt='(A,I0)') 'cpt:ivar=', ivar
! - search for tag -
  ifail = 7
  Do
     Read (Unit=iin, Fmt='(A)', Err=1, End=3) cline
     i = Index(cline, Trim(cpt_tag))
     If (i > 0) Exit
  End Do
!
! Identify number of data sources
  i = Index(cline, 'cpt:nsource=')
  If (i == 0) Then
     ifail = 8
     GoTo 3
  End If
  i = i + Len('cpt:nsource=')
  Read (Unit=cline(i:), Fmt=*) nsource
!
! Prompt for source
  If (nsource > 1) Then
     If (isource == 0) Then
        Allocate (csrcs(nsource))
        Allocate (isrcs(nsource,2))
! - identify sources -
        ifail = 9
        Do isrc = 1, nsource
           Do
              Read (Unit=iin, Fmt='(A)', Err=1, End=3) cline
              i = Index(cline, 'cpt:isource=')
              If (i > 0) Then
                 i = i + Len('cpt:isource=')
                 Read (Unit=cline(i:), Fmt=*) isrcs(isrc,1)
                 i = Index(cline, 'cpt:ioldsource=') ! - old source file identifier -
                 If (i > 0) Then
                    i = i + Len('cpt:ioldsource=')
                    Read (Unit=cline(i:), Fmt=*) isrcs(isrc,2)
                 Else
                    isrcs(isrc,2) = isrcs(isrc,1)
                 End If
                 i = Index(cline, 'cpt:csource=')
                 If (i == 0) GoTo 3
                 i = i + Len('cpt:csource=')
                 Read (Unit=cline(i:), Fmt=*) csrcs(isrc)
                 Read (Unit=iin, Fmt=*, Err=1, End=2) lpwdr ! - password protected? -
                 If (iddt == id_gcm) Then
                    Read (Unit=iin, Fmt=*)
                    Read (Unit=iin, Fmt=*, Err=1, End=2) lnmmem ! - NMME? -
                    If (lnmmem) csrcs(isrc) = Trim(csrcs(isrc))//' (NMME)'
                 End If
                 If (lpwdr) csrcs(isrc) = Trim(csrcs(isrc))//' (password)'
                 Exit
              End If
           End Do
        End Do
! - prompt for source -
        isrc = 1
        Call init_win (title='Data sources')
        Call prompt_list ('Data sources', .false., csrcs, nsource, isrc, &
             ns=isrcs(:,1), altns=isrcs(:,2))
        iw = win_prompt(.true.)
        If (iw == 1 .and. isrc /= 0) Then
           isource = isrcs(isrc,1)
        Else
           isource = 0
        End If
        Deallocate (isrcs)
        Deallocate (csrcs)
! - rewind to relevant data source -
        If (isource > 0) Then
           Rewind (Unit=iin)
           GoTo 4
        Else
           ifail = -1
           GoTo 5
        End If
! - search for requested dataset -
     Else
        Write (Unit=cpt_tag, Fmt='(A,I0)') 'cpt:isource=', isource
        Do
           Read (Unit=iin, Fmt='(A)', Err=1, End=2) cline
           i = Index(cline, Trim(cpt_tag))
           If (i > 0) Exit
        End Do
     End If
  Else
     Read (Unit=iin, Fmt='(A)', Err=1, End=2) cline
  End If
!
! Identify dataset
  i = Index(cline, 'cpt:csource=')
  If (i == 0) Then
     ifail = 9
     GoTo 3
  End If
  i = i + Len('cpt:csource=')
  Read (Unit=cline(i:), Fmt=*) cg_data                           ! - dataset -
!
! Determine structure
  Read (Unit=iin, Fmt=*, Err=1, End=2) lpwd                      ! - password protected? -
  Read (Unit=iin, Fmt=*, Err=1, End=2) cds                       ! - dataset structure -
  If (iddt == id_gcm) Read (Unit=iin, Fmt=*, Err=1, End=2) lnmme ! - first date of second dataset -
  Read (Unit=iin, Fmt=*, Err=1, End=2) ndats                     ! - number of datasets -
  Read (Unit=iin, Fmt=*, Err=1, End=2) iymn                      ! - first year of data -
  iymx = dnow%iyr
  If (ndats > 1) Read (Unit=iin, Fmt=*, Err=1, End=2) ds2        ! - first date of second dataset -
!
! Read appropriate structure
! - index data -
  Select Case (cds)
   Case ('index')
     Read (Unit=iin, Fmt=*, Err=1, End=2) nvar                       ! - number of indices -
     nplv = 0
! - grid data -
   Case ('grid')
     Read (Unit=iin, Fmt=*, Err=1, End=2) lsinv                      ! - date is in real-time variable URL? -
     If (iddt == 3) Read (Unit=iin, Fmt=*, Err=1, End=2) nem         ! - number of ensemble members -
     Read (Unit=iin, Fmt=*, Err=1, End=2) nplv                       ! - number of pressure levels -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%nlt              ! - number of latitudes -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%nlg              ! - number of longitudes -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%region%alim%rltn ! - northern limit -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%region%alim%rlts ! - southern limit -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%region%alim%rlgw ! - western limit -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%region%alim%rlge ! - eastern limit -
     Read (Unit=iin, Fmt=*, Err=1, End=2) dfield(1)%ln2s             ! - north-south orientation? -
     Read (Unit=iin, Fmt=*, Err=1, End=2) lyrange                    ! - Y-range? -
     nvar = 1
  End Select
!
! Allocate memory
  Call download_mem_alloc (nvar, ndats, nplv)
!
! Read variable controls
  Select Case (cds)
   Case ('grid ')
! - summed variables? -
     Do k = 1, nvar
        Read (Unit=iin, Fmt=*, Err=1, End=2) llsum(k)
     End Do
! - unit conversions -
     Do k = 1, nvar
        Read (Unit=iin, Fmt='(A)', Err=1, End=2) cline
        i = Index(cline, ' ')
        Read (Unit=cline(:i), Fmt='(A)', Err=1, End=2) ccvt(k)
     End Do
! - variables -
   Case ('index')
     Do k = 1, nvar
        Read (Unit=iin, Fmt=*, Err=1, End=2) cvars(k)
     End Do
  End Select
!
! Read URLs
! - dataset portal -
  Read (Unit=iin, Fmt='(A)', Err=1, End=2) cline
  i = Index(cline, ' ')
  Read (Unit=cline(:i), Fmt='(A)', Err=1, End=2) url_dl
! - source URL -
  Read (Unit=iin, Fmt='(A)', Err=1, End=2) cline
  i = Index(cline, ' ')
  Read (Unit=cline(:i), Fmt='(A)', Err=1, End=2) url_src
! - variables URLs -
  Do j = 1, ndats
     Do k = 1, nvar
        Read (Unit=iin, Fmt='(A)', Err=1, End=2) cline
        i = Index(cline, ' ')
        Read (Unit=cline(:i), Fmt='(A)', Err=1, End=2) url_vars(k,j)
     End Do
  End Do
!
! Read pressure levels
  Select Case (cds)
   Case ('grid')
     If (nplv > 0) Read (Unit=iin, Fmt=*, Err=1, End=2) cpgrd, (cplvs(i), i=1,nplv)
  End Select
  Close (Unit=iin)
!
  ifail = 0
  Return
!
! Errors
1 ifail = 4
  GoTo 3
!
2 ifail = 5
!
3 Call cpt_error ('get_dataset', ifail, .false.)
!
! Close dataset file
5 Close (Unit=iin)
  Return
!
  End Subroutine get_dataset
!
!
!
  Subroutine download_mem_alloc (nvar, ndats, nplv)
!
! Allocates memory for reading download parameters
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nvar  ! - number of variables -
  Integer, Intent(In) :: ndats ! - number of datasets -
  Integer, Intent(In) :: nplv  ! - number of pressure levels -
!
! Executable Statements
!
! Allocate memory
  Allocate (ivars(nvar))
  Allocate (url_vars(nvar,ndats))
  Allocate (ccvt(nvar))
  Allocate (cvars(nvar))
  Allocate (llsum(nvar))
  If (nplv > 0) Then
     Allocate (iplvs(nplv))
     iplvs(:) = 0
     Allocate (cplvs(nplv))
  End If
!
  Return
  End Subroutine download_mem_alloc
!
!
!
  Subroutine download_mem_dealloc ()
!
! Frees memory for reading download parameters
!
! Functions and Subroutines
!
! Instrinsic functions
  Intrinsic Allocated
!
! Executable Statements
!
! Free memory
  If (Allocated(cplvs))    Deallocate (cplvs)
  If (Allocated(iplvs))    Deallocate (iplvs)
  If (Allocated(llsum))    Deallocate (llsum)
  If (Allocated(cvars))    Deallocate (cvars)
  If (Allocated(ccvt))     Deallocate (ccvt)
  If (Allocated(url_vars)) Deallocate (url_vars)
  If (Allocated(ivars))    Deallocate (ivars)
!
  Return
  End Subroutine download_mem_dealloc
!
!
!
  Subroutine encode_url (url)
!
! Replaces non-ASCII characters in URL
!
! Modules
  Use data_iridl, Only: c_invalid, c_replace, nic
!
! Arguments
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: url ! - URL -
!
! Locals
!
! Local scalars
  Integer :: i ! - invalid character index -
  Integer :: j ! - invalid character location -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index, Trim
!
! Executable Statements
!
! Replace left brackets
  Do i = 1, nic
     search_ic: Do
        j = Index(Trim(url), c_invalid(i))
        If (j == 0) Exit search_ic
        url = url(1:j-1)//c_replace(i)//Trim(url(j+1:))
     End Do search_ic
  End Do
!
  Return
  End Subroutine encode_url
 End Function download_data
!
!
!
 Subroutine get_download_url (cg_data, iymn, iymx, url_src, url_vars, ivar, ls0, lfv10, lsinv, ltsv, lyrange, lgrid, ifail, &
            llsum, ccvt, cbox, nvar, cvars, nplv, iplvs, cplvs, cpgrd, nem, ds2, dfield)
!
! Constructs URL for the download file
!
! Modules
  Use data_numbers,        Only: rp, r360
  Use data_time_constants, Only: nmn
  Use data_iridl,          Only: dext, url_dcpt, url_grid, url_index, url_scpt10, url_sources, url_tcpt10, url_yrange
  Use labels,              Only: l_download
  Use time,                Only: pdate, pprd, &
                                 Operator(>=), Operator(<), Operator(>), &
                                 dnow, &
                                 add_to_month
  Use gui,                 Only: iw, &
                                 box_close, box_open, boxes_open, init_win, print_text, prompt_list, prompt_select, win_prompt
  Use errors,              Only: cpt_error
  Use settings,            Only: set_options
  Use fields,              Only: field
  Use season,              Only: cssel, &
                                 get_cssel
  Use missing,             Only: xmiss
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iymn ! - first year available -
  Integer, Intent(In) :: iymx ! - last year available -
!
  Character(Len=*), Intent(In) :: cg_data ! - data -
  Character(Len=*), Intent(In) :: url_src ! - source URL -
!
  Logical, Intent(In) :: lfv10   ! - format in expert mode for v10? -
  Logical, Intent(In) :: ls0     ! - prompt for start date? -
  Logical, Intent(In) :: lsinv   ! - date is in real-time variable URL? -
  Logical, Intent(In) :: ltsv    ! - point to .tsv file?
  Logical, Intent(In) :: lyrange ! - Y-range?
  Logical, Intent(In) :: lgrid   ! - gridded data?
!
! - optional input scalars -
  Integer, Intent(In), Optional :: nvar ! - number of variables -
  Integer, Intent(In), Optional :: nplv ! - number of pressure levels -
  Integer, Intent(In), Optional :: nem  ! - number of ensemble members -
!
  Character(Len=*), Intent(In), Optional :: cbox  ! - box heading for variable prompt -
  Character(Len=*), Intent(In), Optional :: cpgrd ! - pressure level grid code -
!
! Input/output scalars
  Integer, Intent(InOut) :: ivar ! - variable index -
!
! - optional input/output scalars -
  Type(pdate), Intent(InOut), Optional :: ds2 ! - first date of second dataset -
!
! Input arrays
! - optional input arrays -
  Character(Len=*), Dimension(:), Intent(In), Optional :: ccvt  ! - unit conversions -
  Character(Len=*), Dimension(:), Intent(In), Optional :: cvars ! - variable -
  Character(Len=*), Dimension(:), Intent(In), Optional :: cplvs ! - pressure levels -
!
  Logical, Dimension(:), Intent(In), Optional :: llsum ! - sum rather than average? -
!
! Input/output scalars
  Character(Len=*), Dimension(:,:), Intent(InOut) :: url_vars ! - variable URL -
!
! Input/output arrays
! - optional input/output scalars -
  Integer, Dimension(:), Intent(InOut), Optional :: iplvs ! - pressure level selections -
!
! Output scalars
  Integer, Intent(Out) :: ifail ! - error indicator -
!
! Pointer arrays
! - optional pointer arrays -
  Type(field), Dimension(:), Pointer, Optional :: dfield ! - field -
!
! Locals
!
! Local scalars
  Integer :: i   ! - latitude/longitude index -
  Integer :: ids ! - dataset identifier -
  Integer :: iyd ! - year difference -
  Integer :: np  ! - number of prompts -
!
  Real(Kind=rp) :: d ! - grid resolution -
!
  Logical :: ldom   ! - prompt for domain? -
  Logical :: lplv   ! - prompt for pressure level? -
  Logical :: lsprev ! - previous source? -
  Logical :: lsum   ! - sum rather than average? -
  Logical :: lvar   ! - prompt for variable? -
!
  Type(pdate) :: d1 ! - start date -
  Type( pprd) :: d2 ! - end period -
!
! Functions and Subroutines
!
! Instrinsic functions
  Intrinsic Count, Present, Sum, Trim
!
! Executable Statements
!
! Backup
  ifail = set_options(i1=iyd1, i2=iyd2, i3=imd0, i4=imd1, i5=lsd, i6=nld, i7=ivar)
!
! Identify which prompts to use
  ldom = Present(dfield) ! - domain -
  lvar = Present(nvar)   ! - variable -
  lplv = Present(nplv)   ! - pressure level -
  np = Count((/lvar, lplv/))
!
! Initialise domain
  If (ldom) Call init_domain (dfield(:))
!
! Initialise season
  If (.not.ls0) imd0 = 0
  i = imd1 - 1
  If (imd0 > 0) Then
     cssel = get_cssel(imd0, imd1, lsd, 1, 0)
  Else
     cssel = get_cssel(imd0, imd1, lsd, nld, 0)
  End If
!
! Identify whether source is the same as previous
  lsprev = (Trim(csprev) == cg_data)
!
! Reset defaults if outside of range
  If (iyd1 < iymn .or. .not.lsprev) iyd1 = iymn
  If (iyd2 < iyd1) iyd2 = iyd1
!
! Prompt for download parameters
1 Call init_win (title=Trim(l_download%c)//' - '//Trim(cg_data))
  If (.not.ls0) imd0 = 0
! - variable -
  If (np > 0) Then
     Call boxes_open (2, 1, .true., .false.)
     If (np > 1) Call boxes_open (1, 2, .true., .false.)
  End If
  If (lvar) Then
     Call box_open (.false., .false., &
          title=cbox)
     Call prompt_list ('Variable', .true., cvars, nvar, ivar)
     Call box_close (.false.)
     If (np > 1) Call box_close (.false.)
     Call print_text (' ')
  Else
     ivar = 1
  End If
! - pressure level -
  If (lplv) Then
     Call box_open (.false., .false., &
          title='Pressure levels')
     Call prompt_select ('Pressure levels', nplv, iplvs, cplvs)
     Call box_close (.false.)
     If (np > 1) Call box_close (.false.)
  End If
  If (np > 0) Call box_close (.false.)
! - dates -
  Call prompt_dates (iymn,iymx,ls0)
  If (np > 0) Call box_close (.false.)
  If (ldom) Then
     iw = win_prompt(.true., &
          cextra='Select sub-domain', fextra=prompt_domain)
     ldom = ldset
  Else
     iw = win_prompt(.true.)
  End If
!
! Set option
  If (set_options(iw=iw, i1=iyd1, i2=iyd2, i3=imd0, i4=imd1, i5=lsd, i6=nld, i7=ivar, lrc=.false.) >= 3) Then
     ifail = -1
     Return
  End If
!
! Check number of selections
  If (lplv) Then
     If (Sum(iplvs(1:nplv)) == 0) Then
         If (nplv == 1) Then
            iplvs(:) = 1
         Else
            ifail = 1
            Call cpt_error ('get_download_url', ifail, .false.)
            GoTo 1
        End If
     End If
  End If
!
! Set date limits
  If (.not.ls0) imd0 = imd1
  d1%iyr = iyd1
  d1%imn = imd0
  d1%idy = 1
  d2%sdate%iyr = iyd2
  d2%sdate%imn = imd0
  d2%sdate%idy = 1
  d2%edate = d2%sdate
  If (.not.ls0) Then
     d2%edate%imn = add_to_month(d2%edate%imn, lsd*nld, &
                    iyr=d2%edate%iyr)
  End If
!
! Check for valid settings
! - last date is before first date -
  If (d2%sdate < d1) Then
     iyd2 = iyd1
     ifail = 2
     Call cpt_error ('get_download_url', ifail, .false.)
     GoTo 1
! - first date is in the future -
  Else If (d1 > dnow) Then
     If (iyd1 == dnow%iyr) Then
        iyd1 = dnow%iyr - 1
     Else
        iyd1 = dnow%iyr
     End If
     d1 = d2%sdate
     ifail = 3
     Call cpt_error ('get_download_url', ifail, .false.)
     GoTo 1
! - second date is in the future -
  Else If (d2%edate > dnow) Then
     If (d2%edate%iyr == dnow%iyr) Then
        d2%edate%iyr = dnow%iyr - 1
        d2%sdate%iyr = dnow%iyr - 1
     Else
        iyd = d2%edate%iyr - dnow%iyr
        d2%edate%iyr = dnow%iyr
        d2%sdate%iyr = d2%sdate%iyr - iyd
        If (d2%edate > dnow) Then
           d2%edate%iyr = d2%edate%iyr - 1
           d2%sdate%iyr = d2%sdate%iyr - 1
        End If
     End If
     iyd2 = d2%sdate%iyr
     If (d1 > d2%sdate) d1 = d2%sdate
     ifail = 4
     Call cpt_error ('get_download_url', ifail, .false.)
     GoTo 1
  End If
  iyd2 = d2%edate%iyr
!
! Check for need for two datasets
  If (Present(ds2)) Then
     If (ds2%imn < 1) Then ! - determine relative date -
        ds2%iyr = dnow%iyr
        ds2%imn = dnow%imn + ds2%imn
        ds2%idy = 1
        Do
           If (ds2%imn > 0) Exit
           ds2%imn = ds2%imn + nmn
           ds2%iyr = ds2%iyr - 1
        End Do
     End If
     If (d2%sdate >= ds2) Then ! - second dataset is needed -
        If (lsinv) Call put_s_in_var_url (d2%sdate, url_vars(ivar,2))
        If (d1 < ds2) Then ! - first and second datasets are needed -
           ids = 4
           If (ds2%imn < 1) ids=ids+1 ! - construct for ECHAM's incompatible time grids -
        Else ! - only the second dataset is needed -
           ids = 2
        End If
     Else ! - only the first dataset is needed -
        ids = 1
     End If
  Else
     ids = 1
  End If
!
! Determine whether to sum units
  If (Present(llsum)) Then
     lsum = llsum(ivar)
  Else
     lsum = .false.
  End If
!
! Construct URL
  url_dfile = Trim(url_dl)//url_sources//url_src
  If (ls0) Then
     Select Case (ids)
      Case (1) ! - first dataset only -
        url_dfile = Trim(url_dfile)//Trim(url_vars(ivar,ids))
        url_dfile = Trim(url_dfile)//get_url_s(iyd1,iyd2)                               ! - select start dates -
      Case (2) ! - second dataset only -
        url_dfile = Trim(url_dfile)//Trim(url_vars(ivar,ids))
        If (lsinv) url_dfile=Trim(url_dfile)//change_time_grids()                       ! - create compatible time grids -
      Case (4) ! - first and second datasets -
        url_dfile = Trim(url_dfile)//Trim(url_vars(ivar,1))
        url_dfile = Trim(url_dfile)//get_url_s(iyd1, ds2%iyr-1)                         ! - select start dates -
        url_dfile = Trim(url_dfile)//url_sources//Trim(url_src)//Trim(url_vars(ivar,2)) ! - add second dataset -
        url_dfile = Trim(url_dfile)//get_url_s(ds2%iyr, iyd2)                           ! - select start dates -
        url_dfile = Trim(url_dfile)//get_url_append()                                   ! - append datasets -
      Case (5) ! - first and second datasets -
        url_dfile = Trim(url_dfile)//Trim(url_vars(ivar,1))
        url_dfile = Trim(url_dfile)//get_url_s(iyd1, iyd2-1)                            ! - select start dates -
        url_dfile = Trim(url_dfile)//url_sources//Trim(url_src)//Trim(url_vars(ivar,2)) ! - add second dataset -
        url_dfile = Trim(url_dfile)//equate_y_grid()                                    ! - equate Y grids -
        If (lsinv) url_dfile = Trim(url_dfile)//change_time_grids()                     ! - create compatible time grids -
        url_dfile = Trim(url_dfile)//get_url_l(1, 1, 7, .false.)                        ! - select target season by lag -
        url_dfile = Trim(url_dfile)//get_url_append()                                   ! - append datasets -
     End Select
     url_dfile = Trim(url_dfile)//get_url_l(imd0, imd1, lsd, .true.)                    ! - select target season by lag -
  Else
     url_dfile = Trim(url_dfile)//Trim(url_vars(ivar,1))                                ! - download variable -
     url_dfile = Trim(url_dfile)//get_url_t(lsum)                                       ! - select target dates -
  End If
! - add pressure level -
  If (lplv) Then
     If (nplv > 1) url_dfile = Trim(url_dfile)//get_url_p(cpgrd, nplv, iplvs, cplvs)
  End If
! - calculate ensemble mean -
  If (Present(nem)) url_dfile = Trim(url_dfile)//get_url_m(nem)
! - add domain -
  If (ldom) url_dfile = Trim(url_dfile)//get_url_xy()
! - convert units -
  If (Present(ccvt)) Then
     If (Trim(ccvt(ivar)) /= 'N/A') url_dfile = Trim(url_dfile)//get_url_convert(Trim(ccvt(ivar)))
  End If
! - set missing value flag -
  If (Present(dfield)) Then
     url_dfile = Trim(url_dfile)//get_url_missing(dfield(1)%rmiss)
  Else
     url_dfile = Trim(url_dfile)//get_url_missing(xmiss)
  End If
! - add Y-range -
  If (lyrange) url_dfile = Trim(url_dfile)//url_yrange
! - request CPT10 format -
  If ((lfv10 .or. ltsv) .and. lgrid) Then
     url_dfile = Trim(url_dfile)//url_grid
     If (ls0) Then
        url_dfile = Trim(url_dfile)//url_scpt10
     Else
        url_dfile = Trim(url_dfile)//url_tcpt10
     End If
     If (.not.ltsv) Then
        url_dfile = Trim(url_dfile)//'/'
     Else
        url_dfile = Trim(url_dfile)//dext
     End If
  Else
     If (.not.ltsv) Then
        url_dfile = Trim(url_dfile)//url_dcpt
     Else
        url_dfile = Trim(url_dfile)//url_index//url_tcpt10//dext
     End If
  End If
!
! Reset previous data source
  csprev = cg_data
  ifail = 0
!
  Return
!
 Contains
!
!
  Subroutine prompt_dates (iymn, iymx, ls0)
!
! Prompt for date download parameters
!
! Modules
#if GUI == 1
  Use data_text,           Only: ilang
  Use data_time_constants, Only: cma, cmon
#endif
  Use data_time_constants, Only: nmn
  Use labels,              Only: ca_seqs, cg_cselection
  Use gui,                 Only: box_close, box_open, boxes_open, print_cvalue, print_text, prompt_integer
  Use season,              Only: cssel, lssel
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iymn ! - minimum year -
  Integer, Intent(In) :: iymx ! - maximum year -
!
  Logical, Intent(In) :: ls0 ! - prompt for start date? -
!
! Functions and Subroutines
!
! Instrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Prompt for download date parameters
  Call box_open (.false., .false., &
       title='Date limits')
  Call boxes_open (1, 2, .true., .false.)
! - prompt for year limits -
  Call box_open (.false., .false., &
       title=Trim(ca_seqs(1)))
  Call prompt_integer ('First year of data to download', iyd1, &
       itab=60, ilow=iymn, ihgh=iymx, lfocus=.true.)
  Call prompt_integer ('Last year of data to download', iyd2, &
       itab=60, ilow=iymn, ihgh=iymx, lfocus=.true.)
  Call box_close (.false.)
  Call box_close (.false.)
! - prompt for months limits -
  Call box_open (.false., .false., &
       title=Trim(ca_seqs(2))//' / '//Trim(ca_seqs(3)))
  If (ls0) Then
     Call prompt_integer ('First month of target period', imd1, &
          itab=60, ilow=1, ihgh=nmn, lfocus=.false., fcb1=update_dssel)
  Else
     Call prompt_integer ('First month of data to download', imd1, &
          itab=60, ilow=1, ihgh=nmn, lfocus=.false., fcb1=update_dssel)
  End If
  Call prompt_integer ('Length of season to average', lsd, &
       itab=60, ilow=1, ihgh=nmn, lfocus=.false., fcb1=update_dssel)
  If (ls0) Then
     Call prompt_integer ('Month forecasts were initialised', imd0, &
          itab=60, ilow=1, ihgh=nmn, lfocus=.false., fcb1=update_dssel)
  Else
     Call prompt_integer ('Number of seasons', nld, &
          itab=60, ilow=1, ihgh=nmn, lfocus=.false., fcb1=update_dssel)
  End If
  Call print_text (' ')
  Call print_cvalue (Trim(cg_cselection), cssel, &
       iwid=lssel)
#if GUI == 1
  Call print_text (' ')
  Call print_text ('Examples:')
  If (ls0) Then
     Call print_text ('   For '//cmon(6,ilang)//' from '//cmon(5,ilang)//' data enter 6, 1, 5')
     Call print_text ('   For '//cma(ilang)(6:9)//' from '//cmon(5,ilang)//' data enter 6, 4, 5')
  Else
     Call print_text ('   For '//cmon(4,ilang)//' data enter 4, 1, 1')
     Call print_text ('   For '//cma(ilang)(4:6)//' data enter 4, 3, 1')
     Call print_text ('   For '//cmon(4,ilang)//', '//cmon(5,ilang)//', and '//cmon(6,ilang)//' data enter 4, 1, 3')
     Call print_text ('   For '//cma(ilang)(4:6)//', and '//cma(ilang)(7:9)//' data enter 4, 3, 2')
  End If
#endif
  Call box_close (.false.)
  Call box_close (.false.)
  Call box_close (.false.)
!
  Return
  End Subroutine prompt_dates
!
!
!
  Subroutine init_domain (dfield)
!
! Initialises domain download parameters prompt
!
! Modules
  Use space, Only: dlats=>rlat1, dlngs=>rlng1, cdstn=>cstn1, dlatsd=>rlatd1, dlngsd=>rlngd1, cdstnd=>cstnd1, idom=>idom1, &
                   reset_grids
!
! Arguments
!
! Input arrays
  Type(field), Dimension(:), Intent(In) :: dfield ! - field settings -
!
! Functions and Subroutines
!
! Instrinsic functions
  Intrinsic Abs, Real
!
! Executable Statements
!
! Initialise domain arrays
  ldset = .false.
  Call reset_grids (1, 1, dfield(:)%nlt, dlats, dlngs, cdstn, dlatsd, dlngsd, cdstnd, idom, ifail, &
                    nlg=dfield(:)%nlg)
!
! Set latitudes
  dlats(1,1) = dfield(1)%region%alim%rltn
  If (dfield(1)%nlt > 2) Then
     d = Abs(dfield(1)%region%alim%rltn - dfield(1)%region%alim%rlts)/Real(dfield(1)%nlt-1, Kind=rp)
     Do i = 2, dfield(1)%nlt - 1
        dlats(i,1) = dlats(i-1,1) - d
     End Do
  End If
  dlats(dfield(1)%nlt,1) = dfield(1)%region%alim%rlts
!
! Set longitudes
  dlngs(1,1) = dfield(1)%region%alim%rlgw
  If (dfield(1)%nlg > 2) Then
     If (dfield(1)%region%alim%rlge > dfield(1)%region%alim%rlgw) Then
        d = (dfield(1)%region%alim%rlge - dfield(1)%region%alim%rlgw)/Real(dfield(1)%nlg-1, Kind=rp)
     Else
        d = (r360+dfield(1)%region%alim%rlge-dfield(1)%region%alim%rlgw)/Real(dfield(1)%nlg-1, Kind=rp)
     End If
     Do i = 2, dfield(1)%nlg-1
        dlngs(i,1) = dlngs(i-1,1) + d
     End Do
  End If
  dlngs(dfield(1)%nlg,1) = dfield(1)%region%alim%rlge
!
  Return
  End Subroutine init_domain
 End Subroutine get_download_url
!
!
!
 Function update_dssel()
!
! Updates download season selection indicator
!
! Modules
  Use screen, Only: window_update
  Use season, Only: cssel, &
                    get_cssel
!
! Function type
  Integer :: update_dssel
!
! Executable Statements
!
! Update season selection indicator
  If (imd0 > 0) Then
     cssel = get_cssel(imd0, imd1, lsd, 1, 0)
  Else
     cssel = get_cssel(imd0, imd1, lsd, nld, 0)
  End If
  Call window_update (cssel)
  update_dssel = 2
!
  Return
 End Function update_dssel
!
!
!
 Function prompt_domain()
!
! Prompts for domain settings
!
! Modules
  Use space,  Only: dlats=>rlat1, dlngs=>rlng1, cdstn=>cstn1, dlatsd=>rlatd1, dlngsd=>rlngd1, cdstnd=>cstnd1, idom=>idom1
  Use fields, Only: dfield=>tfield, &
                    get_area
!
! Function type
  Integer :: prompt_domain
!
! Locals
!
! Local scalars
  Integer :: nv    ! - number of variables -
  Integer :: ifail ! - error indicator -
!
! Executable Statements
!
! Prompt for domain
  Call get_area ('1', 1, 1, 1, 1, 1, .false., dfield, dlats, dlngs, cdstn, nv, dlatsd, dlngsd, cdstnd, idom, ifail)
  Select Case (ifail)
   Case (0)
     ldset = .true.
     prompt_domain = 1
   Case Default
     prompt_domain = 2
  End Select
!
  Return
 End Function prompt_domain
!
!
!
 Subroutine put_s_in_var_url (d, url)
!
! Inserts date into URL of variable for real-time data
!
! Modules
  Use data_time_constants, Only: cmon_l
  Use time,                Only: pdate
!
! Arguments
!
! Input scalars
  Type(pdate), Intent(In) :: d ! - date -
!
! Input/output scalars
  Character(Len=*), Intent(InOut) :: url ! - URL -
!
! Locals
!
! Local scalars
  Integer :: i ! - location of date -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Index
!
! Executable Statements
!
! Locate date in URL
  i= Index(url, 'mmmyyyy')
  If (i == 0) Return
!
! Insert date
  Write (Unit=url(i:i+6), Fmt='(A3,I4)') cmon_l(d%imn), d%iyr
!
  Return
 End Subroutine put_s_in_var_url
!
!
!
 Function change_time_grids() &
          Result (url)
!
! Changes time grids of real-time forecasts to those of the hindcasts
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Add conversion
! - reset S grid
  url = 'S/(months since 1960-01-01)streamgridunitconvert/'
! - change T grid to L grid -
  url = Trim(url)//'T//L/(months)ordered/-1.5/1/6.5/NewEvenGRID/replaceGRID/'
!
  Return
 End Function change_time_grids
!
!
!
 Function equate_y_grid() &
          Result (url)
!
! Equates Y grids (in case of rounding errors on gaussian grid)
!
! Function type
  Character(Len=32) :: url ! - URL -
!
! Executable Statements
!
! Equate Y grids
  url = 'Y/2/index/.Y/replaceGRID/'
!
  Return
 End Function equate_y_grid
!
!
!
 Function get_url_append() &
          Result (url)
!
! Constructs unit conversion component of URL
!
! Modules
  Use data_iridl, Only: url_append
!
! Function type
  Character(Len=16) :: url ! - URL -
!
! Executable Statements
!
! Convert units
  url = url_append
!
  Return
 End Function get_url_append
!
!
!
 Function get_url_convert(cvt) &
          Result (url)
!
! Constructs unit conversion component of URL
!
! Modules
  Use data_iridl, Only: url_convert
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Character(Len=*), Intent(In) :: cvt ! - conversion units -
!
! Executable Statements
!
! Convert units
  url = cvt//url_convert
!
  Return
 End Function get_url_convert
!
!
!
 Function get_url_l(imd0, imd1, lsd, lave) &
          Result (url)
!
! Constructs time component of URL using lag
!
! Modules
  Use time,  Only: add_to_month
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: imd0 ! - month of forecast initialization -
  Integer, Intent(In) :: imd1 ! - first month of data to download -
  Integer, Intent(In) :: lsd  ! - length of season to average -
!
  Logical, Intent(In) :: lave  ! - average? -
!
! Locals
!
! Local scalars
  Integer :: l1 ! - first lag -
  Integer :: l2 ! - second lag -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Determine lags
  If (imd0 /= imd1) Then
     l1 = add_to_month(imd1, -imd0)
  Else
     l1 = 0
  End If
  l2 = l1 + lsd - 1
!
! Construct time component of URL 
  Select Case (lsd)
   Case (1) ! - no seasonal averaging -
     Write (Unit=url, Fmt='(A,I0,A)') 'L/(', l1, '.5)VALUES/'
   Case Default ! - seasonal averaging -
     Write (Unit=url, Fmt='(A,I0,A,I0,A)') 'L/', l1, '.5/', l2, '.5/RANGE'
     If (lave) url = Trim(url)//'[L]//keepgrids/average'
     url = Trim(url)//'/'
  End Select
!
  Return
 End Function get_url_l
!
!
!
 Function get_url_m(nem) &
          Result (url)
!
! Constructs ensemble averaging component of URL
!
! Function type
  Character(Len=32) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: nem ! - number of ensemble members -
!
! Executable Statements
!
! Construct time component of URL 
  Write (Unit=url, Fmt='(A,I0,A)') 'M/1/', nem, '/RANGE[M]average/'
!
  Return
 End Function get_url_m
!
!
!
 Function get_url_missing(rmiss) &
          Result (url)
!
! Constructs missing value component of URL
!
! Modules
  Use data_numbers, Only: rp
  Use data_iridl,   Only: url_missing
  Use maths,        Only: get_cnumber
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Real(Kind=rp), Intent(In) :: rmiss ! - missing value flag -
!
! Executable Statements
!
! Construct missing value component of URL
! - set missing value flag -
  url = get_cnumber(rmiss,6)
! - construct URL -
  url = Trim(url)//url_missing
!
  Return
 End Function get_url_missing
!
!
!
 Function get_url_p(cg, n, i, p) &
          Result (url)
!
! Constructs pressure level component of URL
!
! Function type
  Character(Len=256) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: n ! - number of pressure level -
!
  Character(Len=*), Intent(In) :: cg ! - grid code -
!
! Input scalars
  Integer, Dimension(:), Intent(In) :: i ! - selected levels -
!
  Character(Len=*), Dimension(:), Intent(In) :: p ! - pressure levels -
!
! Locals
!
! Local scalars
  Integer :: ilev ! - current level -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Construct pressure level component of URL
  Write (Unit=url, Fmt='(2A)') cg,'/'
! - add pressure levels -
  Do ilev = 1, n
     If (i(ilev) == 1) url=Trim(url)//'('//Trim(p(ilev))//')'
  End Do
  url = Trim(url)//'VALUES/'
!
  Return
 End Function get_url_p
!
!
!
 Function get_url_s(iyd1,iyd2) &
          Result (url)
!
! Constructs time component of URL using start date
!
! Modules
  Use data_time_constants, Only: cmon
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Integer, Intent(In) :: iyd1 ! - first year of data to download -
  Integer, Intent(In) :: iyd2 ! - last year of data to download -
!
! Executable Statements
!
! Construct time component of URL 
  Write (Unit=url, Fmt='(3A,I4,A,I4,A)') 'S/(0000 1 ', cmon(imd0,1), ' ', iyd1, '-', iyd2, ')VALUES/'
!
  Return
 End Function get_url_s
!
!
!
 Function get_url_t(lsum) &
          Result (url)
!
! Constructs time component of URL 
!
! Modules
  Use data_time_constants, Only: cmon
  Use time,                Only: add_to_month
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Arguments
!
! Input scalars
  Logical, Intent(In) :: lsum ! - sum rather than average? -
!
! Locals
!
! Local scalars
  Character(Len=24) :: cfmt ! - format statement -
  Character(Len= 7) :: cmul ! - multiply (for sum) -
!
! Executable Statements
!
! Construct time component of URL 
  Select Case (lsd)
   Case (1) ! - no seasonal averaging -
     Select Case (nld)
      Case (1) ! - no lagged fields -
        Write (Unit=url, Fmt='(3A,I4,A,I4,A)') 'T/(', cmon(imd1,1), ' ', iyd1, '-', iyd2, ')VALUES/'
      Case Default ! - lagged fields -
        imd2 = add_to_month(imd1, nld-1)
        Write (Unit=url, Fmt='(5A,I4,A,I4,A)') 'T/(', cmon(imd1,1), '-', cmon(imd2,1), ' ', iyd1, '-', iyd2, ')VALUES/'
     End Select
   Case Default ! - seasonal averaging -
     Select Case (nld)
      Case (1) ! - no lagged fields -
        imd2 = add_to_month(imd1, lsd-1)
        If (lsd < 10) Then
           cfmt = '(3A,I4,3A,I4,A,I1,A)'
        Else
           cfmt = '(3A,I4,3A,I4,A,I2,A)'
        End If
        Write (Unit=url, Fmt=cfmt) &
          'T/(', cmon(imd1,1), ' ', iyd1, ')(', cmon(imd2,1), ' ', iyd2, ')RANGEEDGES/T/', lsd, '/boxAverage/T/12/STEP/'
      Case Default ! - lagged fields -
        imd2 = add_to_month(imd1, lsd*nld-1)
        If (lsd < 10) Then
           cfmt = '(3A,I4,3A,I4,A,I1,6A)'
        Else
           cfmt = '(3A,I4,3A,I4,A,I2,6A)'
        End If
        Write (Unit=url, Fmt=cfmt) &
          'T/(', cmon(imd1,1), ' ', iyd1, ')(', cmon(imd2,1), ' ', iyd2, ')RANGEEDGES/T/', lsd, '/boxAverage/'//&
         &'T/(', cmon(imd1,1), '-', cmon(imd2,1), ')VALUES/'
     End Select
  End Select
!
! Convert to sum if desired
  If (lsum .and. lsd > 1) Then
     If (lsd < 10) Then
        cfmt = '(I1,A)'
     Else
        cfmt = '(I2,A)'
     End If
     Write (Unit=cmul, Fmt=cfmt) lsd, '/mul/'
     url=Trim(url)//Trim(cmul)
  End If
!
  Return
 End Function get_url_t
!
!
!
 Function get_url_xy() &
          Result (url)
!
! Constructs space component of URL 
!
! Modules
  Use data_numbers, Only: r360
  Use maths,        Only: get_cnumber
  Use fields,       Only: dfield=>tfield
!
! Function type
  Character(Len=128) :: url ! - URL -
!
! Locals
!
! Local scalars
  Character(Len=32) :: clatn ! - northernmost latitude -
  Character(Len=32) :: clats ! - southernmost latitude -
  Character(Len=32) :: clngw ! - westernmost longitude -
  Character(Len=32) :: clnge ! - easternmost longitude -
!
! Functions and Subroutines
!
! Intrinsic functions
  Intrinsic Trim
!
! Executable Statements
!
! Construct latitude/longitude specifiers
  clatn = get_cnumber(dfield(1)%region%alim%rltn, 5)
  clats = get_cnumber(dfield(1)%region%alim%rlts, 5)
  clngw = get_cnumber(dfield(1)%region%alim%rlgw, 5)
  If (dfield(1)%region%alim%rlgw < dfield(1)%region%alim%rlge) Then
     clnge = get_cnumber(dfield(1)%region%alim%rlge, 5)
  Else
     clnge = get_cnumber(dfield(1)%region%alim%rlge+r360, 5)
  End If
!
! Construct space component of URL 
  Write (Unit=url, Fmt='(9A)') &
    'Y/(', Trim(clatn), ')(', Trim(clats), ')RANGEEDGES/X/(', Trim(clngw), ')(', Trim(clnge), ')RANGEEDGES/'
!
  Return
 End Function get_url_xy
End Module downloads
