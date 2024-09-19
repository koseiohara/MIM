!#########################################################
!
!  Namelistlist  
!
!#########################################################
module namelist

    use params, only : rkp, kp, FILENAME_MAX

    implicit none
  
    integer, parameter :: km_max=512    ! maximum km
    integer, parameter :: ko_max=512    ! maximum ko
    integer, parameter :: jm_max=2048   ! maximum jm

    character(FILENAME_MAX) :: INPUT_UVT_FILENAME
    character(FILENAME_MAX) :: INPUT_U_FILENAME                 ! Zonal Wind [m/s]
    character(FILENAME_MAX) :: INPUT_V_FILENAME                 ! Meridional WInd [m/s]
    character(FILENAME_MAX) :: INPUT_T_FILENAME                 ! Temperature [K]
    character(FILENAME_MAX) :: INPUT_PS_FILENAME                ! 
    character(FILENAME_MAX) :: INPUT_MSL_FILENAME
    character(FILENAME_MAX) :: INPUT_TS_FILENAME
    character(FILENAME_MAX) :: INPUT_Z_FILENAME                 ! Geopotential [m] or Geopotential Height [m^2/s^2]
    character(FILENAME_MAX) :: INPUT_OMEGA_FILENAME             ! Vertical Velocity [Pa/s]
    character(FILENAME_MAX) :: INPUT_TOPO_FILENAME
    character(FILENAME_MAX) :: INPUT_Q_FILENAME
    character(FILENAME_MAX) :: INPUT_SHORTWAVE_FILENAME         ! Q by short wave radiation
    character(FILENAME_MAX) :: INPUT_LONGWAVE_FILENAME          ! Q by long wave radiation
    character(FILENAME_MAX) :: INPUT_LHR_LARGE_FILENAME         ! Q by latent heat release due to large scale precipitation
    character(FILENAME_MAX) :: INPUT_LHR_CONV_FILENAME          ! Q by latent heat release due to convective heating
    character(FILENAME_MAX) :: INPUT_DIFFUSION_FILENAME         ! Q by vertical diffusion
  
    character(32) :: INPUT_UNIT_Z
    character(32) :: INPUT_UNIT_PS
    character(32) :: INPUT_UNIT_MSL
    character(32) :: INPUT_UNIT_TOPO
    character(32) :: INPUT_UNIT_Q
    
    real(rkp) :: INPUT_UNDEF_DEFAULT  ! below derivatives
    real(rkp) :: INPUT_UNDEF_UVT
    real(rkp) :: INPUT_UNDEF_U
    real(rkp) :: INPUT_UNDEF_V
    real(rkp) :: INPUT_UNDEF_T
    real(rkp) :: INPUT_UNDEF_PS
    real(rkp) :: INPUT_UNDEF_MSL
    real(rkp) :: INPUT_UNDEF_TS
    real(rkp) :: INPUT_UNDEF_Z
    real(rkp) :: INPUT_UNDEF_OMEGA
    real(rkp) :: INPUT_UNDEF_Q
    real(rkp) :: INPUT_UNDEF_SHORTWAVE
    real(rkp) :: INPUT_UNDEF_LONGWAVE
    real(rkp) :: INPUT_UNDEF_LHR_LARGE
    real(rkp) :: INPUT_UNDEF_LHR_CONV
    real(rkp) :: INPUT_UNDEF_DIFFUSION
  
    character(16) :: INPUT_ENDIAN_DEFAULT
    character(16) :: INPUT_ENDIAN_UVT
    character(16) :: INPUT_ENDIAN_U
    character(16) :: INPUT_ENDIAN_V
    character(16) :: INPUT_ENDIAN_T
    character(16) :: INPUT_ENDIAN_PS
    character(16) :: INPUT_ENDIAN_MSL
    character(16) :: INPUT_ENDIAN_TS
    character(16) :: INPUT_ENDIAN_Z
    character(16) :: INPUT_ENDIAN_OMEGA
    character(16) :: INPUT_ENDIAN_Q
    character(16) :: INPUT_ENDIAN_SHORTWAVE
    character(16) :: INPUT_ENDIAN_LONGWAVE
    character(16) :: INPUT_ENDIAN_LHR_LARGE
    character(16) :: INPUT_ENDIAN_LHR_CONV
    character(16) :: INPUT_ENDIAN_DIFFUSION
    character(16) :: INPUT_ENDIAN_TOPO
  
    integer :: INPUT_XDEF_NUM
  
    character(16) :: INPUT_YDEF_TYPE
    integer       :: INPUT_YDEF_NUM
    real(kp)      :: INPUT_YDEF_LEVEL(jm_max)
    real(kp)      :: INPUT_YDEF_SOUTH
    real(kp)      :: INPUT_YDEF_NORTH
    logical       :: INPUT_YDEF_YREV_DEFAULT
    logical       :: INPUT_YDEF_YREV_TOPO
  
    integer  :: INPUT_ZDEF_NUM
    !integer  :: INPUT_ZDEF_NUM_OMEGA
    real(kp) :: INPUT_ZDEF_LEVEL(km_max)
    logical  :: INPUT_ZDEF_ZREV
  
    character(16) :: INPUT_TDEF_TYPE
    integer       :: INPUT_TDEF_DAYNUM
    integer       :: INPUT_TDEF_TSTEP
    integer       :: INPUT_TDEF_YEAR
    integer       :: INPUT_TDEF_MONTH
    integer       :: INPUT_TDEF_365DAY
    real(kp)      :: INPUT_TDEF_DT  ! derivative
  
    integer :: WAVE_MAX_NUMBER
  
    integer  :: OUTPUT_ZDEF_NUM
    real(kp) :: OUTPUT_ZDEF_LEVEL(ko_max)
  
    character(FILENAME_MAX) :: OUTPUT_ZONAL_FILENAME
    character(FILENAME_MAX) :: OUTPUT_VINT_FILENAME
    character(FILENAME_MAX) :: OUTPUT_GMEAN_FILENAME
    character(FILENAME_MAX) :: OUTPUT_WAVE_FILENAME
    character(FILENAME_MAX) :: OUTPUT_ERROR_FILENAME
    character(FILENAME_MAX) :: OUTPUT_WARN_FILENAME

    logical :: Q_EXIST
    logical :: Q_COMPS_EXIST
  

    contains


    ! Read Namelist from standard input
    subroutine namelist_init()
      
        integer, parameter :: nml_input = 5
    
        !***** declare *****!
        namelist / INPUT  / INPUT_UVT_FILENAME      , &
                          & INPUT_U_FILENAME        , &
                          & INPUT_V_FILENAME        , &
                          & INPUT_T_FILENAME        , &
                          & INPUT_PS_FILENAME       , &
                          & INPUT_MSL_FILENAME      , &
                          & INPUT_TS_FILENAME       , &
                          & INPUT_Z_FILENAME        , &
                          & INPUT_OMEGA_FILENAME    , &
                          & INPUT_TOPO_FILENAME     , &
                          & INPUT_Q_FILENAME        , &
                          & INPUT_SHORTWAVE_FILENAME, &
                          & INPUT_LONGWAVE_FILENAME , &
                          & INPUT_LHR_LARGE_FILENAME, &
                          & INPUT_LHR_CONV_FILENAME , &
                          & INPUT_DIFFUSION_FILENAME
    
        namelist / INPUT_UNIT / INPUT_UNIT_Z   , &
                              & INPUT_UNIT_PS  , &
                              & INPUT_UNIT_MSL , &
                              & INPUT_UNIT_TOPO, &
                              & INPUT_UNIT_Q
    
        namelist / INPUT_UNDEF / INPUT_UNDEF_DEFAULT  , &
                               & INPUT_UNDEF_UVT      , &
                               & INPUT_UNDEF_U        , &
                               & INPUT_UNDEF_V        , &
                               & INPUT_UNDEF_T        , &
                               & INPUT_UNDEF_PS       , &
                               & INPUT_UNDEF_MSL      , &
                               & INPUT_UNDEF_TS       , &
                               & INPUT_UNDEF_Z        , &
                               & INPUT_UNDEF_OMEGA    , &
                               & INPUT_UNDEF_Q        , &
                               & INPUT_UNDEF_SHORTWAVE, &
                               & INPUT_UNDEF_LONGWAVE , &
                               & INPUT_UNDEF_LHR_LARGE, &
                               & INPUT_UNDEF_LHR_CONV , &
                               & INPUT_UNDEF_DIFFUSION
    
        namelist / INPUT_ENDIAN / INPUT_ENDIAN_DEFAULT  , &
                                & INPUT_ENDIAN_UVT      , &
                                & INPUT_ENDIAN_U        , &
                                & INPUT_ENDIAN_V        , &
                                & INPUT_ENDIAN_T        , &
                                & INPUT_ENDIAN_PS       , &
                                & INPUT_ENDIAN_MSL      , &
                                & INPUT_ENDIAN_TS       , &
                                & INPUT_ENDIAN_Z        , &
                                & INPUT_ENDIAN_OMEGA    , &
                                & INPUT_ENDIAN_Q        , &
                                & INPUT_ENDIAN_SHORTWAVE, &
                                & INPUT_ENDIAN_LONGWAVE , &
                                & INPUT_ENDIAN_LHR_LARGE, &
                                & INPUT_ENDIAN_LHR_CONV , &
                                & INPUT_ENDIAN_DIFFUSION, &
                                & INPUT_ENDIAN_TOPO
    
    
        namelist / INPUT_XDEF / INPUT_XDEF_NUM
    
        namelist / INPUT_YDEF / INPUT_YDEF_TYPE        , &
                              & INPUT_YDEF_NUM         , &
                              & INPUT_YDEF_LEVEL       , &
                              & INPUT_YDEF_SOUTH       , &
                              & INPUT_YDEF_NORTH       , &
                              & INPUT_YDEF_YREV_DEFAULT, &
                              & INPUT_YDEF_YREV_TOPO

        namelist / INPUT_ZDEF / INPUT_ZDEF_NUM      , &
                              & INPUT_ZDEF_LEVEL    , &
                              & INPUT_ZDEF_ZREV
                              !& INPUT_ZDEF_NUM_OMEGA, &
    
        namelist / INPUT_TDEF / INPUT_TDEF_TYPE  , &
                              & INPUT_TDEF_DAYNUM, &
                              & INPUT_TDEF_TSTEP , &
                              & INPUT_TDEF_YEAR  , &
                              & INPUT_TDEF_MONTH , &
                              & INPUT_TDEF_365DAY
    
        namelist / WAVE / WAVE_MAX_NUMBER
    
        namelist / OUTPUT / OUTPUT_ZONAL_FILENAME, &
                          & OUTPUT_VINT_FILENAME , &
                          & OUTPUT_GMEAN_FILENAME, &
                          & OUTPUT_WAVE_FILENAME , &
                          & OUTPUT_ERROR_FILENAME, &
                          & OUTPUT_WARN_FILENAME
    
        namelist / OUTPUT_ZDEF / OUTPUT_ZDEF_NUM  , &
                               & OUTPUT_ZDEF_LEVEL
    
    
        !***** default values *****!
        INPUT_UVT_FILENAME       = ''
        INPUT_U_FILENAME         = ''
        INPUT_V_FILENAME         = ''
        INPUT_T_FILENAME         = ''
        INPUT_PS_FILENAME        = ''
        INPUT_MSL_FILENAME       = ''
        INPUT_TS_FILENAME        = ''
        INPUT_Z_FILENAME         = ''
        INPUT_OMEGA_FILENAME     = ''
        INPUT_TOPO_FILENAME      = ''
        INPUT_Q_FILENAME         = ''
        INPUT_SHORTWAVE_FILENAME = ''
        INPUT_LONGWAVE_FILENAME  = ''
        INPUT_LHR_LARGE_FILENAME = ''
        INPUT_LHR_CONV_FILENAME  = ''
        INPUT_DIFFUSION_FILENAME = ''
    
        INPUT_UNIT_Z    = 'm'
        INPUT_UNIT_PS   = 'hPa'
        INPUT_UNIT_MSL  = 'hPa'
        INPUT_UNIT_TOPO = 'm'
        INPUT_UNIT_Q    = 'K/s'
    
        INPUT_UNDEF_DEFAULT   = 9.999E+20_rkp
        INPUT_UNDEF_UVT       = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_U         = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_V         = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_T         = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_PS        = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_MSL       = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_TS        = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_Z         = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_OMEGA     = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_Q         = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_SHORTWAVE = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_LONGWAVE  = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_LHR_LARGE = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_LHR_CONV  = INPUT_UNDEF_DEFAULT 
        INPUT_UNDEF_DIFFUSION = INPUT_UNDEF_DEFAULT 
    
        INPUT_ENDIAN_DEFAULT   = 'NATIVE'
        INPUT_ENDIAN_UVT       = ''
        INPUT_ENDIAN_U         = ''
        INPUT_ENDIAN_V         = ''
        INPUT_ENDIAN_T         = ''
        INPUT_ENDIAN_PS        = ''
        INPUT_ENDIAN_MSL       = ''
        INPUT_ENDIAN_TS        = ''
        INPUT_ENDIAN_Z         = ''
        INPUT_ENDIAN_OMEGA     = ''
        INPUT_ENDIAN_TOPO      = ''
        INPUT_ENDIAN_Q         = ''
        INPUT_ENDIAN_SHORTWAVE = ''
        INPUT_ENDIAN_LONGWAVE  = ''
        INPUT_ENDIAN_LHR_LARGE = ''
        INPUT_ENDIAN_LHR_CONV  = ''
        INPUT_ENDIAN_DIFFUSION = ''

        INPUT_XDEF_NUM = 0
    
        INPUT_YDEF_NUM          = 0
        INPUT_YDEF_SOUTH        = -90._kp
        INPUT_YDEF_NORTH        =  90._kp
        INPUT_YDEF_YREV_DEFAULT = .False.
        INPUT_YDEF_YREV_TOPO    = .False.
    
        !INPUT_ZDEF_NUM_OMEGA = 0
        INPUT_ZDEF_NUM  = 0
        INPUT_ZDEF_ZREV = .False.
   
        INPUT_TDEF_TYPE   = ''
        INPUT_TDEF_DAYNUM = 0
        INPUT_TDEF_365DAY = 0
    
        WAVE_MAX_NUMBER = 0
    
        OUTPUT_ZONAL_FILENAME = ''
        OUTPUT_VINT_FILENAME  = ''
        OUTPUT_GMEAN_FILENAME = ''
        OUTPUT_WAVE_FILENAME  = ''
        OUTPUT_ERROR_FILENAME = ''
        OUTPUT_WARN_FILENAME  = ''
    
    
        read(nml_input, nml=INPUT       )
        read(nml_input, nml=INPUT_ENDIAN)
        read(nml_input, nml=INPUT_UNDEF )
        read(nml_input, nml=INPUT_UNIT  )
        read(nml_input, nml=INPUT_TDEF  )
        read(nml_input, nml=INPUT_XDEF  )
        read(nml_input, nml=INPUT_YDEF  )
        read(nml_input, nml=INPUT_ZDEF  )
        read(nml_input, nml=OUTPUT      )
        read(nml_input, nml=OUTPUT_ZDEF )
        read(nml_input, nml=WAVE        )
    

        if ( trim(INPUT_ENDIAN_DEFAULT) /= ''              .AND. &
           & trim(INPUT_ENDIAN_DEFAULT) /= 'little_endian' .AND. &
           & trim(INPUT_ENDIAN_DEFAULT) /= 'big_endian'          ) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Unexpected endian name : ' // trim(INPUT_ENDIAN_DEFAULT)
            ERROR STOP
        endif
    
        call endian(INPUT_ENDIAN_UVT      , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_U        , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_V        , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_T        , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_PS       , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_MSL      , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_TS       , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_Z        , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_OMEGA    , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_Q        , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_SHORTWAVE, &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_LONGWAVE , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_LHR_LARGE, &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_LHR_CONV , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_DIFFUSION, &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
        call endian(INPUT_ENDIAN_TOPO     , &  !! INOUT
                  & INPUT_ENDIAN_DEFAULT    )  !! IN
    
        !if (INPUT_YDEF_YREV_TOPO == -1) then
        !    INPUT_YDEF_YREV_TOPO = INPUT_YDEF_YREV_DEFAULT
        !endif
    
        !if (INPUT_ZDEF_NUM_OMEGA == 0) then
        !    INPUT_ZDEF_NUM_OMEGA = INPUT_ZDEF_NUM
        !endif
        
        INPUT_TDEF_DT = real(24*60*60, kind=kp) / real(INPUT_TDEF_DAYNUM, kind=kp)

        Q_EXIST       = (trim(INPUT_Q_FILENAME) /= '')
        Q_COMPS_EXIST = (trim(INPUT_SHORTWAVE_FILENAME) /= '' .AND. &
                       & trim(INPUT_LONGWAVE_FILENAME ) /= '' .AND. &
                       & trim(INPUT_LHR_LARGE_FILENAME) /= '' .AND. &
                       & trim(INPUT_LHR_CONV_FILENAME ) /= '' .AND. &
                       & trim(INPUT_DIFFUSION_FILENAME) /= ''       )
    
        !***** check *****!
        call namelist_check()
        
    end subroutine namelist_init
  
  
    subroutine endian(endian_out, endian_default)
      character(*), intent(inout) :: endian_out
      character(*), intent(in)    :: endian_default
  
      if (trim(endian_out) == '') then
         endian_out = endian_default
      else if (trim(endian_out) /= 'big_endian' .AND. endian_out /= 'little_endian') then
         write(0,'(A)') 'ERROR STOP'
         write(0,'(A)') 'Unexpected endian name : ' // trim(endian_out)
         ERROR STOP
      endif
  
    end subroutine endian
  
  
    !
    ! check namelist
    !
    subroutine namelist_check()

        if (trim(INPUT_UVT_FILENAME) == '' .AND. &
          & (trim(INPUT_U_FILENAME)  == '' .OR.  &
          &  trim(INPUT_V_FILENAME)  == '' .OR.  &
          &  trim(INPUT_T_FILENAME)  == ''     ) ) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_UVT_FILENAME is not specified'
            write(0,'(A)') 'One of INPUT_U_FILENAME, INPUT_V_FILENAME, or INPUT_T_FILENAME is not specified'
            write(0,'(A)') 'Define only INPUT_UVT_FILENAME, or all three : INPUT_U_FILENAME, INPUT_V_FILENAME, and INPUT_T_FILENAME'
            ERROR STOP
        endif

        if (trim(INPUT_PS_FILENAME)   == '' .AND. &
          & (trim(INPUT_MSL_FILENAME) == '' .OR.  &
          &  trim(INPUT_TS_FILENAME)  == ''     ) ) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_PS_FILENAME is not specified'
            write(0,'(A)') 'INPUT_MSL_FILENAME or INPUT_TS_FILENAME is not specified'
            write(0,'(A)') 'Define only INPUT_PS_FILENAME, or both INPUT_MSL_FILENAME and INPUT_TS_FILENAME'
            ERROR STOP
        endif

        if (trim(INPUT_Z_FILENAME) == '') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_Z_FILENAME is not specified'
            ERROR STOP
        endif

        if (trim(INPUT_OMEGA_FILENAME) == '') then
            write(*,*)
            write(*,'(A)') 'INPUT_OMEGA_FILENAME is not specified'
            write(*,'(A)') 'Vertical velocity will be estimated by the continuity equation'
        endif

        if (trim(INPUT_TOPO_FILENAME) == '') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_TOPO_FILENAME is not specified'
            ERROR STOP
        endif

        if ((.NOT. Q_EXIST) .AND. (.NOT. Q_COMPS_EXIST)) then
            write(*,*)
            write(*,'(A)')              'INPUT_Q_FILENAME is not specified'
            write(*,'(A)',ADVANCE='NO') 'One of INPUT_SHORTWAVE_FILENAME, INPUT_LONGWAVE_FILENAME, INPUT_LHR_LARGE_FILENAME, '
            write(*,'(A)')              'INPUT_LHR_CONV_FILENAME, INPUT_DIFFUSION_FILENAME is not specified'
            write(*,'(A)') 'Diabatic heating will be estimated from the temporal derivative of the potential temperature'
        endif
  
        if (trim(INPUT_UNIT_Z) /= 'm' .AND. trim(INPUT_UNIT_Z) /= 'm^2/s^2') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_UNIT_Z : ' // trim(INPUT_UNIT_Z)
            write(0,'(A)') '"m" or "m^2/s^2" is valid as input'
            ERROR STOP
        endif
    
        if (trim(INPUT_UNIT_PS) /= 'Pa' .AND. trim(INPUT_UNIT_PS) /= 'hPa') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_UNIT_PS : ' // trim(INPUT_UNIT_PS)
            write(0,'(A)') '"Pa" or "hPa" is valid as input'
            ERROR STOP
        endif
    
        if (trim(INPUT_UNIT_MSL) /= 'Pa' .AND. trim(INPUT_UNIT_MSL) /= 'hPa') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_UNIT_MSL : ' // trim(INPUT_UNIT_MSL)
            write(0,'(A)') '"Pa" or "hPa" is valid as input'
            ERROR STOP
        end if
    
        if (trim(INPUT_UNIT_TOPO) /= 'm' .AND. trim(INPUT_UNIT_TOPO) /= 'm^2/s^2') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_UNIT_TOPO : ' // trim(INPUT_UNIT_TOPO)
            write(0,'(A)') '"m" or "m^2/s^2" is valid as input'
            ERROR STOP
        endif

        if (trim(INPUT_UNIT_Q) /= 'K/s' .AND. trim(INPUT_UNIT_Q) /= 'K/day') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_UNIT_Q : ' // trim(INPUT_UNIT_Q)
            write(0,'(A)') '"K/s" or "K/day" is valid as input'
            ERROR STOP
        endif
    
        if (INPUT_XDEF_NUM < 1) then
            write(0,'(A)')    'ERROR STOP'
            write(0,'(A,I0)') 'Invalid INPUT_XDEF_NUM : ', INPUT_XDEF_NUM
            write(0,'(A)')    'INPUT_XDEF_NUM must be more than 1'
            ERROR STOP
        endif
    
        if (trim(INPUT_YDEF_TYPE) /= 'lat_degree' .AND. &
          & trim(INPUT_YDEF_TYPE) /= 'lat_radian' .AND. &
          & trim(INPUT_YDEF_TYPE) /= 'linear'           ) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_YDEF_TYPE : ' // trim(INPUT_YDEF_TYPE)
            write(0,'(A)') 'INPUT_YDEF_TYPE must be one of "lat_degree", "lat_radian", or "linear"'
            write(0,'(A)') ''
            write(0,'(A)') 'If INPUT_YDEF_TYPE=lat_degree, the latitude is defined by INPUT_YDEF_LEVEL (degree)'
            write(0,'(A)') 'If INPUT_YDEF_TYPE=lat_radian, the latitude is defined by INPUT_YDEF_LEVEL (radian)'
            write(0,'(A)',ADVANCE='NO') 'If INPUT_YDEF_TYPE=linear, the latitude is defined from INPUT_YDEF_NORTH, '
            write(0,'(A)') 'INPUT_YDEF_SOUTH, and INPUT_YDEF_NUM'
            ERROR STOP
        endif
    
        if (INPUT_YDEF_NUM < 1 .OR. INPUT_YDEF_NUM > jm_max) then
            write(0,'(A)')    'ERROR STOP'
            write(0,'(A,I0)') 'Invalid INPUT_YDEF_NUM : ', INPUT_YDEF_NUM
            write(0,'(A,I0)') 'INPUT_YDEF_NUM should be between 1 and ', jm_max
            ERROR STOP
        endif

        if (INPUT_YDEF_SOUTH < -90._kp .OR. INPUT_YDEF_SOUTH >= 90._kp) then
            write(0,'(A)')       'ERROR STOP'
            write(0,'(A,F0.3)') 'Invalid INPUT_YDEF_SOUTH : ', INPUT_YDEF_SOUTH
            write(0,'(A)')       'INPUT_YDEF_SOUTH must be between -90 and 90'
            ERROR STOP
        endif
    
        if (INPUT_YDEF_NORTH <= -90._kp .OR. INPUT_YDEF_NORTH > 90._kp) then
            write(0,'(A)')       'ERROR STOP'
            write(0,'(A,F0.3)') 'Invalid INPUT_YDEF_NORTH : ', INPUT_YDEF_NORTH
            write(0,'(A)')       'INPUT_YDEF_NORTH must be between -90 and 90'
            ERROR STOP
        endif
    
        if (INPUT_ZDEF_NUM < 1 .OR. INPUT_ZDEF_NUM > km_max) then
            write(0,'(A)')    'ERROR STOP'
            write(0,'(A,I0)') 'Invalid INPUT_ZDEF_NUM : ', INPUT_ZDEF_NUM
            write(0,'(A,I0)') 'INPUT_ZDEF_NUM should be between 1 and ', km_max
            ERROR STOP
        endif
    
        !if (INPUT_ZDEF_NUM_OMEGA > INPUT_ZDEF_NUM) then
        !    write(*,*) 'error: INPUT_ZDEF_NUM_OMEGA(', INPUT_ZDEF_NUM_OMEGA, &
        !        &     ') > INPUT_ZDEF_NUM(', INPUT_ZDEF_NUM, ')'
        !    ERROR STOP
        !endif
    
        if (trim(INPUT_TDEF_TYPE) /= 'tstep'   .AND. &
          & trim(INPUT_TDEF_TYPE) /= 'monthly' .AND. &
          & trim(INPUT_TDEF_TYPE) /= 'annual'        ) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_TDEF_TYPE : ' // trim(INPUT_TDEF_TYPE)
            write(0,'(A)') 'INPUT_TDEF_TYPE must be one of "tstep", "monthly", or "annual"'
            write(0,'(A)',ADVANCE='NO') 'If INPUT_TDEF_TYPE=tstep, '
            write(0,'(A)')              'specify INPUT_TDEF_TSTEP and its value is defined as the number of time steps'
            write(0,'(A)',ADVANCE='NO') 'If INPUT_TDEF_TYPE=monthlly, specify INPUT_TDEF_YEAR, INPUT_TDEF_MONTH, '
            write(0,'(A)')              'and INPUT_TDEF_365DAY and MIM computes for one month'
            write(0,'(A)',ADVANCE='NO') 'If INPUT_TDEF_TYPE=annual, specify INPUT_TDEF_YEAR and INPUT_TDEF_365DAY '
            write(0,'(A)')              'and MIM computes for the year'
            ERROR STOP
        endif

        if (INPUT_TDEF_DAYNUM <= 0) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A,I0)') 'Invalid INPUT_TDEF_DAYNUM : ', INPUT_TDEF_DAYNUM
            write(0,'(A)')    'INPUT_TDEF_DAYNUM must be more than 0'
            ERROR STOP
        endif

        if (WAVE_MAX_NUMBER < 0) then
            write(0,'(A)')    'ERROR STOP'
            write(0,'(A,I0)') 'Invalid WAVE_MAX_NUMBER : ', WAVE_MAX_NUMBER
            write(0,'(A)')    'WAVE_MAX_NUMBER must be equal or more than 0'
            ERROR STOP
        endif
    
        if (OUTPUT_ZDEF_NUM < 1 .OR. OUTPUT_ZDEF_NUM > km_max) then
            write(0,'(A)')    'ERROR STOP'
            write(0,'(A,I0)') 'Invalid OUTPUT_ZDEF_NUM : ', OUTPUT_ZDEF_NUM
            write(0,'(A,I0)') 'OUTPUT_ZDEF_NUM must be between 1 and', km_max
            ERROR STOP
        endif
    
        if (OUTPUT_WAVE_FILENAME /= '') then
            if (WAVE_MAX_NUMBER > INPUT_XDEF_NUM) then
                write(0,'(A)')    'ERROR STOP'
                write(0,'(A)')    'WAVE_MAX_NUMBER must not be bigger than INPUT_XDEF_NUM'
                write(0,'(A,I0)') 'INPUT_XDEF_NUM  : ', INPUT_XDEF_NUM
                write(0,'(A,I0)') 'WAVE_MAX_NUMBER : ', WAVE_MAX_NUMBER
                ERROR STOP
            else if (WAVE_MAX_NUMBER < 1) then
                write(0,'(A)') 'ERROR STOP'
                write(0,'(A)') 'OUTPUT_WAVE_FILENAME is specified, but WAVE_MAX_NUMBER is invalid'
                write(0,'(A)') 'WAVE_MAX_NUMBER must be more than 1'
                write(0,'(A,I0)') 'WAVE_MAX_NUMBER : ', WAVE_MAX_NUMBER
                ERROR STOP
            endif
        endif

        if (trim(OUTPUT_ERROR_FILENAME) == '') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'OUTPUT_ERROR_FILENAME is not specified'
            ERROR STOP
        endif

        if (trim(OUTPUT_WARN_FILENAME) == '') then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'OUTPUT_WARN_FILENAME is not specified'
            ERROR STOP
        endif
  
    end subroutine namelist_check


end module namelist

