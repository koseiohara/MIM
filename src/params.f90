!#########################################################
!
!  Parameter  
!
!#########################################################
module params
    implicit none
  
    integer, parameter :: rkp = 4  ! Kind parameter for the input files
    integer, parameter :: wkp = 4  ! Kind parameter for the output files
    integer, parameter :: kp  = 4  ! Kind parameter for computation
    integer, parameter :: ckp = 4   ! Kind parameter for fft
  
    integer, parameter :: FILENAME_MAX = 256

    real(kp), parameter :: rkappa  = 0.286_kp                                   ! r/cp
    real(kp), parameter :: grav    = 9.81_kp                                    ! m/sec**2
    real(kp), parameter :: pi      = 3.141592653589793238462643383279_kp
    real(kp), parameter :: radian  = pi / 180._kp                               ! const. for angle
    real(kp), parameter :: radius  = 6.371E+6_kp                                ! earth radius
    real(kp), parameter :: h0      = 7000._kp                                   ! scale height
    real(kp), parameter :: cp      = 1003._kp                                   ! heat capacity (j/(kg*k))
    real(kp), parameter :: gasr    = rkappa * cp                                !
    real(kp), parameter :: twomg   = 1.45849E-4_kp                              !
    real(kp), parameter :: gamma   = 6.5E-3_kp                                  ! lapse rate
  
    ! possible range of each variable ( e.g. used in check_range() )
  !  real(4),parameter :: t_min     = 100.0      ! Temperature [K]
    real(kp), parameter :: t_min     =  50._kp                                  ! Temperature [K]
    real(kp), parameter :: t_max     =  10000._kp                               ! Temperature [K]
    real(kp), parameter :: pt_min    =  100._kp                                 ! Potential Temperature [K]
    real(kp), parameter :: pt_max    =  10000._kp                               ! Potential Temperature [K]
    real(kp), parameter :: wind_min  = -500._kp                                 ! Horizontal Wind [m/s]
    real(kp), parameter :: wind_max  =  500._kp                                 ! Horizontal Wind [m/s]
    real(kp), parameter :: omega_min = -1.0E+10_kp                              ! Omega-velocity [Pa/s]
    real(kp), parameter :: omega_max =  1.0E+10_kp                              ! Omega-velocity [Pa/s]
    real(kp), parameter :: w_min     = -100._kp                                 ! Vertical Wind [m/s]
    real(kp), parameter :: w_max     =  100._kp                                 ! Vertical Wind [m/s]
    real(kp), parameter :: p_min     =  1.0E-10_kp                              ! Pressure [hPa]
    real(kp), parameter :: p_max     =  1200._kp                                ! Pressure [hPa]
    real(kp), parameter :: z_min     = -3000._kp                                ! Height [m] (-1000 is better)
    real(kp), parameter :: z_max     =  1.0E+6_kp                               ! Height [m]
    real(kp), parameter :: alt_min   = -1000._kp                                ! Altitude [m]
    real(kp), parameter :: alt_max   =  10000._kp                               ! Altitude [m]
    real(kp), parameter :: st_min    = -1.0E+12_kp                              ! Streamfunction [kg/s]
    real(kp), parameter :: st_max    =  1.0E+12_kp                              ! Streamfunction [kg/s]
    real(kp), parameter :: divf_min  = -0.1_kp                                  ! EP Flux divergence [m/s^2]
    real(kp), parameter :: divf_max  =  0.1_kp                                  ! EP Flux divergence [m/s^2]
    real(kp), parameter :: econv_min = -10._kp                                  ! Energy Conversion [W/kg]
    real(kp), parameter :: econv_max =  10._kp                                  ! Energy Conversion [W/kg]

end module params

