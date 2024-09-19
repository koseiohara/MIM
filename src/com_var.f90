!############################################################
!
!  Common variables throughout the MIM program
!
!    Don't change variables out of com_var.f90
!
!############################################################
module com_var

    use params  , only : kp, grav, pi, h0
    use namelist, only : INPUT_XDEF_NUM                                                                       , &
                       & INPUT_YDEF_TYPE, INPUT_YDEF_NUM, INPUT_YDEF_LEVEL, INPUT_YDEF_SOUTH, INPUT_YDEF_NORTH, &
                       & INPUT_ZDEF_NUM, INPUT_ZDEF_LEVEL                                                     , &
                       & OUTPUT_ZDEF_NUM, OUTPUT_ZDEF_LEVEL                                                   , &
                       & WAVE_MAX_NUMBER

    implicit none

    private
    public :: im, jm, km, ko, wmax                                                  , &
            & pin, pout, zd, rho, alat, sintbl, costbl, tantbl                      , &
            & com_var_ini, com_var_end, get_pin, get_pout, get_zd, get_rho, get_alat, &
            & get_sintbl, get_costbl, get_tantbl
  
    integer :: im      ! x (East -> West) direction grid
    integer :: jm      ! y (North -> South) direction grid
    integer :: km     ! input data z direction levels (p levels)
    integer :: ko     ! output data z direction levels (p+ levels, p++ levels)
    integer :: wmax    ! maximum wave number to analyze
  
    real(kp), allocatable :: pin(:)     ! standard pressure levels [hPa]
    real(kp), allocatable :: pout(:)    ! standard p+/p++ levels [hPa]
    real(kp), allocatable :: zd(:)      ! z_dagger corresponding to pout [m]
    real(kp), allocatable :: rho(:)     ! standard density corresponding to pout [kg/m^3]
    real(kp), allocatable :: alat(:)    ! latitude [rad]
    real(kp), allocatable :: sintbl(:)  ! sin(latitude)
    real(kp), allocatable :: costbl(:)  ! cos(latitude)
    real(kp), allocatable :: tantbl(:)  ! tan(latitude)
  
  
    contains


    subroutine com_var_ini()
  
        im   = INPUT_XDEF_NUM
        jm   = INPUT_YDEF_NUM
        km   = INPUT_ZDEF_NUM
        ko   = OUTPUT_ZDEF_NUM
        wmax = WAVE_MAX_NUMBER

        allocate(pin(km)  )
        allocate(pout(ko) )
        allocate(zd(ko)   )
        allocate(rho(ko)  )
        allocate(alat(jm)  )
        allocate(sintbl(jm))  
        allocate(costbl(jm))
        allocate(tantbl(jm))

    end subroutine com_var_ini
  
  
    subroutine com_var_end()

        deallocate(pin   )
        deallocate(pout  )
        deallocate(zd    )
        deallocate(rho   )
        deallocate(alat  )
        deallocate(sintbl)  
        deallocate(costbl)
        deallocate(tantbl)
        
    end subroutine com_var_end
  
  
    !*********************************!
    !                                 !
    !          Set variables          !
    !                                 !
    !*********************************!
  
    ! pin [hPa]
    subroutine get_pin()

        if (INPUT_ZDEF_LEVEL(1) < INPUT_ZDEF_LEVEL(2)) then    ! Upper -> Lower
            pin(1:km) = INPUT_ZDEF_LEVEL(1:km)
        else                                                  ! Lower -> Upper
            pin(1:km) = INPUT_ZDEF_LEVEL(km:1:-1)
        endif

    end subroutine get_pin
  
  
    ! pout [hPa]
    subroutine get_pout()

        if (OUTPUT_ZDEF_LEVEL(1) < OUTPUT_ZDEF_LEVEL(2)) then  ! Upper -> Lower
           pout(1:ko) = OUTPUT_ZDEF_LEVEL(1:ko)
        else                                                  ! Lower -> Upper
           pout(1:ko) = OUTPUT_ZDEF_LEVEL(ko:1:-1)
        endif

    end subroutine get_pout
  
  
    ! zd [m]
    subroutine get_zd()

        zd(1:ko) = -h0 * log(pout(1:ko) * 1.E-3_kp)

    end subroutine get_zd
  
  
    ! rho [kg/m^3]
    subroutine get_rho()

        rho(1:ko) = pout(1:ko) * 100._kp / (h0*grav)
        
    end subroutine get_rho
  
  
    ! alat [rad]
    subroutine get_alat()
        real(kp) :: rjm
        integer  :: j
        
        ! - all variables should be YREV (north->south) 
        !   and ZREV (upper->lower) in mim.f90.
        if (trim(INPUT_YDEF_TYPE) == 'lat_radian') then 
  
            if (INPUT_YDEF_LEVEL(1) < INPUT_YDEF_LEVEL(2)) then
                alat(1:jm) = INPUT_YDEF_LEVEL(jm:1:-1)
            else
                alat(1:jm) = INPUT_YDEF_LEVEL(1:jm)
            endif
           
        else if (trim(INPUT_YDEF_TYPE) == 'lat_degree') then
  
            if (INPUT_YDEF_LEVEL(1) < INPUT_YDEF_LEVEL(2)) then
                alat(1:jm) = INPUT_YDEF_LEVEL(jm:1:-1) * pi / 180._kp
            else
                alat(1:jm) = INPUT_YDEF_LEVEL(1:jm) * pi / 180._kp
            endif
           
        else if (trim(INPUT_YDEF_TYPE) == 'linear') then
  
            rjm = ( INPUT_YDEF_NORTH - INPUT_YDEF_SOUTH ) / real(jm-1, kind=kp)

            do j = 1, jm
                alat(j) = ( INPUT_YDEF_NORTH - real(j-1, kind=kp) * rjm ) * pi / 180._kp
            enddo
  
        endif
  
    end subroutine get_alat
  
  
    ! sine table
    subroutine get_sintbl()

        sintbl(1:jm) = sin(alat(1:jm))

    end subroutine get_sintbl
  
  
    ! cosine table
    subroutine get_costbl()

        costbl(1:jm) = cos(alat(1:jm))

    end subroutine get_costbl
  
  
    ! tangent table
    subroutine get_tantbl()

        tantbl(1:jm) = tan(alat(1:jm))
       
    end subroutine get_tantbl
  
  
end module com_var

