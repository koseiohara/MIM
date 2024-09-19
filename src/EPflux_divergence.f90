module EPflux_divergence

    use params       , only : kp, radius, divf_min, divf_max
    use com_var      , only : jm, ko, zd, rho, alat, costbl
    use status_output, only : warn_write
    use derivative   , only : derivative_z

    implicit none

    private
    public :: epflux_div_y, epflux_div_z

    contains


    !
    ! EP flux divergence (meridional component)
    !
    subroutine epflux_div_y(y_component, y_divergence)
        real(kp), intent(in)  :: y_component(jm,ko)
        real(kp), intent(out) :: y_divergence(jm,ko)
    
        integer :: k
    
        ! OPENMP
        do k = 1, ko
    
            y_divergence(1,k) = 0._kp

            y_divergence(2:jm-1,k) = (y_component(3:jm,k)*costbl(3:jm) - y_component(1:jm-2,k)*costbl(1:jm-2)) &
                                 & / ((alat(3:jm)-alat(1:jm-2)) * rho(k) * (radius*costbl(2:jm-1))**2)
   
            y_divergence(jm,k) = 0._kp
    
        enddo
    
        call warn_write(1              , &  !! IN
                      & jm             , &  !! IN
                      & ko             , &  !! IN
                      & y_divergence   , &  !! IN
                      & divf_min       , &  !! IN
                      & divf_max       , &  !! IN
                      & 'y_divergence' , &  !! IN
                      & 'epflux_dv_y()'  )  !! IN
    
    end subroutine epflux_div_y
    
    
    
    !
    ! EP flux divergence (vertical  component)
    !
    subroutine epflux_div_z(z_component, z_divergence)
        real(kp), intent(in)  :: z_component(jm,ko)
        real(kp), intent(out) :: z_divergence(jm,ko)
    
        call derivative_z(1                     , &  !! IN
                        & jm                    , &  !! IN
                        & ko                    , &  !! IN
                        & zd(1:ko)              , &  !! IN
                        & z_component(1:jm,1:ko), &  !! IN
                        & z_divergence(1:jm,1:ko) )  !! OUT

        z_divergence(1:jm,1:ko) = z_divergence(1:jm,1:ko) / (spread(rho(1:ko), 1, jm) * spread(costbl(1:jm)*radius, 2, ko))
    
        z_divergence( 1,1:ko) = 0._kp
        z_divergence(jm,1:ko) = 0._kp
    
        call warn_write(1                         , &  !! IN
                      & jm                        , &  !! IN
                      & ko                        , &  !! IN
                      & z_divergence              , &  !! IN
                      & divf_min                  , &  !! IN
                      & divf_max                  , &  !! IN
                      & '(one of the)z_divergence', &  !! IN
                      & 'epflux_div_z()'            )  !! IN
    
    end subroutine epflux_div_z


end module EPflux_divergence

