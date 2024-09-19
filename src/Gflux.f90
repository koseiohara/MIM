module Gflux

    use params            , only : kp, radius
    use com_var           , only : jm, ko, zd, rho, alat, costbl
    use mim_var           , only : v_zm, pt_zm, v_v_x_zm, v_pt_dot_x_zm, gy, dgy, gz, dgz
    use derivative        , only : derivative_y, derivative_z
    use w_advection_expand, only : w_advection_components

    implicit none

    private
    public :: gflux_y, gflux_div_y, gflux_z, gflux_div_z

    contains


    !
    ! meridional component of G flux 
    ! (G flux: eddy transport in meridional momentum equation)
    !
    subroutine gflux_y()

        gy(1:jm,1:ko) = -spread(rho(1:ko)*radius, 1, jm) * v_v_x_zm(1:jm,1:ko)

    end subroutine gflux_y


    !
    ! G flux divergence (meridional component)
    !
    subroutine gflux_div_y()
        real(kp) :: work_gy_cos_filter(jm,ko)

        work_gy_cos_filter(1:jm,1:ko) = gy(1:jm,1:ko) * spread(costbl(1:jm), 2, ko)

        call derivative_y(1                            , &  !! IN
                        & jm                           , &  !! IN
                        & ko                           , &  !! IN
                        & alat(1:jm)                   , &  !! IN
                        & work_gy_cos_filter(1:jm,1:ko), &  !! IN
                        & dgy(1:jm,1:ko)                 )  !! OUT

        dgy(1:jm,1:ko) = dgy(1:jm,1:ko) / (radius*radius * spread(rho(1:ko), 1, jm) * spread(costbl(1:jm), 2, ko))
        dgy(1,1:ko) = 0._kp
        dgy(jm,1:ko) = 0._kp
      
    end subroutine gflux_div_y


    !
    ! vertical component of G flux 
    ! (G flux: eddy transport in meridional momentum equation)
    !
    subroutine gflux_z()
        real(kp) :: dummy_uv(jm,ko)
        real(kp) :: dummy_ut(jm,ko)

        call w_advection_components(v_v_x_zm(1:jm,1:ko)     , &  !! IN
                                  & v_pt_dot_x_zm(1:jm,1:ko), &  !! IN
                                  & dummy_uv(1:jm,1:ko)     , &  !! OUT
                                  & dummy_ut(1:jm,1:ko)     , &  !! OUT
                                  & gz(1:jm,1:ko)             )  !! OUT

        gz(1:jm,1:ko) = -gz(1:jm,1:ko) * spread(rho(1:ko)*radius, 1, jm)

    end subroutine gflux_z


    !
    ! G flux divergence (vertical component)
    !
    subroutine gflux_div_z()

        call derivative_z(1            , &  !! IN
                        & jm           , &  !! IN
                        & ko           , &  !! IN
                        & zd(1:ko)     , &  !! IN
                        & gz(1:jm,1:ko), &  !! IN
                        & dgz(1:jm,1:ko) )  !! OUT

        dgz(1:jm,1:ko) = dgz(1:jm,1:ko) * spread(1._kp/(rho(1:ko)*radius), 1, jm)

    end subroutine gflux_div_z


end module Gflux

