module hight_derivative_products

    use params         , only : kp
    use namelist       , only : dt=>INPUT_TDEF_DT
    use com_var        , only : im, jm, km, ko, pout, alat
    use mim_var        , only : u, v, z, p_pd, dz_dlat_zm, v_dz_dlat_zm, u_dz_dlon_zm, phi_dagger, &
                              & p_dphi_dt, z_pd, z_pd_past, p_pd_past
    use derivative     , only : derivative_x, derivative_y
    use isentropic_mean, only : zonalMean

    implicit none

    private
    public :: z_yderiv, z_xderiv, z_tderiv

    contains


    ! get zonal mean of dz/dlat and zonal mean of v*(dz/dlat)
    ! output : dz_dlat_zm and v_dz_dx_zm
    subroutine z_yderiv()
        real(kp) :: dz_dlat(im,jm,km)

        call derivative_y(im                     , &  !! IN
                        & jm                     , &  !! IN
                        & km                     , &  !! IN
                        & alat(1:jm)             , &  !! IN
                        & z(1:im,1:jm,1:km)      , &  !! IN
                        & dz_dlat(1:im,1:jm,1:km)  )  !! OUT

        call zonalMean(dz_dlat(1:im,1:jm,1:km), &  !! IN
                     & dz_dlat_zm(1:jm,1:ko)    )  !! OUT

        dz_dlat(1:im,1:jm,1:km) = dz_dlat(1:im,1:jm,1:km) * v(1:im,1:jm,1:km)

        call zonalMean(dz_dlat(1:im,1:jm,1:km), &  !! IN
                     & v_dz_dlat_zm(1:jm,1:ko)  )  !! OUT

    end subroutine z_yderiv


    ! get zonal mean of u*(dz/dlon)
    ! output : u_dz_dx_zm
    subroutine z_xderiv()
        real(kp) :: dz_dlon(im,jm,km)

        call derivative_x(im                     , &  !! IN
                        & jm                     , &  !! IN
                        & km                     , &  !! IN
                        & z(1:im,1:jm,1:km)      , &  !! IN
                        & dz_dlon(1:im,1:jm,1:km)  )  !! OUT

        dz_dlon(1:im,1:jm,1:km) = dz_dlon(1:im,1:jm,1:km) * u(1:im,1:jm,1:km)

        call zonalMean(dz_dlon(1:im,1:jm,1:km), &  !! IN
                     & u_dz_dlon_zm(1:jm,1:ko)  )  !! OUT

    end subroutine z_xderiv


    subroutine z_tderiv(p_dz_dt, p_dz_dt_zm)
        real(kp), intent(out) :: p_dz_dt(im,jm,ko)
        real(kp), intent(out) :: p_dz_dt_zm(jm,ko)

        p_dz_dt(1:im,1:jm,1:ko) = (z_pd(1:im,1:jm,1:ko) - z_pd_past(1:im,1:jm,1:ko)) &
                              & * (p_pd(1:im,1:jm,1:ko) + p_pd_past(1:im,1:jm,1:ko)) &
                              & / (2._kp * 100._kp * dt)

        p_dz_dt_zm(1:jm,1:ko) = sum(p_dz_dt(1:im,1:jm,1:ko), dim=1) / real(im, kind=kp)

    end subroutine z_tderiv


end module hight_derivative_products

