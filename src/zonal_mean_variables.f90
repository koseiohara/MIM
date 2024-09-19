module zonal_mean_variables

    use params         , only : kp, grav, pi, radius, wind_min, wind_max, st_min, st_max
    use com_var        , only : im, jm, km, ko, pout, costbl
    use mim_var        , only : u, v, pt_dot, u_zm, v_zm, st_zm, pt_dot_zm, p_pds
    use status_output  , only : warn_write
    use derivative     , only : derivative_p
    use isentropic_mean, only : zonalMean, p2pd_integral

    implicit none

    private
    public :: zonal_mean_u, get_v_zm_st, zonal_mean_pt_dot, zonal_mean_square, zonal_mean_cube

    contains


    subroutine zonal_mean_u()

        call zonalMean(u(1:im,1:jm,1:km), &  !! IN
                     & u_zm(1:jm,1:ko)    )  !! OUT

        call warn_write(1              , &  !! IN
                      & jm             , &  !! IN
                      & ko             , &  !! IN
                      & u_zm(1:jm,1:ko), &  !! IN
                      & wind_min       , &  !! IN
                      & wind_max       , &  !! IN
                      & 'u_zm'         , &  !! IN
                      & 'zonal_mean_u()' )  !! IN
        
    end subroutine zonal_mean_u


    subroutine get_v_zm_st()
        real(kp), parameter :: coef = 2._kp * pi * radius * 100._kp / grav
        real(kp) :: v_zm_vint(jm,ko)
        real(kp) :: dummy(im,jm,ko)

        integer :: k

        call p2pd_integral(v(1:im,1:jm,1:km)    , &  !! IN
                         & dummy(1:im,1:jm,1:ko), &  !! OUT
                         & v_zm_vint(1:jm,1:ko)   )  !! OUT

        call derivative_p(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & pout(1:ko)          , &  !! IN
                        & p_pds(1:jm)         , &  !! IN
                        & v_zm_vint(1:jm,1:ko), &  !! IN
                        & v_zm(1:jm,1:ko)       )  !! OUT

        call warn_write(1              , &  !! IN
                      & jm             , &  !! IN
                      & ko             , &  !! IN
                      & v_zm(1:jm,1:ko), &  !! IN
                      & wind_min       , &  !! IN
                      & wind_max       , &  !! IN
                      & 'v_zm'         , &  !! IN
                      & 'zonal_mean_v()' )  !! IN

        do k = 1, ko
            st_zm(1:jm,k) = v_zm_vint(1:jm,k) * coef * costbl(1:jm)
        enddo

        call warn_write(1               , &  !! IN
                      & jm              , &  !! IN
                      & ko              , &  !! IN
                      & st_zm(1:jm,1:ko), &  !! IN
                      & st_min          , &  !! IN
                      & st_max          , &  !! IN
                      & 'streamfunction', &  !! IN
                      & 'get_v_zm_st'     )  !! IN
        
    end subroutine get_v_zm_st


    subroutine zonal_mean_pt_dot()

        call zonalMean(pt_dot(1:im,1:jm,1:km), &  !! IN
                      & pt_dot_zm(1:jm,1:ko)    )  !! OUT

    end subroutine zonal_mean_pt_dot


    subroutine zonal_mean_square(var1, var2, var1_zm, var2_zm, prod_zm, eddy_prod_zm)
        real(kp), intent(in)  :: var1(im,jm,km)
        real(kp), intent(in)  :: var2(im,jm,km)
        real(kp), intent(in)  :: var1_zm(jm,ko)
        real(kp), intent(in)  :: var2_zm(jm,ko)
        real(kp), intent(out) :: prod_zm(jm,ko)
        real(kp), intent(out) :: eddy_prod_zm(jm,ko)

        real(kp) :: prod(im,jm,km)

        prod(1:im,1:jm,1:km) = var1(1:im,1:jm,1:km) * var2(1:im,1:jm,1:km)

        call zonalMean(prod(1:im,1:jm,1:km), &  !! IN
                     & prod_zm(1:jm,1:ko)    )  !! OUT

        eddy_prod_zm(1:jm,1:ko) = prod_zm(1:jm,1:ko) - var1_zm(1:jm,1:ko) * var2_zm(1:jm,1:ko)

    end subroutine zonal_mean_square


    subroutine zonal_mean_cube(var1, var2, var3, prod_zm)
        real(kp), intent(in)  :: var1(im,jm,km)
        real(kp), intent(in)  :: var2(im,jm,km)
        real(kp), intent(in)  :: var3(im,jm,km)
        real(kp), intent(out) :: prod_zm(jm,ko)

        real(kp) :: prod(im,jm,km)

        prod(1:im,1:jm,1:km) = var1(1:im,1:jm,1:km) * var2(1:im,1:jm,1:km) * var3(1:im,1:jm,1:km)

        call zonalMean(prod(1:im,1:jm,1:km), &  !! IN
                     & prod_zm(1:jm,1:ko)    )  !! OUT

    end subroutine zonal_mean_cube


end module zonal_mean_variables

