module estimate_w

    use params       , only : kp, pi, radian, radius, omega_min, omega_max
    use namelist     , only : INPUT_UNDEF_OMEGA
    use com_var      , only : im, jm, km, pin, alat
    use mim_var      , only : u, v, omega, p_sfc
    use status_output, only : warn_write

    implicit none

    private
    public :: get_omega

    contains


    !
    ! Function
    !   estimate vertical velocity using continuity equation
    !
    ! Arguments (in)
    !   u        : x-wind
    !   v        : y-wind
    !   p_sfc    : surface pressure [hPa]
    !
    ! Arguments (out)
    !   omega    : vertical velocity (Pa/s)
    !
    ! Note
    !   - This subroutine is used, if input does not have vertical velocity
    !
    subroutine get_omega()
        real(kp) :: div(im,jm,km)
        real(kp) :: p_d(im,jm)
        real(kp) :: p_u(im,jm)
        real(kp) :: var_d(im,jm)
        real(kp) :: var_u(im,jm)
        real(kp) :: dp(im,jm)
        integer  :: k

        ! calculate horizontal divergence
        call get_div(div(1:im,1:jm,1:km))  !! OUT


        ! calculate omega from vertical integration of horizonal divergence
        omega(1:im,1:jm,1) = 0._kp

        do k = 2, km
            where (pin(k) < p_sfc(1:im,1:jm))
                p_d(1:im,1:jm) = pin(k)
            else where
                p_d(1:im,1:jm) = p_sfc(1:im,1:jm)
            endwhere

            var_d(1:im,1:jm) = div(1:im,1:jm,k)

            where (pin(k-1) < p_sfc(1:im,1:jm))
                p_u(1:im,1:jm) = pin(k-1)
            else where
                p_u(1:im,1:jm) = p_sfc(1:im,1:jm)
            endwhere

            var_u(1:im,1:jm) = div(1:im,1:jm,k-1)

            dp(1:im,1:jm) = max(p_d(1:im,1:jm)-p_u(1:im,1:jm), 0._kp)

            omega(1:im,1:jm,k) = omega(1:im,1:jm,k-1) - 0.5_kp * (var_d(1:im,1:jm) + var_u(1:im,1:jm)) * dp(1:im,1:jm)*100._kp
        enddo

        call warn_write(im                   , &  !! IN
                      & jm                   , &  !! IN
                      & km                   , &  !! IN
                      & omega(1:im,1:jm,1:km), &  !! IN
                      & omega_min            , &  !! IN
                      & omega_max            , &  !! IN
                      & 'omega'              , &  !! IN
                      & 'get_omega()'          )  !! IN

    end subroutine get_omega


    subroutine get_div(div)
        real(kp), intent(out) :: div(im,jm,km)

        real(kp) :: phi(jm-1)
        real(kp) :: stg_costbl(jm-1)
        real(kp) :: stg_sintbl(jm-1)
        real(kp) :: uft(0:im+1,1:jm,1:km)
        real(kp) :: vft(im,jm,km)
        real(kp) :: dlon
        real(kp) :: dlat(jm)
        real(kp) :: fn(1:im,2:jm-1)
        real(kp) :: fs(1:im,2:jm-1)
        real(kp) :: fe(1:im,2:jm-1)
        real(kp) :: fw(1:im,2:jm-1)
        real(kp) :: integral(1:im,2:jm-1)
        real(kp) :: ds(2:jm-1)
        real(kp) :: dsn
        real(kp) :: dss

        integer :: k

        dlon = 360._kp * radian / real(im, kind=kp)     ! [radian]

        !--- get staggard sin&cos
        phi(1:jm-1) = 0.5_kp * (alat(1:jm-1) + alat(2:jm))

        stg_sintbl(1:jm-1) = sin(phi(1:jm-1))
        stg_costbl(1:jm-1) = cos(phi(1:jm-1))

        !--- check lower boundary & substitute 0 below surface
        !
        uft(1:im,1:jm,1:km) = u(1:im,1:jm,1:km)
        vft(1:im,1:jm,1:km) = v(1:im,1:jm,1:km)

        do k = 1, km
            where (pin(k) > p_sfc(1:im,1:jm))
                uft(1:im,1:jm,k) = 0._kp
                vft(1:im,1:jm,k) = 0._kp
            endwhere
        enddo
        uft(   0,1:jm,1:km) = uft(im,1:jm,1:km)
        uft(im+1,1:jm,1:km) = uft(1 ,1:jm,1:km)

        !---- calculate divergence
        dlat(2:jm-1) = 0.5_kp * (alat(1:jm-2) - alat(3:jm))
        ds(2:jm-1) = radius*radius * (stg_sintbl(1:jm-2) - stg_sintbl(2:jm-1)) * dlon
        do k = 1, km

            fn(1:im,2:jm-1) = 0.5_kp * (vft(1:im,2:jm-1,k)   + vft(1:im,1:jm-2,k)  )
            fs(1:im,2:jm-1) = 0.5_kp * (vft(1:im,3:jm,k)     + vft(1:im,2:jm-1,k)  )
            fe(1:im,2:jm-1) = 0.5_kp * (uft(1:im,2:jm-1,k)   + uft(2:im+1,2:jm-1,k))
            fw(1:im,2:jm-1) = 0.5_kp * (uft(0:im-1,2:jm-1,k) + uft(1:im,2:jm-1,k)  )

            integral(1:im,2:jm-1) = radius * (  spread(dlat(2:jm-1), 1, im) * (fe(1:im,2:jm-1) - fw(1:im,2:jm-1)) &
                                            & - spread(stg_costbl(2:jm-1), 1, im) * dlon * fs(1:im,2:jm-1)        &
                                            & + spread(stg_costbl(1:jm-2), 1, im) * dlon * fn(1:im,2:jm-1)        )
            
            where (pin(k) < p_sfc(1:im,2:jm-1))
                div(1:im,2:jm-1,k) = integral(1:im,2:jm-1) / spread(ds(2:jm-1), 1, im)
            else where
                div(1:im,2:jm-1,k) = 0._kp
            endwhere

            div(1, 1,k) = -sum(vft(1:im,   1:2 ,k)) * 0.5_kp * radius * stg_costbl(1)    * dlon
            div(1,jm,k) =  sum(vft(1:im,jm-1:jm,k)) * 0.5_kp * radius * stg_costbl(jm-1) * dlon 


            dsn = 2._kp * pi * radius*radius * ( 1._kp-stg_sintbl(1)   )
            dss = 2._kp * pi * radius*radius * (-1._kp-stg_sintbl(jm-1))
            !
            if (pin(k) < p_sfc(1,1)) then
                div(1:im,1,k) = div(1,1,k) / dsn
            else
                div(1:im,1,k) = 0._kp
            endif


            if (pin(k) < p_sfc(1,jm)) then
                div(1:im,jm,k) = div(1,jm,k) / dss
            else
                div(1:im,jm,k) = 0._kp
            endif

        enddo

    end subroutine get_div


end module estimate_w

