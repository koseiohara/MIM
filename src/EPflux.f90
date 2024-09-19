module EPflux

    use params            , only : kp, ckp, grav, pi, radius
    use fft_utility       , only : fft
    use com_var           , only : im, jm, ko, wmax, pout, rho, alat, costbl
    use mim_var           , only : p_sfc, pt_zm, u_v_x_zm, u_pt_dot_x_zm    , &
                                 & epy, epz_form, epz_uv, epz_ut, epz_uw    , &
                                 & p_pd, z_pd                               , &
                                 & epz_wave, z_pt_wave, p_pt_wave
    use w_advection_expand, only : w_advection_components

    implicit none

    private
    public :: epflux_y, epflux_z_uw, epflux_z_form_wrap, epflux_z_form_wave

    contains


    !
    ! meridional component of EP flux
    !
    ! epy is computed
    !
    subroutine epflux_y()

        epy(1:jm,1:ko) = -u_v_x_zm(1:jm,1:ko) * spread(radius*rho(1:ko), 1, jm) * spread(costbl(1:jm), 2, ko)
        
    end subroutine epflux_y


    !
    ! one of the vertical component of EP flux
    !
    ! epz_uv, epz_ut, and epz_uw are computed
    !
    subroutine epflux_z_uw()

        call w_advection_components(u_v_x_zm(1:jm,1:ko)     , &  !! IN
                                  & u_pt_dot_x_zm(1:jm,1:ko), &  !! IN
                                  & epz_uv(1:jm,1:ko)       , &  !! OUT
                                  & epz_ut(1:jm,1:ko)       , &  !! OUT
                                  & epz_uw(1:jm,1:ko)         )  !! OUT

        epz_uv(1:jm,1:ko) = -spread(rho(1:ko)*radius, 1, jm) * spread(costbl(1:jm), 2, ko) * epz_uv(1:jm,1:ko)
        epz_ut(1:jm,1:ko) = -spread(rho(1:ko)*radius, 1, jm) * spread(costbl(1:jm), 2, ko) * epz_ut(1:jm,1:ko)
        epz_uw(1:jm,1:ko) = -spread(rho(1:ko)*radius, 1, jm) * spread(costbl(1:jm), 2, ko) * epz_uw(1:jm,1:ko)

    end subroutine epflux_z_uw


    subroutine epflux_z_form_wrap()

        call epflux_z_form(p_pd(1:im,1:jm,1:ko), &  !! IN
                         & z_pd(1:im,1:jm,1:ko), &  !! IN
                         & epz_form(1:jm,1:ko)   )  !! OUT
        
    end subroutine epflux_z_form_wrap


    !
    ! Function
    !   get form drag (one of the vertical component of EP flux)
    !   
    ! Arguements (in)
    !   p_pd   : pressure at the p+ surface (if j,k=fixed -> pt surface) [hPa]
    !   z_pt   : geopotential height z(x,y,p+) [m]
    !   p_sfc  : surface pressure [hPa]
    !
    ! Arguements (out)
    !   epz_form : form drag Fz(y,p+) [kg/s^2]
    !
    ! Renamed Variable
    ! p_pd -> isentrop_p
    ! z_pt -> z_pt (not changed)
    ! p_sfc -> (removed)
    ! epz_form -> form_drag
    subroutine epflux_z_form(isentrop_p, isentrop_z, form_drag)
        real(kp), intent(in)  :: isentrop_p(im,jm,ko)
        real(kp), intent(in)  :: isentrop_z(im,jm,ko)
        real(kp), intent(out) :: form_drag(jm,ko)
        real(kp) :: dz
        real(kp) :: pdz
        real(kp) :: p_pd2(im,jm,ko)

        integer  :: i
        integer  :: j
        integer  :: k

        p_pd2(1:im,1:jm,1:ko) = isentrop_p(1:im,1:jm,1:ko)
        do k = 1, ko
            where (p_pd2(1:im,1:jm,k) >= p_sfc(1:im,1:jm))
                p_pd2(1:im,1:jm,k) = p_sfc(1:im,1:jm)
            endwhere
        enddo
        
        ! p_pd2 [hPa] -> p_pd2 [hPa]
        p_pd2(1:im,1:jm,1:ko) = p_pd2(1:im,1:jm,1:ko) * 100._kp

        do k = 1, ko            ! Changed the roop order
            do j = 1, jm
               
                pdz = 0._kp 
                
                !*** calculate p * 2dz  (center difference)

                ! i=1
                dz = isentrop_z(1,j,k) - isentrop_z(im-1,j,k)
                pdz = pdz + p_pd2(im,j,k) * dz
                
                do i = 2, im-1
                    dz = isentrop_z(i+1,j,k) - isentrop_z(i-1,j,k) 
                    pdz = pdz + p_pd2(i,j,k) * dz
                enddo
                
                ! i=im
                dz = isentrop_z(2,j,k) - isentrop_z(im,j,k)
                pdz = pdz + p_pd2(1,j,k) * dz

                ! Iwasaki(1989) e.q.(D.3) M' -> Phi, p' -> p
                ! F_z = a cos(phi) (1/N) sum_{i=1}^N p_i (dz_i/dx)
                !     = a cos(phi) 1/(Ndx) sum_{i=1}^N p_i dz_i
                !     = 1/(4 pai) sum_{i=1}^N p_i 2dz_i
                form_drag(j,k) = pdz * (1._kp / (4._kp * pi))
               
            enddo
        enddo

    end subroutine epflux_z_form


    !
    ! Function
    !   wave number deconposition of the form drag
    !   
    ! Arguements (in)
    !   p_pd   : pressure at the p+ surface (if j,k=fixed -> pt surface) [hPa]
    !   z_pd   : geopotential height z(x,y,p+) [m]
    !   p_sfc  : surface pressure [hPa]
    !   wmax   : maximum wave number to decompose
    !
    ! Arguements (out)
    !   p_pd_wave     : #w component of p
    !   z_pd_wave     : #w component of z
    !   epz_form_wave : #w component of form drag
    !
    ! Renamed Variable
    ! p_pd -> (removed)
    ! z_pd -> (removed)
    ! p_sfc -> (removed)
    ! wmax -> (removed)
    ! p_pd_wave -> p_pt_wave -> (removed)
    ! z_pd_wave -> z_pt_wave -> (removed)
    ! epz_form_wave -> epz_wave -> (removed)
    subroutine epflux_z_form_wave()
        complex(ckp) :: wave(0:im-1)
        complex(ckp) :: kwave(0:im-1)

        integer :: w
        integer :: j
        integer :: k


        !***** wavenumber decomposition *****!
        do k = 1, ko
            do j = 1, jm

                wave(0:im-1) = cmplx(p_pd(1:im,j,k), 0._ckp, kind=ckp)

                call fft(im        , &  !! IN
                       & 1         , &  !! IN
                       & wave(0:im-1))  !! INOUT


                do w = 1, wmax

                    kwave(0:im-1) = cmplx(0._ckp, 0._ckp, kind=ckp)
                    kwave(w) = wave(w)
                    kwave(im-w) = wave(im-w)

                    call fft(im         , &  !! IN
                           & -1         , &  !! IN
                           & kwave(0:im-1))  !! INOUT

                    kwave(0:im-1) = kwave(0:im-1) / real(im, kind=ckp)
                    
                    p_pt_wave(w,1:im,j,k) = real(kwave(0:im-1), kind=kp)

                enddo

                wave(0:im-1) = cmplx(z_pd(1:im,j,k), 0._ckp, kind=ckp)

                call fft(im        , &  !! IN
                       & 1         , &  !! IN
                       & wave(0:im-1))  !! INOUT

                do w = 1, wmax

                    kwave(0:im-1) = cmplx(0._ckp, 0._ckp, kind=ckp)
                    kwave(w) = wave(w)
                    kwave(im-w) = wave(im-w)

                    call fft(im         , &  !! IN
                           & -1         , &  !! IN
                           & kwave(0:im-1))  !! INOUT

                    kwave(0:im-1) = kwave(0:im-1) / real(im, kind=ckp)

                    z_pt_wave(w,1:im,j,k) = real(kwave(0:im-1), kind=kp)

                enddo

            enddo
        enddo

        !***** epz_form *****!
        do w = 1, wmax
            ! +0 is cruicial to force p_pt_wave and z_pt_wave to be processed as passed-by-value like manner
            call epflux_z_form(p_pt_wave(w,1:im,1:jm,1:ko)+0._kp, &  !! IN
                             & z_pt_wave(w,1:im,1:jm,1:ko)+0._kp, &  !! IN
                             & epz_wave(1:jm,1:ko,w)              )  !! OUT
        enddo

    end subroutine epflux_z_form_wave


end module epflux

