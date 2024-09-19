module geopotential

    use params       , only : kp, rkappa, grav, cp, h0, gasr, z_min, z_max
    use com_var      , only : im, jm, km, ko, pin, pout
    use mim_var      , only : alt, t, z, t_dagger, z_zm, p_pd, pt_pds, phi_dagger, z_pd, p_pds
    use status_output, only : warn_write

    implicit none

    private
    public :: geopotential_zm, geopotential_zonal_mean_state

    contains


    !
    ! get_z_pt() - get geopotential height on p+ levels
    !
    !   t      : 3D temperature on p levels [K]
    !   z      : 3D geopotential height on p levels [m]
    !   p_pt   : 3D pressure on p+ levels [hPa]
    !   p_pds  : p+ at the surface [hPa]
    !   alt    : altitude [m] at the surface
    !   z_pt   : 3D geopotential height on p+ levels [m]
    !   z_zm   : zonal mean geopotential height z(y,p+) [m]
    !            for EP Flux calculation
    !
    ! Renamed Variables
    ! alt -> (removed)
    ! t -> (removed)
    ! z -> (removed)
    ! p_pt -> p_pd -> (removed)
    ! p_pds -> (removed)
    ! z_pt -> z_pd -> (removed)
    ! z_zm -> (removed)
    subroutine geopotential_zm()
        real(kp), parameter :: const = rkappa*cp/grav
        real(kp) :: t_zdiff(im,jm,km)
        real(kp) :: hwork1
        real(kp) :: hwork2
        real(kp) :: grad
        integer  :: l
        integer  :: i
        integer  :: k
        integer  :: j
        integer  :: kl
        integer  :: ku
        integer  :: sw  ! interpolation method

        t_zdiff(1:im,1:jm,1)    = (t(1:im,1:jm,2)    - t(1:im,1:jm,1)     ) / (z(1:im,1:jm,2)    - z(1:im,1:jm,1)     )
        t_zdiff(1:im,1:jm,2:km) = (t(1:im,1:jm,2:km) - t(1:im,1:jm,1:km-1)) / (z(1:im,1:jm,2:km) - z(1:im,1:jm,1:km-1))


        do k = 1, ko
            do j = 1, jm
                do i = 1, im

                    ! get kl & ku
                    do l = 2, km
                        if (pin(l-1) < p_pd(i,j,k) .AND. p_pd(i,j,k) <= pin(l)) then
                            kl = l - 1
                            ku = l
                            sw = 1
                            exit
                        endif
                    enddo

                    if (p_pd(i,j,k) <= pin(1)) then
                        kl = 1
                        sw = 2
                    else if (p_pd(i,j,k) > pin(km)) then
                        kl = km
                        sw = 2
                    endif
                    
                    ! interpolation
                    if (sw == 1) then
                        hwork1 = -const * (log(p_pd(i,j,k) / pin(ku))) * t(i,j,ku) &
                             & + 0.5_kp * const * h0 * t_zdiff(i,j,ku) * (log(p_pd(i,j,k) / pin(ku)))**2
                        
                        hwork2 = -const * (log( pin(kl) / pin(ku))) * t(i,j,ku) &
                             & + 0.5_kp * const * h0 * t_zdiff(i,j,ku) * (log(pin(kl) / pin(ku)))**2 
                        
                        z_pd(i,j,k) = z(i,j,ku) + (z(i,j,kl) - z(i,j,ku)) * hwork1 / hwork2
                       
                    else if (sw == 2) then
                        z_pd(i,j,k) = z(i,j,kl) &
                                  & - const * (log(p_pd(i,j,k) / pin(kl))) * t(i,j,kl) &
                                  & + 0.5_kp * const * h0 * t_zdiff(i,j,kl) * (log(p_pd(i,j,k) / pin(kl)))**2 
                        
                    endif

                    ! modify
                    if (z_pd(i,j,k) < 0._kp) then
                        z_pd(i,j,k) = 0._kp
                    endif
                    if (z_pd(i,j,k) <= alt(i,j)) then
                        z_pd(i,j,k) = alt(i,j)
                    endif
                    if (k /= 1) then
                        if (z_pd(i,j,k) > z_pd(i,j,k-1)) then
                            z_pd(i,j,k) = z_pd(i,j,k-1)
                        endif
                    endif
                   
                enddo
              
            enddo
        enddo


        call warn_write(im                  , &  !! IN
                      & jm                  , &  !! IN
                      & ko                  , &  !! IN
                      & z_pd(1:im,1:jm,1:ko), &  !! IN
                      & z_min               , &  !! IN
                      & z_max               , &  !! IN
                      & 'z_pd'              , &  !! IN
                      & 'geopotential_zm()'   )  !! IN

        
        ! z_zm
        z_zm(1:jm,1:ko) = sum(z_pd(1:im,1:jm,1:ko), dim=1) / real(im, kind=kp)

        ! modify
        do j = 1, jm

            grad = 0._kp
            do k = 2, ko

                if (pout(k-1) < p_pds(j) .AND. p_pds(j) <= pout(k)) then
                    grad = (z_zm(j,k-1) - z_zm(j,k-2)) / log(pout(k-1) / pout(k-2))
                endif

                if (p_pds(j) < pout(k)) then
                    z_zm(j,k) = z_zm(j,k-1) - grad * log(pout(k-1) / pout(k))
                endif

            enddo

        enddo

        call warn_write(1                 , &  !! IN
                      & jm                , &  !! IN
                      & ko                , &  !! IN
                      & z_zm(1:jm,1:ko)   , &  !! IN
                      & z_min             , &  !! IN
                      & z_max             , &  !! IN
                      & 'z_zm'            , &  !! IN
                      & 'geopotential_zm()' )  !! IN

    end subroutine geopotential_zm


    !
    ! Function
    !   get geopotential by integrating temperature
    !
    ! Arguements (in)
    !   alt      : altitude
    !   p_pds    : p+s
    !   pt_pds   : pt at the surface in p+ system
    !   t_dagger : T+ (temperature in p+ system)
    !
    ! Arguements (out)
    !   phi_dagger : Phi+ (geopotential in p+ system)
    !
    ! Note
    !   -phi_dagger and z_zm is different:
    !      phi_dagger is obtained from temperature in the zonal mean state.
    !      z_zm is obtained directly from input geopotential data (i.e. z_zm is the zonal mean hight).
    !
    subroutine geopotential_zonal_mean_state()
        real(kp) :: alt_zm(jm)
        real(kp) :: t_dagger_pds(jm)
        real(kp) :: grad
        integer  :: j
        integer  :: k

        alt_zm(1:jm) = sum(alt(1:im,1:jm), dim=1) / real(im, kind=kp)
        t_dagger_pds(1:jm) = pt_pds(1:jm) * (p_pds(1:jm) * 0.001_kp)**rkappa

        ! integrate : R int[ T_dagger ] dlog(p+)
        do j = 1, jm
            do k = ko, 1, -1
               
                if (pout(k) >= p_pds(j))then
                    phi_dagger(j,k) = alt_zm(j) * grav
                else if (k == ko .OR. pout(k+1) > p_pds(j)) then
                !else if (k == ko) then
                    phi_dagger(j,k) = alt_zm(j) * grav + 0.5_kp*gasr * (t_dagger_pds(j) + t_dagger(j,k)) * log(p_pds(j)/pout(k))
                !else if(pout(k+1) > p_pds(j)) then
                !    phi_dagger(j,k) = alt_zm(j) * grav + 0.5_kp*gasr * (t_dagger_pds(j) + t_dagger(j,k)) * log(p_pds(j)/pout(k))
                else 
                    phi_dagger(j,k) = phi_dagger(j,k+1) + 0.5_kp*gasr * (t_dagger(j,k+1) + t_dagger(j,k)) * log(pout(k+1)/pout(k))
                endif
               
            enddo
        enddo


        ! modify
        do j = 1, jm
            grad = 0._kp
            do k = 3, ko

               if (pout(k) >= p_pds(j) .AND. pout(k-1) < p_pds(j)) then
                   grad = (phi_dagger(j,k-1) - phi_dagger(j,k-2)) / log(pout(k-1)/pout(k-2))
               endif

               if (pout(k) > p_pds(j)) then
                   phi_dagger(j,k) = phi_dagger(j,k-1) - grad * log(pout(k-1)/pout(k))
               endif

            enddo
        enddo

    end subroutine geopotential_zonal_mean_state


end module geopotential

