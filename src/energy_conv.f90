module energy_conv

    use params       , only : kp, rkappa, grav, radius, h0, cp, gasr, econv_min, econv_max
    use com_var      , only : im, jm, ko, pout, alat, costbl, tantbl
    use mim_var      , only : alt, u_zm, v_zm, pt_zm, u_u_x_zm, depy, depz_form, depz_uw, dgy, dgz , &
                            & c_az_kz, c_kz_ae, c_kz_ae_u, c_kz_ae_v, c_ae_ke, c_ae_ke_u, c_ae_ke_v, &
                            & c_kz_ke, c_kz_ke_uy, c_kz_ke_uz, c_kz_ke_vy, c_kz_ke_vz, c_kz_ke_tan , &
                            & c_kz_w                                                               , &
                            & pt_sfc, dz_dlat_zm, v_dz_dlat_zm, u_dz_dlon_zm, phi_dagger, p_pds
    use status_output, only : warn_write

    implicit none

    private
    public :: energy_conversion

    contains


    subroutine energy_conversion()

        call energy_conv_az_kz()

        call energy_conv_kz_ae()

        call energy_conv_ae_ke()

        call energy_conv_kz_ke()

        c_kz_w(1:jm,1:ko) = c_kz_ke(1:jm,1:ko) + c_kz_ae(1:jm,1:ko)

    end subroutine energy_conversion


    subroutine energy_conv_az_kz()

        ! Mountain modify (if necessary)
        call mount_modify(c_az_kz(1:jm,1:ko))  !! OUT

        call warn_write(1                    , &  !! IN
                      & jm                   , &  !! IN
                      & ko                   , &  !! IN
                      & c_az_kz(1:jm,1:ko)   , &  !! IN
                      & econv_min            , &  !! IN
                      & econv_max            , &  !! IN
                      & 'c_az_kz'            , &  !! IN
                      & 'energy_conv_az_kz()'  )  !! IN

    end subroutine energy_conv_az_kz


    subroutine energy_conv_kz_ae()

        c_kz_ae_u(1:jm,1:ko) = -u_zm(1:jm,1:ko) * depz_form(1:jm,1:ko)

        c_kz_ae_v(1:jm,1:ko) = v_zm(1:jm,1:ko) * dz_dlat_zm(1:jm,1:ko) * grav / radius + c_az_kz(1:jm,1:ko)

        c_kz_ae(1:jm,1:ko) = c_kz_ae_u(1:jm,1:ko) + c_kz_ae_v(1:jm,1:ko)

        call warn_write(1                    , &  !! IN
                      & jm                   , &  !! IN
                      & ko                   , &  !! IN
                      & c_kz_ae(1:jm,1:ko)   , &  !! IN
                      & econv_min            , &  !! IN
                      & econv_max            , &  !! IN
                      & 'c_kz_ae'            , &  !! IN
                      & 'energy_conv_kz_ae()'  )  !! IN

    end subroutine energy_conv_kz_ae


    subroutine energy_conv_ae_ke()
        integer :: k

        do k = 1, ko
            c_ae_ke_u(     1,k) = 0._kp
            c_ae_ke_u(2:jm-1,k) = -u_dz_dlon_zm(2:jm-1,k) * (grav / radius) / costbl(2:jm-1) + c_kz_ae_u(2:jm-1,k)
            c_ae_ke_u(    jm,k) = 0._kp
        enddo
        
        c_ae_ke_v(1:jm,1:ko) = -(v_dz_dlat_zm(1:jm,1:ko) - v_zm(1:jm,1:ko) * dz_dlat_zm(1:jm,1:ko)) * grav / radius
        c_ae_ke(1:jm,1:ko)   = c_ae_ke_u(1:jm,1:ko) + c_ae_ke_v(1:jm,1:ko)

        call warn_write(1                    , &  !! IN
                      & jm                   , &  !! IN
                      & ko                   , &  !! IN
                      & c_ae_ke(1:jm,1:ko)   , &  !! IN
                      & econv_min            , &  !! IN
                      & econv_max            , &  !! IN
                      & 'c_ae_ke'            , &  !! IN
                      & 'energy_conv_ae_ke()'  )  !! IN

    end subroutine energy_conv_ae_ke


    subroutine energy_conv_kz_ke()

        c_kz_ke_uy(1:jm,1:ko)  = -u_zm(1:jm,1:ko) * depy(1:jm,1:ko)
        c_kz_ke_uz(1:jm,1:ko)  = -u_zm(1:jm,1:ko) * depz_uw(1:jm,1:ko)
        c_kz_ke_vy(1:jm,1:ko)  = -v_zm(1:jm,1:ko) * dgy(1:jm,1:ko)
        c_kz_ke_vz(1:jm,1:ko)  = -v_zm(1:jm,1:ko) * dgz(1:jm,1:ko)
        c_kz_ke_tan(1:jm,1:ko) = v_zm(1:jm,1:ko) * u_u_x_zm(1:jm,1:ko) * spread(tantbl(1:jm), 2, ko) / radius
        c_kz_ke_tan(   1,1:ko) = 0._kp
        c_kz_ke_tan(  jm,1:ko) = 0._kp

        c_kz_ke(1:jm,1:ko) =  c_kz_ke_uy(1:jm,1:ko) + c_kz_ke_uz(1:jm,1:ko) &
                         & +  c_kz_ke_vy(1:jm,1:ko) + c_kz_ke_vz(1:jm,1:ko) &
                         & + c_kz_ke_tan(1:jm,1:ko)

        call warn_write(1                    , &  !! IN
                      & jm                   , &  !! IN
                      & ko                   , &  !! IN
                      & c_kz_ke(1:jm,1:ko)   , &  !! IN
                      & econv_min            , &  !! IN
                      & econv_max            , &  !! IN
                      & 'c_kz_ke'            , &  !! IN
                      & 'energy_conv_kz_ke()'  )  !! IN

    end subroutine energy_conv_kz_ke


    subroutine mount_modify(c_az_kz_modified)
        real(kp), intent(out) :: c_az_kz_modified(jm,ko)

        real(kp) :: phi_dagger_modify(jm, ko)
        real(kp) :: pout_modify(jm, ko)
        real(kp) :: v_zm_modify(jm, ko)
        real(kp) :: pt_zm_modify(jm, ko)

        real(kp) :: work1_p_sfc(im, jm)
        real(kp) :: work2_p_sfc(im, jm)
        real(kp) :: work2_p_sfc_zm(jm)
        real(kp) :: p_sfc_modify(im, jm)
        real(kp) :: p_sfc_modify_max(jm)
        
        real(kp) :: surface_pt(jm)
        real(kp) :: alt_min(jm)
        real(kp) :: grad
        real(kp) :: normalized_surface_pressure_zm(jm)
        
        real(kp) :: d
        integer  :: i
        integer  :: j
        integer  :: k
        integer  :: l
        integer  :: ju
        integer  :: jl
        integer  :: ku
        integer  :: kl


        ! work1_p_sfc : 2-dimensional surface pressure normalized by p+s
        !               estimated from altitude
        normalized_surface_pressure_zm(1:jm) = sum(exp(-alt(1:im,1:jm)*(1._kp/h0)), dim=1) / real(im, kind=kp)
        do j = 1, jm
            if (normalized_surface_pressure_zm(j) /= 0._kp) then
                work1_p_sfc(1:im,j) = p_pds(j) * exp(-alt(1:im,j)*(1._kp/h0)) / normalized_surface_pressure_zm(j)
            else
                work1_p_sfc(1:im,j) = 0._kp
            endif
        enddo

        
        ! p_sfc_modify : 2-dimensional surface pressure normalized by p+s
        !                estimated from phi_dagger
        do j = 1, jm
            do i = 1, im
                do k = 2, ko

                    if (work1_p_sfc(i,j) <= pout(k)) then
                       work2_p_sfc(i,j) = pout(k-1) &
                                      & + (grav * alt(i,j) - phi_dagger(j,k-1)) &
                                      & * (pout(k) - pout(k-1)) / (phi_dagger(j,k) - phi_dagger(j,k-1))
                       exit
                    else if (k == ko) then
                       work2_p_sfc(i,j) = pout(ko) &
                                      & + (grav * alt(i,j) - phi_dagger(j,ko)) &
                                      & * (pout(ko) - pout(ko-1)) / (phi_dagger(j,ko) - phi_dagger(j,ko-1))
                       exit
                    endif

                enddo
            enddo
        enddo

        
        work2_p_sfc_zm(1:jm) = sum(work2_p_sfc(1:im,1:jm), dim=1) / real(im, kind=kp)
        
        ! OPENMP
        do j = 1, jm
            p_sfc_modify(1:im,j) = work2_p_sfc(1:im,j) * p_pds(j) / work2_p_sfc_zm(j)
        enddo
        

        ! pout_modify : actual p+
        !               pout_modify != pout only near the surface
        do k = 1, ko
            pout_modify(1:jm,k) = 0._kp
            do i = 1, im
                pout_modify(1:jm,k) = pout_modify(1:jm,k) + min(pout(k), p_sfc_modify(i,1:jm))
            enddo
            pout_modify(1:jm,k) = pout_modify(1:jm,k) / real(im, kind=kp)
        enddo
        

        ! interpolate (pout levels -> pout_modify levels)

        ! interpolate theta & v_zm  
        do j = 1, jm

            ! k=1
            v_zm_modify(j,1)  = v_zm(j,1)     
            pt_zm_modify(j,1) = pt_zm(j,1)       

            do k = 2, ko
                do l = ko, 2, -1

                    if (pout_modify(j,k) == pout(l)) then
                        v_zm_modify(j,k)  =  v_zm(j,l)
                        pt_zm_modify(j,k) = pt_zm(j,l)

                    else if (pout_modify(j,k) < pout(l)   .AND. &
                           & pout_modify(j,k) > pout(l-1)       ) then

                        d = (pout_modify(j,k) - pout(l-1)) / (pout(l) - pout(l-1))
                        v_zm_modify(j,k)  =  v_zm(j,l-1) * (1._kp-d) +  v_zm(j,l) * d
                        pt_zm_modify(j,k) = pt_zm(j,l-1) * (1._kp-d) + pt_zm(j,l) * d

                    endif

                enddo
            enddo

        enddo
        
        alt_min(1:jm) = minval(alt(1:im,1:jm), dim=1)
        p_sfc_modify_max(1:jm) = maxval(p_sfc_modify(1:im,1:jm), dim=1)
        surface_pt(1:jm) = minval(pt_sfc(1:im,1:jm), dim=1)

        !end naiso
        

        ! modify phi_dagger

        !pt_new -> pt_zm  &  pa_sea -->  pasmax   
        do k = ko, 1, -1
            do j = 1, jm

               if (pout(k) >= p_sfc_modify_max(j)) then
                   phi_dagger_modify(j,k) = alt_min(j) * grav

               else if (k == ko) then
                   phi_dagger_modify(j,k) = alt_min(j) * grav                                         &
                                        & + gasr                                                      &
                                        & * (surface_pt(j) * (p_sfc_modify_max(j)*1.E-3_kp)**rkappa   &
                                        & + pt_zm_modify(j,k) * (pout(k)*1.E-3_kp)**rkappa ) * 0.5_kp &
                                        & * log(p_sfc_modify_max(j) / pout(k))

               else if (pout(k+1) >= p_sfc_modify_max(j)) then
                   phi_dagger_modify(j,k) = alt_min(j) * grav                                         &
                                        & + gasr                                                      &
                                        & * (surface_pt(j) * (p_sfc_modify_max(j)*1.E-3_kp)**rkappa   &
                                        & + pt_zm_modify(j,k) * (pout(k)*1.E-3_kp)**rkappa ) * 0.5_kp &
                                        & * log(p_sfc_modify_max(j) / pout(k))

               else 
                   phi_dagger_modify(j,k) = phi_dagger_modify(j,k+1)                                  &
                                        & + gasr                                                      &
                                        & * (pt_zm_modify(j,k+1) * (pout(k+1)*1.E-3_kp)**rkappa       &
                                        & + pt_zm_modify(j,k) * (pout(k)*1.E-3_kp)**rkappa ) * 0.5_kp &
                                        & * log(pout(k+1) / pout(k))
               endif

            enddo
        enddo
        
        
        do j = 1, jm
            grad = 0._kp
            do k = 3, ko

                if (pout(k) >= p_sfc_modify_max(j) .AND. pout(k-1) < p_sfc_modify_max(j)) then
                    grad = ( phi_dagger_modify(j,k-1) - phi_dagger_modify(j,k-2) ) &
                       & / log(pout(k-1) / pout(k-2))
                endif

                if (pout(k) > p_sfc_modify_max(j)) then
                    phi_dagger_modify(j,k) = phi_dagger_modify(j,k-1) - grad * log(pout(k-1)/pout(k))
                endif

            enddo
        enddo


        ! modify C(Az,Kz)
        ! Probably, d(pout_modify)/d(pout) is multipled 
        ! for the accurate vertical integration.
        ! c_az_kz_modified may be on the pout_modify levels, not on the pout levels
        
        !v_new->v_zm
        do k = 1, ko
            kl = k - 1
            ku = k + 1
            if (k == 1) then
                kl = 1
            endif
            if (k == ko) then
                ku = ko
            endif
            
            do j = 1, jm
                jl = j - 1
                ju = j + 1
                if (j == 1) then
                    jl = 1
                endif
                if (j == jm) then
                    ju = jm
                endif
                
                c_az_kz_modified(j,k) = -v_zm_modify(j,k)                                        &
                                    & * (phi_dagger_modify(ju,k) - phi_dagger_modify(jl,k))      &
                                    & * (pout_modify(j,kl) - pout_modify(j,ku))                  &
                                    & / ((alat(ju) - alat(jl)) * (pout(kl) - pout(ku)) * radius)
               
            enddo
        enddo

    end subroutine mount_modify


end module energy_conv

