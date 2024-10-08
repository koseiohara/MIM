module mim

    use params                    , only : kp, rkappa, t_min, t_max
    use namelist                  , only : INPUT_OMEGA_FILENAME, WAVE_MAX_NUMBER, OUTPUT_WAVE_FILENAME                      , &
                                         & OUTPUT_ZONAL_FILENAME, OUTPUT_VINT_FILENAME, OUTPUT_GMEAN_FILENAME               , &
                                         & Q_EXIST, Q_COMPS_EXIST
    use com_var                   , only : im, jm, km, ko
    use mim_var                   , only : u, v, omega, pt, pt_dot, u_zm, v_zm, pt_dot_zm                                   , &
                                         & u_u_zm, u_u_x_zm, v_v_zm, v_v_x_zm, u_v_zm, u_v_x_zm                             , &
                                         & u_u_v_zm, v_v_v_zm, u_pt_dot_zm, u_pt_dot_x_zm, v_pt_dot_zm, v_pt_dot_x_zm       , &
                                         & u_u_pt_dot_zm, v_v_pt_dot_zm                                                     , &
                                         & epy, depy, epz_form, depz_form, epz_uv, depz_uv, epz_ut, depz_ut, epz_uw, depz_uw, &
                                         & epz, depz, divf, p_pd, pt_sfc, pt_pds, pt_pdds                                   , &
                                         & t_dagger, phi_dagger, z_pd                                                       , &
                                         & u_past, v_past, omega_past, pt_past, z_pd_past, p_pd_past
    use io_main                   , only : read_alt, read_uvt, read_p_sfc, read_z, read_omega, read_q                       , &
                                         & write_zonal, write_vint, write_gmean, write_wave
    use surface_pressure          , only : get_surface_pressure
    use potential_temperature_3d  , only : get_pt_3d
    use intpl                     , only : intpl_pd_p, intpl_pdd_pd
    use get_levels                , only : getpt_global, getpt_y
    use estimate_w                , only : get_omega
    use estimate_pt_diabatic      , only : estimate_diabatic_heating, get_pt_dot_q
    use zonal_mean_variables      , only : zonal_mean_u, get_v_zm_st, zonal_mean_pt_dot, zonal_mean_square, zonal_mean_cube
    use get_w_zm                  , only : st2w
    use geopotential              , only : geopotential_zm, geopotential_zonal_mean_state
    use get_t_dagger              , only : get_t_zonal_mean_state
    use EPflux                    , only : epflux_y, epflux_z_uw, epflux_z_form_wrap, epflux_z_form_wave
    use EPflux_divergence         , only : epflux_div_y, epflux_div_z
    use Gflux                     , only : gflux_y, gflux_div_y, gflux_z, gflux_div_z
    use height_derivative_products, only : z_yderiv, z_xderiv
    use energy_az                 , only : get_az
    use energy_ae                 , only : get_ae
    use energy_pz                 , only : get_pz
    use energy_k                  , only : get_kz, get_ke
    use energy_tendency           , only : kz_advection, ke_advection, ke_flux_div
    use energy_conv               , only : energy_conversion
    use diabatic                  , only : diabaticHeating_exec

    implicit none
    
    private
    public :: mim_exec

    contains


    subroutine mim_exec(nt)
        integer, intent(in) :: nt
        integer :: tt

        character(16) :: log_format

        call icount_digit(nt, log_format)

        call read_alt()

        write(*,*)

        do tt = 1, nt

            write(*,log_format) 'ICOUNT = ', tt, ' / ', nt

            call read_uvt()

            call read_p_sfc()

            call read_z()

            call read_omega()

            call read_q()

            ! p_dagger and p_dagger_dagger at the surface are computed from p_sfc
            call get_surface_pressure()

            ! 3d potential temperature (pt) and surface potential temperature (pt_sfc)
            call get_pt_3d()

            ! pt_pds  : Potential temperature at the surface in the p_dagger coordinate
            ! pt_pdds : Potential temperature at the surface in the p_dagger_dagger coordinate
            pt_pds(1:jm) = minval(pt_sfc(1:im,1:jm), dim=1)
            pt_pdds(1)   = minval(pt_pds(1:jm)     , dim=1)

            ! get 5 parameters :
            !   dlev (mass weight)
            !   nlev (grid label)
            !   p_pd (pressure on the isentropic surfaces)
            !   p_zm (zonal mean pressure)
            !   pt_zm (zonal mean potential temperature)
            call getpt_global(tt)  !! IN

            if (tt == 1) then
                p_pd_past(1:im,1:jm,1:ko)  = p_pd(1:im,1:jm,1:ko)
                pt_past(1:im,1:jm,1:km)    = pt(1:im,1:jm,1:km)
                omega_past(1:im,1:jm,1:km) = omega(1:im,1:jm,1:km)
                v_past(1:im,1:jm,1:km)     = v(1:im,1:jm,1:km)
                u_past(1:im,1:jm,1:km)     = u(1:im,1:jm,1:km)
            endif

            ! get t_dagger
            !   t_dagger = pt_zm * (p_dagger/1000)^rkappa
            call get_t_zonal_mean_state()

            ! get pd_p (p_dagger on the input pressure surfaces)
            call intpl_pd_p()

            ! get 5 parameters :
            !   dlev_y (mass weight)
            !   nlev_y (grid label)
            !   pd_pdd (p_dagger on the isentropic surfaces)
            !   pd_ym (global mean pressure on the isentropic surfaces)
            !   pt_ym (potential temperature on the p_dagger_dagger surfaces)
            call getpt_y(tt)  !! IN

            ! get pdd_pd (p_dagger_dagger on the p_dagger surfaces)
            call intpl_pdd_pd()

            ! get zonal mean zonal wind
            call zonal_mean_u()  !! IN

           
            ! get zonal mean meridional wind and the streamfunction
            call get_v_zm_st()

            ! get zonal mean vertical velocity (m/s) from the streamfunction
            call st2w()

            if ((.NOT. Q_EXIST) .AND. (.NOT. Q_COMPS_EXIST)) then
                
                if (trim(INPUT_OMEGA_FILENAME) == '') then
                    ! estimate omega by the continuity equation
                    call get_omega()

                    if (tt == 1) then
                        omega_past(1:im,1:jm,1:km) = omega(1:im,1:jm,1:km)
                    endif
                endif

                ! diabatic heating is estimated from the temporal derivative of the potential temperature
                call estimate_diabatic_heating()

            else
                ! temporal derivative of the potential temperature is estmated from the diabatic heating
                call get_pt_dot_q()
            endif

            ! get zonal mean temporal derivative potential temperature
            call zonal_mean_pt_dot()

            ! get zonal mean of u^2 and (u')^2
            call zonal_mean_square(u(1:im,1:jm,1:km)  , &  !! IN
                                 & u(1:im,1:jm,1:km)  , &  !! IN
                                 & u_zm(1:jm,1:ko)    , &  !! IN
                                 & u_zm(1:jm,1:ko)    , &  !! IN
                                 & u_u_zm(1:jm,1:ko)  , &  !! OUT
                                 & u_u_x_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of uv and (u'v')^2
            call zonal_mean_square(u(1:im,1:jm,1:km)  , &  !! IN
                                 & v(1:im,1:jm,1:km)  , &  !! IN
                                 & u_zm(1:jm,1:ko)    , &  !! IN
                                 & v_zm(1:jm,1:ko)    , &  !! IN
                                 & u_v_zm(1:jm,1:ko)  , &  !! OUT
                                 & u_v_x_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of v^2 and (v')^2
            call zonal_mean_square(v(1:im,1:jm,1:km)  , &  !! IN
                                 & v(1:im,1:jm,1:km)  , &  !! IN
                                 & v_zm(1:jm,1:ko)    , &  !! IN
                                 & v_zm(1:jm,1:ko)    , &  !! IN
                                 & v_v_zm(1:jm,1:ko)  , &  !! OUT
                                 & v_v_x_zm(1:jm,1:ko)  )  !! OUT
            
            ! get zonal mean of u^2v
            call zonal_mean_cube(u(1:im,1:jm,1:km)  , &  !! IN
                               & u(1:im,1:jm,1:km)  , &  !! IN
                               & v(1:im,1:jm,1:km)  , &  !! IN
                               & u_u_v_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of v^3
            call zonal_mean_cube(v(1:im,1:jm,1:km)  , &  !! IN
                               & v(1:im,1:jm,1:km)  , &  !! IN
                               & v(1:im,1:jm,1:km)  , &  !! IN
                               & v_v_v_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of u(Dpt/Dt) and u'(DptDt)'
            call zonal_mean_square(u(1:im,1:jm,1:km)       , &  !! IN
                                 & pt_dot(1:im,1:jm,1:km)  , &  !! IN
                                 & u_zm(1:jm,1:ko)         , &  !! IN
                                 & pt_dot_zm(1:jm,1:ko)    , &  !! IN
                                 & u_pt_dot_zm(1:jm,1:ko)  , &  !! OUT
                                 & u_pt_dot_x_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of v(Dpt/Dt) and v'(DptDt)'
            call zonal_mean_square(v(1:im,1:jm,1:km)       , &  !! IN
                                 & pt_dot(1:im,1:jm,1:km)  , &  !! IN
                                 & v_zm(1:jm,1:ko)         , &  !! IN
                                 & pt_dot_zm(1:jm,1:ko)    , &  !! IN
                                 & v_pt_dot_zm(1:jm,1:ko)  , &  !! OUT
                                 & v_pt_dot_x_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of u^2(Dpt/Dt)
            call zonal_mean_cube(u(1:im,1:jm,1:km)       , &  !! IN
                               & u(1:im,1:jm,1:km)       , &  !! IN
                               & pt_dot(1:im,1:jm,1:km)  , &  !! IN
                               & u_u_pt_dot_zm(1:jm,1:ko)  )  !! OUT

            ! get zonal mean of v^2(Dpt/Dt)
            call zonal_mean_cube(v(1:im,1:jm,1:km)       , &  !! IN
                               & v(1:im,1:jm,1:km)       , &  !! IN
                               & pt_dot(1:im,1:jm,1:km)  , &  !! IN
                               & v_v_pt_dot_zm(1:jm,1:ko)  )  !! OUT


            ! get 3d geopotential on the isentropic surfaces and its zonal mean
            call geopotential_zm()
            if (tt == 1) then
                z_pd_past(1:im,1:jm,1:ko) = z_pd(1:im,1:jm,1:ko)
            endif

            ! get y-component of the EP flux
            call epflux_y()

            ! get meridional divergence of the EP flux
            call epflux_div_y(epy(1:jm,1:ko) , &  !! IN
                            & depy(1:jm,1:ko)  )  !! OUT

            ! get form drag component of the EP flux
            call epflux_z_form_wrap()

            ! get vertical divergence of the form drag
            call epflux_div_z(epz_form(1:jm,1:ko),  &  !! IN
                            & depz_form(1:jm,1:ko)  )  !! OUT

            ! wavenumber decomposition of the form drag
            if (trim(OUTPUT_WAVE_FILENAME) /= '' .AND. WAVE_MAX_NUMBER > 0) then
                call epflux_z_form_wave()
            endif

            ! get u'w' component of the EP flux. the term can be decomposed into u'v' and u'(Dpt/Dt)' terms
            ! all of these three terms are computed in this routine
            call epflux_z_uw()

            ! get divergence of u'w' term and the other two terms
            call epflux_div_z(epz_uw(1:jm,1:ko) , &  !! IN
                            & depz_uw(1:jm,1:ko)  )  !! OUT

            call epflux_div_z(epz_ut(1:jm,1:ko) , &  !! IN
                            & depz_ut(1:jm,1:ko)  )  !! OUT

            call epflux_div_z(epz_uv(1:jm,1:ko) , &  !! IN
                            & depz_uv(1:jm,1:ko)  )  !! OUT

            ! get vertical component of the total EP flux and its vertical divergence
            epz(1:jm,1:ko) = epz_form(1:jm,1:ko) + epz_uw(1:jm,1:ko)
            !call epflux_div_z(epz(1:jm,1:ko) , &  !! IN
            !                & depz(1:jm,1:ko)  )  !! OUT
            depz(1:jm,1:ko) = depz_form(1:jm,1:ko) + depz_uw(1:jm,1:ko)

            divf(1:jm,1:ko) = depy(1:jm,1:ko) + depz(1:jm,1:ko)

            ! get y-component of the G flux
            call gflux_y()

            ! get meridional divergence of the G flux
            call gflux_div_y()

            ! get vertical component of the G flux
            call gflux_z()

            ! get vertical divergence of the G flux
            call gflux_div_z()

            ! get Phi_dagger
            call geopotential_zonal_mean_state()

            !if (tt == 1) then
            !    phi_dagger_past(1:jm,1:ko) = phi_dagger(1:jm,1:ko)
            !endif

            ! get zonal mean of dz/dy and zonal mean of v*(dz/dt)
            call z_yderiv()

            ! get zonal mean of u*(dz/dx)
            call z_xderiv()


            !----- ENERGY -----!
            ! get zonal mean kinetic energy
            call get_kz()

            ! get eddy kinetic energy
            call get_ke()

            ! get zonal potential energy (not the available potential energy)
            call get_pz()

            ! get zonal available potential energy
            call get_az()

            ! get eddy available potential energy
            call get_ae()

            ! get C(Az,Kz), C(Kz,ke), C(kz,Ae), C(Kz,W), and C(Ae,Ke)
            call energy_conversion()

            ! get Kz advection by v_zm and w_zm
            call kz_advection()

            ! get Ke advection by v and w
            call ke_advection()

            ! get Ke flux by the EP flux
            call ke_flux_div()

            ! get all variables related to the diabatic heating
            call diabaticHeating_exec()


            !----- OUTPUT -----!
            if (trim(OUTPUT_ZONAL_FILENAME) /= '') then
                call write_zonal(tt)  !! IN
            endif

            if (trim(OUTPUT_VINT_FILENAME) /= '') then
                call write_vint(tt)   !! IN
            endif

            if (trim(OUTPUT_GMEAN_FILENAME) /= '') then
                call write_gmean(tt)  !! IN
            endif

            if (trim(OUTPUT_WAVE_FILENAME) /= '') then
                call write_wave(tt)  !! IN
            endif

            pt_past(1:im,1:jm,1:km)    = pt(1:im,1:jm,1:km)
            u_past(1:im,1:jm,1:km)     = u(1:im,1:jm,1:km)
            v_past(1:im,1:jm,1:km)     = v(1:im,1:jm,1:km)
            omega_past(1:im,1:jm,1:km) = omega(1:im,1:jm,1:km)
            z_pd_past(1:im,1:jm,1:km)  = z_pd(1:im,1:jm,1:km)
            p_pd_past(1:im,1:jm,1:km)  = p_pd(1:im,1:jm,1:km)

        enddo

    end subroutine mim_exec


    subroutine icount_digit(nt, log_format)
        integer     , intent(in)  :: nt
        character(*), intent(out) :: log_format

        integer :: digit
        integer :: i
        
        character(8) :: digit_char

        digit = 50

        do i = 0, 10
            if (nt/(10**i) == 0) then
                digit = i
                exit
            endif
        enddo

        write(digit_char,'(I0)') digit

        log_format = '(2(A,I' // trim(digit_char) // '))'

    end subroutine icount_digit


end module mim

