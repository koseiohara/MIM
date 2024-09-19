module diabatic

    use params         , only : kp, rkappa, cp
    use namelist       , only : INPUT_Q_FILENAME, INPUT_SHORTWAVE_FILENAME, INPUT_LONGWAVE_FILENAME        , &
                              & INPUT_LHR_LARGE_FILENAME, INPUT_LHR_CONV_FILENAME, INPUT_DIFFUSION_FILENAME, &
                              & Q_EXIST, Q_COMPS_EXIST
    use com_var        , only : im, jm, km, ko, pin, pout
    use mim_var        , only : pd_p, pdd_pd                                                                            , &
                              & q_3d, q_shortwave_3d, q_longwave_3d, q_lhr_large_3d, q_lhr_conv_3d, q_diffusion_3d      , &
                              & q_zm, q_shortwave_zm, q_longwave_zm, q_lhr_large_zm, q_lhr_conv_zm, q_diffusion_zm      , &
                              & qgz_zm, qgz_shortwave_zm, qgz_longwave_zm, qgz_lhr_large_zm, qgz_lhr_conv_zm            , &
                              & qgz_diffusion_zm                                                                        , &
                              & qz_vint, qz_shortwave_vint, qz_longwave_vint, qz_lhr_large_vint, qz_lhr_conv_vint       , &
                              & qz_diffusion_vint                                                                       , &
                              & qz_gmean, qz_shortwave_gmean, qz_longwave_gmean, qz_lhr_large_gmean, qz_lhr_conv_gmean  , &
                              & qz_diffusion_gmean                                                                      , &
                              & qe_zm, qe_shortwave_zm, qe_longwave_zm, qe_lhr_large_zm, qe_lhr_conv_zm, qe_diffusion_zm, &
                              & p_pds, p_pdds
    use integral       , only : integral_p
    use isentropic_mean, only : zonalMean, globalMean

    implicit none

    private
    public :: diabaticHeating_exec

    contains


    subroutine diabaticHeating_exec()

        if (trim(INPUT_SHORTWAVE_FILENAME) /= '') then
            call diabaticHeating(q_shortwave_3d(1:im,1:jm,1:km), &  !! IN
                               & q_shortwave_zm(1:jm,1:ko)     , &  !! OUT
                               & qgz_shortwave_zm(1:jm,1:ko)   , &  !! OUT
                               & qe_shortwave_zm(1:jm,1:ko)    , &  !! OUT
                               & qz_shortwave_vint(1:jm)       , &  !! OUT
                               & qz_shortwave_gmean(1:1)         )  !! OUT
        endif

        if (trim(INPUT_LONGWAVE_FILENAME) /= '') then
            call diabaticHeating(q_longwave_3d(1:im,1:jm,1:km), &  !! IN
                               & q_longwave_zm(1:jm,1:ko)     , &  !! OUT
                               & qgz_longwave_zm(1:jm,1:ko)   , &  !! OUT
                               & qe_longwave_zm(1:jm,1:ko)    , &  !! OUT
                               & qz_longwave_vint(1:jm)       , &  !! OUT
                               & qz_longwave_gmean(1:1)         )  !! OUT
        endif

        if (trim(INPUT_LHR_LARGE_FILENAME) /= '') then
            call diabaticHeating(q_lhr_large_3d(1:im,1:jm,1:km), &  !! IN
                               & q_lhr_large_zm(1:jm,1:ko)     , &  !! OUT
                               & qgz_lhr_large_zm(1:jm,1:ko)   , &  !! OUT
                               & qe_lhr_large_zm(1:jm,1:ko)    , &  !! OUT
                               & qz_lhr_large_vint(1:jm)       , &  !! OUT
                               & qz_lhr_large_gmean(1:1)         )  !! OUT
        endif

        if (trim(INPUT_LHR_CONV_FILENAME) /= '') then
            call diabaticHeating(q_lhr_conv_3d(1:im,1:jm,1:km), &  !! IN
                               & q_lhr_conv_zm(1:jm,1:ko)     , &  !! OUT
                               & qgz_lhr_conv_zm(1:jm,1:ko)   , &  !! OUT
                               & qe_lhr_conv_zm(1:jm,1:ko)    , &  !! OUT
                               & qz_lhr_conv_vint(1:jm)       , &  !! OUT
                               & qz_lhr_conv_gmean(1:1)         )  !! OUT
        endif

        if (trim(INPUT_DIFFUSION_FILENAME) /= '') then
            call diabaticHeating(q_diffusion_3d(1:im,1:jm,1:km), &  !! IN
                               & q_diffusion_zm(1:jm,1:ko)     , &  !! OUT
                               & qgz_diffusion_zm(1:jm,1:ko)   , &  !! OUT
                               & qe_diffusion_zm(1:jm,1:ko)    , &  !! OUT
                               & qz_diffusion_vint(1:jm)       , &  !! OUT
                               & qz_diffusion_gmean(1:1)         )  !! OUT
        endif

        if ((.NOT. Q_EXIST) .AND. Q_COMPS_EXIST) then
            q_zm(1:jm,1:ko)   = q_shortwave_zm(1:jm,1:ko) + q_longwave_zm(1:jm,1:ko) + q_lhr_large_zm(1:jm,1:ko) &
                            & + q_lhr_conv_zm(1:jm,1:ko)  + q_diffusion_zm(1:jm,1:ko)

            qgz_zm(1:jm,1:ko) = qgz_shortwave_zm(1:jm,1:ko) + qgz_longwave_zm(1:jm,1:ko) + qgz_lhr_large_zm(1:jm,1:ko) &
                            & + qgz_lhr_conv_zm(1:jm,1:ko)  + qgz_diffusion_zm(1:jm,1:ko)

            qe_zm(1:jm,1:ko)  = qe_shortwave_zm(1:jm,1:ko) + qe_longwave_zm(1:jm,1:ko) + qe_lhr_large_zm(1:jm,1:ko) &
                            & + qe_lhr_conv_zm(1:jm,1:ko)  + qe_diffusion_zm(1:jm,1:ko)

            qz_vint(1:jm)     = qz_shortwave_vint(1:jm) + qz_longwave_vint(1:jm) + qz_lhr_large_vint(1:jm) &
                            & + qz_lhr_conv_vint(1:jm)  + qz_diffusion_vint(1:jm)

            qz_gmean(1)       = qz_shortwave_gmean(1) + qz_longwave_gmean(1) + qz_lhr_large_gmean(1) &
                            & + qz_lhr_conv_gmean(1)  + qz_diffusion_gmean(1)
        else
            call diabaticHeating(q_3d(1:im,1:jm,1:km), &  !! IN
                               & q_zm(1:jm,1:ko)     , &  !! OUT
                               & qgz_zm(1:jm,1:ko)   , &  !! OUT
                               & qe_zm(1:jm,1:ko)    , &  !! OUT
                               & qz_vint(1:jm)       , &  !! OUT
                               & qz_gmean(1:1)         )  !! OUT
        endif

    end subroutine diabaticHeating_exec


    ! Arguments :
    !     1. 3d distribution of diabatic heating
    !     2. zonal mean diabatic heating
    !     3. diabatic heating to the zonal mean state
    !     4. Generation rate of eddy available potential energy
    !     5. Generation rate of zonal available potential energy
    subroutine diabaticHeating(heating_3d, heating_zm, heating_gz_zm, eddy_gen, zonal_gen_vint, zonal_gen_gmean)
        real(kp), intent(in)  :: heating_3d(im,jm,km)
        real(kp), intent(out) :: heating_zm(jm,ko)
        real(kp), intent(out) :: heating_gz_zm(jm,ko)
        real(kp), intent(out) :: eddy_gen(jm,ko)
        real(kp), intent(out) :: zonal_gen_vint(jm)
        real(kp), intent(out) :: zonal_gen_gmean(1)

        real(kp) :: heating_exner_zm(jm,ko)
        real(kp) :: zonal_gen_zm(jm,ko)
        real(kp) :: zonal_gen(ko)

        ! Zonal Mean Diabatic Heanting
        call zonalMean(heating_3d(1:im,1:jm,1:km), &  !! IN
                     & heating_zm(1:jm,1:ko)       )  !! OUT

        call heating_per_Exner(heating_3d(1:im,1:jm,1:km), &  !! IN
                             & heating_exner_zm(1:jm,1:ko) )  !! OUT

        call heating_ZonalMeanState(heating_exner_zm(1:jm,1:ko), &  !! IN
                                  & heating_gz_zm(1:jm,1:ko)     )  !! OUT

        call heating_Eddy_highPrecision(heating_3d(1:im,1:jm,1:km), &  !! IN
                                      & eddy_gen(1:jm,1:ko)         )  !! OUT

        call heating_Zonal_highPrecision(heating_exner_zm(1:jm,1:ko), &  !! IN
                                       & zonal_gen_zm(1:jm,1:ko)    , &  !! OUT
                                       & zonal_gen(1:ko)              )  !! OUT

        call integral_p(jm                     , &  !! IN
                      & ko                     , &  !! IN
                      & p_pds(1:jm)            , &  !! IN
                      & zonal_gen_zm(1:jm,1:ko), &  !! IN
                      & zonal_gen_vint(1:jm)     )  !! OUT

        call integral_p(1                 , &  !! IN   size in lat-direction
                      & ko                , &  !! IN   size in p-direction
                      & p_pdds(1)         , &  !! IN   p_dagger_dagger at the surface
                      & zonal_gen(1:ko)   , &  !! IN   vertical profile of parameter
                      & zonal_gen_gmean(1)  )  !! OUT  vertically integrated parameter

    end subroutine diabaticHeating


    ! Zonal mean diabatic heating divided by Exner function bar{[Q/Exner(p)]^*}
    subroutine heating_per_Exner(heating_3d, heating_exner_zm)
        real(kp), intent(in)  :: heating_3d(im,jm,km)
        !real(4), intent(out) :: heating_gz_zm(jm,ko)
        real(kp), intent(out) :: heating_exner_zm(jm,ko)

        real(kp) :: heating_exner_3d(im,jm,km)

        integer :: k

        do k = 1, km
            heating_exner_3d(1:im,1:jm,k) = heating_3d(1:im,1:jm,k) / (cp * (pin(k)*1.E-3_kp)**rkappa)
        enddo

        ! Zonal Mean
        call zonalMean(heating_exner_3d(1:im,1:jm,1:km), &  !! IN
                     & heating_exner_zm(1:jm,1:ko)       )  !! OUT

    end subroutine heating_per_Exner


    ! Diabatic heating to the zonal mean state bar{[Q/Exner(p)]^*}Exner(p_dagger)
    subroutine heating_ZonalMeanState(heating_exner_zm, heating_gz_zm)
        real(kp), intent(in)  :: heating_exner_zm(jm,ko)
        real(kp), intent(out) :: heating_gz_zm(jm,ko)

        integer :: k

        do k = 1, ko
            heating_gz_zm(1:jm,k) = heating_exner_zm(1:jm,k) * cp * (pout(k)*1.E-3_kp)**rkappa
        enddo

    end subroutine heating_ZonalMeanState


    ! Generation rate of eddy available potential energy bar{Q^*} - bar{[Q/Exner(p)]^*}Exner(p_dagger)
    ! This should not be used beacause of bad precision
    subroutine heating_Eddy_simple(heating_zm, heating_gz_zm, eddy_generation)
        real(kp), intent(in)  :: heating_zm(jm,ko)
        real(kp), intent(in)  :: heating_gz_zm(jm,ko)
        real(kp), intent(out) :: eddy_generation(jm,ko)

        eddy_generation(1:jm,1:ko) = heating_zm(1:jm,1:ko) - heating_gz_zm(1:jm,1:ko)

    end subroutine heating_Eddy_simple


    ! Generation rate of eddy available potential energy bar{Q^*} - bar{[Q/Exner(p)]^*}Exner(p_dagger)
    ! This should not be used beacause of bad precision
    subroutine heating_Eddy_highPrecision(heating_3d, eddy_generation)
        real(kp), intent(in)  :: heating_3d(im,jm,km)
        real(kp), intent(out) :: eddy_generation(jm,ko)

        real(kp) :: work_eddy_generation(im,jm,km)

        integer :: k

        do k = 1, km
            work_eddy_generation(1:im,1:jm,k) = heating_3d(1:im,1:jm,k) * (1._kp - (pd_p(1:im,1:jm,k) / pin(k))**rkappa)
        enddo

        call zonalMean(work_eddy_generation(1:im,1:jm,1:km), &  !! IN
                     & eddy_generation(1:jm,1:ko)            )  !! OUT

    end subroutine heating_Eddy_highPrecision


    ! Diabatic heating to the ground state bar{bar{[Q/Exner(o)]partial p / partial p_dagger_dagger}}Exner(p_dagger_dagger)
    subroutine heating_GroundState(heating_exner_zm, heating_pdd)
        real(kp), intent(in)  :: heating_exner_zm(jm,ko)
        real(kp), intent(out) :: heating_pdd(ko)

        real(kp) :: work_heating_pdd(jm,ko)

        work_heating_pdd(1:jm,1:ko) = heating_exner_zm(1:jm,1:ko) * cp * (pdd_pd(1:jm,1:ko)*1.E-3_kp)**rkappa

        call globalMean(work_heating_pdd(1:jm,1:ko), &  !! IN
                      & heating_pdd(1:ko)            )  !! OUT

    end subroutine heating_GroundState


    ! Generation rate of zonal mean available potential energy
    ! bar{[Q/Exner(p)]^*}Exner(p_dagger) - bar{bar{[Q/Exner(o)]partial p / partial p_dagger_dagger}}Exner(p_dagger_dagger)
    ! This should not be used beacause bad precision
    subroutine heating_Zonal_simple(heating_gz_zm, heating_pdd, zonal_generation)
        real(kp), intent(in)  :: heating_gz_zm(jm,ko)
        real(kp), intent(in)  :: heating_pdd(ko)
        real(kp), intent(out) :: zonal_generation(ko)

        call globalMean(heating_gz_zm(1:jm,1:ko), &  !! IN
                      & zonal_generation(1:ko)    )  !! OUT

        zonal_generation(1:ko) = zonal_generation(1:ko) - heating_pdd(1:ko)

    end subroutine heating_Zonal_simple


    subroutine heating_Zonal_highPrecision(heating_exner_zm, zonal_generation_zm, zonal_generation)
        real(kp), intent(in)  :: heating_exner_zm(jm,ko)
        real(kp), intent(out) :: zonal_generation_zm(jm,ko)
        real(kp), intent(out) :: zonal_generation(ko)

        !real(kp) :: work_zonal_generation(jm,ko)

        integer :: k

        do k = 1, ko
            zonal_generation_zm(1:jm,k) = heating_exner_zm(1:jm,k) * &
                                        & cp * ((pout(k)*1.E-3_kp)**rkappa - (pdd_pd(1:jm,k)*1.E-3_kp)**rkappa)
        enddo

        call globalMean(zonal_generation_zm(1:jm,1:ko), &  !! IN
                      & zonal_generation(1:ko)          )  !! OUT

    end subroutine heating_Zonal_highPrecision


end module diabatic

