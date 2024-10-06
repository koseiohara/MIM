module energy_ae

    use params         , only : kp, rkappa, grav, cp
    use com_var        , only : im, jm, km, ko, pin, pout
    use mim_var        , only : p_sfc, p_zm, pt_zm, ae_total_zm, ae_zm_vint, p_pd, pt_pds, p_pds
    use derivative     , only : derivative_p
    use integral       , only : integral_pt
    use isentropic_mean, only : zonalMean

    implicit none

    private
    public :: get_ae


    real(kp), parameter :: coeff = cp / (1.E+5_kp**rkappa * (1._kp+rkappa)*grav)


    contains


    subroutine get_ae()

        call energy_ae_total()

        call energy_ae_vint()

    end subroutine get_ae


    !
    ! 2-dimensional A_E
    !   = P - Pz
    !   = const * pt d/d(p+) [ p^(kappa+1) - p+^(kappa+1) ]
    !
    ! Note:
    !   It is not used for estimate of vertically integrated A_E 
    !     since mountain effects have not been neglected here.
    !
    subroutine energy_ae_total()
        real(kp) :: p_diff(jm,ko)

        p_diff(1:jm,1:ko) = sum((p_pd(1:im,1:jm,1:ko)*100._kp)**(rkappa+1._kp), dim=1) / real(im, kind=kp) &
                        & - (p_zm(1:jm,1:ko)*100._kp)**(rkappa+1._kp)

        call derivative_p(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & pout(1:ko)*100._kp  , &  !! IN
                        & p_pds(1:jm)*100._kp , &  !! IN
                        & p_diff(1:jm,1:ko)   , &  !! IN
                        & ae_total_zm(1:jm,1:ko))  !! OUT

        ae_total_zm(1:jm,1:ko) = ae_total_zm(1:jm,1:ko) * pt_zm(1:jm,1:ko) * coeff

    end subroutine energy_ae_total


    !
    ! Vettically integrated A_E
    !
    ! ae_zm_vint : proportional to int[ (p^(rkappa+1))_zm - pd^(rkappa+1) ] d(pt)
    !
    subroutine energy_ae_vint()
        real(kp) :: p_kappa_p1_zm(jm,ko)
        real(kp) :: p_diff(jm,ko)
        real(kp) :: temp1(im,jm)
        real(kp) :: ae_modify(jm)


        ! Note: even if Taylor expansion is used, results are similar to the below
        p_kappa_p1_zm(1:jm,1:ko) = sum((p_pd(1:im,1:jm,1:ko)*100._kp)**(1._kp+rkappa), dim=1) / real(im, kind=kp)
        p_diff(1:jm,1:ko) = coeff * (p_kappa_p1_zm(1:jm,1:ko) - (p_zm(1:jm,1:ko)*100._kp)**(1._kp+rkappa))

        ! integrate with pt
        call integral_pt(jm               , &  !! IN
                       & ko               , &  !! IN
                       & p_diff(1:jm,1:ko), &  !! IN
                       & ae_zm_vint(1:jm)   )  !! OUT

        ! lower boundary modification
        !   it is proportional to pt_ymin * ( ps^{kappa+1} - p+s^{kappa+1} ).
        temp1(1:im,1:jm) = coeff * spread(pt_pds(1:jm), 1, im) * (p_sfc(1:im,1:jm)**(rkappa+1._kp) &
                       & - spread(p_pds(1:jm), 1, im)**(rkappa+1._kp))

        ae_modify(1:jm) = sum(temp1(1:im,1:jm), dim=1) / real(im, kind=kp)

        ae_zm_vint(1:jm) = ae_zm_vint(1:jm) + ae_modify(1:jm)

    end subroutine energy_ae_vint


end module energy_ae

