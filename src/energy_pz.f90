module energy_pz

    use params , only : kp, cp, gasr
    use com_var, only : jm, ko
    use mim_var, only : pz_zm, t_dagger, phi_dagger

    implicit none

    private
    public :: get_pz

    contains


    subroutine get_pz()

        pz_zm(1:jm,1:ko) = (cp-gasr) * t_dagger(1:jm,1:ko) + phi_dagger(1:jm,1:ko)

    end subroutine get_pz


end module energy_pz

