module energy_k

    use params , only : kp
    use com_var, only : jm, ko
    use mim_var, only : u_zm, v_zm, u_u_x_zm, v_v_x_zm, kz_zm, ke_zm

    implicit none

    private
    public :: get_kz, get_ke

    contains


    subroutine get_kz()
        
        kz_zm(1:jm,1:ko) = 0.5_kp * (u_zm(1:jm,1:ko)*u_zm(1:jm,1:ko) + v_zm(1:jm,1:ko)*v_zm(1:jm,1:ko))

    end subroutine get_kz


    subroutine get_ke()

        ke_zm(1:jm,1:ko) = 0.5_kp * (u_u_x_zm(1:jm,1:ko) + v_v_x_zm(1:jm,1:ko))

    end subroutine get_ke


end module energy_k

