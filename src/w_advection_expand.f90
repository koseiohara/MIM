module w_advection_expand

    use params , only : kp, grav, radius
    use com_var, only : jm, ko, pout, rho, alat, costbl
    use mim_var, only : pt_zm

    implicit none

    private
    public :: w_advection_components

    contains

    !
    ! Input
    !   var_v_x_zm : zonal mean of (A'v')
    !   var_pt_dot_x_zm : zonal mean of (A'(dtheta/dt)')
    !                           where A is any, v is the eddy meridional wind, theta is potential temperature, t is time, 
    !                                 ( )' is eddy component.
    ! Output
    !   output_v : (1/(rho*a*g)) * var_v_x_zm * (dp_dagger/dphi)_theta
    !   output_t : (1/(rho*g)) * var_pt_dot_x_zm * (dp_dagger/dtheta)_phi
    !   output_w : output_v + output_t
    !                           where rho is the reference density, a is the earth radius, g is the gravity constant, 
    !                                 p_dagger is zonal mean pressure, phi is latitude, theta is potential temperature
    !
    ! ret_var_w : (var' w_dagger')_zm
    ! ret_var_v : var' v' component
    ! ret_var_t : var' pt_dot' component
    !
    subroutine w_advection_components(var_v_x_zm, var_pt_dot_x_zm, output_v, output_t, output_w)
        real(kp), intent(in)  :: var_v_x_zm(jm,ko)
        real(kp), intent(in)  :: var_pt_dot_x_zm(jm,ko)
        real(kp), intent(out) :: output_v(jm,ko)
        real(kp), intent(out) :: output_t(jm,ko)
        real(kp), intent(out) :: output_w(jm,ko)
        
        !***** output_v *****!
        call get_var_v_x_zm(var_v_x_zm(1:jm,1:ko), &            !! IN
                          & output_v(1:jm,1:ko)    )            !! OUT

        !***** output_t *****!
        call get_var_pt_dot_x_zm(var_pt_dot_x_zm(1:jm,1:ko), &  !! IN
                               & output_t(1:jm,1:ko)         )  !! OUT

        !***** output_w *****!
        output_w(1:jm,1:ko) = output_v(1:jm,1:ko) + output_t(1:jm,1:ko)

    end subroutine w_advection_components


    subroutine get_var_v_x_zm(input, output)
        real(kp), intent(in)  :: input(jm,ko)
        real(kp), intent(out) :: output(jm,ko)

        integer :: k
        integer :: kl
        integer :: ku

        k  = 1
        kl = 1
        ku = 2
        output(1,k) = 0._kp
        where (abs(pt_zm(2:jm-1,ku)-pt_zm(2:jm-1,kl)) <= 1.E-2_kp)
            output(2:jm-1,k) = 0._kp
        elsewhere
            ! 100 is multiplied to convert pout [hPa] to pout [Pa]
            output(2:jm-1,k) = input(2:jm-1,k) * (100._kp / (rho(k)*grav*radius)) &
                           & * (pt_zm(3:jm,k) - pt_zm(1:jm-2,k)) / (alat(3:jm) - alat(1:jm-2)) &
                           & * (pout(ku) - pout(kl)) / (pt_zm(2:jm-1,ku) - pt_zm(2:jm-1,kl))
        endwhere
        output(jm,k) = 0._kp

        do k = 2, ko-1
            output(1,k) = 0._kp
            
            kl = k - 1
            ku = k + 1

            where (abs(pt_zm(2:jm-1,ku)-pt_zm(2:jm-1,kl)) <= 1.E-2_kp)
                output(2:jm-1,k) = 0._kp
            elsewhere
                ! 100 is multiplied to convert pout [hPa] to pout [Pa]
                output(2:jm-1,k) = input(2:jm-1,k) * (100._kp / (rho(k)*grav*radius)) &
                               & * (pt_zm(3:jm,k) - pt_zm(1:jm-2,k)) / (alat(3:jm) - alat(1:jm-2)) &
                               & * (pout(ku) - pout(kl)) / (pt_zm(2:jm-1,ku) - pt_zm(2:jm-1,kl))
            endwhere

            output(jm,k) = 0._kp
        enddo

        k  = ko
        kl = ko-1
        ku = ko
        output(1,k) = 0._kp
        where (abs(pt_zm(2:jm-1,ku)-pt_zm(2:jm-1,kl)) <= 1.E-2_kp)
            output(2:jm-1,k) = 0._kp
        elsewhere
            ! 100 is multiplied to convert pout [hPa] to pout [Pa]
            output(2:jm-1,k) = input(2:jm-1,k) * (100._kp / (rho(k)*grav*radius)) &
                           & * (pt_zm(3:jm,k) - pt_zm(1:jm-2,k)) / (alat(3:jm) - alat(1:jm-2)) &
                           & * (pout(ku) - pout(kl)) / (pt_zm(2:jm-1,ku) - pt_zm(2:jm-1,kl))
        endwhere
        output(jm,k) = 0._kp

    end subroutine get_var_v_x_zm


    subroutine get_var_pt_dot_x_zm(input, output)
        real(kp), intent(in)  :: input(jm,ko)
        real(kp), intent(out) :: output(jm,ko)

        integer :: k
        integer :: kl
        integer :: ku

        output(1:jm,1) = 0._kp

        do k = 2, ko-1

            kl = k-1
            ku = k+1
            where (abs(pt_zm(1:jm,ku)-pt_zm(1:jm,kl)) < 1.E-2_kp)
                output(1:jm,k) = 0._kp
            elsewhere
                ! 100 is multiplied to convert pout [hPa] to pout [Pa]
                output(1:jm,k) = -input(1:jm,k) * (100._kp / (rho(k) * grav)) &
                             & * (pout(ku) - pout(kl)) / (pt_zm(1:jm,ku) - pt_zm(1:jm,kl))
            endwhere

        enddo

        k  = ko
        kl = ko-1
        ku = ko
        where (abs(pt_zm(1:jm,ku)-pt_zm(1:jm,kl)) < 1.E-2_kp)
            output(1:jm,k) = 0._kp
        elsewhere
            ! 100 is multiplied to convert pout [hPa] to pout [Pa]
            output(1:jm,k) = -input(1:jm,k) * (100._kp / (rho(k) * grav)) &
                         & * (pout(ku) - pout(kl)) / (pt_zm(1:jm,ku) - pt_zm(1:jm,kl))
        endwhere

    end subroutine get_var_pt_dot_x_zm


end module w_advection_expand

