module energy_az

    use params  , only : kp, rkappa, grav, cp
    use com_var , only : jm, ko, pout, alat, costbl
    use mim_var , only : az_zm, az_zm_vint, az_gmean, pd_pdd, pd_ym, pt_ym, pt_pdds, p_pds, p_pdds
    use integral, only : integral_meridional, integral_pt_ym

    implicit none

    private
    public :: get_az


    real(kp), parameter :: coeff = cp / (1.E+5_kp**rkappa * (1._kp + rkappa) * grav)


    contains


    subroutine get_az()

        call get_az_zm_latprof()

        call get_az_vint()

        call get_az_gmean()

    end subroutine get_az


    !
    ! A_Z
    !
    ! Note:
    !   It is not used for global mean A_Z (because of the accuracy)
    !
    ! pd_ym : global mean p+ at the p++ levels
    !         pd_ym must be almost equal to standard p++ levels 
    !         except under the ground.
    ! 
    ! Contribution of each latitude to the az can be seen with this method.
    ! Taylor expansion is used (refer Lorenz(1995))
    ! 
    !!!!! This routine is highly recommended to compute az_zm !!!!!
    subroutine get_az_zm_latprof()
        real(kp) :: global_mean_pressure(ko)

        call integral_meridional(1                         , &  !! IN
                               & jm                        , &  !! IN
                               & ko                        , &  !! IN
                               & pd_pdd(1:jm,1:ko)         , &  !! IN
                               & global_mean_pressure(1:ko)  )  !! OUT

     
        ! get integrand
        call energy_az_decompose(global_mean_pressure(1:ko))  !! IN

    end subroutine get_az_zm_latprof


    !
    ! Vertically integrated A_Z
    !
    ! Note:
    !   It is not used for global mean A_Z (because of the accuracy)
    !
    ! Simple method to compute the zonal mean az, but the contribution of each latitude to the az cannot be seen in this method
    !
    subroutine get_az_zm_simple()
        real(kp) :: pd_ym(ko)
        integer :: k

        ! pd_ym : global mean p+ at the p++ levels
        !         pd_ym must be almost equal to standard p++ levels 
        !         except under the ground.
        call integral_meridional(1                , &  !! IN
                               & jm               , &  !! IN
                               & ko               , &  !! IN
                               & pd_pdd(1:jm,1:ko), &  !! IN
                               & pd_ym(1:ko)        )  !! OUT

        ! get integrand
        do k = 1, ko
            az_zm(1:jm,k) = coeff * ((pd_pdd(1:jm,k)*100._kp)**(1._kp+rkappa) - (pd_ym(k)*100._kp)**(1._kp+rkappa))
        enddo

    end subroutine get_az_zm_simple


    subroutine energy_az_decompose(global_mean_pressure)
        real(kp), intent(in) :: global_mean_pressure(ko)

        real(kp), parameter :: C2 = (rkappa+1._kp)*rkappa/2._kp
        real(kp), parameter :: C3 = (rkappa+1._kp)*rkappa*(rkappa-1._kp)/6._kp
        real(kp), parameter :: C4 = (rkappa+1._kp)*rkappa*(rkappa-1._kp)*(rkappa-2._kp)/24._kp
        real(kp), parameter :: C5 = (rkappa+1._kp)*rkappa*(rkappa-1._kp)*(rkappa-2._kp)*(rkappa-3._kp)/120._kp
        real(kp), parameter :: C6 = (rkappa+1._kp)*rkappa*(rkappa-1._kp)*(rkappa-2._kp)*(rkappa-3._kp)*(rkappa-4._kp)/720._kp

        real(kp) :: ratio(jm)

        integer :: k

        do k = 1, ko

            ratio(1:jm) = pd_pdd(1:jm,k) / global_mean_pressure(k)

            az_zm(1:jm,k) = C6 * ratio(1:jm)**6
            az_zm(1:jm,k) = az_zm(1:jm,k) + (                                          C5 -  6._kp*C6) * ratio(1:jm)**5
            az_zm(1:jm,k) = az_zm(1:jm,k) + (                              C4 -  5._kp*C5 + 15._kp*C6) * ratio(1:jm)**4
            az_zm(1:jm,k) = az_zm(1:jm,k) + (                   C3 - 4._kp*C4 + 10._kp*C5 - 20._kp*C6) * ratio(1:jm)**3
            az_zm(1:jm,k) = az_zm(1:jm,k) + (        C2 - 3._kp*C3 + 6._kp*C4 - 10._kp*C5 + 15._kp*C6) * ratio(1:jm)**2
            az_zm(1:jm,k) = az_zm(1:jm,k) + (- 2._kp*C2 + 3._kp*C3 - 4._kp*C4 +  5._kp*C5 -  6._kp*C6) * ratio(1:jm)
            az_zm(1:jm,k) = az_zm(1:jm,k) + (        C2 -       C3 +       C4 -        C5 +        C6)
            
            az_zm(1:jm,k) = az_zm(1:jm,k) * coeff * (global_mean_pressure(k)*100._kp)**(1._kp+rkappa)

        enddo

    end subroutine energy_az_decompose


    subroutine get_az_vint()
        real(kp) :: az_modify(jm)

        ! integrate with pt
        call integral_pt_ym(jm              , &  !! IN
                          & ko              , &  !! IN
                          & az_zm(1:jm,1:ko), &  !! IN
                          & az_zm_vint(1:jm)  )  !! OUT

        ! lower boundary modification
        !   it is proportional to pt_ymin * ( p+s^{kappa+1} - p++s^{kappa+1} ).
        az_modify(1:jm) = coeff * pt_pdds(1) * (p_pds(1:jm)**(rkappa+1._kp) - p_pdds(1)**(rkappa+1._kp))

        az_zm_vint(1:jm) = az_zm_vint(1:jm) + az_modify(1:jm)

    end subroutine get_az_vint    



    !
    ! Global mean A_Z
    !
    ! Note:
    !   order of the integration is different from energy_az_vint()
    !
    subroutine get_az_gmean()
        real(kp) :: integ(jm,ko)
        real(kp) :: integ_temp(ko)
        real(kp) :: temp1(jm)
        real(kp) :: az_modify(1)
        integer :: k
        
        ! get zonal mean az
        do k = 1, ko
            integ(1:jm,k) = coeff * ((pd_pdd(1:jm,k)*100._kp)**(1._kp+rkappa) - (pd_ym(k)*100._kp)**(1._kp+rkappa))
        enddo
        
        ! meridional mean
        call integral_meridional(1               , &  !! IN
                               & jm              , &  !! IN
                               & ko              , &  !! IN
                               & integ(1:jm,1:ko), &  !! IN
                               & integ_temp(1:ko)  )  !! OUT


        ! integrate with pt
        call integral_pt_ym(1               , &  !! IN
                          & ko              , &  !! IN
                          & integ_temp(1:ko), &  !! IN
                          & az_gmean(1)       )  !! OUT

        ! lower boundary modification
        !   it is proportional to pt_ymin * ( p+s^{kappa+1} - p++s^{kappa+1} ).
        temp1(1:jm) = coeff * pt_pdds(1) * (p_pds(1:jm)**(rkappa+1._kp) - p_pdds(1)**(rkappa+1._kp))

        call integral_meridional(1          , &  !! IN
                               & jm         , &  !! IN
                               & 1          , &  !! IN
                               & temp1(1:jm), &  !! IN
                               & az_modify(1) )  !! OUT
           
        az_gmean(1) = az_gmean(1) + az_modify(1)
      
    end subroutine get_az_gmean


end module energy_az

