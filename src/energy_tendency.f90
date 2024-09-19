module energy_tendency

    use params            , only : kp, radius, econv_min, econv_max
    use com_var           , only : jm, ko, zd, rho, alat, costbl
    use mim_var           , only : u_zm, v_zm, w_zm, pt_dot_zm, v_v_zm, u_v_zm, u_u_v_zm, v_v_v_zm                 , &
                                 & u_pt_dot_zm, v_pt_dot_zm, u_u_pt_dot_zm, v_v_pt_dot_zm, epy, epz_uw, epz, gy, gz, &
                                 & kz_zm, ke_zm, dkzdt_vkz, dkzdt_wkz, dkedt_uy, dkedt_vy, dkedt_uz, dkedt_vz      , &
                                 & dkedt_vke, dkedt_wke
    use status_output     , only : warn_write
    use derivative        , only : derivative_y, derivative_z
    use w_advection_expand, only : w_advection_components

    implicit none

    private
    public :: kz_advection, ke_advection, ke_flux_div

    contains


    subroutine kz_advection
        
        call kz_advection_v()

        call kz_advection_w()

    end subroutine kz_advection


    subroutine ke_advection()
        real(kp) :: v_ke_zm(jm,ko)

        v_ke_zm(1:jm,1:ko) = 0.5_kp * (u_u_v_zm(1:jm,1:ko) - u_v_zm(1:jm,1:ko) * (u_zm(1:jm,1:ko)+u_zm(1:jm,1:ko))  &
                                    & + u_zm(1:jm,1:ko) * u_zm(1:jm,1:ko) * v_zm(1:jm,1:ko)                         &
                                    & + v_v_v_zm(1:jm,1:ko) - v_v_zm(1:jm,1:ko) * (v_zm(1:jm,1:ko)+v_zm(1:jm,1:ko)) &
                                    & + v_zm(1:jm,1:ko) * v_zm(1:jm,1:ko) * v_zm(1:jm,1:ko)                         )

        call ke_advection_v(v_ke_zm(1:jm,1:ko))  !! IN

        call ke_advection_w(v_ke_zm(1:jm,1:ko))  !! IN

    end subroutine ke_advection


    subroutine ke_flux_div()

        call ke_flux_u_y()

        call ke_flux_v_y()

        call ke_flux_u_z()      ! Only u'w' component of the EP flux

        call ke_flux_v_z()

    end subroutine ke_flux_div


    ! renamed from energy_tendency_dkzdt_vkz
    ! advection by v
    ! -1 / [ a cos(phi) ] * d/d(phi) [ Kz v cos(phi) ]
    subroutine kz_advection_v()
        real(kp) :: work_v_kz(jm,ko)
    
        work_v_kz(1:jm,1:ko) = v_zm(1:jm,1:ko) * kz_zm(1:jm,1:ko) * spread(costbl(1:jm), 2, ko)
    
        call derivative_y(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & alat(1:jm)          , &  !! IN
                        & work_v_kz(1:jm,1:ko), &  !! IN
                        & dkzdt_vkz(1:jm,1:ko)  )  !! OUT

        dkzdt_vkz(1:jm,1:ko) = -dkzdt_vkz(1:jm,1:ko) * spread(1._kp/(costbl(1:jm)*radius), 2, ko)
        dkzdt_vkz( 1,1:ko) = 0._kp                                                          ! suppress sivergence
        dkzdt_vkz(jm,1:ko) = 0._kp                                                          ! suppress sivergence
    
        call warn_write(1                   , &  !! IN
                      & jm                  , &  !! IN
                      & ko                  , &  !! IN
                      & dkzdt_vkz(1:jm,1:ko), &  !! IN
                      & econv_min           , &  !! IN
                      & econv_max           , &  !! IN
                      & 'dkzdt_vkz'         , &  !! IN
                      & 'kz_advection_v()'    )  !! IN
    
    end subroutine kz_advection_v
    
   
    ! 
    ! renamed from energy_tendency_dkzdt_wkz
    ! advection by w+
    ! -1/rho * d/d(z+) ( Kz w )
    ! 
    subroutine kz_advection_w()
        real(kp) :: work_w_kz(jm,ko)
    
        work_w_kz(1:jm,1:ko) = w_zm(1:jm,1:ko) * kz_zm(1:jm,1:ko) * spread(rho(1:ko), 1, jm) 
    
        call derivative_z(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & zd(1:ko)            , &  !! IN
                        & work_w_kz(1:jm,1:ko), &  !! IN
                        & dkzdt_wkz(1:jm,1:ko)  )  !! OUT
    
        dkzdt_wkz(1:jm,1:ko) = -dkzdt_wkz(1:jm,1:ko) * spread(1._kp/rho(1:ko), 1, jm)
    
        call warn_write(1                   , &  !! IN
                      & jm                  , &  !! IN
                      & ko                  , &  !! IN
                      & dkzdt_wkz(1:jm,1:ko), &  !! IN
                      & econv_min           , &  !! IN
                      & econv_max           , &  !! IN
                      & 'dkzdt_wkz'         , &  !! IN
                      & 'kz_advection_w()'    )  !! IN
    
    end subroutine kz_advection_w
    

    !
    ! renamed from energy_tendency_dkedt_vke
    ! advection by v
    ! -1/ [ a cos(phi) ] * d/d(phi) [ Ke v cos(phi) ]
    !
    subroutine ke_advection_v(v_ke_zm)
        real(kp), intent(in) :: v_ke_zm(jm,ko)

        real(kp) :: work_v_ke(jm,ko)
    
        work_v_ke(1:jm,1:ko) = v_ke_zm(1:jm,1:ko) * spread(costbl(1:jm), 2, ko)
    
        call derivative_y(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & alat(1:jm)          , &  !! IN
                        & work_v_ke(1:jm,1:ko), &  !! IN
                        & dkedt_vke(1:jm,1:ko)  )  !! OUT

        dkedt_vke(1:jm,1:ko) = -dkedt_vke(1:jm,1:ko) * spread(1._kp/(costbl(1:jm)*radius), 2, ko)
        dkedt_vke(   1,1:ko) = 0._kp                                                  ! suppress sivergence
        dkedt_vke(  jm,1:ko) = 0._kp                                                  ! suppress sivergence
    
        call warn_write(1                   , &  !! IN
                      & jm                  , &  !! IN
                      & ko                  , &  !! IN
                      & dkedt_vke(1:jm,1:ko), &  !! IN
                      & econv_min           , &  !! IN
                      & econv_max           , &  !! IN
                      & 'dkedt_vke'         , &  !! IN
                      & 'ke_advection_v()'    )  !! IN
    
    end subroutine ke_advection_v
    
    
    !
    ! renamed from energy_tendency_dkedt_wke
    ! advection by w+
    ! -1/rho * d/d(z+) ( Ke w )
    !
    subroutine ke_advection_w(v_ke_zm)
        real(kp), intent(in) :: v_ke_zm(jm,ko)

        real(kp) :: pt_dot_ke_zm(jm,ko)
        real(kp) :: w_ke_zm(jm,ko)
        real(kp) :: dummy1(jm,ko)
        real(kp) :: dummy2(jm,ko)

        pt_dot_ke_zm(1:jm,1:ko) = 0.5_kp * (u_u_pt_dot_zm(1:jm,1:ko) - u_pt_dot_zm(1:jm,1:ko) * (u_zm(1:jm,1:ko)+u_zm(1:jm,1:ko))  &
                                         & + u_zm(1:jm,1:ko)*u_zm(1:jm,1:ko) * pt_dot_zm(1:jm,1:ko)                                &
                                         & + v_v_pt_dot_zm(1:jm,1:ko) - v_pt_dot_zm(1:jm,1:ko) * (v_zm(1:jm,1:ko)+v_zm(1:jm,1:ko)) &
                                         & + v_zm(1:jm,1:ko)*v_zm(1:jm,1:ko) * pt_dot_zm(1:jm,1:ko)                                )
    
        call w_advection_components(v_ke_zm(1:jm,1:ko)     , &  !! IN
                                  & pt_dot_ke_zm(1:jm,1:ko), &  !! IN
                                  & dummy1(1:jm,1:ko)      , &  !! OUT
                                  & dummy2(1:jm,1:ko)      , &  !! OUT
                                  & w_ke_zm(1:jm,1:ko)       )  !! OUT
                                  
        w_ke_zm(1:jm,1:ko) = w_ke_zm(1:jm,1:ko) * spread(rho(1:ko), 1, jm)
    
        call derivative_z(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & zd(1:ko)            , &  !! IN
                        & w_ke_zm(1:jm,1:ko)  , &  !! IN
                        & dkedt_wke(1:jm,1:ko)  )  !! OUT
    
        dkedt_wke(1:jm,1:ko) = -dkedt_wke(1:jm,1:ko) * spread(1._kp/rho(1:ko), 1, jm)
        dkedt_wke( 1,1:ko)   = 0._kp
        dkedt_wke(jm,1:ko)   = 0._kp
    
        call warn_write(1                   , &  !! IN
                      & jm                  , &  !! IN
                      & ko                  , &  !! IN
                      & dkedt_wke(1:jm,1:ko), &  !! IN
                      & econv_min           , &  !! IN
                      & econv_max           , &  !! IN
                      & 'dkedt_wke'         , &  !! IN
                      & 'ke/advection_w()'    )  !! IN
    
    end subroutine ke_advection_w
    
    
    !
    ! renamed from energy_tendency_dkedt_uy
    ! 1 / [ rho a^2 cos(phi) ] * d/d(phi) ( u * Fy cos(phi) )
    !
    subroutine ke_flux_u_y()
        real(kp) :: work_u_epy(jm,ko)
    
        work_u_epy(1:jm,1:ko) = u_zm(1:jm,1:ko) * epy(1:jm,1:ko)
    
        call derivative_y(1                    , &  !! IN
                        & jm                   , &  !! IN
                        & ko                   , &  !! IN
                        & alat(1:jm)           , &  !! IN
                        & work_u_epy(1:jm,1:ko), &  !! IN
                        & dkedt_uy(1:jm,1:ko)    )  !! OUT
    
        dkedt_uy(1:jm,1:ko) = dkedt_uy(1:jm,1:ko) / (spread(rho(1:ko)*radius*radius, 1, jm) * spread(costbl(1:jm), 2, ko))
        dkedt_uy(1,1:ko)    = 0._kp   ! suppress sivergence
        dkedt_uy(jm,1:ko)   = 0._kp  ! suppress sivergence
    
        call warn_write(1                  , &  !! IN
                      & jm                 , &  !! IN
                      & ko                 , &  !! IN
                      & dkedt_uy(1:jm,1:ko), &  !! IN
                      & econv_min          , &  !! IN
                      & econv_max          , &  !! IN
                      & 'dkedt_uy'         , &  !! IN
                      & 'ke_flux_u_y()'      )  !! IN
    
    end subroutine ke_flux_u_y
    
    
    !
    ! renamed from energy_tendency_dkedt_vy
    ! 1 / [ rho a^2 cos(phi) ] d/d(phi) ( v * Gy cos(phi) )
    !
    subroutine ke_flux_v_y()
        real(kp) :: work_v_gy(jm,ko)
    
        work_v_gy(1:jm,1:ko) = v_zm(1:jm,1:ko) * gy(1:jm,1:ko) * spread(costbl(1:jm), 2, ko)
    
        call derivative_y(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & alat(1:jm)          , &  !! IN
                        & work_v_gy(1:jm,1:ko), &  !! IN
                        & dkedt_vy(1:jm,1:ko)   )  !! OUT
    
        dkedt_vy(1:jm,1:ko) = dkedt_vy(1:jm,1:ko) / (spread(rho(1:ko)*radius*radius, 1, jm) * spread(costbl(1:jm), 2, ko))
        dkedt_vy( 1,1:ko) = 0._kp                                                 ! suppress sivergence
        dkedt_vy(jm,1:ko) = 0._kp                                                 ! suppress sivergence
    
        call warn_write(1                  , &  !! IN
                      & jm                 , &  !! IN
                      & ko                 , &  !! IN
                      & dkedt_vy(1:jm,1:ko), &  !! IN
                      & econv_min          , &  !! IN
                      & econv_max          , &  !! IN
                      & 'dkedt_vy'         , &  !! IN
                      & 'ke_flux_v_y()'      )  !! IN
    
    end subroutine ke_flux_v_y
    
    
    !
    ! renamed from energy_tendency_dkedt_uz
    ! 1 / [ rho a cos(phi) ] d/dz ( u * (-a rho cos(phi) bar{u'w'}) )
    !
    subroutine ke_flux_u_z()
        real(kp) :: work_u_epz_uw(jm,ko)
    
        work_u_epz_uw(1:jm,1:ko) = u_zm(1:jm,1:ko) * epz_uw(1:jm,1:ko)
    
        call derivative_z(1                       , &  !! IN
                        & jm                      , &  !! IN
                        & ko                      , &  !! IN
                        & zd(1:ko)                , &  !! IN
                        & work_u_epz_uw(1:jm,1:ko), &  !! IN
                        & dkedt_uz(1:jm,1:ko)       )  !! OUT
    
        dkedt_uz(1:jm,1:ko) = dkedt_uz(1:jm,1:ko) / (spread(rho(1:ko)*radius, 1, jm) * spread(costbl(1:jm), 2, ko))

        dkedt_uz( 1,1:ko) = 0._kp
        dkedt_uz(jm,1:ko) = 0._kp
    
        call warn_write(1                  , &  !! IN
                      & jm                 , &  !! IN
                      & ko                 , &  !! IN
                      & dkedt_uz(1:jm,1:ko), &  !! IN
                      & econv_min          , &  !! IN
                      & econv_max          , &  !! IN
                      & 'dkedt_uz'         , &  !! IN
                      & 'ke_flux_u_z()'      )  !! IN
    
    end subroutine ke_flux_u_z
    
    
    !
    ! renamed from energy_tendency_dkedt_vz
    ! 1 / [ rho a ] d/dz (v * Gz)
    !
    subroutine ke_flux_v_z()
        real(kp) :: work_v_gz(jm,ko)
    
        work_v_gz(1:jm,1:ko) = v_zm(1:jm,1:ko) * gz(1:jm,1:ko)
    
        call derivative_z(1                   , &  !! IN
                        & jm                  , &  !! IN
                        & ko                  , &  !! IN
                        & zd(1:ko)            , &  !! IN
                        & work_v_gz(1:jm,1:ko), &  !! IN
                        & dkedt_vz(1:jm,1:ko)   )  !! OUT
    
        dkedt_vz(1:jm,1:ko) = dkedt_vz(1:jm,1:ko) * spread(1._kp/(rho(1:ko)*radius), 1, jm)
    
        dkedt_vz( 1,1:ko) = 0._kp
        dkedt_vz(jm,1:ko) = 0._kp
    
        call warn_write(1                  , &  !! IN
                      & jm                 , &  !! IN
                      & ko                 , &  !! IN
                      & dkedt_vz(1:jm,1:ko), &  !! IN
                      & econv_min          , &  !! IN
                      & econv_max          , &  !! IN
                      & 'dkedt_vz'         , &  !! IN
                      & 'ke_flux_v_z()'      )  !! IN
    
    end subroutine ke_flux_v_z


end module energy_tendency

