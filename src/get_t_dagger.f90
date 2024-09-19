module get_t_dagger

    use params       , only : kp, rkappa, t_min, t_max
    use com_var      , only : jm, ko, pout
    use mim_var      , only : pt_zm, t_dagger
    use status_output, only : warn_write

    implicit none

    private
    public :: get_t_zonal_mean_state

    contains


    subroutine get_t_zonal_mean_state()
        
        t_dagger(1:jm,1:ko) = pt_zm(1:jm,1:ko) * spread((pout(1:ko)*1.E-3_kp)**rkappa, 1, jm)

        call warn_write(1                         , &  !! IN
                      & jm                        , &  !! IN
                      & ko                        , &  !! IN
                      & t_dagger(1:jm,1:ko)       , &  !! IN
                      & t_min                     , &  !! IN
                      & t_max                     , &  !! IN
                      & 't_dagger'                , &  !! IN
                      & 'get_t_zonal_mean_state()'  )  !! IN
        
    end subroutine get_t_zonal_mean_state

 
end module get_t_dagger

