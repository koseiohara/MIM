module get_w_zm

    use params       , only : kp, pi, grav, radius, h0, w_min, w_max
    use com_var      , only : jm, ko, sintbl, pout
    use mim_var      , only : p_pds, st_zm, w_zm
    use status_output, only : warn_write

    implicit none

    private
    public :: st2w

    contains


    !
    ! get vertical velocity from mass streamfunction
    ! 
    ! w (m/s) is computed
    !
    subroutine st2w()
        real(kp), parameter :: coef = h0*grav / (2._kp*pi*radius*radius)
        real(kp) :: work_st
        integer  :: k


        do k = 1, ko
            work_st = (st_zm(1,k) + st_zm(1,k)) - (st_zm(2,k) + st_zm(2,k))
            w_zm(1,k) = coef * work_st / ((sintbl(1) - sintbl(2))*2._kp * 100._kp*pout(k))

            w_zm(2:jm-1,k) = coef * (st_zm(1:jm-2,k) - st_zm(3:jm,k)) / ((sintbl(1:jm-2) - sintbl(3:jm)) * 100._kp*pout(k))

            work_st = (st_zm(jm-1,k) + st_zm(jm-1,k)) - (st_zm(jm,k) + st_zm(jm,k))
            w_zm(jm,k) = coef * work_st / ((sintbl(jm-1) - sintbl(jm))*2._kp * 100._kp*pout(k))

            where (p_pds(1:jm) <= pout(k))
                w_zm(1:jm,k) = 0._kp
            endwhere
        enddo

        call warn_write(1              , &  !! IN
                      & jm             , &  !! IN
                      & ko             , &  !! IN
                      & w_zm(1:jm,1:ko), &  !! IN
                      & w_min          , &  !! IN
                      & w_max          , &  !! IN
                      & 'w_zm'         , &  !! IN
                      & 'st2w()'         )  !! IN

    end subroutine st2w


end module get_w_zm

