module intpl

    use params , only : kp
    use com_var, only : im, jm, km, ko
    use mim_var, only : pt, p_zm, pt_zm, pd_p, pd_ym, pt_ym, pdd_pd

    implicit none

    private
    public :: intpl_pd_p, intpl_pdd_pd

    contains


    !
    ! p_zm : p+ (zonal mean pressure) at the p+ levels 
    ! pd_p : 3-dimensional p+ at the pressure levels
    !
    subroutine intpl_pd_p()
        real(kp) :: w
        integer  :: i
        integer  :: j
        integer  :: k
        integer  :: h
        integer  :: kl
        integer  :: ku

        ku = 0
        kl = 0

        do k = 1, km            !!! Changed the roop order
            do j = 1, jm
                do i = 1, im

                    if (pt(i,j,k) > pt_zm(j,1)) then
                        kl = 1
                        ku = 2
                    else if (pt(i,j,k) < pt_zm(j,ko)) then
                        kl = ko - 1
                        ku = ko
                    else
                        do h = 1, ko-1
                            if ( pt_zm(j,h+1) <= pt(i,j,k) .AND. &
                               & pt(i,j,k)    <= pt_zm(j,h)      ) then
                                kl = h
                                ku = h + 1
                                exit
                            endif
                        enddo
                    endif

                    ! linear
                    w = ( pt(i,j,k) - pt_zm(j,kl) ) / ( pt_zm(j,ku) - pt_zm(j,kl) )
                    pd_p(i,j,k) = w * p_zm(j,ku) + (1._kp-w) * p_zm(j,kl)

                    ! log(p)
                    if (pd_p(i,j,k) <= 0._kp) then
                        w = (pt(i,j,k) - pt_zm(j,kl)) / (pt_zm(j,ku) - pt_zm(j,kl))
                        pd_p(i,j,k) = p_zm(j,kl) * (p_zm(j,ku) / p_zm(j,kl))**w
                    endif

                enddo
            enddo
        enddo

    end subroutine intpl_pd_p


    !
    ! pd_ym : p++ (global mean pressure) at the p++ levels 
    ! pdd_pd : 2-dimensional p++ at the p+ levels
    !
    subroutine intpl_pdd_pd()
        real(kp) :: w
        integer  :: j
        integer  :: k
        integer  :: h
        integer  :: kl
        integer  :: ku

        ku = 0
        kl = 0

        do k = 1, ko              !!! Changed the roop order
            do j = 1, jm
               
                if (pt_zm(j,k) > pt_ym(1)) then
                    kl = 1
                    ku = 2
                else if (pt_zm(j,k) < pt_ym(ko)) then
                    kl = ko - 1
                    ku = ko
                else
                    do h = 1, ko-1
                       if ( pt_ym(h+1) <= pt_zm(j,k) .AND. &
                          & pt_zm(j,k) <= pt_ym(h)         ) then
                           kl = h
                           ku = h + 1
                           exit
                       endif
                    enddo
                endif

                ! linear interpolation
                w = (pt_zm(j,k) - pt_ym(kl)) / (pt_ym(ku) - pt_ym(kl))
                pdd_pd(j,k) = w * pd_ym(ku) + (1._kp-w) * pd_ym(kl)

                if (pdd_pd(j,k) <= 0._kp) then
                    ! log(p) interpolation
                    w = (pt_zm(j,k) - pt_ym(kl)) / (pt_ym(ku) - pt_ym(kl))
                    pdd_pd(j,k) = pd_ym(kl) * (pd_ym(ku) / pd_ym(kl))**w
                endif

            enddo
        enddo

    end subroutine intpl_pdd_pd


end module intpl

