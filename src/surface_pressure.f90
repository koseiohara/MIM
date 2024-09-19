module surface_pressure

    use params  , only : kp
    use com_var , only : im, jm
    use mim_var , only : p_sfc, p_pds, p_pdds
    use integral, only : integral_meridional

    implicit none

    private
    public :: get_surface_pressure

    contains

    ! p_pds and p_pdds are computed from p_sfc
    subroutine get_surface_pressure()

        p_pds(1:jm) = sum(p_sfc(1:im,1:jm), dim=1) / real(im, kind=kp)

        call integral_meridional(1          , &  !! IN
                               & jm         , &  !! IN
                               & 1          , &  !! IN
                               & p_pds(1:jm), &  !! IN
                               & p_pdds(1)    )  !! OUT

    end subroutine get_surface_pressure


end module surface_pressure

