module potential_temperature_3d

    use params       , only : kp, rkappa, pt_max, pt_min
    use com_var      , only : im, jm, km, pin
    use mim_var      , only : t, p_sfc, pt, pt_sfc
    use status_output, only : warn_write

    implicit none

    private
    public :: get_pt_3d

    contains


    !
    ! Function
    !   get potential temperature at the pressure levels
    !
    ! Arguements (in)
    !   t      : temperature [K]
    !   p_sfcs : surface pressure [hPa]
    !
    ! Arguements (out)
    !   pt     : potential temperature [K]
    !   pt_sfc : surface potential temperature [K]
    !
    ! Note
    !   -If the atmosphere is unstable, potential temperature profile is adjusted
    !    (not so appropriate until now)
    !   -Potential temperature under the ground is set to its surface value
    !
    subroutine get_pt_3d()
        real(kp) :: d
        integer  :: i
        integer  :: j
        integer  :: k
        integer  :: l

        ! temperature -> potential temperature
        do k = 1, km
            pt(1:im,1:jm,k) = t(1:im,1:jm,k) * (1.E3_kp / pin(k))**rkappa
        enddo
        
        ! if atmosphere is instable, modify pt by extrapolation
        ! Note: energy conservation is NOT satisfied until now
        !       (it should be modified later...)
        do k = 2, km
            do j = 1, jm
                do i = 1, im
                    if (k > 1 .AND. pt(i,j,k) >= pt(i,j,k-1)) then
                        pt(i,j,k) = pt(i,j,k-1) - 1._kp
                    endif
                enddo
            enddo
        enddo
        
        ! get pt_sfc
        do i = 1, im
            do j = 1, jm

                ! search for the level which is nearest to the ground
                l = 1
                do while (pin(l+1) < p_sfc(i,j) .AND. l <= km-2)
                    l = l + 1
                enddo
                
                ! interpolate to ground level
                d = (p_sfc(i,j) - pin(l)) / (pin(l+1) - pin(l))
                pt_sfc(i,j) = pt(i,j,l) * (1._kp-d) + pt(i,j,l+1) * d
                
                ! if (pt_sfc(i,j) < 0._kp)then
                !     write(*,'(A)') 'WANING'
                !     write(*,'(A)') 'Surface potential temperature is less than 0'
                !     write(*,'("pt_sfc(",I0,",",I0,")="ES0.4)') i, j, pt_sfc(i,j)
                !     write(*,'("pt(",I0,",",I0,","I0,")="ES0.4)') i, j, l, pt(i,j,l)
                !     write(*,'("pt(",I0,",",I0,","I0,"+1)="ES0.4)') i, j, l, pt(i,j,l+1)
                ! endif

                ! pt = pt_sfc under the ground
                do k = 1, km
                    if (pin(k) > p_sfc(i,j)) then
                        pt(i,j,k) = pt_sfc(i,j)
                    endif
                enddo
               
            enddo
        enddo
        
        ! check value
        call warn_write(im                , &  !! IN
                      & jm                , &  !! IN
                      & km                , &  !! IN
                      & pt(1:im,1:jm,1:km), &  !! IN
                      & pt_min            , &  !! IN
                      & pt_max            , &  !! IN
                      & 'pt'              , &  !! IN
                      & 'getpt_3d()'        )  !! IN

        call warn_write(im                , &  !! IN
                      & jm                , &  !! IN
                      & 1                 , &  !! IN
                      & pt_sfc(1:im,1:jm) , &  !! IN
                      & pt_min            , &  !! IN
                      & pt_max            , &  !! IN
                      & 'pt_sfc'          , &  !! IN
                      & 'getpt_3d()'        )  !! IN

    end subroutine get_pt_3d


end module potential_temperature_3d

