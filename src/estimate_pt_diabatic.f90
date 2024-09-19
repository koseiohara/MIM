module estimate_pt_diabatic
  
    use params  , only : kp, rkappa, pi, radius, cp
    use namelist, only : dt=>INPUT_TDEF_DT
    use com_var , only : im, jm, km, pin, alat, costbl
    use mim_var , only : u, v, omega, pt, pt_dot, u_past, v_past, omega_past, pt_past, q_3d

    implicit none

    private
    public :: estimate_diabatic_heating, get_pt_dot_q

    contains


    subroutine estimate_diabatic_heating()

        call get_pt_dot_omega_vector()
        !call get_pt_dot_omega_roop()

        call get_pt_dot_q_inv()

    end subroutine estimate_diabatic_heating


    !
    ! Function
    !   estimate D(pt)/Dt (i.e. diabatic heating) from total differential
    !
    ! Arguements (in)
    !   dt        : timestep of the input date
    !   u         : x-wind
    !   v         : y-wind
    !   omega     : omega velocity
    !   pt        : potential temperature at the pressure levels
    !   pt_before : pt 1-step before
    !
    ! Arguements (out)
    !   pt_dot    : D(pt)/Dt
    !
    ! Note
    !   -If omega=0 (no input), vertical advection of pt will be neglected.
    !     - This is aborted (2020/8/21)
    !     - omega is estimated in subroutine get_omega
    !   -If pt=pt_before (e.g. first time step), local pt change will be neglected.
    !   -Use data in one step before for advection. This reduces errors arose from
    !    miss matching in time of tenddency and advection terms (2020/8/21).
    !
    ! Using this routine is not recommended bevause this is not checked
    !
    subroutine get_pt_dot_omega_roop()
        integer, parameter :: ntm = 1      ! estimave advection at (t+t-1)/2
        real(kp) :: tint
        real(kp) :: pt_east
        real(kp) :: pt_west
        real(kp) :: pt_nth
        real(kp) :: pt_sth
        real(kp) :: pt_d
        real(kp) :: pt_u
        real(kp) :: u_int
        real(kp) :: v_int
        real(kp) :: omega_int
        real(kp) :: dphi
        real(kp) :: dp

        integer :: i
        integer :: j
        integer :: k
        integer :: il
        integer :: iu
        integer :: jl
        integer :: ju
        integer :: kd
        integer :: ku ! lower and upper indices
        integer :: nt

        pt_dot(1:im,1:jm,1:km) = 0._kp

        ! tendency of PT
        pt_dot(1:im,1:jm,1:km) = (pt(1:im,1:jm,1:km) - pt_past(1:im,1:jm,1:km)) / dt

        do nt=1, ntm
            ! tint (0 < tint < 1) time interpolation parameter
            tint = (nt - 0.5_kp) / real(ntm, kind=kp)

            ! check tint
            !if (tint <= 0._kp .OR. tint >= 1._kp) then
            !    write(0,'(A)') 'ERROR STOP from get_pt_dot_omega()'
            !    write(0,'(A,I0)') "tint is out of range : ", tint
            !    ERROR STOP
            !endif

            do k = 1, km

                kd = k + 1
                ku = k - 1
                if (kd == km+1) then
                    kd = k
                endif
               
                do j = 1, jm

                    jl = j - 1
                    ju = j + 1
                    if (jl == 0) then
                        jl = 1
                    endif
                    if (ju == jm+1) then
                        ju = jm
                    endif

                    do i = 1, im

                        ! u d(pt)/dx
                        il = i - 1
                        iu = i + 1
                        if (il == 0) then
                            il = im
                        endif
                        if (iu == im+1) then
                            iu = 1
                        endif

                        pt_east = pt_past(iu,j,k) * (1._kp-tint) + pt(iu,j,k) * tint
                        pt_west = pt_past(il,j,k) * (1._kp-tint) + pt(il,j,k) * tint

                        u_int   = u_past(i,j,k) * (1._kp-tint) + u(i,j,k) * tint

                        pt_dot(i,j,k) = pt_dot(i,j,k) &
                                    & + u_int * (pt_east - pt_west) &
                                    & * real(im, kind=kp) / (4._kp*pi*radius * costbl(j) * real(ntm, kind=kp))


                        ! v d(pt)/dy
                        pt_nth = pt_past(i,ju,k) * (1._kp-tint) + pt(i,ju,k) * tint
                        pt_sth = pt_past(i,jl,k) * (1._kp-tint) + pt(i,jl,k) * tint
                        dphi = (alat(ju) - alat(jl))

                        v_int  = v_past(i,j,k) * (1._kp-tint) + v(i,j,k) * tint

                        pt_dot(i,j,k) = pt_dot(i,j,k) + v_int * ( pt_nth - pt_sth ) / (dphi * radius * real(ntm, kind=kp))

                        ! omega d(pt)/dp
                        if (k == 1) then
                            continue
                        else

                            pt_d = pt_past(i,j,kd) * (1._kp-tint) + pt(i,j,kd) * tint
                            pt_u = pt_past(i,j,ku) * (1._kp-tint) + pt(i,j,ku) * tint
                            dp = (pin(kd) - pin(ku)) * 100._kp

                            omega_int = omega_past(i,j,k) * (1._kp-tint) + omega(i,j,k) * tint

                            pt_dot(i,j,k) = pt_dot(i,j,k) + omega_int * (pt_d - pt_u) / (dp * real(ntm, kind=kp))
                        endif


                    enddo
                enddo
            enddo

        enddo

    end subroutine get_pt_dot_omega_roop


    subroutine get_pt_dot_omega_vector()
        integer, parameter :: ntm = 1      ! estimave advection at (t+t-1)/2
        real(kp) :: pt_east(im,jm,km)
        real(kp) :: pt_west(im,jm,km)
        real(kp) :: pt_nth(im,jm,km)
        real(kp) :: pt_sth(im,jm,km)
        real(kp) :: pt_up(im,jm,km)
        real(kp) :: pt_lo(im,jm,km)
        real(kp) :: lat_dif(jm)
        real(kp) :: p_dif(km)
        real(kp) :: tint

        integer :: j
        integer :: k
        integer :: nt

        !pt_dot(1:im,1:jm,1:km) = 0._kp

        ! tendency of PT
        pt_dot(1:im,1:jm,1:km) = (pt(1:im,1:jm,1:km) - pt_past(1:im,1:jm,1:km)) / dt

        lat_dif(     1) = alat(   2) - alat(     1)
        lat_dif(2:jm-1) = alat(3:jm) - alat(1:jm-2)
        lat_dif(    jm) = alat(  jm) - alat(  jm-1)

        p_dif(2:km-1) = (pin(3:km) - pin(1:km-2)) * 100._kp
        p_dif(    km) = (pin(  km) - pin(  km-1)) * 100._kp

        do nt = 1, ntm
            ! tint (0 < tint < 1) time interpolation parameter
            tint = (real(nt, kind=kp) - 0.5_kp) / real(ntm, kind=kp)

            ! check tint
            !if (tint <= 0._kp .OR. tint >= 1._kp) then
            !    write(0,'(A)') 'ERROR STOP from get_pt_dot_omega()'
            !    write(0,'(A,I0)') "tint is out of range : ", tint
            !    ERROR STOP
            !endif

            pt_east(1:im-1,1:jm,1:km) = pt_past(2:im,1:jm,1:km) * (1._kp-tint) + pt(2:im,1:jm,1:km) * tint
            pt_east(    im,1:jm,1:km) = pt_past(   1,1:jm,1:km) * (1._kp-tint) + pt(   1,1:jm,1:km) * tint

            pt_west(   1,1:jm,1:km) = pt_past(    im,1:jm,1:km) * (1._kp-tint) + pt(    im,1:jm,1:km) * tint
            pt_west(2:im,1:jm,1:km) = pt_past(1:im-1,1:jm,1:km) * (1._kp-tint) + pt(1:im-1,1:jm,1:km) * tint

            do j = 1, jm
                pt_dot(1:im,j,1:km) = pt_dot(1:im,j,1:km)                                                    &
                                  & + (u_past(1:im,j,1:km) * (1._kp-tint) + u(1:im,j,1:km) * tint)           &
                                  & * (pt_east(1:im,j,1:km) - pt_west(1:im,j,1:km))                          &
                                  & * real(im, kind=kp) / (4._kp*pi*radius * costbl(j) * real(ntm, kind=kp))
            enddo


            pt_nth(1:im,1:jm-1,1:km) = pt_past(1:im,2:jm,1:km) * (1._kp-tint) + pt(1:im,2:jm,1:km) * tint
            pt_nth(1:im,    jm,1:km) = pt_past(1:im,  jm,1:km) * (1._kp-tint) + pt(1:im,  jm,1:km) * tint

            pt_sth(1:im,   1,1:km) = pt_past(1:im,     1,1:km) * (1._kp-tint) + pt(1:im,     1,1:km) * tint
            pt_sth(1:im,2:jm,1:km) = pt_past(1:im,1:jm-1,1:km) * (1._kp-tint) + pt(1:im,1:jm-1,1:km) * tint

            do j = 1, jm
                pt_dot(1:im,j,1:km) = pt_dot(1:im,j,1:km) &
                                  & + (v_past(1:im,j,1:km) * (1._kp-tint) + v(1:im,j,1:km) * tint) &
                                  & * (pt_nth(1:im,j,1:km) - pt_sth(1:im,j,1:km))                  &
                                  & / (lat_dif(j) * radius * real(ntm, kind=kp))
            enddo


            pt_up(1:im,1:jm,2:km) = pt_past(1:im,1:jm,1:km-1) * (1._kp-tint) + pt(1:im,1:jm,1:km-1) * tint

            pt_lo(1:im,1:jm,2:km-1) = pt_past(1:im,1:jm,3:km) * (1._kp-tint) + pt(1:im,1:jm,3:km) * tint
            pt_lo(1:im,1:jm,    km) = pt_past(1:im,1:jm,  km) * (1._kp-tint) + pt(1:im,1:jm,  km) * tint

            do k = 2, km
                pt_dot(1:im,1:jm,k) = pt_dot(1:im,1:jm,k)                                                         &
                                  & + (omega_past(1:im,1:jm,k) * (1._kp-tint) + omega(1:im,1:jm,k) * tint)        &
                                  & * (pt_lo(1:im,1:jm,k) - pt_up(1:im,1:jm,k)) / (p_dif(k) * real(ntm, kind=kp))
            enddo

        enddo

    end subroutine get_pt_dot_omega_vector


    !
    ! Function
    !   estimate D(pt)/Dt (i.e. diabatic heating) from diabatic heating input
    !
    ! Arguements (in)
    !   q_3d      : input diabatic heating data [J/(kg s)]
    !
    ! Arguements (out)
    !   pt_dot    : D(pt)/Dt
    !
    subroutine get_pt_dot_q_inv()
        integer :: k

        do k = 1, km
            q_3d(1:im,1:jm,k) = pt_dot(1:im,1:jm,k) * cp * (pin(k) * 1.E-3_kp)**rkappa
        enddo

    end subroutine get_pt_dot_q_inv


    subroutine get_pt_dot_q()
        integer :: k

        do k = 1, km
            pt_dot(1:im,1:jm,k) = q_3d(1:im,1:jm,k) / (cp * (pin(k) * 1.E-3_kp)**rkappa)
        enddo

    end subroutine get_pt_dot_q


end module estimate_pt_diabatic

