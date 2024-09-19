module fft_utility

    use params, only : ckp, pi

    implicit none

    private
    public :: fft, dft, fft_quick

    contains


    !
    ! Wrapper routine fo fft_in
    !
    !  N     : Data length
    !  sign  : 1 or -1 (the sign of the exponent in exp)
    !  a     : input and output data
    !
    subroutine fft(N, sign, a)
        integer    , intent(in)    :: N
        integer    , intent(in)    :: sign      
        complex(ckp), intent(inout) :: a(0:N-1)

        complex(ckp) :: t(0:N-1)
        real(ckp)    :: theta

        theta = sign * 2._ckp * pi / real(N, kind=ckp)

        call fft_in(N       , &  !! IN
                  & theta   , &  !! IN
                  & a(0:N-1), &  !! INOUT
                  & t(0:N-1)  )  !! INOUT

    end subroutine fft


    !
    ! FFT
    !
    !  N     : Data length
    !  theta : +-2 * pi / N (where the sign specify transform/inverse)
    !  a     : input and output data
    !  t     : work space
    ! 
    ! A_k = sum_{j=0}^{N-1} a_j exp(i theta j k / N)
    ! 2,3,5... can be used for tha base
    !
    recursive subroutine fft_in(N, theta, a, t)
        integer     , intent(in)    :: N
        real(ckp)   , intent(in)    :: theta
        complex(ckp), intent(inout) :: a(0:N-1)
        complex(ckp), intent(inout) :: t(0:N-1)

        integer :: base
        integer :: Nb
        integer :: j
        integer :: k
        integer :: alpha
        integer :: p

        if (N == 1) then
            return
        endif

        if      (mod(N,  2) == 0) then
            base = 2
        else if (mod(N,  3) == 0) then
            base = 3
        else if (mod(N,  5) == 0) then
            base = 5
        else if (mod(N,  7) == 0) then
            base = 7
        else if (mod(N, 11) == 0) then
            base = 11
        else 
            base = N
        end if
           
        Nb = N / base

        ! Make some clusters for the next step
        do j = 0, Nb-1
            do alpha = 0, base-1

                t(j+Nb*alpha) = cmplx(0._ckp, 0._ckp)
                do p = 0, base-1
                    t(j+Nb*alpha) = t(j+Nb*alpha) + a(j+p*Nb) * exp(cmplx(0._ckp, theta*real((j + p*Nb)*alpha, kind=ckp)))
                enddo

            enddo
        enddo

        ! recurse
        do alpha = 0, base-1
            call fft_in(Nb                         , &  !! IN
                      & real(base, kind=ckp)*theta , &  !! IN
                      & t(Nb*alpha:Nb*(alpha+1)-1) , &  !! INOUT
                      & a(0:N-1)                     )  !! INOUT
        enddo

        ! save the result
        do k = 0, Nb-1
            do p = 0, base-1
                a(p+base*k) = t(k+Nb*p)
            enddo
        enddo

    end subroutine fft_in



    ! N : Data length
    ! thetain = +-2 * pi (where the sign specify transform/inverse)
    ! A_k = sum_{j=0}^{N-1} a_j exp(i theta j k / N)
    subroutine dft(N, theta, fr, fi, ar, ai)
        integer , intent(in)  :: N
        real(ckp), intent(in)  :: theta
        real(ckp), intent(in)  :: fr(0:N-1), fi(0:N-1)
        real(ckp), intent(out) :: ar(0:N-1), ai(0:N-1)

        real(ckp) :: temp
        integer  :: j
        integer  :: k

        ar(0:N-1) = 0._ckp
        ai(0:N-1) = 0._ckp
        
        do k = 0, N-1
            do j = 0, N-1
                temp  = theta * real(k * j, kind=ckp) / real(N, kind=ckp)
                ar(k) = ar(k) + fr(j) * cos(temp) - fi(j) * sin(temp)
                ai(k) = ai(k) + fi(j) * cos(temp) + fr(j) * sin(temp)
            enddo
        enddo

    end subroutine dft




    ! FFT¤Î¥µ¥Ö¥ë¡¼¥Á¥ó
    ! ¹âÂ®¡¢¹â¸úÎ¨¡£¤ä¤äÊ¬¤«¤ê¤Ë¤¯¤¤¡£
    ! thetain = +-2 * pi (Éä¹æ¤ÏÊÑ´¹¡¢µÕÊÑ´¹¤òÄêµÁ)
    ! A_k = ¦²_{j=0}^{N-1} a_j exp(i theta j k / N)
    subroutine fft_quick(N, thetain, ar, ai)
        integer , intent(in)    :: N
        real(ckp), intent(in)    :: thetain
        real(ckp), intent(inout) :: ar(0:N-1)
        real(ckp), intent(inout) :: ai(0:N-1)

        integer :: m
        integer :: mh
        integer :: i
        integer :: j
        integer :: k
        integer :: L
        integer :: NLOG
        real(ckp) :: wr
        real(ckp) :: wi
        real(ckp) :: xr
        real(ckp) :: xi
        real(ckp) :: tmp
        real(ckp) :: theta


        theta = thetain

        ! srambler
        i = 0
        do j = 1, N-1
            ! k = N / 2
            k = ISHFT(N, -1)

            do while( k <= i )
                i = i - k
                ! k = k / 2
                k = ISHFT(k, -1)
            enddo

            i = i + k
            
            if (j < i) then
                tmp   = ar(j)
                ar(j) = ar(i)
                ar(i) = tmp
                tmp   = ai(j)
                ai(j) = ai(i)
                ai(i) = tmp
            endif
        enddo

        NLOG = int( log(real(N+1, kind=ckp)) / log(2._ckp) ) - 1 
        ! write(*,*) "N=",N,"NLOG=",NLOG

        do L = 0, NLOG
            mh = 2**L
            m  = 2**(L+1)
            theta = theta * 0.5_ckp
            ! write(*,*) "m=",m,"mh=",mh,"theta=",theta

            do i = 0, mh-1
                wr = cos(theta*real(i, kind=ckp))
                wi = sin(theta*real(i, kind=ckp))
                do j = i, N-1, M
                    k = j + mh
                    xr = wr * ar(k) - wi * ai(k)
                    xi = wr * ai(k) + wi * ar(k)
                    ar(k) = ar(j) - xr
                    ai(k) = ai(j) - xi
                    ar(j) = ar(j) + xr
                    ai(j) = ai(j) + xi
                enddo
            enddo
        enddo

    end subroutine fft_quick


end module fft_utility

