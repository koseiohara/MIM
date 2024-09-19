module undef

    use params , only : rkp
    use com_var, only : pin

    implicit none

    private
    public :: undef_fill

    contains


    !
    ! Function
    !   interpolate/extrapolate undef data (pressure level only)
    !
    ! Arguements (in)
    !   im    : number of x-direction grid
    !   jm    : number of y-direction grid
    !   km    : number of z-direction grid
    !   undef : undef value
    !   pin   : pressure levels
    !
    ! Arguements (inout)
    !   var   : variable
    !
    subroutine undef_fill(nx, ny, nz, undef, var)
        integer , intent(in)    :: nx
        integer , intent(in)    :: ny
        integer , intent(in)    :: nz
        real(rkp), intent(in)    :: undef
        real(rkp), intent(inout) :: var(nx,ny,nz)

        real(rkp) :: new_value
        integer   :: i
        integer   :: j
        integer   :: k

        if (any(var(1:nx,1:ny,1:nz) == undef)) then

            do k = 3, nz
                do j = 1, ny
                    do i = 1, nx

                        if (var(i,j,k) == undef) then
                            call undef_interp(real(pin(k), kind=rkp)  , &  !! IN
                                            & real(pin(k-1), kind=rkp), &  !! IN
                                            & var(i,j,k-1)            , &  !! IN
                                            & real(pin(k-2), kind=rkp), &  !! IN
                                            & var(i,j,k-2)            , &  !! IN
                                            & new_value                 )  !! OUT
                            var(i,j,k) = new_value
                        endif

                    enddo
                enddo
            enddo

        endif

    end subroutine undef_fill


    !
    ! Function
    !   interpolate/extrapolate undef data
    !
    ! Arguements (in)
    !   p    : interpolated/extrapolated level
    !   p1   : level (1) to be used in interpolation/extrapolation
    !   A1   : value at p=p1
    !   p2   : level (2) to be used in interpolation/extrapolation
    !   A2   : value at p=p2
    !
    ! Arguements (inout)
    !   ret  : value at p
    !
    ! Note
    !   -log(p) linear interpolation/extrapolation is used.
    !
    subroutine undef_interp(p, p1, A1, p2, A2, output)
        real(rkp), intent(in)  :: p
        real(rkp), intent(in)  :: p1
        real(rkp), intent(in)  :: A1
        real(rkp), intent(in)  :: p2
        real(rkp), intent(in)  :: A2
        real(rkp), intent(out) :: output
  
        output = (A2*log(p/p1) + A1*log(p2/p)) / log(p2/p1)
        !output = ( (A2-A1)*log(p) + A1*log(p2) - A2*log(p1) ) / (log(p2)-log(p1))

        !  output = ( (A2-A1)*p + A1*p2 - A2*p1 ) / (p2-p1)

    end subroutine undef_interp


end module undef

