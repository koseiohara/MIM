module NaNchecker

    use, intrinsic :: IEEE_ARITHMETIC, only : IEEE_IS_NAN

    implicit none

    private
    public :: isNaN

    interface isNaN
        module procedure    &
          & isNaN_scalar_s, &
          & isNaN_scalar_d, &
          & isNaN_1dim_s  , &
          & isNaN_1dim_d  , &
          & isNaN_2dim_s  , &
          & isNaN_2dim_d  , &
          & isNaN_3dim_s  , &
          & isNaN_3dim_d 
    end interface isNaN


    contains


    function isNaN_scalar_s(input) result(output)
        real(4), intent(in) :: input
        logical :: output

        output = IEEE_IS_NAN(input)

    end function isNaN_scalar_s


    function isNaN_scalar_d(input) result(output)
        real(8), intent(in) :: input
        logical :: output

        output = IEEE_IS_NAN(input)

    end function isNaN_scalar_d


    function isNaN_1dim_s(input) result(output)
        real(4), intent(in) :: input(:)
        logical :: output

        output = any(IEEE_IS_NAN(input(:)))

    end function isNaN_1dim_s


    function isNaN_1dim_d(input) result(output)
        real(8), intent(in) :: input(:)
        logical :: output

        output = any(IEEE_IS_NAN(input(:)))

    end function isNaN_1dim_d


    function isNaN_2dim_s(input) result(output)
        real(4), intent(in) :: input(:,:)
        logical :: output

        output = any(IEEE_IS_NAN(input(:,:)))

    end function isNaN_2dim_s


    function isNaN_2dim_d(input) result(output)
        real(8), intent(in) :: input(:,:)
        logical :: output

        output = any(IEEE_IS_NAN(input(:,:)))

    end function isNaN_2dim_d


    function isNaN_3dim_s(input) result(output)
        real(4), intent(in) :: input(:,:,:)
        logical :: output

        output = any(IEEE_IS_NAN(input(:,:,:)))

    end function isNaN_3dim_s


    function isNaN_3dim_d(input) result(output)
        real(8), intent(in) :: input(:,:,:)
        logical :: output

        output = any(IEEE_IS_NAN(input(:,:,:)))

    end function isNaN_3dim_d


end module NaNchecker

