module status_output

    use NaNchecker, only : isNaN
    use params    , only : kp
    use namelist  , only : OUTPUT_WARN_FILENAME

    implicit none

    private
    public :: warn_open, warn_close, nan_write, warn_write

    ! If unit number is -1, the variable will be overwritten by the NEWUNIT argument
    ! -1 or 6 is strongly recommended to avoid bugs
    ! 6 is used for standard output
    integer, save :: warn_unit = -1

    contains


    subroutine warn_open()
        integer :: iostat
        character(128) :: iomsg

        if (warn_unit == 6) then
            return
        endif
        
        if (warn_unit == -1) then
            open(NEWUNIT=warn_unit           , &
               & FILE   =OUTPUT_WARN_FILENAME, &
               & ACTION ='WRITE'             , &
               & IOSTAT =iostat              , &
               & IOMSG  =iomsg                 )
        else
            open(UNIT  =warn_unit           , &
               & FILE  =OUTPUT_WARN_FILENAME, &
               & ACTION='WRITE'             , &
               & IOSTAT=iostat              , &
               & IOMSG =iomsg                 )
        endif

        if (iostat /= 0) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Failed to open ' // trim(OUTPUT_WARN_FILENAME)
            write(0,'(A)') trim(iomsg)
            ERROR STOP
        endif

    end subroutine warn_open


    subroutine warn_close()

        close(warn_unit)

    end subroutine warn_close


    subroutine nan_write(tt, nx, ny, nz, var, varname)
        integer     , intent(in) :: tt
        integer     , intent(in) :: nx
        integer     , intent(in) :: ny
        integer     , intent(in) :: nz
        real(kp)    , intent(in) :: var(nx,ny,nz)
        character(*), intent(in) :: varname

        integer :: x
        integer :: y
        integer :: z

        if (isNaN(var(1:nx,1:ny,1:nz))) then
            do z = 1, nz
                do y = 1, ny
                    do x = 1, nx
                        if (isNaN(var(x,y,z))) then
                            write(warn_unit,'(A)') 'WARNING'
                            write(warn_unit,'(A,I0)') 'T = ', tt
                            write(warn_unit,'(A,I0,A,I0,A,I0,A)') trim(varname) // '(', x, ',', y, ',', z, ') is NaN'
                        endif
                    enddo
                enddo
            enddo
        endif

    end subroutine nan_write


    subroutine warn_write(nx, ny, nz, var, min, max, varname, func)
        integer     , intent(in) :: nx
        integer     , intent(in) :: ny
        integer     , intent(in) :: nz
        real(kp)    , intent(in) :: var(nx,ny,nz)
        real(kp)    , intent(in) :: min
        real(kp)    , intent(in) :: max
        character(*), intent(in) :: varname
        character(*), intent(in) :: func

        integer :: x
        integer :: y
        integer :: z

        if (any(var(1:nx,1:ny,1:nz) < min .OR. var(1:nx,1:ny,1:nz) > max)) then
            do z = 1, nz
                do y = 1, ny
                    do x = 1, nx
                        if (var(x,y,z) < min .OR. var(x,y,z) > max) then
                            write(warn_unit,'(A)') 'WARNING in ' // trim(func)
                            write(warn_unit,'(A,ES15.3,A,ES15.3)') trim(varname) // ' should be between ', min, ' and ', max
                            write(warn_unit,'(A,I0,A,I0,A,I0,A,ES15.3)') trim(varname) // '(', x, ',', y, ',', z, ')=',var(x,y,z)
                        endif
                    enddo
                enddo
            enddo
        endif
    
    end subroutine warn_write


end module status_output

