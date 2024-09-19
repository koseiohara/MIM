module derivative

    use params, only : kp, pi
    
    implicit none

    private
    public :: derivative_x, derivative_y, derivative_p, derivative_z

    contains


    !
    ! d(var)/dx
    !
    subroutine derivative_x(nx, ny, nz, input, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: input(nx,ny,nz)
        real(kp), intent(out) :: output(nx,ny,nz)
      
        ! d(output) / d(lon) using central derivative.
        ! d(lon) = 2pi/nx -> 1/d(lon) = nx/2pi.
        ! another /2 in this program comes from the coefficient of central derivative.
        output(1     ,1:ny,1:nz) = (input(   2,1:ny,1:nz) - input(    nx,1:ny,1:nz))
        output(2:nx-1,1:ny,1:nz) = (input(3:nx,1:ny,1:nz) - input(1:nx-2,1:ny,1:nz))
        output(nx    ,1:ny,1:nz) = (input(   1,1:ny,1:nz) - input(  nx-1,1:ny,1:nz))

        output(1:nx,1:ny,1:nz) = output(1:nx,1:ny,1:nz) * real(nx, kind=kp) * (1._kp/(4._kp*pi))

    end subroutine derivative_x
    
    
    !
    ! d(var)/dy
    !
    subroutine derivative_y(nx, ny, nz, lat, input, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: lat(ny)
        real(kp), intent(in)  :: input(nx,ny,nz)
        real(kp), intent(out) :: output(nx,ny,nz)

        integer :: y
      
        ! Single precision derivative for both edges (=pole)
        ! Central derivative for the other grids
        output(1:nx,1,1:nz) = (input(1:nx,2,1:nz) - input(1:nx,1,1:nz)) / (lat(2) - lat(1))
        do y = 2, ny-1
            output(1:nx,y,1:nz) = (input(1:nx,y+1,1:nz) - input(1:nx,y-1,1:nz)) / (lat(y+1) - lat(y-1))
        enddo
        output(1:nx,ny,1:nz) = (input(1:nx,ny,1:nz) - input(1:nx,ny-1,1:nz)) / (lat(ny) - lat(ny-1))
      
    end subroutine derivative_y
    
    
    !
    ! d(var)/dp
    !
    ! Note: if the unit of plev is [hPa], then the result will be [*/hPa]
    !
    subroutine derivative_p(nx, ny, nz, plev, ps, input, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: plev(nz)
        real(kp), intent(in)  :: ps(nx,ny)
        real(kp), intent(in)  :: input(nx,ny,nz)
        real(kp), intent(out) :: output(nx,ny,nz)
    
        real(kp) :: x_hlf(nx,ny,nz+1)
        real(kp) :: p_hlf(nx,ny,nz+1)
        real(kp) :: d
        integer  :: x
        integer  :: y
        integer  :: z
    
        !! Check the order of x-y roop
        ! the original is do-x{do-y}, but chenged to do-y{do-x}
        ! Can the z-roop move to the most outside?
        do y = 1, ny
            do x = 1, nx
                do z = 2, nz
                   
                    if (plev(z) <= ps(x,y)) then
                       
                        x_hlf(x,y,z) = ( input(x,y,z) - input(x,y,z-1) ) / ( plev(z) - plev(z-1) )
                        p_hlf(x,y,z) = 0.5_kp * ( plev(z-1) + plev(z) )      
                       
                    else if (plev(z-1) < ps(x,y)) then
    
                        if (abs(ps(x,y) - plev(z-1)) / plev(z-1) < 0.01_kp) then
                            x_hlf(x,y,z) = ( input(x,y,z-1) - input(x,y,z-2) ) / ( plev(z-1) - plev(z-2) )
                        else
                            x_hlf(x,y,z) = ( input(x,y,z) - input(x,y,z-1) ) / ( ps(x,y) - plev(z-1) )
                        endif
    
                        p_hlf(x,y,z) = 0.5_kp * ( plev(z-1) + ps(x,y) )
    
                    else
                       
                        x_hlf(x,y,z) = ( input(x,y,z) - input(x,y,z-1) ) / ( plev(z) - plev(z-1) )
                        p_hlf(x,y,z) = 0.5_kp * ( plev(z-1) + plev(z) )     
                       
                    endif
                   
                enddo
                
                x_hlf(x,y,1) = x_hlf(x,y,2)
                x_hlf(x,y,nz+1) = x_hlf(x,y,nz)
                
                p_hlf(x,y,1) = max( 0.5_kp*plev(1), -0.5_kp*plev(2)+1.5_kp*plev(1) )
                p_hlf(x,y,nz+1) = -0.5_kp * plev(nz-1) + 1.5_kp * plev(nz)
                
                !intpl   
                do z = 2, nz-1
                    d = ( plev(z) - p_hlf(x,y,z) ) / ( p_hlf(x,y,z+1) - p_hlf(x,y,z) )
                    output(x,y,z) = x_hlf(x,y,z) * (1._kp-d) + x_hlf(x,y,z+1) * d
                enddo
                
                output(x,y,nz) = x_hlf(x,y,nz)
                output(x,y,1)    = x_hlf(x,y,1)
                
            enddo
        enddo
    
    end subroutine derivative_p
    
    
    !
    ! d(var)/dp
    !
    subroutine derivative_p_nops(nx, ny, nz, plev, input, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: plev(nz)
        real(kp), intent(in)  :: input(nz,ny,nz)
        real(kp), intent(out) :: output(nx,ny,nz)
    
        real(kp) :: ps(nx,ny)
    
        ps(1:nx,1:ny) = 1E+10_kp  ! dummy
    
        call derivative_p(nz, ny, nz, plev, ps, input, output)
    
    end subroutine derivative_p_nops
    
    
    !
    ! d(var)/dz
    !
    subroutine derivative_z(nx, ny, nz, zlev, input, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: zlev(nz)
        real(kp), intent(in)  :: input(nx,ny,nz)
        real(kp), intent(out) :: output(nx,ny,nz)
    
        real(kp) :: x_hlf(nx,ny,nz+1)
        real(kp) :: z_hlf(nx,ny,nz+1)
        real(kp) :: d
        integer  :: x
        integer  :: y
        integer  :: z
    
        !! Check the order of x-y roop
        ! the original is do-x{do-y}, but chenged to do-y{do-x}
        ! Can the z-roop move to the most outside?
        do x = 1, nx
            do y = 1, ny
    
                do z = 2, nz
                    x_hlf(x,y,z) = ( input(x,y,z) - input(x,y,z-1) ) / ( zlev(z) - zlev(z-1) )
                    z_hlf(x,y,z) = 0.5_kp * ( zlev(z-1) + zlev(z) )      
                enddo
                
                x_hlf(x,y,1)    = x_hlf(x,y,2)
                x_hlf(x,y,nz+1) = x_hlf(x,y,nz)
                
                z_hlf(x,y,1)      = max( 0.5_kp*zlev(1), -0.5_kp*zlev(2)+1.5_kp*zlev(1) )
                z_hlf(x,y,nz+1) = -0.5_kp * zlev(nz-1) + 1.5_kp * zlev(nz)
                
                ! interpolate
                do z = 2, nz-1
                    d = ( zlev(z) - z_hlf(x,y,z) ) / ( z_hlf(x,y,z+1) - z_hlf(x,y,z) )
                    output(x,y,z) = x_hlf(x,y,z) * (1._kp-d) + x_hlf(x,y,z+1) * d
                enddo
                
                output(x,y,nz) = x_hlf(x,y,nz)
                output(x,y,1)  = x_hlf(x,y,1)
               
            enddo
        enddo
    
    end subroutine derivative_z


end module derivative

