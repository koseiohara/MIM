module integral

    use params , only : kp, grav, pi
    use com_var, only : pout, alat
    use mim_var, only : p_pds, p_pdds, pt_zm, pt_ym, pt_pds, pt_pdds

    implicit none

    private
    public :: integral_meridional, integral_p, integral_pt, integral_pt_ym

    contains

    
    !
    ! 1/2 int[ var cos(phi) ] d(phi)
    !
    subroutine integral_meridional(nx, ny, nz, input, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: input(nx,ny,nz)
        real(kp), intent(out) :: output(nx,nz)
    
        real(kp), parameter :: pi2 = 0.5_kp*pi
        integer :: j

        
        output(1:nx,1:nz) = input(1:nx,1,1:nz) * cos(alat(1)) * (pi2 - alat(1)) * 0.25_kp

        do j = 2, ny
            output(1:nx,1:nz) = output(1:nx,1:nz) &
                             & + (input(1:nx,j-1,1:nz)*cos(alat(j-1)) + input(1:nx,j,1:nz)*cos(alat(j))) &
                             & * (alat(j-1) - alat(j)) * 0.25_kp
        enddo

        output(1:nx,1:nz) = output(1:nx,1:nz) + input(1:nx,ny,1:nz) * cos(alat(ny)) * (alat(ny) + pi2) * 0.25_kp
         
    end subroutine integral_meridional
    
    
    !
    ! xint = 1/g int[ x_zm ] dp
    !
    !   p: p+ or p++
    !
    !  if p=p+  -> ps_zm = p_pds
    !  if p=p++ -> ps_zm = p_pdds
    !
    subroutine integral_p(ny, nz, ps_zm, x, xint)
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: ps_zm(ny)
        real(kp), intent(in)  :: x(ny,nz)
        real(kp), intent(out) :: xint(ny)
    
        real(kp) :: x_tmp(ny)
        integer  :: j
        integer  :: k
        
        ! trapeziodal integration
        do j = 1, ny
           
            xint(j) = x(j,1) * pout(1)
            
            do k = 2, nz
               
                if (pout(k) > ps_zm(j) .AND. pout(k-1) < ps_zm(j)) then

                    x_tmp(j) = (x(j,k) + x(j,k-1)) * (ps_zm(j) - pout(k-1)) / ((pout(k) - pout(k-1)) * 2._kp)
                    xint(j) = 0.5_kp * (x_tmp(j) + x(j,k-1)) * (ps_zm(j) - pout(k-1)) + xint(j)
                    
                else if (pout(k) <= ps_zm(j)) then
                   
                    xint(j) = 0.5_kp * (x(j,k) + x(j,k-1)) * (pout(k) - pout(k-1)) + xint(j)
                   
                else
                    xint(j) = xint(j)
                endif
               
            enddo
            
            xint(j) = xint(j) * 100._kp / grav
           
        enddo
        
    end subroutine integral_p
    
    
    !
    ! x_vint = int[ x_zm ] d(pt)
    !   [ pt_min : pt_top ]
    !
    subroutine integral_pt(ny, nz, x_zm, x_vint)
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: x_zm(ny,nz)
        real(kp), intent(out) :: x_vint(ny)

        integer :: j
        integer :: k
        
        ! trapeziodal integration
        do j = 1, ny
    
            x_vint(j) = 0._kp  ! above top -> 0
    
            do k = 2, nz
    
                if (pout(k) <= p_pds(j)) then
    
                    x_vint(j) = x_vint(j) + (pt_zm(j,k-1) - pt_zm(j,k)) * (x_zm(j,k-1) + x_zm(j,k)) * 0.5_kp
    
                    if (k == nz .AND. pt_zm(j,nz) >= pt_pds(j)) then
                        x_vint(j) = x_vint(j) + (pt_zm(j,nz) - pt_pds(j)) * x_zm(j,nz)
                    endif
    
                else if (pout(k)      >  p_pds(j)  .AND. &
                       & pout(k-1)    <= p_pds(j)  .AND. &
                       & pt_zm(j,k-1) >= pt_pds(j)       ) then
    
                    x_vint(j) = x_vint(j) + (pt_zm(j,k-1) - pt_pds(j)) * (x_zm(j,k-1) + x_zm(j,k)) * 0.5_kp
    
                endif
    
            enddo
        enddo
      
    end subroutine integral_pt
    
    
    !
    ! x_vint = int[ x_pdd ] d(pt)
    !   [ pt_ymin : pt_top ]
    !
    subroutine integral_pt_ym(ny, nz, x_pdd, x_vint)
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: x_pdd(ny,nz)
        real(kp), intent(out) :: x_vint(ny)
        integer :: j
        integer :: k
        
        ! trapeziodal integration
        do j = 1, ny
      
            x_vint(j) = 0._kp  ! above top -> 0
      
            do k = 2, nz
      
                if (pout(k) <= p_pdds(1)) then
      
                    x_vint(j) = x_vint(j) + (pt_ym(k-1) - pt_ym(k)) * (x_pdd(j,k-1) +  x_pdd(j,k)) * 0.5_kp
      
                    if (k == nz .AND. pt_ym(nz) >= pt_pdds(1)) then
                        x_vint(j) = x_vint(j) + (pt_ym(nz) - pt_pdds(1)) * x_pdd(j,nz)
                    endif
      
                else if (pout(k)    >  p_pdds(1)  .AND. &
                       & pout(k-1)  <= p_pdds(1)  .AND. &
                       & pt_ym(k-1) >= pt_pdds(1)       ) then
      
                    x_vint(j) = x_vint(j) + (pt_ym(k-1) - pt_pdds(1)) * (x_pdd(j,k-1) + x_pdd(j,k)) * 0.5_kp
      
                endif
      
            enddo
        enddo
      
    end subroutine integral_pt_ym


end module integral


