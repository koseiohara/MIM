module isentropic_mean

    use params    , only : kp
    use com_var   , only : im, jm, km, ko, pout, alat
    use mim_var   , only : p_pd, p_sfc, p_pds, p_pdds, pd_pdd, nlev, dlev, nlev_y, dlev_y
    use derivative, only : derivative_p
    use integral  , only : integral_meridional

    implicit none

    private
    public :: zonalMean, p2pd_integral, globalMean

    contains


    !
    ! Function
    !   call p2pd_integral() and derivative_p()
    !
    ! Npte
    !   -Normally, it is not necessary to call p2pd_integral()
    !    and derivative_p(), independently.
    ! 
    subroutine zonalMean(x, x_zm)
        real(kp), intent(in)  :: x(im,jm,km)
        real(kp), intent(out) :: x_zm(jm,ko)

        real(kp) :: x_pt(im,jm,ko)
        real(kp) :: xint_zm(jm,ko)
  
        call p2pd_integral(x(1:im,1:jm,1:km)   , &  !! IN
                         & x_pt(1:im,1:jm,1:ko), &  !! OUT
                         & xint_zm(1:jm,1:ko)    )  !! OUT

        call derivative_p(1                 , &  !! IN  : nx
                        & jm                , &  !! IN  : ny
                        & ko                , &  !! IN  : nz
                        & pout(1:ko)        , &  !! IN  : plev
                        & p_pds(1:jm)       , &  !! IN  : ps
                        & xint_zm(1:jm,1:ko), &  !! IN  : var
                        & x_zm(1:jm,1:ko)     )  !! OUT

    end subroutine zonalMean
  
  
  
    !
    ! Function
    !   integrate vertically for coordinate transformation from p to p+
    !
    ! Arguements (in)
    !   x       : 3-dimensional pressure coordinate date
    !
    ! Arguements (out)
    !   x_pd    : 3-dimensional p+ coordinate data
    !   xint_zm : vertically integrated data for the next step (derivative_p)
    !
    ! Note
    !   -Normally, derivative_p() is called just after this subroutine.
    !
    subroutine p2pd_integral(x, x_pd, xint_zm)
        real(kp), intent(in)  :: x(im,jm,km)
        real(kp), intent(out) :: x_pd(im,jm,ko)
        real(kp), intent(out) :: xint_zm(jm,ko)
        
        real(kp) :: xint(im,jm,ko)
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        
        do j = 1, jm
            do i = 1, im
                
                ! 3D pressure coordinate -> 3D p+ coordinate
                do k = 1, ko
                    n = nlev(i,j,k)
                    !d = dlev(i,j,k)
                    x_pd(i,j,k) = x(i,j,n) * (1._kp-dlev(i,j,k)) + x(i,j,n+1) * dlev(i,j,k)
                enddo
    
                
                ! integrate with respect to p
                xint(i,j,1) = x_pd(i,j,1) * p_pd(i,j,1)
                
                do k = 2, ko  ! upper -> lower
    
                    if (p_pd(i,j,k) > p_sfc(i,j)) then  ! near (or under) the ground
                        xint(i,j,k) = 0.5_kp * (x_pd(i,j,k) + x_pd(i,j,k-1)) * (p_sfc(i,j) - p_pd(i,j,k-1)) + xint(i,j,k-1)
                    else
                        xint(i,j,k) = 0.5_kp * (x_pd(i,j,k) + x_pd(i,j,k-1)) * (p_pd(i,j,k) - p_pd(i,j,k-1)) + xint(i,j,k-1)
                    endif
     
                enddo
                
            enddo
        enddo
        
        ! zonal mean
        xint_zm = sum(xint, dim=1) / real(im, kind=kp)
      
    end subroutine p2pd_integral
  
  
  
    !
    ! Function
    !   differentiate with respect to p+ in order to transform coordinate
    !
    ! Arguement (in)
    !   xint_zm : result of p2pd_integral()
    !
    ! Arguement (out)
    !   x_zm    : 2-dimensional p+ coordinate data
    !
    ! Note
    !   -Normally, p2pd_integral() is run before this subroutine.
    !
    subroutine p2pd_derivative(xint_zm, x_zm)
        real(kp), intent(in)  :: xint_zm(jm,ko)
        real(kp), intent(out) :: x_zm(jm,ko)
        
        real(kp) :: x_hlf(jm, ko+1)
        real(kp) :: p_hlf(jm, ko+1)
        real(kp) :: d
        integer  :: k
        integer  :: j
        integer  :: ncalc
    
        ! switch for derivation
        ncalc = 1   ! nacalc=1 ---> standard
        !ncalc = 2   ! for comparision with mochi-version
        
        if (ncalc == 1) then
            call derivative_p(1, jm, ko, pout, p_pds, xint_zm, &
                 &             x_zm )
    
        else
           
           !========intpl ---->> bibun=============================
           !   if vertical grid is uniform ,calc.2 will be better than calc.1
           !-----------interpolate xint_zm(k) on half level-------------
           
           !set new half level
           do j=1, jm
              p_hlf(j,1) = max( 0.5*pout(1), -0.5*pout(2)+1.5*pout(1) )
              
              do k=2, ko+1
                 
                 p_hlf(j,k) = 2 * pout(k-1) - p_hlf(j,k-1)
                 
                 if( p_hlf(j,k) >= pout(k) .and. k /= ko+1 )then
                    
                    write(6,*) &
                         &  "warning: vertical grid is not appropriate. ", &
                         &  "p must be pout(k+1)>p_hlf(k)>pout(k))", &
                         &  k, p_hlf(j,k), pout(k)
                    
                 end if
                 
              end do
              
              d = ( p_hlf(j,1) - pout(1) ) / ( pout(2) - pout(1) )
              x_hlf(j,1) = xint_zm(j,1) * (1.0-d) + xint_zm(j,2) * d !k=1 extpl d>0
              
              do k=2, ko
                 
                 d = ( p_hlf(j,k) - pout(k-1) ) / ( pout(k) - pout(k-1) ) !intpl
                 x_hlf(j,k) = xint_zm(j,k-1) * (1.0-d) + xint_zm(j,k) * d
                 
              end do
              
              d = ( p_hlf(j,ko+1) - pout(ko-1) ) / ( pout(ko) - pout(ko-1) )
              !k=ko+1 extpl d>1
              x_hlf(j,ko+1) = xint_zm(j,ko-1) * (1.0-d) + xint_zm(j,ko) * d 
              
              !----------bibun -------------------------
              
              do k = 1, ko
                 
                 x_zm(j,k) = ( x_hlf(j,k+1) - x_hlf(j,k) ) / &
                      &      ( p_hlf(j,k+1) - p_hlf(j,k) )
                 
              enddo
           enddo
           
        endif
      
    end subroutine p2pd_derivative

  
    !
    ! Function
    !   call pd2pdd_integral() and derivative_p()
    !
    ! Npte
    !   -Normally, it is not necessary to call pd2pdd_integral()
    !    and derivative_p(), independently.
    ! 
    subroutine globalMean(x_zm, x_ym)
        real(kp), intent(in)  :: x_zm(jm,ko)
        real(kp), intent(out) :: x_ym(ko)
        real(kp) :: x_pt(jm,ko)
        real(kp) :: xint_ym(ko)

        
        call pd2pdd_integral(x_zm(1:jm,1:ko), &  !! IN
                           & x_pt(1:jm,1:ko), &  !! OUT
                           & xint_ym(1:ko)    )  !! OUT

        call derivative_p(1            , &  !! IN  : nx
                        & 1            , &  !! IN  : ny
                        & ko           , &  !! IN  : nz 
                        & pout(1:ko)   , &  !! IN
                        & p_pdds(1)    , &  !! IN
                        & xint_ym(1:ko), &  !! IN
                        & x_ym(1:ko)     )  !! OUT

    end subroutine globalMean
  
  
    !
    ! Function
    !   integrate vertically for coordinate transformation from p+ to p++
    !
    ! Arguements (in)
    !   x_zm    : 2-dimensional p+ coordinate date
    !
    ! Arguements (out)
    !   x_pdd   : 2-dimensional p++ coordinate data
    !   xint_ym : vertically integrated data for the next step (derivative_p)
    !
    ! Note
    !   -Normally, derivative_p() is called just after this subroutine.
    !
    subroutine pd2pdd_integral(x_zm, x_pdd, xint_ym)
        real(kp), intent(in)  :: x_zm(jm,ko)
        real(kp), intent(out) :: x_pdd(jm,ko)
        real(kp), intent(out) :: xint_ym(ko)
        
        real(kp) :: xint(jm,ko)
        integer  :: j
        integer  :: k
        integer  :: n
        
        do j = 1, jm
    
            ! 2D p+ coordinate -> 2D p++ coordinate
            do k = 1, ko
                n = nlev_y(j,k)
                x_pdd(j,k) = x_zm(j,n) * (1._kp-dlev_y(j,k)) + x_zm(j,n+1) * dlev_y(j,k)
            enddo
               
    
            ! integrate with respect to p+
            xint(j,1) = x_pdd(j,1) * pd_pdd(j,1)
               
            do k = 2, ko  ! upper -> lower
                   
                if (pd_pdd(j,k) > p_pds(j)) then
                    xint(j,k) = 0.5_kp * (x_pdd(j,k) + x_pdd(j,k-1)) * (p_pds(j) - pd_pdd(j,k-1)) + xint(j,k-1)
                else
                    xint(j,k) = 0.5_kp * (x_pdd(j,k) + x_pdd(j,k-1)) * (pd_pdd(j,k) - pd_pdd(j,k-1)) + xint(j,k-1)
                endif
            enddo
           
        enddo
    
        ! meridional mean
        call integral_meridional(1              , &  !! IN  : nx
                               & jm             , &  !! IN  : ny
                               & ko             , &  !! IN  : nz
                               & xint(1:jm,1:ko), &  !! IN
                               & xint_ym(1:ko)    )  !! OUT
          
    end subroutine pd2pdd_integral
  
  
end module isentropic_mean

