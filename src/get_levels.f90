module get_levels

    use params       , only : kp, rkappa, p_min, p_max, pt_min, pt_max
    use com_var      , only : im, jm, km, ko, pin, pout
    use mim_var      , only : p_sfc, pt, p_zm, pt_zm, dlev, nlev, p_pd, pt_sfc, pt_pds, &
                            & nlev_y, dlev_y, pd_pdd, pd_ym, pt_ym, pt_pdds, p_pds, p_pdds
    use integral     , only : integral_meridional
    use status_output, only : warn_write

    implicit none

    private
    public :: getpt_global, getpt_y

    contains 


    subroutine getpt_global(tt)
        integer, intent(in) :: tt
        real(kp) :: potential_temp(im,km)
        real(kp) :: weight(im,ko)
        integer  :: grid_label(im,ko)
        real(kp) :: isentrop_p(im,ko)
        real(kp) :: zonal_mean_p(ko)
        real(kp) :: zonal_mean_pt(ko)

        integer :: j

        do j = 1, jm
            ! work space
            potential_temp(1:im,1:km) = pt(1:im,j,1:km)

            call getpt(tt                       , &  !! IN
                     & p_sfc(1:im,j)            , &  !! IN
                     & p_pds(j)                 , &  !! IN
                     & potential_temp(1:im,1:km), &  !! IN
                     & pt_sfc(1:im,j)           , &  !! IN
                     & pt_pds(j)                , &  !! IN
                     & weight(1:im,1:ko)        , &  !! OUT
                     & grid_label(1:im,1:ko)    , &  !! OUT
                     & isentrop_p(1:im,1:ko)    , &  !! OUT
                     & zonal_mean_p(1:ko)       , &  !! OUT
                     & zonal_mean_pt(1:ko)        )  !! OUT

            ! work space
            dlev(1:im,j,1:ko) = weight(1:im,1:ko)
            nlev(1:im,j,1:ko) = grid_label(1:im,1:ko)
            p_pd(1:im,j,1:ko) = isentrop_p(1:im,1:ko)
            p_zm(j,1:ko)      = zonal_mean_p(1:ko)
            pt_zm(j,1:ko)     = zonal_mean_pt(1:ko)
        enddo

    end subroutine getpt_global


    !
    ! Function
    !   prepare for p -> p+ coordinate transformation
    !
    ! Arguements (in)
    !   im     : number of input data grid point in x-direction
    !   km     : number of input data grid point in z-direction
    !   ko     : number of output data grid point in z-direction
    !   icount : current time step
    !   pin    : pressure levels
    !   pout   : p+ levels
    !   p_sfc  : surface pressure
    !   p_pds  : p+s
    !   pt     : potential temperature
    !   pt_sfc : potential temperature at the surface
    !   pt_pds : potential temperature at the surface in the p+ coordinate
    !
    ! Arguements (out)
    !   dlev   : interpolation parameter (weight)
    !   nlev   : interpolation parameter (grid number)
    !   p_pd   : pressure at the p+ levels
    !   p_zm   : zonal mean pressure at the p+ levels (it should approximate pout)
    !   pt_zm  : (zonal mean) potential temperature at the p+ levels
    !
    ! Note
    !   p_zm will be close to pout after running this subroutine.
    !
    ! Renamed Variables
    ! im -> (removed)
    ! km -> (removed)
    ! ko -> (removed)
    ! icount -> icount (not changed)
    ! pin -> (removed)
    ! pout -> (removed)
    ! p_sfc -> surface_p
    ! p_pds -> surface_p_pd
    ! pt -> potential_temp
    ! pt_sfc -> surface_pt
    ! pt_pds -> surface_pt_pd
    ! dlev -> weight
    ! nlev -> grid_label
    ! p_pd -> isentrop_p        ! pressure on isentropic surfaces
    ! p_zm -> pressure_zm
    ! pt_zm -> potential_temp_zm
    subroutine getpt(icount, surface_p, surface_p_pd, potential_temp, surface_pt, surface_pt_pd, &  !! IN
                   & weight, grid_label, isentrop_p, pressure_zm, potential_temp_zm            )  !! OUT
        integer , intent(in)  :: icount
        real(kp), intent(in)  :: surface_p(im)
        real(kp), intent(in)  :: surface_p_pd
        real(kp), intent(in)  :: potential_temp(im,km)
        real(kp), intent(in)  :: surface_pt(im)
        real(kp), intent(in)  :: surface_pt_pd
        real(kp), intent(out) :: weight(im,ko)
        integer , intent(out) :: grid_label(im,ko)
        real(kp), intent(out) :: isentrop_p(im,ko)
        real(kp), intent(out) :: pressure_zm(ko)
        real(kp), intent(out) :: potential_temp_zm(ko)

        real(kp) :: pt_pd(im,ko)  ! potential temperature at the p+ levels
        real(kp) :: potential_temp_zm_old(ko)
        real(kp) :: dr
        
        integer, parameter :: itmax = 10  ! maximum iteration number
        integer :: k
        integer :: it

        real(kp) :: dr_old
        real(kp) :: weight_bst(im,ko)
        real(kp) :: isentrop_p_bst(im,ko)
        real(kp) :: pressure_zm_bst(ko)
        real(kp) :: potential_temp_zm_bst(ko)
        real(kp) :: dr_all(ko)
        integer  :: grid_label_bst(im,ko)

        dr = 0._kp

        ! get 1st approximation of pt_zm
        call getpt_pt1(im                       , &  !! IN
                     & 1                        , &  !! IN
                     & km                       , &  !! IN
                     & ko                       , &  !! IN
                     & pin(1:km)                , &  !! IN
                     & pout(1:ko)               , &  !! IN
                     & potential_temp(1:im,1:km), &  !! IN
                     & surface_p(1:im)          , &  !! IN
                     & pt_pd(1:im,1:ko)         , &  !! OUT
                     & grid_label(1:im,1:ko)    , &  !! OUT
                     & weight(1:im,1:ko)          )  !! OUT

        potential_temp_zm(1:ko) = sum(pt_pd(1:im,1:ko), dim=1) / real(im, kind=kp)

        ! get p_pd ( pressure at the standard p+ levels ) and p_zm ( = p+ )
        call getpt_p(im                       , &  !! IN
                   & 1                        , &  !! IN
                   & km                       , &  !! IN
                   & ko                       , &  !! IN
                   & grid_label(1:im,1:ko)    , &  !! IN
                   & weight(1:im,1:ko)        , &  !! IN
                   & pin(1:km)                , &  !! IN
                   & potential_temp(1:im,1:km), &  !! IN
                   & potential_temp_zm(1:ko)  , &  !! IN
                   & isentrop_p(1:im,1:ko)      )  !! OUT

        pressure_zm(1:ko) = sum(isentrop_p(1:im,1:ko), dim=1) / real(im, kind=kp)         ! zonal mean
        where (pout > spread(surface_p_pd, 1, ko))
            potential_temp_zm = spread(surface_pt_pd, 1, ko)     ! if underground, pt_zm = pt_pds
        endwhere

        !***** iteration (itmax or less times) *****!
        do it = 1, itmax
             ! write(*,*) it

            ! unstable -> stable
            do k = 2, ko-1
                if (potential_temp_zm(k) < potential_temp_zm(ko)) then
                    potential_temp_zm(k) = (potential_temp_zm(k-1) + potential_temp_zm(k+1)) * 0.5_kp
                endif
            enddo

            ! previous value
            potential_temp_zm_old(1:ko) = potential_temp_zm(1:ko)

            ! interpolate pt_zm using pt_zm_old in order to make pressure_zm close to pout
            call getpt_ptiter(ko                         , &  !! IN
                            & pout(1:ko)                 , &  !! IN
                            & potential_temp_zm_old(1:ko), &  !! IN
                            & surface_p_pd               , &  !! IN
                            & surface_pt_pd              , &  !! IN
                            & pressure_zm(1:ko)          , &  !! IN
                            & potential_temp_zm(1:ko)      )  !! OUT

            ! get interpolation parameter nlev & dlev
            call getpt_lev(im                       , &  !! IN
                         & km                       , &  !! IN
                         & ko                       , &  !! IN
                         & pin(1:km)                , &  !! IN
                         & potential_temp(1:im,1:km), &  !! IN
                         & surface_pt(1:im)         , &  !! IN
                         & potential_temp_zm(1:ko)  , &  !! IN
                         & surface_p(1:im)          , &  !! IN
                         & grid_label(1:im,1:ko)    , &  !! OUT
                         & weight(1:im,1:ko)          )  !! OUT

            ! get isentrop_p (pressure on standard p+ levels) and pressure_zm ( = p+ )
            call getpt_p(im                       , &  !! IN
                       & 1                        , &  !! IN
                       & km                       , &  !! IN
                       & ko                       , &  !! IN
                       & grid_label(1:im,1:ko)    , &  !! IN
                       & weight(1:im,1:ko)        , &  !! IN
                       & pin(1:km)                , &  !! IN
                       & potential_temp(1:im,1:km), &  !! IN
                       & potential_temp_zm(1:ko)  , &  !! IN
                       & isentrop_p(1:im,1:ko)      )  !! OUT

            pressure_zm(1:ko) = sum(isentrop_p(1:im,1:ko), dim=1) / real(im, kind=kp)         ! zonal mean


            where (pout > spread(surface_p_pd, 1, ko))
                potential_temp_zm(1:ko) = spread(surface_pt_pd, 1, ko)     ! if underground, pt_zm = pt_pds
            endwhere

            ! check convergence condition and finish if appropriate
            dr_all(1:ko) = abs(pressure_zm(1:ko) / pout(1:ko) - 1._kp)

            ! If the results are betther than previous one, then store it
            if (it == 1 .OR. sum(dr_all(1:ko)) < dr_old) then
                weight_bst(1:im,1:ko)        = weight(1:im,1:ko)
                grid_label_bst(1:im,1:ko)    = grid_label(1:im,1:ko)
                pressure_zm_bst(1:ko)        = pressure_zm(1:ko)
                potential_temp_zm_bst(1:ko)  = potential_temp_zm(1:ko)
                isentrop_p_bst(1:im,1:ko)    = isentrop_p(1:im,1:ko)
                dr_old = sum(dr_all(1:ko))
            endif

            if (maxval(dr_all(1:ko)) <= 1.E-3_kp) then
                exit
            endif

        enddo

        ! If above iteration is not converged, then use best values
        if (it == itmax + 1) then
           weight(1:im,1:ko)        = weight_bst(1:im,1:ko)
           grid_label(1:im,1:ko)    = grid_label_bst(1:im,1:ko)
           pressure_zm(1:ko)        = pressure_zm_bst(1:ko)
           potential_temp_zm(1:ko)  = potential_temp_zm_bst(1:ko)
           isentrop_p(1:im,1:ko)    = isentrop_p_bst(1:im,1:ko)
        endif

        ! check pressure_zm and warn
        do k = 1, ko
            dr = abs(pressure_zm(k) / pout(k) - 1._kp)
            if (surface_p_pd > pout(k) .AND. dr > 1.E-2_kp) then
                write(*,'(A,I0,A,I3,4(A,ES10.3))') &            !!!!! BUG : Unit 65 does not exist (fixed)
                    & "icount = ", icount, ",   k = ", k, &
                    & ",   pout = ", pout(k), ",   p_zm = ", pressure_zm(k), ",   err = ", dr, ",   p_pds = ", surface_p_pd
            endif
        enddo

        ! check whether the atmosphere is stable or not
        call getpt_stable(1                      , &  !! IN
                        & ko                     , &  !! IN
                        & potential_temp_zm(1:ko), &  !! IN
                        & [surface_p_pd]         , &  !! IN
                        & [surface_pt]             )  !! IN

        call getp_stable(1                       , &  !! IN
                       & ko                      , &  !! IN
                       & pressure_zm(1:ko)       , &  !! IN
                       & [surface_p_pd]            )  !! IN

        call warn_write(im                     , &  !! IN
                      & 1                      , &  !! IN
                      & ko                     , &  !! IN
                      & isentrop_p(1:im,1:ko)  , &  !! IN
                      & p_min                  , &  !! IN
                      & p_max                  , &  !! IN
                      & 'p_pd'                 , &  !! IN
                      & 'getpt()'                )  !! IN

        call warn_write(1                      , &  !! IN
                      & 1                      , &  !! IN
                      & ko                     , &  !! IN
                      & potential_temp_zm(1:ko), &  !! IN
                      & pt_min                 , &  !! IN
                      & pt_max                 , &  !! IN
                      & 'pt_zm'                , &  !! IN
                      & 'getpt()'                )  !! IN

    end subroutine getpt


    !
    ! Function
    !   prepare for p+ -> p++ coordinate transformation
    !
    ! Arguements (in)
    !   jm     : number of input data grid point in y-direction
    !   km     : number of input (p+) data grid point in z-direction
    !   ko     : number of output (p++) data grid point in z-direction
    !   icount : current time step
    !   pin    : p+ levels (NOT equals to pin in com_var.f90)
    !   pout   : p++ levels
    !   alat   : latitude in radian
    !   p_pds  : p+s
    !   p_pdds : p++s
    !   pt_zm  : (zonal mean) potential temperature at the p+ levels
    !   pt_pds : potential temperature at the surface in the p+ coordinate
    !   pt_pdds: potential temperature at the surface in the p++ coordinate
    !
    ! Arguements (out)
    !   dlev_y : interpolation parameter (weight)
    !   nlev_y : interpolation parameter (grid number)
    !   pd_pdd : p+ at the p++ levels
    !   pd_ym  : global mean pressure at the p++ levels
    !            (it should approximate pout)
    !   pt_ym  : (global mean) potential temperature at the p++ levels
    !
    ! Note
    !   pd_ym will be close to pout after running this subroutine.
    !
    ! Renamed Variables
    ! jm -> (removed)
    ! km -> (removed)
    ! ko -> (removed)
    ! icount -> icount (not changed)
    ! pin -> (removed)
    ! pout -> (removed)
    ! alat -> (removed)
    ! p_pds -> (removed)
    ! p_pdds -> (removed)
    ! pt_zm -> (removed)
    ! pt_pds -> (removed)
    ! pt_pdds -> (removed)
    ! dlev_y -> (removed)
    ! nlev_y -> (removed)
    ! pd_pdd -> (removed)
    ! pd_ym -> (removed)
    ! pt_ym -> (removed)
    subroutine getpt_y(icount)

        integer, intent(in)  :: icount

        integer, parameter :: itmax = 10            ! maximum iteration number
        real(kp) :: pt_pdd(jm,ko)                  ! potential temperature at the p++ levels
        real(kp) :: pt_ym_old(ko)
        real(kp) :: dr
        integer  :: k
        integer  :: it
        !
        real(kp) :: dr_old
        integer  :: nlev_y_bst(jm,ko)
        real(kp) :: dlev_y_bst(jm,ko)
        real(kp) :: pd_pdd_bst(jm,ko)
        real(kp) :: pd_ym_bst(ko)
        real(kp) :: pt_ym_bst(ko)
        real(kp) :: dr_all(ko)

        ! get 1st approximation of pt_ym
        call getpt_pt1(1                , &  !! IN  : nx
                     & jm               , &  !! IN  : ny
                     & ko               , &  !! IN  : nzi
                     & ko               , &  !! IN  : nzo
                     & pout(1:ko)       , &  !! IN
                     & pout(1:ko)       , &  !! IN
                     & pt_zm(1:jm,1:ko) , &  !! IN
                     & p_pds(1:jm)      , &  !! IN
                     & pt_pdd(1:jm,1:ko), &  !! OUT
                     & nlev_y(1:jm,1:ko), &  !! OUT
                     & dlev_y(1:jm,1:ko)  )  !! OUT

        call integral_meridional(1                , &  !! IN
                               & jm               , &  !! IN
                               & ko               , &  !! IN
                               & pt_pdd(1:jm,1:ko), &  !! IN
                               & pt_ym(1:ko)        )  !! OUT

        ! get pd_pdd ( p+ at the standard p++ levels ) and pd_ym ( = p++ )
        call getpt_p(1                , &  !! IN
                   & jm               , &  !! IN
                   & ko               , &  !! IN
                   & ko               , &  !! IN
                   & nlev_y(1:jm,1:ko), &  !! IN
                   & dlev_y(1:jm,1:ko), &  !! IN
                   & pout(1:ko)       , &  !! IN
                   & pt_zm(1:jm,1:ko) , &  !! IN
                   & pt_ym(1:ko)      , &  !! IN            !!!!! BUG : The dimension of this argument must be (jm,ko)
                   & pd_pdd(1:jm,1:ko)  )  !! OUT

        call integral_meridional(1                , &  !! IN
                               & jm               , &  !! IN
                               & ko               , &  !! IN
                               & pd_pdd(1:jm,1:ko), &  !! IN
                               & pd_ym(1:ko)        )  !! OUT

        where (pout > spread(p_pdds(1), 1, ko))
            pt_ym(1:ko) = spread(pt_pdds(1), 1, ko)    ! if underground, pt_ym = pt_pdds
        endwhere


        !***** iteration (itmax times or less) *****!
        do it = 1, itmax

           ! unstable -> stable
           do k = 2, ko-1
               if (pt_ym(k) < pt_ym(ko)) then
                   pt_ym(k) = (pt_ym(k-1) + pt_ym(k+1)) * 0.5_kp
               endif
           enddo


           ! previous value
           pt_ym_old(1:ko) = pt_ym(1:ko)

           ! interpolate pt_ym using pt_ym_old in order to make pd_ym close to pout
           call getpt_ptiter(ko             , &  !! IN 
                           & pout(1:ko)     , &  !! IN 
                           & pt_ym_old(1:ko), &  !! IN 
                           & p_pdds(1)      , &  !! IN 
                           & pt_pdds(1)     , &  !! IN 
                           & pd_ym(1:ko)    , &  !! IN 
                           & pt_ym(1:ko)      )  !! OUT

           ! get interpolation parameter nlev_y & dlev_y
           call getpt_lev(jm               , &  !! IN
                        & ko               , &  !! IN
                        & ko               , &  !! IN
                        & pout(1:ko)       , &  !! IN
                        & pt_zm(1:jm,1:ko) , &  !! IN
                        & pt_pds(1:jm)     , &  !! IN
                        & pt_ym(1:ko)      , &  !! IN
                        & p_pds(1:jm)      , &  !! IN
                        & nlev_y(1:jm,1:ko), &  !! OUT
                        & dlev_y(1:jm,1:ko)  )  !! OUT

           ! get pd_pdd ( p+ on standard p++ levels ) and p_ym ( = p++ )
           call getpt_p(1                , &  !! IN
                      & jm               , &  !! IN
                      & ko               , &  !! IN
                      & ko               , &  !! IN
                      & nlev_y(1:jm,1:ko), &  !! IN
                      & dlev_y(1:jm,1:ko), &  !! IN
                      & pout(1:ko)       , &  !! IN
                      & pt_zm(1:jm,1:ko) , &  !! IN
                      & pt_ym(1:ko)      , &  !! IN            !!!!! BUG : The dimension of this argument must be (jm,ko)
                      & pd_pdd(1:jm,1:ko)  )  !! OUT

           call integral_meridional(1                , &  !! IN
                                  & jm               , &  !! IN
                                  & ko               , &  !! IN
                                  & pd_pdd(1:jm,1:ko), &  !! IN
                                  & pd_ym(1:ko)        )  !! OUT

           where (pout > spread(p_pdds(1), 1, ko))
               pt_ym(1:ko) = spread(pt_pdds(1), 1, ko)     ! if underground, pt_ym = pt_pdds
           endwhere

           ! check convergence condition and finish if appropriate
           dr_all(1:ko) = abs(pd_ym(1:ko)/pout(1:ko) - 1._kp)

           ! If the results are betther than previous one, then store it
           if (it == 1 .OR. sum(dr_all(1:ko)) < dr_old) then
               dlev_y_bst(1:jm,1:ko) = dlev_y(1:jm,1:ko)
               nlev_y_bst(1:jm,1:ko) = nlev_y(1:jm,1:ko)
               pd_ym_bst(1:ko)  = pd_ym(1:ko)
               pt_ym_bst(1:ko) = pt_ym(1:ko)
               pd_pdd_bst(1:jm,1:ko)   = pd_pdd(1:jm,1:ko)
               dr_old = sum(dr_all(1:ko))
           endif

           if (maxval(dr_all(1:ko)) <= 0.001_kp) then
               exit
           endif

        enddo

        ! If above iteration is not converged, then use best values
        if (it == itmax + 1) then
            dlev_y(1:jm,1:ko) = dlev_y_bst(1:jm,1:ko)
            nlev_y(1:jm,1:ko) = nlev_y_bst(1:jm,1:ko)
            pd_ym(1:ko)  = pd_ym_bst(1:ko)
            pt_ym(1:ko) = pt_ym_bst(1:ko)
            pd_pdd(1:jm,1:ko)   = pd_pdd_bst(1:jm,1:ko)
        endif

        ! check pd_ym and warn
        do k = 1, ko
            dr = abs(pd_ym(k) / pout(k) - 1._kp)
            if (p_pdds(1) > pout(k) .AND. dr > 1.E-2_kp) then
                write(*,'(A,I0,A,I3,4(A,ES10.3))') &            !!!!! BUG : Unit 65 does not exist (fixed)
                    & "icount = ", icount, ",   k = ", k, &
                    & ",   pout = ", pout(k), ",   pd_ym = ", pd_ym(k), ",   err = ", dr, ",   p_pdds = ", p_pdds
            endif
        enddo

        ! check whether the atmosphere is stable or not
        call getpt_stable(1          , &  !! IN
                        & ko         , &  !! IN
                        & pt_ym(1:ko), &  !! IN
                        & p_pdds(1:1), &  !! IN
                        & pt_pdds(1:1) )  !! IN

        call getp_stable(1           , &  !! IN
                       & ko          , &  !! IN
                       & pd_ym(1:ko) , &  !! IN
                       & p_pdds(1:1)   )  !! IN

        call warn_write(1                , &  !! IN
                      & jm               , &  !! IN
                      & ko               , &  !! IN
                      & pd_pdd(1:jm,1:ko), &  !! IN
                      & p_min            , &  !! IN
                      & p_max            , &  !! IN
                      & 'pd_pdd'         , &  !! IN
                      & 'getpt_y()'        )  !! IN

        call warn_write(1               , &  !! IN
                      & 1               , &  !! IN
                      & ko              , &  !! IN
                      & pt_ym(1:ko)     , &  !! IN
                      & pt_min          , &  !! IN
                      & pt_max          , &  !! IN
                      & 'pt_ym'         , &  !! IN
                      & 'getpt_y()'       )  !! IN

    end subroutine getpt_y


    !
    ! Function
    !   get 1st approximation of potential temperature
    !
    ! Arguements (in)
    !   im     : number of input data grid point in x-direction
    !   jm     : number of input data grid point in y-direction
    !   km     : number of input data grid point in z-direction
    !   ko     : number of output data grid point in z-direction
    !   pin    : input pressure levels
    !   pout   : output p+ levels
    !   pt     : potential temperature at the pressure levels [K]
    !   p_sfc  : surface pressure [hPa]
    !
    ! Arguements (out)
    !   pt_pd  : potential temperature approximately at the p+ levels
    !   nlev   : interpolation parameter
    !   dlev   : interpolation parameter
    !
    ! Note
    !   -Exactly speaking, pt_zm must be potential temperature at the p+ levels.
    !    Instead, potential temperature at the pressure levels is used
    !    as a first approximation. Not only pt_zm but also nlev and dlev
    !    will be modified in the other subroutines.
    !
    ! Renamed Variables
    !   im -> nx
    !   jm -> ny
    !   km -> nzi
    !   ko -> nzo
    !   pin -> plev_input
    !   pout -> plev_output
    !   pt_pd -> potential_temp_lev
    !   nlev -> grid_label
    !   dlev -> weight
    !
    subroutine getpt_pt1(nx, ny, nzi, nzo, plev_input, plev_output, potential_temp, surface_p, &  !! IN
                       & potential_temp_lev, grid_label, weight                                )  !! OUT
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nzi
        integer , intent(in)  :: nzo
        real(kp), intent(in)  :: plev_input(nzi)
        real(kp), intent(in)  :: plev_output(nzo)
        real(kp), intent(in)  :: potential_temp(nx,ny,nzi)
        real(kp), intent(in)  :: surface_p(nx,ny)
        real(kp), intent(out) :: potential_temp_lev(nx,ny,nzo)  ! pt on p+ levels
        integer , intent(out) :: grid_label(nx,ny,nzo)
        real(kp), intent(out) :: weight(nx,ny,nzo)

        integer  :: i
        integer  :: j
        integer  :: k
        integer  :: n
        real(kp) :: d

        grid_label(1:nx,1:ny,1:nzo) = 0
        weight(1:nx,1:ny,1:nzo) = 0._kp

        do k = 1, nzo
            do j = 1, ny
                do i = 1, nx

                   ! determine n which is one level lower than or equal to k-level
                   ! if pin == pout -> k=n
                   n = 1
                   if (k > 1) then
                       n = max(1, grid_label(i,j,k-1))
                   endif

                   do while (plev_input(n+1) < plev_output(k) .AND. n <= nzi-2)
                       n = n + 1
                   enddo

                   ! Prepare for interpolation
                   !   pin(n) < pout(k) < pin(n+1)
                   if (surface_p(i,j) > plev_input(n+1)) then  ! in the atmosphere

                       d = min(1._kp, ( plev_output(k)-plev_input(n) ) / ( plev_input(n+1)-plev_input(n) ) )
                       grid_label(i,j,k) = n
                       weight(i,j,k) = d
                       potential_temp_lev(i,j,k) = potential_temp(i,j,n) * (1._kp-d) + potential_temp(i,j,n+1) * d

                   else  ! underground

                       n = grid_label(i,j,k-1)
                       d = ( surface_p(i,j)-plev_input(n) ) / ( plev_input(n+1)-plev_input(n) )
                       grid_label(i,j,k) = n
                       weight(i,j,k) = d
                       potential_temp_lev(i,j,k) = potential_temp(i,j,n) * (1._kp-d) + potential_temp(i,j,n+1) * d

                   endif

                enddo
            enddo
        enddo

    end subroutine getpt_pt1


    !
    ! getpt_lev() - get interpolation parameter
    !
    !   pt      pressure  variable
    !
    !  pt(n)     pin(n)    A(n)         <- standard pressure levels
    !     |        |
    !  pt_zm     (ppp)      A           <- standard p+ levels (interpolation point)
    !     |        |
    !  pt(n+1)   pin(n+1)  A(n+1)       <- standard pressure levels
    !
    ! interpolation eq. : A = d A(n+1) + (1-d) A(n)
    !   nlev <-> n
    !   dlev <-> d (interpolation weight)
    !
    ! if p+ -> p++ interpolation, i is replaced by j
    ! (e.g. im -> jm, pt_zm -> pt_ym)
    !
    ! Renamed Variables
    !   im -> nx
    !   km -> nzi
    !   ko -> nzo
    !   pin -> plev_input
    !   pt -> potential_temp
    !   pt_sfc -> potential_temp_sfc
    !   pt_zm -> potential_temp_mean
    !   p_sfc -> surface_p
    !   nlev -> grid_label
    !   dlev -> weight
    subroutine getpt_lev(nx, nzi, nzo, plev_input, potential_temp, potential_temp_sfc, potential_temp_mean, surface_p, &  !! IN
                       & grid_label, weight                                                                            )  !! OUT
        integer , intent(in)  :: nx
        integer , intent(in)  :: nzi
        integer , intent(in)  :: nzo
        real(kp), intent(in)  :: plev_input(nzi)
        real(kp), intent(in)  :: potential_temp(nx,nzi)
        real(kp), intent(in)  :: potential_temp_sfc(nx)
        real(kp), intent(in)  :: potential_temp_mean(nzo)
        real(kp), intent(in)  :: surface_p(nx)
        integer , intent(out) :: grid_label(nx,nzo)
        real(kp), intent(out) :: weight(nx,nzo)

        real(kp) :: ppp
        real(kp) :: dppp
        integer  :: i
        integer  :: k
        integer  :: n

        do k = 1, nzo
            do i = 1, nx

                !********** get pressure ppp corresponding to pt_zm **********!

                if (potential_temp_mean(k) > potential_temp_sfc(i)) then  ! if pt_zm(k) is in the atmosphere

                    !***** find n which satisfies pt(n) > pt_zm(k) > pt(n+1)
                    n = 1
                    do while( potential_temp(i,n+1) > potential_temp_mean(k) .AND. n <= nzo-2 )
                        n = n + 1
                    enddo

                    !***** interpolate to ppp *****!

                    ! if pt_zm(k) is below the lowermost grid point,
                    ! linear interpolation using pt(i,km) and pt_sfc(i)
                    if (potential_temp_mean(k) <= potential_temp(i,nzi)) then
                        ppp = plev_input(nzi) &
                          & + (potential_temp_mean(k) - potential_temp(i,nzi)) * (surface_p(i) - plev_input(nzi)) &
                          & / (potential_temp_sfc(i) - potential_temp(i,nzi))

                    ! if pin(n+1) is under the ground,
                    ! linear interpolation using pt(i,n) and pt_sfc(i)
                    else if (plev_input(n+1) > surface_p(i)) then
                        ppp = plev_input(n) + (potential_temp_mean(k) - potential_temp(i,n)) * (surface_p(i) - plev_input(n)) &
                                      & / (potential_temp_sfc(i) - potential_temp(i,n))

                    !  if pt_zm(k) is above the uppermost grid point,
                    !  use definition of potential temperature and assume T=const.
                    else if (potential_temp_mean(k) > potential_temp(i,1)) then
                        ppp = ( potential_temp(i,1) / potential_temp_mean(k) )**(1._kp/rkappa) * plev_input(1)

                    else ! log(p) interpolation
                        dppp = log(plev_input(n)) + (potential_temp_mean(k) &
                           & - potential_temp(i,n)) * log(plev_input(n+1)/plev_input(n)) &
                           & / (potential_temp(i,n+1) - potential_temp(i,n))
                        ppp = exp(dppp)

                    endif

                else ! if pt_zm(k) is under the ground

                    !***** find n which satisfies pin(n) < p_sfc(k) < pin(n+1)
                    n = 1
                    do while (plev_input(n+1) < surface_p(i) .AND. n <= nzo-2)
                        n = n + 1
                    enddo

                    ! use surface value (from the definition)
                    ppp = surface_p(i)

                endif


                grid_label(i,k) = n

                !********** calculate interpolation weight **********!
                ! log(p) interpolation
                weight(i,k) = log(ppp/plev_input(n)) / log(plev_input(n+1)/plev_input(n))

                ! if pt_zm(k) is above the uppermost grid point,
                ! use linear interpolation to avoid too small dlev.
                if (weight(i,k) < 0._kp .AND. ppp < plev_input(1)) then
                    weight(i,k) = (ppp - plev_input(n)) / (plev_input(n+1) - plev_input(n))
                endif

            enddo
        enddo

    end subroutine getpt_lev


    !
    ! Function
    !   get pressure at the p+ levels
    !     or get p+ at the p++ levels
    !
    ! if im != 1,
    !   pin   : pressure levels
    !   pt    : potential temperature at the pressure levels
    !   pt_zm : potential temperature at the p+ levels
    !   p_pd  : pressure at the p+ levels
    !
    ! if im == 1, be careful...
    !   pin   : p+ levels
    !   pt    : potential temperature at the p+ levels
    !   pt_zm : potential temperature at the p++ levels
    !   p_pd  : p+ at the p++ levels
    !
    ! Renamed Variables
    !     im -> nx
    !     jm -> ny
    !     km -> nzi
    !     ko -> nzo
    !     nlev -> grid_lebel
    !     dlev -> weight
    !     pin -> plev_input
    !     pt -> potential_temp
    !     pt_zm -> potential_temp_mean
    !     p_pd -> output
    subroutine getpt_p(nx, ny, nzi, nzo, grid_label, weight, plev_input, potential_temp, potential_temp_mean, output)
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nzi
        integer , intent(in)  :: nzo
        integer , intent(in)  :: grid_label(nx,ny,nzo)
        real(kp), intent(in)  :: weight(nx,ny,nzo)
        real(kp), intent(in)  :: plev_input(nzi)
        real(kp), intent(in)  :: potential_temp(nx,ny,nzi)
        real(kp), intent(in)  :: potential_temp_mean(ny,nzo)
        real(kp), intent(out) :: output(nx,ny,nzo)

        real(kp) :: d
        integer  :: i
        integer  :: j
        integer  :: k
        integer  :: l


        !! Use OPENMP
        do k = 1, nzo
            do j = 1, ny
                do i = 1, nx

                    ! linear interpolation/extrapolation
                    !!!!!!! p_pd(i,j,k) = (1.0-d) * pin(l) + d * pin(l+1)        !!! Not Used
                    !! nlev and dlev is evaluated with log(p)-interpolation,
                    !! so below is better than above.
                    output(i,j,k) = exp((1._kp-weight(i,j,k)) * log(plev_input(grid_label(i,j,k))) &
                                & + weight(i,j,k) * log(plev_input(grid_label(i,j,k)+1)))

                    ! if p<=0 -> extrapolate again with log(p) interpolation
                    if (output(i,j,k) <= 0._kp) then
                        l = grid_label(i,j,k)
                        d = (potential_temp_mean(j,k) - potential_temp(i,j,l)) / (potential_temp(i,j,l+1) - potential_temp(i,j,l))
                        output(i,j,k) = plev_input(l) * ( plev_input(l+1) / plev_input(l) )**d
                    endif

                enddo
            enddo
        enddo

    end subroutine getpt_p


    !
    ! Function
    !   interpolate pt_zm using pt_zm_old in order to make p_zm close to pout
    !
    ! pt_zm_old -> pt_zm
    !
    ! ko -> nz
    ! pout -> plev_output
    ! pt_zm_old -> pt_zm_old
    ! p_zm -> pressure_zm
    ! p_pds -> surface_p
    ! pt_pds -> surface_pt
    ! pt_zm -> pt_zm_new
    subroutine getpt_ptiter(nz, plev_output, pt_zm_old, surface_p, surface_pt, pressure_zm, pt_zm_new)
        integer , intent(in)  :: nz
        real(kp), intent(in)  :: plev_output(nz)
        real(kp), intent(in)  :: pt_zm_old(nz)
        real(kp), intent(in)  :: surface_p
        real(kp), intent(in)  :: surface_pt
        real(kp), intent(in)  :: pressure_zm(nz)
        real(kp), intent(out) :: pt_zm_new(nz)

        !-- pressure/potential temperature arrays including surface values
        real(kp) :: p_zm_itr(nz+1)
        real(kp) :: pt_zm_itr(nz+1)

        real(kp) :: pl
        real(kp) :: pu
        real(kp) :: ptl
        real(kp) :: ptu
        real(kp) :: a

        integer :: k
        integer :: kk
        
        !
        p_zm_itr(1:nz) = pressure_zm(1:nz)
        p_zm_itr(nz+1) = surface_p
        pt_zm_itr(1:nz) = pt_zm_old(1:nz)
        pt_zm_itr(nz+1) = surface_pt
        do k = 1, nz+1
            if (p_zm_itr(k) > surface_p) then
                p_zm_itr(k) = surface_p
            endif
            if (pt_zm_itr(k) < surface_pt) then
                pt_zm_itr(k) = surface_pt
            endif
        enddo


        pt_zm_new(1:nz) = -999._kp             ! substitute unrealistic values for checking


        kk = 1
        kloop: do k = 2, nz+1                  ! loop for p_zm_itr

            do
                !-- find layers that sandwitch pout(kk)
                a = ( p_zm_itr(k-1) - plev_output(kk) ) * ( p_zm_itr(k) - plev_output(kk) )
                if (plev_output(kk) < p_zm_itr(1) .OR. a <= 0._kp) then

                    if (plev_output(kk) < p_zm_itr(1)) then ! use k=1 and 2
                        pu = p_zm_itr(1)
                        pl = p_zm_itr(2)
                        ptu = pt_zm_itr(1)
                        ptl = pt_zm_itr(2)
                    else                    ! a <= 0.0
                        pu = p_zm_itr(k-1)
                        pl = p_zm_itr(k)
                        ptu = pt_zm_itr(k-1)
                        ptl = pt_zm_itr(k)
                    endif

                    pt_zm_new(kk) = ptl + (ptl - ptu) * (plev_output(kk) - pl) / (pl - pu)

                    kk = kk + 1

                    if (kk == nz+1) then     ! check whether all pt_zm are updated or not
                        exit kloop
                    endif

                    cycle

                endif

                exit

            enddo

        enddo kloop

        !-- surface check
        do kk = 1, nz
            if (pt_zm_new(kk) < 0._kp) then
               if (plev_output(kk) < surface_p) then
                   pt_zm_new(kk) = pt_zm_old(kk) + (pt_zm_old(kk) &
                                  & - pt_zm_old(kk-1)) * (plev_output(kk) - pressure_zm(kk)) / (pressure_zm(kk) - pressure_zm(kk-1))
               else
                   pt_zm_new(kk) = surface_pt
               endif
            endif
        enddo

    end subroutine getpt_ptiter


    !
    ! getpt_stable() - check whether the atmosphere is stable or not
    !
    subroutine getpt_stable(ny, nz, pt, surface_p, surface_pt)
        integer , intent(in) :: ny
        integer , intent(in) :: nz
        real(kp), intent(in) :: pt(ny,nz)
        real(kp), intent(in) :: surface_p(ny)
        real(kp), intent(in) :: surface_pt(ny)

        integer :: j, k

        do k = 2, nz
            do j = 1, ny

                if (pt(j,k) > pt(j,k-1)) then
                    write(0,'(A)') "ERROR STOP"
                    write(0,'(A)') "Potential temperature is instable"
                    write(0,'("pt(",I0,",",I0,")  =",ES15.3)') j, k, pt(j,k)
                    write(0,'("pt(",I0,",",I0,"-1)=",ES15.3)') j, k, pt(j,k-1)
                    write(0,'(A,ES15.3)') "Surface pressure :", surface_p(j)
                    write(0,'(A,ES15.3)') "Surface potential temperature :", surface_pt(j)
                    ERROR STOP
                endif

            enddo
        enddo

    end subroutine getpt_stable


    !
    ! getp_stable() - check whether the atmosphere is stable or not
    !
    subroutine getp_stable(ny, nz, p, surface_p)
        integer , intent(in) :: ny
        integer , intent(in) :: nz
        real(kp), intent(in) :: p(ny,nz)
        real(kp), intent(in) :: surface_p(ny)
  
        integer :: j, k

        do k = 2, nz
            do j = 1, ny
  
                if (p(j,k) < p(j,k-1)) then
                    write(0,'(A)') "ERROR STOP"
                    write(0,'(A)') "Pressure is instable"
                    write(0,'("p(",I0,",",I0,")  =",ES15.3)') j, k, p(j,k)
                    write(0,'("p(",I0,",",I0,"-1)=",ES15.3)') j, k, p(j,k-1)
                    write(0,'(A,ES15.3)') "Surface pressure :", surface_p(j)
                    ERROR STOP
                endif
  
            enddo
        enddo

    end subroutine getp_stable


end module get_levels

