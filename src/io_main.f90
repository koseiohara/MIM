module io_main

    use params       , only : rkp, kp, grav, cp, gasr, gamma                                                             , &
                            & t_min, t_max, wind_min, wind_max, omega_min, omega_max, p_min, p_max, z_min, z_max         , &
                            & alt_min, alt_max
    use namelist     , only : INPUT_UVT_FILENAME, INPUT_U_FILENAME, INPUT_V_FILENAME, INPUT_T_FILENAME, INPUT_PS_FILENAME, &
                            & INPUT_MSL_FILENAME, INPUT_TS_FILENAME, INPUT_Z_FILENAME, INPUT_OMEGA_FILENAME              , &
                            & INPUT_TOPO_FILENAME, INPUT_Q_FILENAME, INPUT_SHORTWAVE_FILENAME, INPUT_LONGWAVE_FILENAME   , &
                            & INPUT_LHR_LARGE_FILENAME, INPUT_LHR_CONV_FILENAME, INPUT_DIFFUSION_FILENAME                , &
                            & INPUT_UNIT_Z, INPUT_UNIT_PS, INPUT_UNIT_MSL, INPUT_UNIT_TOPO, INPUT_UNIT_Q                 , &
                            & INPUT_UNDEF_UVT, INPUT_UNDEF_U, INPUT_UNDEF_V, INPUT_UNDEF_T                               , &
                            & INPUT_UNDEF_T, INPUT_UNDEF_PS, INPUT_UNDEF_MSL, INPUT_UNDEF_TS, INPUT_UNDEF_Z              , &
                            & INPUT_UNDEF_OMEGA, INPUT_UNDEF_Q, INPUT_UNDEF_SHORTWAVE, INPUT_UNDEF_LONGWAVE              , &
                            & INPUT_UNDEF_LHR_LARGE, INPUT_UNDEF_LHR_CONV, INPUT_UNDEF_DIFFUSION                         , &
                            & INPUT_ENDIAN_UVT, INPUT_ENDIAN_U, INPUT_ENDIAN_V, INPUT_ENDIAN_T, INPUT_ENDIAN_PS          , &
                            & INPUT_ENDIAN_MSL, INPUT_ENDIAN_TS, INPUT_ENDIAN_Z, INPUT_ENDIAN_OMEGA, INPUT_ENDIAN_Q      , &
                            & INPUT_ENDIAN_SHORTWAVE, INPUT_ENDIAN_LONGWAVE, INPUT_ENDIAN_LHR_LARGE                      , &
                            & INPUT_ENDIAN_LHR_CONV, INPUT_ENDIAN_DIFFUSION, INPUT_ENDIAN_TOPO                           , &
                            & INPUT_YDEF_YREV_DEFAULT, INPUT_YDEF_YREV_TOPO                                              , &
                            & INPUT_ZDEF_ZREV                                                                            , &
                            & WAVE_MAX_NUMBER                                                                            , &
                            & OUTPUT_ZONAL_FILENAME, OUTPUT_VINT_FILENAME, OUTPUT_GMEAN_FILENAME, OUTPUT_WAVE_FILENAME   , &
                            & Q_EXIST, Q_COMPS_EXIST
                            !& OUTPUT_ERROR_FILENAME                                                                      , &
    use com_var      , only : im, jm, km, ko, wmax
    use mim_var      , only : u, v, t, z, omega, p_sfc, alt                                                              , &
                            & u_zm, v_zm, pt_zm, t_dagger, st_zm, w_zm, z_zm, u_u_x_zm, epy, depy, epz_form, depz_form   , &
                            & epz_uv, depz_uv, epz_ut, depz_ut, epz_uw, depz_uw, epz, depz, divf, gy, dgy, gz, dgz       , &
                            & kz_zm, ke_zm, pz_zm, ae_zm_vint, ae_total_zm, az_zm, az_zm_vint, az_gmean, c_az_kz, c_kz_ae, &
                            & c_kz_ae_u, c_kz_ae_v, c_ae_ke, c_ae_ke_u, c_ae_ke_v, c_kz_ke, c_kz_ke_uy, c_kz_ke_uz       , &
                            & c_kz_ke_vy, c_kz_ke_vz, c_kz_ke_tan , c_kz_w                                               , &
                            & dkzdt_vkz, dkzdt_wkz, dkedt_uy, dkedt_vy                                                   , &
                            & dkedt_uz, dkedt_vz, dkedt_vke, dkedt_wke                                                   , &
                            & epz_wave                                                                                   , &
                            & q_3d, q_shortwave_3d, q_longwave_3d, q_lhr_large_3d, q_lhr_conv_3d, q_diffusion_3d         , &
                            & q_zm, q_shortwave_zm, q_longwave_zm, q_lhr_large_zm, q_lhr_conv_zm, q_diffusion_zm         , &
                            & qgz_zm, qgz_shortwave_zm, qgz_longwave_zm, qgz_lhr_large_zm, qgz_lhr_conv_zm               , &
                            & qgz_diffusion_zm, qz_vint, qz_shortwave_vint, qz_longwave_vint, qz_lhr_large_vint          , &
                            & qz_lhr_conv_vint, qz_diffusion_vint, qz_gmean, qz_shortwave_gmean, qz_longwave_gmean       , &
                            & qz_lhr_large_gmean, qz_lhr_conv_gmean, qz_diffusion_gmean, qe_zm, qe_shortwave_zm          , &
                            & qe_longwave_zm, qe_lhr_large_zm, qe_lhr_conv_zm, qe_diffusion_zm                           , &
                            & p_pds
                            !& divz_tzm, divphi_t, dwdt, d_u_epz
    use status_output, only : warn_write
    use BInIO        , only : finfo
    use undef        , only : undef_fill
    use integral     , only : integral_p

    implicit none

    private
    public :: files_open, files_close, &
            & read_alt, read_uvt, read_p_sfc, read_z, read_omega, read_q, &
            & write_zonal, write_vint, write_gmean, write_wave

    ! Input Files
    type(finfo) :: uvt_file
    type(finfo) :: u_file
    type(finfo) :: v_file
    type(finfo) :: t_file
    type(finfo) :: z_file
    type(finfo) :: omega_file
    type(finfo) :: ps_file
    type(finfo) :: msl_file
    type(finfo) :: ts_file
    type(finfo) :: q_file
    type(finfo) :: shortwave_file
    type(finfo) :: longwave_file
    type(finfo) :: lhr_large_file
    type(finfo) :: lhr_conv_file
    type(finfo) :: diffusion_file
    type(finfo) :: topo_file

    ! Output Files
    type(finfo) :: zonal_file
    type(finfo) :: vint_file
    type(finfo) :: gmean_file
    type(finfo) :: wave_file

    contains


    subroutine files_open()

        ! Input Files
        if (trim(INPUT_UVT_FILENAME) /= '') then

            uvt_file = finfo(INPUT_UVT_FILENAME     , &  !! IN : File Name
                           & 'READ'                 , &  !! IN : Action
                           & im                     , &  !! IN : Xnum
                           & jm                     , &  !! IN : Ynum
                           & km                     , &  !! IN : Znum
                           & .False.                , &  !! IN : Xrev
                           & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                           & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                           & 'UVT'                  , &  !! IN : Variable Name
                           & INPUT_ENDIAN_UVT         )  !! IN : Endian

        else if (trim(INPUT_U_FILENAME) /= '' .AND. &
               & trim(INPUT_V_FILENAME) /= '' .AND. &
               & trim(INPUT_T_FILENAME) /= ''       ) then

            u_file = finfo(INPUT_U_FILENAME       , &  !! IN : File Name
                         & 'READ'                 , &  !! IN : Action
                         & im                     , &  !! IN : Xnum
                         & jm                     , &  !! IN : Ynum
                         & km                     , &  !! IN : Znum
                         & .False.                , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                         & 'U (ZONAL VELOCITY)'   , &  !! IN : Variable Name
                         & INPUT_ENDIAN_U           )  !! IN : Endian
                           
            v_file = finfo(INPUT_V_FILENAME         , &  !! IN : File Name
                         & 'READ'                   , &  !! IN : Action
                         & im                       , &  !! IN : Xnum
                         & jm                       , &  !! IN : Ynum
                         & km                       , &  !! IN : Znum
                         & .False.                  , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT  , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV          , &  !! IN : Zrev
                         & 'V (MERIDIONAL VELOCITY)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_V             )  !! IN : Endian
                           
            t_file = finfo(INPUT_T_FILENAME       , &  !! IN : File Name
                         & 'READ'                 , &  !! IN : Action
                         & im                     , &  !! IN : Xnum
                         & jm                     , &  !! IN : Ynum
                         & km                     , &  !! IN : Znum
                         & .False.                , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                         & 'T (TEMPERATURE)'      , &  !! IN : Variable Name
                         & INPUT_ENDIAN_T           )  !! IN : Endian
                           
        else
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_UVT_FILENAME is empty'
            write(0,'(A)') 'One of INPUT_U_FILENAME, INPUT_V_FILENAME, and INPUT_T_FILENAME is empty'
            write(0,'(A)') 'Specify INPUT_UVT_FILENAME'
            write(0,'(A)') 'Otherwise, fill all of INPUT_U_FILENAME, INPUT_V_FILENAME, and INPUT_T_FILENAME'
            ERROR STOP
        endif

        if (trim(INPUT_Z_FILENAME) /= '') then

            z_file = finfo(INPUT_Z_FILENAME       , &  !! IN : File Name
                         & 'READ'                 , &  !! IN : Action
                         & im                     , &  !! IN : Xnum
                         & jm                     , &  !! IN : Ynum
                         & km                     , &  !! IN : Znum
                         & .False.                , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                         & 'Z (HIGHT)'            , &  !! IN : Variable Name
                         & INPUT_ENDIAN_Z           )  !! IN : Endian

        else
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_Z_FILENAME is empty'
            ERROR STOP
        endif

        if (trim(INPUT_OMEGA_FILENAME) /= '') then

            omega_file = finfo(INPUT_OMEGA_FILENAME   , &  !! IN : File Name
                             & 'READ'                 , &  !! IN : Action
                             & im                     , &  !! IN : Xnum
                             & jm                     , &  !! IN : Ynum
                             & km                     , &  !! IN : Znum
                             & .False.                , &  !! IN : Xrev
                             & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                             & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                             & 'OMEGA (P-VELOCITY)'   , &  !! IN : Variable Name
                             & INPUT_ENDIAN_OMEGA       )  !! IN : Endian
      
        endif

        if (trim(INPUT_PS_FILENAME) /= '') then

            ps_file = finfo(INPUT_PS_FILENAME     , &  !! IN : File Name
                         & 'READ'                 , &  !! IN : Action
                         & im                     , &  !! IN : Xnum
                         & jm                     , &  !! IN : Ynum
                         & 1                      , &  !! IN : Znum
                         & .False.                , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                         & 'PS (SURFACE PRESSURE)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_PS          )  !! IN : Endian

        else if (trim(INPUT_MSL_FILENAME) /= '' .AND. &
               & trim(INPUT_TS_FILENAME ) /= ''       ) then

            msl_file = finfo(INPUT_MSL_FILENAME           , &  !! IN : File Name
                         & 'READ'                         , &  !! IN : Action
                         & im                             , &  !! IN : Xnum
                         & jm                             , &  !! IN : Ynum
                         & 1                              , &  !! IN : Znum
                         & .False.                        , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT        , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV                , &  !! IN : Zrev
                         & 'MSL (MEAN SEA LEVEL PRESSURE)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_MSL                 )  !! IN : Endian

            ts_file = finfo(INPUT_TS_FILENAME        , &  !! IN : File Name
                         & 'READ'                    , &  !! IN : Action
                         & im                        , &  !! IN : Xnum
                         & jm                        , &  !! IN : Ynum
                         & 1                         , &  !! IN : Znum
                         & .False.                   , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT   , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV           , &  !! IN : Zrev
                         & 'TS (SURFACE TEMPERATURE)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_TS             )  !! IN : Endian

        else
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'INPUT_PS_FILENAME is empty'
            write(0,'(A)') 'One of INPUT_MSL_FILENAME and INPUT_TS_FILENAME is empty'
            write(0,'(A)') 'Specify INPUT_PS_FILENAME'
            write(0,'(A)') 'Otherwise, fill both INPUT_MSL_FILENAME and INPUT_TS_FILENAME'
            ERROR STOP
        endif

        if (trim(INPUT_Q_FILENAME) /= '') then

            q_file = finfo(INPUT_Q_FILENAME       , &  !! IN : File Name
                         & 'READ'                 , &  !! IN : Action
                         & im                     , &  !! IN : Xnum
                         & jm                     , &  !! IN : Ynum
                         & km                     , &  !! IN : Znum
                         & .False.                , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT, &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV        , &  !! IN : Zrev
                         & 'Q (DIABATIC HEATING)' , &  !! IN : Variable Name
                         & INPUT_ENDIAN_Q           )  !! IN : Endian

        endif

        if (trim(INPUT_SHORTWAVE_FILENAME) /= '') then

            shortwave_file = finfo(INPUT_SHORTWAVE_FILENAME  , &  !! IN : File Name
                         & 'READ'                            , &  !! IN : Action
                         & im                                , &  !! IN : Xnum
                         & jm                                , &  !! IN : Ynum
                         & km                                , &  !! IN : Znum
                         & .False.                           , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT           , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV                   , &  !! IN : Zrev
                         & 'SHORTWAVE (SHORT WAVE RADIATION)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_SHORTWAVE              )  !! IN : Endian

        endif

        if (trim(INPUT_LONGWAVE_FILENAME) /= '') then

            longwave_file = finfo(INPUT_LONGWAVE_FILENAME  , &  !! IN : File Name
                         & 'READ'                          , &  !! IN : Action
                         & im                              , &  !! IN : Xnum
                         & jm                              , &  !! IN : Ynum
                         & km                              , &  !! IN : Znum
                         & .False.                         , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT         , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV                 , &  !! IN : Zrev
                         & 'LONGWAVE (LONG WAVE RADIATION)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_LONGWAVE             )  !! IN : Endian

        endif

        if (trim(INPUT_LHR_LARGE_FILENAME) /= '') then

            lhr_large_file = finfo(INPUT_LHR_LARGE_FILENAME, &  !! IN : File Name
                         & 'READ'                          , &  !! IN : Action
                         & im                              , &  !! IN : Xnum
                         & jm                              , &  !! IN : Ynum
                         & km                              , &  !! IN : Znum
                         & .False.                         , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT         , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV                 , &  !! IN : Zrev
                         & 'LHR_LARGE (LARGE SCALE LHR)'   , &  !! IN : Variable Name
                         & INPUT_ENDIAN_LHR_LARGE            )  !! IN : Endian

        endif

        if (trim(INPUT_LHR_CONV_FILENAME) /= '') then

            lhr_conv_file = finfo(INPUT_LHR_CONV_FILENAME, &  !! IN : File Name
                         & 'READ'                        , &  !! IN : Action
                         & im                            , &  !! IN : Xnum
                         & jm                            , &  !! IN : Ynum
                         & km                            , &  !! IN : Znum
                         & .False.                       , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT       , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV               , &  !! IN : Zrev
                         & 'LHR_CONV (CONVECTIVE LHR)  ' , &  !! IN : Variable Name
                         & INPUT_ENDIAN_LHR_CONV           )  !! IN : Endian

        endif

        if (trim(INPUT_DIFFUSION_FILENAME) /= '') then

            diffusion_file = finfo(INPUT_DIFFUSION_FILENAME, &  !! IN : File Name
                         & 'READ'                          , &  !! IN : Action
                         & im                              , &  !! IN : Xnum
                         & jm                              , &  !! IN : Ynum
                         & km                              , &  !! IN : Znum
                         & .False.                         , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_DEFAULT         , &  !! IN : Yrev
                         & INPUT_ZDEF_ZREV                 , &  !! IN : Zrev
                         & 'DIFFUSION (VERTICAL DIFFUSION)', &  !! IN : Variable Name
                         & INPUT_ENDIAN_DIFFUSION            )  !! IN : Endian

        endif

        if (trim(INPUT_TOPO_FILENAME) /= '') then

            topo_file = finfo(INPUT_TOPO_FILENAME, &  !! IN : File Name
                         & 'READ'                , &  !! IN : Action
                         & im                    , &  !! IN : Xnum
                         & jm                    , &  !! IN : Ynum
                         & 1                     , &  !! IN : Znum
                         & .False.               , &  !! IN : Xrev
                         & INPUT_YDEF_YREV_TOPO  , &  !! IN : Yrev
                         & .False.               , &  !! IN : Zrev
                         & 'TOPO (TOPOGRAPHY)'   , &  !! IN : Variable Name
                         & INPUT_ENDIAN_TOPO       )  !! IN : Endian

        endif


        ! Output Files
        if (trim(OUTPUT_ZONAL_FILENAME) /= '') then

            zonal_file = finfo(OUTPUT_ZONAL_FILENAME, &  !! IN : File Name
                             & 'WRITE'              , &  !! IN : Xnum
                             & 1                    , &  !! IN : Xnum
                             & jm                   , &  !! IN : Ynum
                             & ko                   , &  !! IN : Znum
                             & .False.              , &  !! IN : Xrev
                             & .False.              , &  !! IN : Yrev
                             & .True.               , &  !! IN : Zrev
                             & 'ZONAL'              , &  !! IN : Variable Name
                             & 'NATIVE'               )  !! IN : Endian

        endif

        if (trim(OUTPUT_VINT_FILENAME) /= '') then

            vint_file = finfo(OUTPUT_VINT_FILENAME, &  !! IN : File Name
                             & 'WRITE'            , &  !! IN : Action
                             & 1                  , &  !! IN : Xnum
                             & jm                 , &  !! IN : Ynum
                             & 1                  , &  !! IN : Znum
                             & .False.            , &  !! IN : Xrev
                             & .False.            , &  !! IN : Yrev
                             & .True.             , &  !! IN : Zrev
                             & 'VINT'             , &  !! IN : Variable Name
                             & 'NATIVE'             )  !! IN : Endian

        endif

        if (trim(OUTPUT_GMEAN_FILENAME) /= '') then

            gmean_file = finfo(OUTPUT_GMEAN_FILENAME, &  !! IN : File Name
                             & 'WRITE'              , &  !! IN : Action
                             & 1                    , &  !! IN : Xnum
                             & 1                    , &  !! IN : Ynum
                             & 1                    , &  !! IN : Znum
                             & .False.              , &  !! IN : Xrev
                             & .False.              , &  !! IN : Yrev
                             & .True.               , &  !! IN : Zrev
                             & 'GMEAN'              , &  !! IN : Variable Name
                             & 'NATIVE'               )  !! IN : Endian

        endif

        if (trim(OUTPUT_WAVE_FILENAME) /= '') then

            wave_file = finfo(OUTPUT_WAVE_FILENAME, &  !! IN : File Name
                             & 'WRITE'            , &  !! IN : Action
                             & jm                 , &  !! IN : Xnum
                             & ko                 , &  !! IN : Ynum
                             & wmax               , &  !! IN : Znum
                             & .False.            , &  !! IN : Xrev
                             & .True.             , &  !! IN : Yrev
                             & .False.            , &  !! IN : Zrev
                             & 'WAVE'             , &  !! IN : Variable Name
                             & 'NATIVE'             )  !! IN : Endian

        endif

    end subroutine files_open


    subroutine files_close()

        call uvt_file%fclose()

        call u_file%fclose()

        call v_file%fclose()

        call t_file%fclose()

        call z_file%fclose()

        call omega_file%fclose()

        call ps_file%fclose()

        call msl_file%fclose()

        call ts_file%fclose()

        call q_file%fclose()

        call shortwave_file%fclose()

        call longwave_file%fclose()

        call lhr_large_file%fclose()

        call lhr_conv_file%fclose()

        call diffusion_file%fclose()

        call topo_file%fclose()

        call zonal_file%fclose()

        call vint_file%fclose()

        call gmean_file%fclose()

        call wave_file%fclose()

    end subroutine files_close


    subroutine read_alt()
        real(rkp) :: reader_alt(im,jm)

        call topo_file%fread(reader_alt(1:im,1:jm))  !! OUT

        alt(1:im,1:jm) = real(reader_alt(1:im,1:jm), kind=kp)

        if (trim(INPUT_UNIT_TOPO) == 'm^2/s^2') then
            alt(1:im,1:jm) = alt(1:im,1:jm) / grav
        else if (trim(INPUT_UNIT_TOPO) == 'm') then
            continue
        else
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Invalid INPUT_UNIT_TOPO : ' // trim(INPUT_UNIT_TOPO)
            ERROR STOP
        endif

        call warn_write(im            , &  !! IN
                      & jm            , &  !! IN
                      & 1             , &  !! IN
                      & alt(1:im,1:jm), &  !! IN
                      & alt_min       , &  !! IN
                      & alt_max       , &  !! IN
                      & 'alt'         , &  !! IN
                      & 'read_alt()'    )  !! IN

        if (sum(alt(1:im,jm/4)) < sum(alt(1:im,jm-jm/4))) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Altitude seems to be wrong'
            write(0,'(A)') 'Probably, INPUT_YDEF_YREV_TOPO in namelist should be changed'
            ERROR STOP
        endif

        if (maxval(alt(1:im,1:jm)) > 1.E+4_kp) then
            write(0,'(A)') 'ERROR STOP'
            write(0,'(A)') 'Ground surface altitude is too high'
            write(0,'(A)') 'Probably, INPUT_UNIT_TOPO is wrong'
            ERROR STOP
        endif

    end subroutine read_alt


    subroutine read_uvt()
        real(rkp) :: reader_u(im,jm,km)
        real(rkp) :: reader_v(im,jm,km)
        real(rkp) :: reader_t(im,jm,km)

        if (trim(INPUT_UVT_FILENAME) /= '') then

            call uvt_file%fread(reader_u(1:im,1:jm,1:km))  !! OUT

            call uvt_file%fread(reader_v(1:im,1:jm,1:km))  !! OUT

            call uvt_file%fread(reader_t(1:im,1:jm,1:km))  !! OUT

        else

            call u_file%fread(reader_u(1:im,1:jm,1:km))  !! OUT

            call v_file%fread(reader_v(1:im,1:jm,1:km))  !! OUT

            call t_file%fread(reader_t(1:im,1:jm,1:km))  !! OUT

        endif

        call undef_fill(im                      , &  !! IN
                      & jm                      , &  !! IN
                      & km                      , &  !! IN
                      & INPUT_UNDEF_U           , &  !! IN
                      & reader_u(1:im,1:jm,1:km)  )  !! INOUT

        call undef_fill(im                      , &  !! IN
                      & jm                      , &  !! IN
                      & km                      , &  !! IN
                      & INPUT_UNDEF_V           , &  !! IN
                      & reader_v(1:im,1:jm,1:km)  )  !! INOUT

        call undef_fill(im                      , &  !! IN
                      & jm                      , &  !! IN
                      & km                      , &  !! IN
                      & INPUT_UNDEF_T           , &  !! IN
                      & reader_t(1:im,1:jm,1:km)  )  !! INOUT

        u(1:im,1:jm,1:km) = real(reader_u(1:im,1:jm,1:km), kind=kp)
        v(1:im,1:jm,1:km) = real(reader_v(1:im,1:jm,1:km), kind=kp)
        t(1:im,1:jm,1:km) = real(reader_t(1:im,1:jm,1:km), kind=kp)

        call warn_write(im               , &  !! IN
                      & jm               , &  !! IN
                      & km               , &  !! IN
                      & u(1:im,1:jm,1:km), &  !! IN
                      & wind_min         , &  !! IN
                      & wind_max         , &  !! IN
                      & 'u'              , &  !! IN
                      & 'read_uvt()'       )  !! IN

        call warn_write(im               , &  !! IN
                      & jm               , &  !! IN
                      & km               , &  !! IN
                      & v(1:im,1:jm,1:km), &  !! IN
                      & wind_min         , &  !! IN
                      & wind_max         , &  !! IN
                      & 'v'              , &  !! IN
                      & 'read_uvt()'       )  !! IN

        call warn_write(im               , &  !! IN
                      & jm               , &  !! IN
                      & km               , &  !! IN
                      & t(1:im,1:jm,1:km), &  !! IN
                      & t_min            , &  !! IN
                      & t_max            , &  !! IN
                      & 't'              , &  !! IN
                      & 'read_uvt()'       )  !! IN

    end subroutine read_uvt


    subroutine read_p_sfc()
        real(rkp) :: reader_p_sfc(im,jm)
        real(rkp) :: reader_msl(im,jm)
        real(rkp) :: reader_t_sfc(im,jm)
        real(kp)  :: msl(im,jm)
        real(kp)  :: t_sfc(im,jm)

        if (trim(INPUT_PS_FILENAME) /= '') then

            call ps_file%fread(reader_p_sfc(1:im,1:jm))  !! OUT

            p_sfc(1:im,1:jm) = real(reader_p_sfc(1:im,1:jm), kind=kp)

            if (trim(INPUT_UNIT_PS) == 'Pa') then
                p_sfc(1:im,1:jm) = p_sfc(1:im,1:jm) * 0.01_kp
            endif

        else

            call msl_file%fread(reader_msl(1:im,1:jm))   !! OUT

            call ts_file%fread(reader_t_sfc(1:im,1:jm))  !! OUT

            msl(1:im,1:jm)   = real(reader_msl(1:im,1:jm)  , kind=kp)
            t_sfc(1:im,1:jm) = real(reader_t_sfc(1:im,1:jm), kind=kp)

            if (trim(INPUT_UNIT_MSL) == 'Pa') then
                msl(1:im,1:jm) = msl(1:im,1:jm) * 0.01_kp
            endif

            p_sfc(1:im,1:jm) = msl(1:im,1:jm) * (1._kp + gamma*alt(1:im,1:jm)/t_sfc(1:im,1:jm))**(-grav / (gasr*gamma))

        endif

        call warn_write(im              , &  !! IN
                      & jm              , &  !! IN
                      & 1               , &  !! IN
                      & p_sfc(1:im,1:jm), &  !! IN
                      & p_min           , &  !! IN
                      & p_max           , &  !! IN
                      & 'p_sfc'         , &  !! IN
                      & 'read_p_sfc()'    )  !! IN

    end subroutine read_p_sfc


    subroutine read_z()
        real(rkp) :: reader_z(im,jm,km)

        call z_file%fread(reader_z(1:im,1:jm,1:km))

        call undef_fill(im                    , &  !! IN
                      & jm                    , &  !! IN
                      & km                    , &  !! IN
                      & INPUT_UNDEF_Z         , &  !! IN
                      & reader_z(1:im,1:jm,1:km)  )  !! INOUT

        z(1:im,1:jm,1:km) = real(reader_z(1:im,1:jm,1:km), kind=kp)

        if (trim(INPUT_UNIT_Z) == 'm^2/s^2') then
            z(1:im,1:jm,1:km) = z(1:im,1:jm,1:km) * (1._kp/grav)
        endif

        call warn_write(im               , &  !! IN
                      & jm               , &  !! IN
                      & km               , &  !! IN
                      & z(1:im,1:jm,1:km), &  !! IN
                      & z_min            , &  !! IN
                      & z_max            , &  !! IN
                      & 'z'              , &  !! IN
                      & 'read_z()'         )  !! IN

    end subroutine read_z


    subroutine read_omega()
        real(rkp) :: reader_omega(im,jm,km)

        if (trim(INPUT_OMEGA_FILENAME) /= '') then

            call omega_file%fread(reader_omega(1:im,1:jm,1:km))  !! OUT

            call undef_fill(im                        , &  !! IN
                          & jm                        , &  !! IN
                          & km                        , &  !! IN
                          & INPUT_UNDEF_OMEGA         , &  !! IN
                          & reader_omega(1:im,1:jm,1:km)  )  !! INOUT

            omega(1:im,1:jm,1:km) = real(reader_omega(1:im,1:jm,1:km), kind=kp)

        else
            omega(1:im,1:jm,1:km) = 0._kp
        endif

        call warn_write(im                   , &  !! IN
                      & jm                   , &  !! IN
                      & km                   , &  !! IN
                      & omega(1:im,1:jm,1:km), &  !! IN
                      & omega_min            , &  !! IN
                      & omega_max            , &  !! IN
                      & 'omega'              , &  !! IN
                      & 'read_omega()'         )  !! IN

    end subroutine read_omega


    subroutine read_q()
        real(rkp) :: read_q_3d(im,jm,km)
        real(rkp) :: read_q_shortwave_3d(im,jm,km)
        real(rkp) :: read_q_longwave_3d(im,jm,km)
        real(rkp) :: read_q_lhr_large_3d(im,jm,km)
        real(rkp) :: read_q_lhr_conv_3d(im,jm,km)
        real(rkp) :: read_q_diffusion_3d(im,jm,km)

        if (trim(INPUT_Q_FILENAME) /= '') then
            
            call q_file%fread(read_q_3d(1:im,1:jm,1:km))  !! OUT
  
            call undef_fill(im                       , &  !! IN
                          & jm                       , &  !! IN
                          & km                       , &  !! IN
                          & INPUT_UNDEF_Q            , &  !! IN
                          & read_q_3d(1:im,1:jm,1:km)  )  !! INOUT

            q_3d(1:im,1:jm,1:km) = real(read_q_3d(1:im,1:jm,1:km), kind=kp)

            if (trim(INPUT_UNIT_Q) == 'K/day') then
                q_3d(1:im,1:jm,1:km) = q_3d(1:im,1:jm,1:km) * (1._kp / (24._kp*60._kp*60._kp))
            endif

        endif

        if (trim(INPUT_SHORTWAVE_FILENAME) /= '') then
            
            call shortwave_file%fread(read_q_shortwave_3d(1:im,1:jm,1:km))  !! OUT

            call undef_fill(im                                 , &  !! IN
                          & jm                                 , &  !! IN
                          & km                                 , &  !! IN
                          & INPUT_UNDEF_SHORTWAVE              , &  !! IN
                          & read_q_shortwave_3d(1:im,1:jm,1:km)  )  !! INOUT

            q_shortwave_3d(1:im,1:jm,1:km) = real(read_q_shortwave_3d(1:im,1:jm,1:km), kind=kp)

            if (trim(INPUT_UNIT_Q) == 'K/day') then
                q_shortwave_3d(1:im,1:jm,1:km) = q_shortwave_3d(1:im,1:jm,1:km) * (1._kp / (24._kp*60._kp*60._kp))
            endif

        endif

        if (trim(INPUT_LONGWAVE_FILENAME) /= '') then
            
            call longwave_file%fread(read_q_longwave_3d(1:im,1:jm,1:km))  !! OUT

            call undef_fill(im                                , &  !! IN
                          & jm                                , &  !! IN
                          & km                                , &  !! IN
                          & INPUT_UNDEF_LONGWAVE              , &  !! IN
                          & read_q_longwave_3d(1:im,1:jm,1:km)  )  !! INOUT

            q_longwave_3d(1:im,1:jm,1:km) = real(read_q_longwave_3d(1:im,1:jm,1:km), kind=kp)

            if (trim(INPUT_UNIT_Q) == 'K/day') then
                q_longwave_3d(1:im,1:jm,1:km) = q_longwave_3d(1:im,1:jm,1:km) * (1._kp / (24._kp*60._kp*60._kp))
            endif

        endif

        if (trim(INPUT_LHR_LARGE_FILENAME) /= '') then
            
            call lhr_large_file%fread(read_q_lhr_large_3d(1:im,1:jm,1:km))  !! OUT

            call undef_fill(im                                 , &  !! IN
                          & jm                                 , &  !! IN
                          & km                                 , &  !! IN
                          & INPUT_UNDEF_LHR_LARGE              , &  !! IN
                          & read_q_lhr_large_3d(1:im,1:jm,1:km)  )  !! INOUT

            q_lhr_large_3d(1:im,1:jm,1:km) = real(read_q_lhr_large_3d(1:im,1:jm,1:km), kind=kp)

            if (trim(INPUT_UNIT_Q) == 'K/day') then
                q_lhr_large_3d(1:im,1:jm,1:km) = q_lhr_large_3d(1:im,1:jm,1:km) * (1._kp / (24._kp*60._kp*60._kp))
            endif

        endif

        if (trim(INPUT_LHR_CONV_FILENAME) /= '') then
            
            call lhr_conv_file%fread(read_q_lhr_conv_3d(1:im,1:jm,1:km))  !! OUT

            call undef_fill(im                                , &  !! IN
                          & jm                                , &  !! IN
                          & km                                , &  !! IN
                          & INPUT_UNDEF_LHR_CONV              , &  !! IN
                          & read_q_lhr_conv_3d(1:im,1:jm,1:km)  )  !! INOUT

            q_lhr_conv_3d(1:im,1:jm,1:km) = real(read_q_lhr_conv_3d(1:im,1:jm,1:km), kind=kp)

            if (trim(INPUT_UNIT_Q) == 'K/day') then
                q_lhr_conv_3d(1:im,1:jm,1:km) = q_lhr_conv_3d(1:im,1:jm,1:km) * (1._kp / (24._kp*60._kp*60._kp))
            endif

        endif

        if (trim(INPUT_DIFFUSION_FILENAME) /= '') then
            
            call diffusion_file%fread(read_q_diffusion_3d(1:im,1:jm,1:km))  !! OUT

            call undef_fill(im                            , &  !! IN
                          & jm                            , &  !! IN
                          & km                            , &  !! IN
                          & INPUT_UNDEF_DIFFUSION         , &  !! IN
                          & read_q_diffusion_3d(1:im,1:jm,1:km)  )  !! INOUT

            q_diffusion_3d(1:im,1:jm,1:km) = real(read_q_diffusion_3d(1:im,1:jm,1:km), kind=kp)

            if (trim(INPUT_UNIT_Q) == 'K/day') then
                q_diffusion_3d(1:im,1:jm,1:km) = q_diffusion_3d(1:im,1:jm,1:km) * (1._kp / (24._kp*60._kp*60._kp))
            endif

        endif

        if ((.NOT. Q_EXIST) .AND. Q_COMPS_EXIST) then
            
            q_3d(1:im,1:jm,1:km) = q_shortwave_3d(1:im,1:jm,1:km) &
                               & + q_longwave_3d(1:im,1:jm,1:km)  &
                               & + q_lhr_large_3d(1:im,1:jm,1:km) &
                               & + q_lhr_conv_3d(1:im,1:jm,1:km)  &
                               & + q_diffusion_3d(1:im,1:jm,1:km)

        endif

        q_3d(1:im,1:jm,1:km)           = q_3d(1:im,1:jm,1:km) * cp
        q_shortwave_3d(1:im,1:jm,1:km) = q_shortwave_3d(1:im,1:jm,1:km) * cp
        q_longwave_3d(1:im,1:jm,1:km)  = q_longwave_3d(1:im,1:jm,1:km)  * cp
        q_lhr_large_3d(1:im,1:jm,1:km) = q_lhr_large_3d(1:im,1:jm,1:km) * cp
        q_lhr_conv_3d(1:im,1:jm,1:km)  = q_lhr_conv_3d(1:im,1:jm,1:km)  * cp
        q_diffusion_3d(1:im,1:jm,1:km) = q_diffusion_3d(1:im,1:jm,1:km) * cp

    end subroutine read_q


    subroutine write_zonal(tt)
        integer, intent(in) :: tt
        
        call zonal_file%fwrite( tt, 'u_zm'             , u_zm )
        call zonal_file%fwrite( tt, 'v_zm'             , v_zm )
        call zonal_file%fwrite( tt, 'pt_zm'            , pt_zm )
        call zonal_file%fwrite( tt, 't_dagger'         , t_dagger )
        call zonal_file%fwrite( tt, 'st_zm'            , st_zm )

        call zonal_file%fwrite( tt, 'w_zm'             , w_zm )
        call zonal_file%fwrite( tt, 'z_zm'             , z_zm )
        call zonal_file%fwrite( tt, 'epy'              , epy )
        call zonal_file%fwrite( tt, 'depy'             , depy )
        call zonal_file%fwrite( tt, 'epz_form'         , epz_form )

        call zonal_file%fwrite( tt, 'depz_form'        , depz_form )
        call zonal_file%fwrite( tt, 'epz_uv'           , epz_uv )
        call zonal_file%fwrite( tt, 'depz_uv'          , depz_uv )
        call zonal_file%fwrite( tt, 'epz_ut'           , epz_ut )
        call zonal_file%fwrite( tt, 'depz_ut'          , depz_ut )

        call zonal_file%fwrite( tt, 'epz'              , epz )
        call zonal_file%fwrite( tt, 'depz'             , depz )
        call zonal_file%fwrite( tt, 'divf'             , divf )
        call zonal_file%fwrite( tt, 'gy'               , gy )
        call zonal_file%fwrite( tt, 'dgy'              , dgy )

        call zonal_file%fwrite( tt, 'gz'               , gz )
        call zonal_file%fwrite( tt, 'dgz'              , dgz )
        call zonal_file%fwrite( tt, 'u_u_x_zm'         , u_u_x_zm )
        call zonal_file%fwrite( tt, 'c_az_kz'          , c_az_kz )
        call zonal_file%fwrite( tt, 'c_kz_ae_u'        , c_kz_ae_u )

        call zonal_file%fwrite( tt, 'c_kz_ae_v'        , c_kz_ae_v )
        call zonal_file%fwrite( tt, 'c_kz_ae'          , c_kz_ae )
        call zonal_file%fwrite( tt, 'c_ae_ke_u'        , c_ae_ke_u )
        call zonal_file%fwrite( tt, 'c_ae_ke_v'        , c_ae_ke_v )
        call zonal_file%fwrite( tt, 'c_ae_ke'          , c_ae_ke )

        call zonal_file%fwrite( tt, 'c_kz_ke_uy'       , c_kz_ke_uy )
        call zonal_file%fwrite( tt, 'c_kz_ke_uz'       , c_kz_ke_uz )
        call zonal_file%fwrite( tt, 'c_kz_ke_vy'       , c_kz_ke_vy )
        call zonal_file%fwrite( tt, 'c_kz_ke_vz'       , c_kz_ke_vz )
        call zonal_file%fwrite( tt, 'c_kz_ke_tan'      , c_kz_ke_tan )

        call zonal_file%fwrite( tt, 'c_kz_ke'          , c_kz_ke )
        call zonal_file%fwrite( tt, 'c_kz_w'           , c_kz_w )
        call zonal_file%fwrite( tt, 'q_zm'             , q_zm )
        call zonal_file%fwrite( tt, 'q_shortwave_zm'   , q_shortwave_zm )
        call zonal_file%fwrite( tt, 'q_longwave_zm'    , q_longwave_zm )

        call zonal_file%fwrite( tt, 'q_lhr_large_zm'   , q_lhr_large_zm )
        call zonal_file%fwrite( tt, 'q_lhr_conv_zm'    , q_lhr_conv_zm )
        call zonal_file%fwrite( tt, 'q_diffusion_zm'   , q_diffusion_zm )
        call zonal_file%fwrite( tt, 'qgz_zm'           , qgz_zm )
        call zonal_file%fwrite( tt, 'qgz_shortwave_zm' , qgz_shortwave_zm )

        call zonal_file%fwrite( tt, 'qgz_longwave_zm'  , qgz_longwave_zm )
        call zonal_file%fwrite( tt, 'qgz_lhr_large_zm' , qgz_lhr_large_zm )
        call zonal_file%fwrite( tt, 'qgz_lhr_conv_zm'  , qgz_lhr_conv_zm )
        call zonal_file%fwrite( tt, 'qgz_diffusion_zm' , qgz_diffusion_zm )
        call zonal_file%fwrite( tt, 'qe_zm'            , qe_zm )

        call zonal_file%fwrite( tt, 'qe_shortwave_zm'  , qe_shortwave_zm )
        call zonal_file%fwrite( tt, 'qe_longwave_zm'   , qe_longwave_zm )
        call zonal_file%fwrite( tt, 'qe_lhr_large_zm'  , qe_lhr_large_zm )
        call zonal_file%fwrite( tt, 'qe_lhr_conv_zm'   , qe_lhr_conv_zm )
        call zonal_file%fwrite( tt, 'qe_diffusion_zm'  , qe_diffusion_zm )

        call zonal_file%fwrite( tt, 'kz_zm'            , kz_zm )
        call zonal_file%fwrite( tt, 'ke_zm'            , ke_zm )
        call zonal_file%fwrite( tt, 'pz_zm'            , pz_zm )
        !call zonal_file%fwrite( tt, 'az_zm'            , az_zm )
        call zonal_file%fwrite( tt, 'ae_total_zm'      , ae_total_zm )

        call zonal_file%fwrite( tt, 'dkzdt_vkz'        , dkzdt_vkz )
        call zonal_file%fwrite( tt, 'dkzdt_wkz'        , dkzdt_wkz )
        call zonal_file%fwrite( tt, 'dkedt_uy'         , dkedt_uy )
        call zonal_file%fwrite( tt, 'dkedt_vy'         , dkedt_vy )
        call zonal_file%fwrite( tt, 'dkedt_uz'         , dkedt_uz )

        call zonal_file%fwrite( tt, 'dkedt_vz'         , dkedt_vz )
        call zonal_file%fwrite( tt, 'dkedt_vke'        , dkedt_vke )
        call zonal_file%fwrite( tt, 'dkedt_wke'        , dkedt_wke )
        !call zonal_file%fwrite( tt, 'd_u_epz'          , d_u_epz )   ! not computed
        !call zonal_file%fwrite( tt, 'divz_tzm'         , divz_tzm )  ! not computed
        !call zonal_file%fwrite( tt, 'divphi_t'         , divphi_t )  ! not computed

        !call zonal_file%fwrite( tt, 'dwdt'             , dwdt )      ! not computed

    end subroutine write_zonal


    subroutine write_vint(tt)
        integer, intent(in) :: tt
        real(kp) :: work_vint(jm)

        call integral_p(jm, ko, p_pds, kz_zm, work_vint)
        call vint_file%fwrite( tt, 'kz_zm_vint', work_vint )

        call integral_p(jm, ko, p_pds, ke_zm, work_vint)
        call vint_file%fwrite(tt, 'ke_zm_vint', work_vint )

        call vint_file%fwrite( tt, 'az_zm_vint', az_zm_vint )

        call vint_file%fwrite( tt, 'ae_zm_vint', ae_zm_vint )

        call integral_p(jm, ko, p_pds, c_az_kz, work_vint)
        call vint_file%fwrite( tt, 'c_az_kz_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ae_u, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ae_u_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ae_v, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ae_v_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ae, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ae_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_ae_ke_u, work_vint)
        call vint_file%fwrite( tt, 'c_ae_ke_u_vint', work_vint )
        ! 10
        call integral_p(jm, ko, p_pds, c_ae_ke_v, work_vint)
        call vint_file%fwrite( tt, 'c_ae_ke_v_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_ae_ke, work_vint)
        call vint_file%fwrite( tt, 'c_ae_ke_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ke_uy, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ke_uy_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ke_uz, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ke_uz_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ke_vy, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ke_vy_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ke_vz, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ke_vz_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ke_tan, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ke_tan_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_ke, work_vint)
        call vint_file%fwrite( tt, 'c_kz_ke_vint', work_vint )

        call integral_p(jm, ko, p_pds, c_kz_w, work_vint)
        call vint_file%fwrite( tt, 'c_kz_w_vint', work_vint )

        call integral_p( jm, ko, p_pds, q_zm, work_vint )
        call vint_file%fwrite( tt, 'q_zm_vint', work_vint )
        ! 20
        call integral_p( jm, ko, p_pds, q_shortwave_zm, work_vint )
        call vint_file%fwrite( tt, 'q_shortwave_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, q_longwave_zm, work_vint )
        call vint_file%fwrite( tt, 'q_longwave_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, q_lhr_large_zm, work_vint )
        call vint_file%fwrite( tt, 'q_lhr_large_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, q_lhr_conv_zm, work_vint )
        call vint_file%fwrite( tt, 'q_lhr_conv_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, q_diffusion_zm, work_vint )
        call vint_file%fwrite( tt, 'q_diffusion_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qgz_zm, work_vint )
        call vint_file%fwrite( tt, 'qgz_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qgz_shortwave_zm, work_vint )
        call vint_file%fwrite( tt, 'qgz_shortwave_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qgz_longwave_zm, work_vint )
        call vint_file%fwrite( tt, 'qgz_longwave_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qgz_lhr_large_zm, work_vint )
        call vint_file%fwrite( tt, 'qgz_lhr_large_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qgz_lhr_conv_zm, work_vint )
        call vint_file%fwrite( tt, 'qgz_lhr_conv_zm_vint', work_vint )
        ! 30
        call integral_p( jm, ko, p_pds, qgz_diffusion_zm, work_vint )
        call vint_file%fwrite( tt, 'qgz_diffusion_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qe_zm, work_vint )
        call vint_file%fwrite( tt, 'qe_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qe_shortwave_zm, work_vint )
        call vint_file%fwrite( tt, 'qe_shortwave_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qe_longwave_zm, work_vint )
        call vint_file%fwrite( tt, 'qe_longwave_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qe_lhr_large_zm, work_vint )
        call vint_file%fwrite( tt, 'qe_lhr_large_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qe_lhr_conv_zm, work_vint )
        call vint_file%fwrite( tt, 'qe_lhr_conv_zm_vint', work_vint )

        call integral_p( jm, ko, p_pds, qe_diffusion_zm, work_vint )
        call vint_file%fwrite( tt, 'qe_diffusion_zm_vint', work_vint )

        call vint_file%fwrite( tt, 'qz_vint', qz_vint )
        call vint_file%fwrite( tt, 'qz_shortwave_vint', qz_shortwave_vint )
        call vint_file%fwrite( tt, 'qz_longwave_vint' , qz_longwave_vint )
        ! 40
        call vint_file%fwrite( tt, 'qz_lhr_large_vint', qz_lhr_large_vint )
        call vint_file%fwrite( tt, 'qz_lhr_conv_vint' , qz_lhr_conv_vint )
        call vint_file%fwrite( tt, 'qz_diffusion_vint', qz_diffusion_vint )

        call integral_p( jm, ko, p_pds, dkzdt_vkz, work_vint )
        call vint_file%fwrite( tt, 'dkzdt_vkz_vint', work_vint )

        call integral_p( jm, ko, p_pds, dkzdt_wkz, work_vint )
        call vint_file%fwrite( tt, 'dkzdt_wkz_vint', work_vint )

        call integral_p( jm, ko, p_pds, dkedt_uy, work_vint )
        call vint_file%fwrite( tt, 'dkedt_uy_vint', work_vint )

        call integral_p( jm, ko, p_pds, dkedt_vy, work_vint )
        call vint_file%fwrite( tt, 'dkedt_vy_vint', work_vint )

        call integral_p( jm, ko, p_pds, dkedt_uz, work_vint )
        call vint_file%fwrite( tt, 'dkedt_uz_vint', work_vint )

        call integral_p( jm, ko, p_pds, dkedt_vz, work_vint )
        call vint_file%fwrite( tt, 'dkedt_vz_vint', work_vint )

        call integral_p( jm, ko, p_pds, dkedt_vke, work_vint )
        call vint_file%fwrite( tt, 'dkedt_vke_vint', work_vint )
        ! 50
        call integral_p( jm, ko, p_pds, dkedt_wke, work_vint )
        call vint_file%fwrite( tt, 'dkedt_wke_vint', work_vint )

        !call integral_p( jm, ko, p_pds, d_u_epz, work_vint )
        !call vint_file%fwrite( tt, 'd_u_epz_vint', work_vint )

        !call integral_p( jm, ko, p_pds, divz_tzm, work_vint )
        !call vint_file%fwrite( tt, 'divz_tzm_vint', work_vint )

        !call integral_p( jm, ko, p_pds, divphi_t, work_vint )
        !call vint_file%fwrite( tt, 'divphi_t_vint', work_vint )

    end subroutine write_vint


    subroutine write_gmean(tt)
        integer, intent(in) :: tt

        call gmean_file%fwrite( tt, 'az_gmean'          , az_gmean )
        call gmean_file%fwrite( tt, 'qz_gmean'          , qz_gmean )
        call gmean_file%fwrite( tt, 'qz_shortwave_gmean', qz_shortwave_gmean )
        call gmean_file%fwrite( tt, 'qz_longwave_gmean' , qz_longwave_gmean )
        call gmean_file%fwrite( tt, 'qz_lhr_large_gmean', qz_lhr_large_gmean )
        call gmean_file%fwrite( tt, 'qz_lhr_conv_gmean' , qz_lhr_conv_gmean )
        call gmean_file%fwrite( tt, 'qz_diffusion_gmean', qz_diffusion_gmean )

    end subroutine write_gmean


    subroutine write_wave(tt)
        integer, intent(in) :: tt
        
        call wave_file%fwrite( tt, 'epz_wave', epz_wave(1:jm,1:ko,1:wmax) )

    end subroutine write_wave


end module io_main

