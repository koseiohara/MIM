!#########################################################
!
!  Variable list for mim.f90
!
!    "use variable" should be specified only in mim.f90
!
!#########################################################
module mim_var

    use params , only : kp
    use com_var, only : im, jm, km, ko, wmax

    implicit none
  
    !********** input **********!
  
    !***** pressure-surface *****!
    real(kp), allocatable :: u(:,:,:)     ! zonal wind [m/s]
    real(kp), allocatable :: v(:,:,:)     ! meridional Wind [m/s]
    real(kp), allocatable :: t(:,:,:)     ! temperature [K]
    real(kp), allocatable :: z(:,:,:)     ! geopotential height [m]
    real(kp), allocatable :: omega(:,:,:) ! p-velocity [Pa/s]
  
    !***** surface *****!
    real(kp), allocatable :: p_sfc(:,:)   ! surface pressure [hPa]
    !real(kp), allocatable :: msl(:,:)     ! mean sea level pressure [hPa]
    !real(kp), allocatable :: t_sfc(:,:)   ! surface temperature [K]
    real(kp), allocatable :: alt(:,:)     ! surface altitude [m]
  
  
    !********** derived **********!
  
    !*** 3D basic variables on the standard pressure levels ***!
    real(kp), allocatable :: pt(:,:,:)     ! potential temperature
    real(kp), allocatable :: pt_dot(:,:,:) ! D(pt)/Dt
  
    !*** 2D zonal mean basic variables on the standard p+ levels ***!
    real(kp), allocatable :: p_zm(:,:)    ! exact p+ [hPa] (approx. equal pout)
    real(kp), allocatable :: u_zm(:,:)    ! zonal wind [m/s]
    real(kp), allocatable :: v_zm(:,:)     ! meridional wind [m/s]
    real(kp), allocatable :: pt_zm(:,:)    ! potential temperature [K]
    real(kp), allocatable :: t_dagger(:,:) ! temperature dagger (from pt_zm) [K]
    ! real(kp), allocatable :: t_zm(:,:)     ! temperature (from t) [K]
    real(kp), allocatable :: st_zm(:,:)    ! mass streamfunction [kg/s]
    real(kp), allocatable :: w_zm (:,:)    ! vertical velocity (from st_zm) [m/s]
    real(kp), allocatable :: z_zm(:,:)     ! geopotential height [m]
    real(kp), allocatable :: pt_dot_zm(:,:)! D(pt)/Dt [K/s]
  
    !*** correlarion ***!
    real(kp), allocatable :: u_u_zm(:,:)      ! (u^2)_zm [m^2/s^2]
    real(kp), allocatable :: u_u_x_zm(:,:)    ! (u'u')_zm [m^2/s^2]
    real(kp), allocatable :: v_v_zm(:,:)
    real(kp), allocatable :: v_v_x_zm(:,:)
    real(kp), allocatable :: u_v_zm(:,:)
    real(kp), allocatable :: u_v_x_zm(:,:)
    real(kp), allocatable :: u_u_v_zm(:,:)     ! (u^2 v)_zm
    real(kp), allocatable :: v_v_v_zm(:,:)     ! (v^3)_zm
    real(kp), allocatable :: u_pt_dot_zm(:,:)
    real(kp), allocatable :: u_pt_dot_x_zm(:,:)
    real(kp), allocatable :: v_pt_dot_zm(:,:)
    real(kp), allocatable :: v_pt_dot_x_zm(:,:)
    real(kp), allocatable :: u_u_pt_dot_zm(:,:)
    real(kp), allocatable :: v_v_pt_dot_zm(:,:)
    !real(kp), allocatable :: v_ke_zm(:,:)      ! (v Ke)_zm
    !real(kp), allocatable :: pt_dot_ke_zm(:,:) ! (pt_dot Ke)_zm
  
    !*** EP flux & G Flux ***!
    real(kp), allocatable :: epy(:,:)      ! Fy [kg/s^2]
    real(kp), allocatable :: depy(:,:)     ! divF due to epy [m/s^2]
    real(kp), allocatable :: epz_form(:,:) ! Fz (Form Drag) [kg/s^2]
    real(kp), allocatable :: depz_form(:,:)! divF due to epz_form [m/s^2]
    real(kp), allocatable :: epz_uv(:,:)   ! part of Fz (u'w') [kg/s^2]
    real(kp), allocatable :: depz_uv(:,:)  ! divF due to epz_uv [m/s^2]
    real(kp), allocatable :: epz_ut(:,:)   ! part of Fz (u'w') [kg/s^2]
    real(kp), allocatable :: depz_ut(:,:)  ! divF due to epz_ut [m/s^2]
    real(kp), allocatable :: epz_uw(:,:)   ! Fz (u'w') = epz_uv + epz_ut [kg/s^2]
    real(kp), allocatable :: depz_uw(:,:)  ! divF due to epz_uw [m/s^2]
    real(kp), allocatable :: epz(:,:)      ! Fz (Total) [kg/s^2]
    real(kp), allocatable :: depz(:,:)     ! divF due to epz [m/s^2]
    real(kp), allocatable :: divf(:,:)     ! EP Flux Divergence [m/s^2]
    real(kp), allocatable :: gy(:,:)       ! Gy [kg/s^2]
    real(kp), allocatable :: dgy(:,:)      ! dGy/dy [m/s^2]
    real(kp), allocatable :: gz(:,:)       ! Gz [kg/s^2]
    real(kp), allocatable :: dgz(:,:)      ! dGz/dy [m/s^2]
  
    !*** energy ***!
    real(kp), allocatable :: kz_zm(:,:)    ! zonal kinetic energy [m^2/s^2]
    real(kp), allocatable :: ke_zm(:,:)    ! eddy kinetic energy  [m^2/s^2]
    real(kp), allocatable :: pz_zm(:,:)    ! potential energy [m^2/s^2]
    real(kp), allocatable :: ae_zm_vint(:) ! eddy available potential energy [J/m^2]
    real(kp), allocatable :: ae_total_zm(:,:) ! (not recommended to use)
    real(kp), allocatable :: ae_s_zm(:)    ! surface correction
    real(kp), allocatable :: az_zm(:,:)
    real(kp), allocatable :: az_zm_vint(:) ! (not recommended to use)
    real(kp)              :: az_gmean(1)   ! zonal available potential energy [J/m^2]
  
    !*** energy conversion (for global mean) ***!
    real(kp), allocatable :: c_az_kz(:,:)   ! C(Az,Kz) [W/m^2]
    real(kp), allocatable :: c_az_kz_modify(:,:)
    real(kp), allocatable :: c_kz_ae(:,:)   ! C(Kz,Ae) [W/m^2]
    real(kp), allocatable :: c_kz_ae_u(:,:)
    real(kp), allocatable :: c_kz_ae_v(:,:)
    real(kp), allocatable :: c_ae_ke(:,:)   ! C(Ae,Ke) [W/m^2]
    real(kp), allocatable :: c_ae_ke_u(:,:)
    real(kp), allocatable :: c_ae_ke_v(:,:)
    real(kp), allocatable :: c_kz_ke(:,:)   ! C(Kz,ke) [W/m^2]
    real(kp), allocatable :: c_kz_ke_uy(:,:)
    real(kp), allocatable :: c_kz_ke_uz(:,:)
    real(kp), allocatable :: c_kz_ke_vy(:,:)
    real(kp), allocatable :: c_kz_ke_vz(:,:)
    real(kp), allocatable :: c_kz_ke_tan(:,:)
    real(kp), allocatable :: c_kz_w(:,:)   ! C(Kz,W) [W/m^2]
  
    !*** for p -> pd ***!
    integer , allocatable :: nlev(:,:,:)
    real(kp), allocatable :: dlev(:,:,:)
    real(kp), allocatable :: p_pd(:,:,:)
    real(kp), allocatable :: x_pd(:,:,:)
    real(kp), allocatable :: xint_zm(:,:)
    real(kp), allocatable :: pt_sfc(:,:)
    real(kp), allocatable :: pt_pds(:)    ! pt_{xmin}
    real(kp), allocatable :: pd_p(:,:,:)  ! pd(x,y,p)
  
    !*** for pd -> pdd ***!
    integer , allocatable :: nlev_y(:,:)
    real(kp), allocatable :: dlev_y(:,:)
    real(kp), allocatable :: pd_pdd(:,:)
    real(kp), allocatable :: pd_ym(:)
    real(kp), allocatable :: pt_ym(:)
    real(kp)              :: pt_pdds(1)   ! pt_pds_{ymin}
    real(kp), allocatable :: pdd_pd(:,:)  ! pdd(y,pd)
  
  
    real(kp), allocatable :: dz_dlat_zm(:,:)
    real(kp), allocatable :: v_dz_dlat_zm(:,:)
    real(kp), allocatable :: u_dz_dlon_zm(:,:)
    real(kp), allocatable :: temp_vint(:) !xxx
    real(kp), allocatable :: phi_dagger(:,:)
    real(kp), allocatable :: phi_dagger_y(:,:)
    real(kp), allocatable :: p_dphi_dt(:,:)
    !real(kp), allocatable :: p_dz_dt_zm(:,:)
    !real(kp), allocatable :: p_dz_dt(:,:,:)
    !real(kp), allocatable :: divz_tzm(:,:)
    !real(kp), allocatable :: divphi_t(:,:)
    !real(kp), allocatable :: dwdt(:,:)
    !real(kp), allocatable :: uuv_tmp(:,:)
    ! real(kp), allocatable :: d_u_epy(:,:)
    !real(kp), allocatable :: d_u_epz(:,:) !local
    ! real(kp), allocatable :: dkedt_uvu_y(:,:)
  
    !*** energy conversion (local) ***!
    real(kp), allocatable :: dkzdt_vkz(:,:)
    real(kp), allocatable :: dkzdt_wkz(:,:)
    real(kp), allocatable :: dkedt_uy(:,:)
    real(kp), allocatable :: dkedt_vy(:,:)
    real(kp), allocatable :: dkedt_uz(:,:)
    real(kp), allocatable :: dkedt_vz(:,:)
    real(kp), allocatable :: dkedt_vke(:,:)
    real(kp), allocatable :: dkedt_wke(:,:)
    !real(kp), allocatable :: dpedt_vt(:,:)
    !real(kp), allocatable :: dpedt_wt(:,:)
  
    real(kp), allocatable :: z_pd(:,:,:)
    real(kp), allocatable :: u_past(:,:,:)
    real(kp), allocatable :: v_past(:,:,:)
    real(kp), allocatable :: omega_past(:,:,:)
    real(kp), allocatable :: pt_past(:,:,:)
    !real(kp), allocatable :: phi_dagger_past(:,:)
    real(kp), allocatable :: z_pd_past(:,:,:)
    real(kp), allocatable :: p_pd_past(:,:,:)
  
    !*** wavenumber decomposition ***!
    real(kp), allocatable :: epz_wave(:,:,:)
    real(kp), allocatable :: z_pt_wave(:,:,:,:)
    real(kp), allocatable :: p_pt_wave(:,:,:,:)
  
    !*** diabatic heating ***!
    real(kp), allocatable :: q_3d(:,:,:)                ! 3D diabatic heating
    real(kp), allocatable :: q_shortwave_3d(:,:,:)      ! 3D diabatic heating by short wave radiation
    real(kp), allocatable :: q_longwave_3d(:,:,:)       ! 3D diabatic heating by long  wave radiation
    real(kp), allocatable :: q_lhr_large_3d(:,:,:)      ! 3D diabatic heating by large scale condensation
    real(kp), allocatable :: q_lhr_conv_3d(:,:,:)       ! 3D diabatic heating by convection
    real(kp), allocatable :: q_diffusion_3d(:,:,:)      ! 3D diabatic heating by vertical diffusion
    real(kp), allocatable :: q_zm(:,:)                  ! zonal mean diabatic heating
    real(kp), allocatable :: q_shortwave_zm(:,:)        ! zonal mean diabatic heating by short wave radiation
    real(kp), allocatable :: q_longwave_zm(:,:)         ! zonal mean diabatic heating by long  wave radiation
    real(kp), allocatable :: q_lhr_large_zm(:,:)        ! zonal mean diabatic heating by large scale condensation
    real(kp), allocatable :: q_lhr_conv_zm(:,:)         ! zonal mean diabatic heating by convection
    real(kp), allocatable :: q_diffusion_zm(:,:)        ! zonal mean diabatic heating by vertical diffusion
    !real(kp), allocatable :: q_ex_3d(:,:,:)
    !real(kp), allocatable :: q_ex_zm(:,:)
    real(kp), allocatable :: qgz_zm(:,:)                ! diabatic heating to the zonal mean state
    real(kp), allocatable :: qgz_shortwave_zm(:,:)      ! diabatic heating to the zonal mean state by short wave radiation
    real(kp), allocatable :: qgz_longwave_zm(:,:)       ! diabatic heating to the zonal mean state by long  wave radiation
    real(kp), allocatable :: qgz_lhr_large_zm(:,:)      ! diabatic heating to the zonal mean state by large scale condensation
    real(kp), allocatable :: qgz_lhr_conv_zm(:,:)       ! diabatic heating to the zonal mean state by convection
    real(kp), allocatable :: qgz_diffusion_zm(:,:)      ! diabatic heating to the zonal mean state by vertical diffusion
    !real(kp), allocatable :: qz_pdd(:)
    real(kp), allocatable :: qz_vint(:)
    real(kp), allocatable :: qz_shortwave_vint(:)
    real(kp), allocatable :: qz_longwave_vint(:)
    real(kp), allocatable :: qz_lhr_large_vint(:)
    real(kp), allocatable :: qz_lhr_conv_vint(:)
    real(kp), allocatable :: qz_diffusion_vint(:)
    real(kp)              :: qz_gmean(1)                ! Az generation by diabatic heating
    real(kp)              :: qz_shortwave_gmean(1)      ! Az generation by short wave radiation
    real(kp)              :: qz_longwave_gmean(1)       ! Az generation by long  wave radiation
    real(kp)              :: qz_lhr_large_gmean(1)      ! Az generation by large scale condensation
    real(kp)              :: qz_lhr_conv_gmean(1)       ! Az generation by convection
    real(kp)              :: qz_diffusion_gmean(1)      ! Az generation by vertical diffusion
    real(kp), allocatable :: qe_zm(:,:)                 ! Ae generation by diabatic heating
    real(kp), allocatable :: qe_shortwave_zm(:,:)       ! Ae generation by short wave radiation
    real(kp), allocatable :: qe_longwave_zm(:,:)        ! Ae generation by long  wave radiation 
    real(kp), allocatable :: qe_lhr_large_zm(:,:)       ! Ae generation by large scale condensation
    real(kp), allocatable :: qe_lhr_conv_zm(:,:)        ! Ae generation by convection
    real(kp), allocatable :: qe_diffusion_zm(:,:)       ! Ae generation by vertical diffusion
  
    !*** others ***!
    real(kp), allocatable :: p_pds(:)    ! p+s: zonal mean surface pressure
    real(kp)              :: p_pdds(1)   ! p++s: global mean surface pressure [hPa]
    !real(kp), allocatable :: work(:,:,:) ! temporal work space (im,jm,km)
  
  
    contains


    subroutine mim_var_ini()
      
        allocate(u(im,jm,km))
        allocate(v(im,jm,km))
        allocate(t(im,jm,km))
        allocate(z(im,jm,km))
        allocate(omega(im,jm,km))
      
        allocate(p_sfc(im,jm))
        !allocate(msl(im,jm))
        !allocate(t_sfc(im,jm))
        allocate(alt(im,jm))
      
        allocate(pt(im,jm,km))
        allocate(pt_dot(im,jm,km))
      
        allocate(p_zm(jm,ko))
        allocate(u_zm(jm,ko))
        allocate(v_zm(jm,ko))
        allocate(pt_zm(jm,ko))
        allocate(t_dagger(jm,ko))
        ! allocate(t_zm(jm,ko))
        allocate(st_zm(jm,ko))
        allocate(w_zm(jm,ko))
        allocate(z_zm(jm,ko))
        allocate(pt_dot_zm(jm,ko))
      
        allocate(u_u_zm(jm,ko))
        allocate(u_u_x_zm(jm,ko))
        allocate(v_v_zm(jm,ko))
        allocate(v_v_x_zm(jm,ko))
        allocate(u_v_zm(jm,ko))
        allocate(u_v_x_zm(jm,ko))
        allocate(u_u_v_zm(jm,ko))
        allocate(v_v_v_zm(jm,ko))
        allocate(u_pt_dot_zm(jm,ko))
        allocate(u_pt_dot_x_zm(jm,ko))
        allocate(v_pt_dot_zm(jm,ko))
        allocate(v_pt_dot_x_zm(jm,ko))
        allocate(u_u_pt_dot_zm(jm,ko))
        allocate(v_v_pt_dot_zm(jm,ko))
        !allocate(v_ke_zm(jm,ko))
        !allocate(pt_dot_ke_zm(jm,ko))
      
        allocate(epy(jm,ko))
        allocate(depy(jm,ko))
        allocate(epz_form(jm,ko))
        allocate(depz_form(jm,ko))
        allocate(epz_uv(jm,ko))
        allocate(depz_uv(jm,ko))
        allocate(epz_ut(jm,ko))
        allocate(depz_ut(jm,ko))
        allocate(epz_uw(jm,ko))
        allocate(depz_uw(jm,ko))
        allocate(epz(jm,ko))
        allocate(depz(jm,ko))
        allocate(divf(jm,ko))
        allocate(gy(jm,ko))
        allocate(dgy(jm,ko))
        allocate(gz(jm,ko))
        allocate(dgz(jm,ko))
      
        allocate(kz_zm(jm,ko))
        allocate(ke_zm(jm,ko))
        allocate(pz_zm(jm,ko))
        allocate(ae_zm_vint(jm))
        allocate(ae_total_zm(jm,ko))
        allocate(ae_s_zm(jm))
        allocate(az_zm(jm,ko))
        allocate(az_zm_vint(jm))
      
        allocate(c_az_kz(jm,ko))
        allocate(c_az_kz_modify(jm,ko))
        allocate(c_kz_ae(jm,ko))
        allocate(c_kz_ae_u(jm,ko))
        allocate(c_kz_ae_v(jm,ko))
        allocate(c_ae_ke(jm,ko))
        allocate(c_ae_ke_u(jm,ko))
        allocate(c_ae_ke_v(jm,ko))
        allocate(c_kz_ke(jm,ko))
        allocate(c_kz_ke_uy(jm,ko))
        allocate(c_kz_ke_uz(jm,ko))
        allocate(c_kz_ke_vy(jm,ko))
        allocate(c_kz_ke_vz(jm,ko))
        allocate(c_kz_ke_tan(jm,ko))
        allocate(c_kz_w(jm,ko))
      
      
        allocate(nlev(im,jm,ko))
        allocate(dlev(im,jm,ko))
        allocate(p_pd(im,jm,ko))
        allocate(x_pd(im,jm,km))
        allocate(xint_zm(jm,ko))
        allocate(pt_sfc(im,jm))
        allocate(pt_pds(jm))
        allocate(pd_p(im,jm,km))
      
        allocate(dlev_y(jm,ko))
        allocate(nlev_y(jm,ko))
        allocate(pd_pdd(jm,ko))
        allocate(pd_ym(ko))
        allocate(pt_ym(ko))
        allocate(pdd_pd(jm,ko))
      
      
      
      
        allocate(dz_dlat_zm(jm,ko))
        allocate(v_dz_dlat_zm(jm,ko))
        allocate(u_dz_dlon_zm(jm,ko))
        allocate(temp_vint(jm))
        allocate(phi_dagger(jm,ko))
        allocate(phi_dagger_y(jm,ko))
        allocate(p_dphi_dt(jm,ko))
        !allocate(p_dz_dt_zm(jm,ko))
        !allocate(p_dz_dt(im,jm,km))
        !allocate(divz_tzm(jm,ko))
        !allocate(divphi_t(jm,ko))
        !allocate(dwdt(jm,ko))
        !allocate(uuv_tmp(jm,ko))
        ! allocate(d_u_epy(jm,ko))
        !allocate(d_u_epz(jm,ko))
        ! allocate(dkedt_uvu_y(jm,ko))
      
        allocate(dkzdt_vkz(jm,ko))
        allocate(dkzdt_wkz(jm,ko))
        allocate(dkedt_uy(jm,ko))
        allocate(dkedt_vy(jm,ko))
        allocate(dkedt_uz(jm,ko))
        allocate(dkedt_vz(jm,ko))
        allocate(dkedt_vke(jm,ko))
        allocate(dkedt_wke(jm,ko))
        !allocate(dpedt_vt(jm,ko))
        !allocate(dpedt_wt(jm,ko))
      
      
        allocate(z_pd(im,jm,ko))
        allocate(pt_past(im,jm,km))
        allocate(u_past(im,jm,km))
        allocate(v_past(im,jm,km))
        allocate(omega_past(im,jm,km))
        !allocate(phi_dagger_past(jm,ko))
        allocate(z_pd_past(im,jm,ko))
        allocate(p_pd_past(im,jm,ko))
      
        if (wmax >= 1) then
            allocate(epz_wave(jm,ko,wmax))
            allocate(z_pt_wave(wmax,im,jm,ko))
            allocate(p_pt_wave(wmax,im,jm,ko))
        endif
    
        allocate(q_3d(im,jm,km))
        allocate(q_shortwave_3d(im,jm,km))
        allocate(q_longwave_3d(im,jm,km))
        allocate(q_lhr_large_3d(im,jm,km))
        allocate(q_lhr_conv_3d(im,jm,km))
        allocate(q_diffusion_3d(im,jm,km))
        allocate(q_zm(jm,ko))
        allocate(q_shortwave_zm(jm,ko))
        allocate(q_longwave_zm(jm,ko))
        allocate(q_lhr_large_zm(jm,ko))
        allocate(q_lhr_conv_zm(jm,ko))
        allocate(q_diffusion_zm(jm,ko))
        !allocate(q_ex_3d(im,jm,km))
        !allocate(q_ex_zm(jm,ko))
        allocate(qgz_zm(jm,ko))
        allocate(qgz_shortwave_zm(jm,ko))
        allocate(qgz_longwave_zm(jm,ko))
        allocate(qgz_lhr_large_zm(jm,ko))
        allocate(qgz_lhr_conv_zm(jm,ko))
        allocate(qgz_diffusion_zm(jm,ko))
        !allocate(qz_pdd(ko))
        allocate(qz_vint(jm))
        allocate(qz_shortwave_vint(jm))
        allocate(qz_longwave_vint(jm))
        allocate(qz_lhr_large_vint(jm))
        allocate(qz_lhr_conv_vint(jm))
        allocate(qz_diffusion_vint(jm))
        allocate(qe_zm(jm,ko))
        allocate(qe_shortwave_zm(jm,ko))
        allocate(qe_longwave_zm(jm,ko))
        allocate(qe_lhr_large_zm(jm,ko))
        allocate(qe_lhr_conv_zm(jm,ko))
        allocate(qe_diffusion_zm(jm,ko))
      
        allocate(p_pds(jm))
        !allocate(work(im,jm,km))
      
  
    end subroutine mim_var_ini
  
  
  
    subroutine mim_var_end()

        deallocate(u)
        deallocate(v)
        deallocate(t)
        deallocate(z)
        deallocate(omega)
  
        deallocate(p_sfc)
        !deallocate(msl)
        !deallocate(t_sfc)
        deallocate(alt)
  
        deallocate(pt)
        deallocate(pt_dot)
  
        deallocate(p_zm)
        deallocate(u_zm)
        deallocate(v_zm)
        deallocate(pt_zm)
        deallocate(t_dagger)
        ! deallocate(t_zm)
        deallocate(st_zm)
        deallocate(w_zm)
        deallocate(z_zm)
        deallocate(pt_dot_zm)
  
        deallocate(u_u_zm)
        deallocate(u_u_x_zm)
        deallocate(v_v_zm)
        deallocate(v_v_x_zm)
        deallocate(u_v_zm)
        deallocate(u_v_x_zm)
        deallocate(u_u_v_zm)
        deallocate(v_v_v_zm)
        deallocate(u_pt_dot_zm)
        deallocate(u_pt_dot_x_zm)
        deallocate(v_pt_dot_zm)
        deallocate(v_pt_dot_x_zm)
        deallocate(u_u_pt_dot_zm)
        deallocate(v_v_pt_dot_zm)
        !deallocate(v_ke_zm)
        !deallocate(pt_dot_ke_zm)
  
        deallocate(epy)
        deallocate(depy)
        deallocate(epz_form)
        deallocate(depz_form)
        deallocate(epz_uv)
        deallocate(depz_uv)
        deallocate(epz_ut)
        deallocate(depz_ut)
        deallocate(epz_uw)
        deallocate(depz_uw)
        deallocate(epz)
        deallocate(depz)
        deallocate(divf)
        deallocate(gy)
        deallocate(dgy)
        deallocate(gz)
        deallocate(dgz)
  
        deallocate(kz_zm)
        deallocate(ke_zm)
        deallocate(pz_zm)
        deallocate(ae_zm_vint)
        deallocate(ae_total_zm)
        deallocate(ae_s_zm)
        deallocate(az_zm)
        deallocate(az_zm_vint)
  
        deallocate(c_az_kz)
        deallocate(c_az_kz_modify)
        deallocate(c_kz_ae)
        deallocate(c_kz_ae_u)
        deallocate(c_kz_ae_v)
        deallocate(c_ae_ke)
        deallocate(c_ae_ke_u)
        deallocate(c_ae_ke_v)
        deallocate(c_kz_ke)
        deallocate(c_kz_ke_uy)
        deallocate(c_kz_ke_uz)
        deallocate(c_kz_ke_vy)
        deallocate(c_kz_ke_vz)
        deallocate(c_kz_ke_tan)
        deallocate(c_kz_w)
  
  
        deallocate(nlev)
        deallocate(dlev)
        deallocate(p_pd)
        deallocate(x_pd)
        deallocate(xint_zm)
        deallocate(pt_sfc)
        deallocate(pt_pds)
        deallocate(pd_p)
  
        deallocate(dlev_y)
        deallocate(nlev_y)
        deallocate(pd_pdd)
        deallocate(pd_ym)
        deallocate(pt_ym)
        deallocate(pdd_pd)
  
  
        deallocate(dz_dlat_zm)
        deallocate(v_dz_dlat_zm)
        deallocate(u_dz_dlon_zm)
        deallocate(temp_vint)
        deallocate(phi_dagger)
        deallocate(phi_dagger_y)
        deallocate(p_dphi_dt)
        !deallocate(p_dz_dt_zm)
        !deallocate(p_dz_dt)
        !deallocate(divz_tzm)
        !deallocate(divphi_t)
        !deallocate(dwdt)
        !deallocate(uuv_tmp)
        ! deallocate(d_u_epy)
        !deallocate(d_u_epz)
        ! deallocate(dkedt_uvu_y)
  
        deallocate(dkzdt_vkz)
        deallocate(dkzdt_wkz)
        deallocate(dkedt_uy)
        deallocate(dkedt_vy)
        deallocate(dkedt_uz)
        deallocate(dkedt_vz)
        deallocate(dkedt_vke)
        deallocate(dkedt_wke)
        !deallocate(dpedt_vt)
        !deallocate(dpedt_wt)
  
  
        deallocate(z_pd)
        deallocate(pt_past)
        deallocate(u_past)
        deallocate(v_past)
        deallocate(omega_past)
        !deallocate(phi_dagger_past)
        deallocate(z_pd_past)
        deallocate(p_pd_past)
  
        if (allocated(epz_wave)) then 
            deallocate(epz_wave)
        endif
        if (allocated(z_pt_wave)) then
            deallocate(z_pt_wave)
        endif
        if (allocated(p_pt_wave)) then
            deallocate(p_pt_wave)
        endif
  
        deallocate(q_3d)
        deallocate(q_shortwave_3d)
        deallocate(q_longwave_3d)
        deallocate(q_lhr_large_3d)
        deallocate(q_lhr_conv_3d)
        deallocate(q_diffusion_3d)
        deallocate(q_zm)
        deallocate(q_shortwave_zm)
        deallocate(q_longwave_zm)
        deallocate(q_lhr_large_zm)
        deallocate(q_lhr_conv_zm)
        deallocate(q_diffusion_zm)
        !deallocate(q_ex_3d)
        !deallocate(q_ex_zm)
        deallocate(qgz_zm)
        deallocate(qgz_shortwave_zm)
        deallocate(qgz_longwave_zm)
        deallocate(qgz_lhr_large_zm)
        deallocate(qgz_lhr_conv_zm)
        deallocate(qgz_diffusion_zm)
        !deallocate(qz_pdd)
        deallocate(qz_vint)
        deallocate(qz_shortwave_vint)
        deallocate(qz_longwave_vint)
        deallocate(qz_lhr_large_vint)
        deallocate(qz_lhr_conv_vint)
        deallocate(qz_diffusion_vint)
        deallocate(qe_zm)
        deallocate(qe_shortwave_zm)
        deallocate(qe_longwave_zm)
        deallocate(qe_lhr_large_zm)
        deallocate(qe_lhr_conv_zm)
        deallocate(qe_diffusion_zm)
  
        deallocate(p_pds)
        !deallocate(work)
  
  
    end subroutine mim_var_end


end module mim_var

