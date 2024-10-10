dset ^JRA3Q_ZONAL.grd
title MIM

undef 9.999e+20

options little_endian yrev

xdef 1 linear 1 1
ydef 145 linear -90.0 1.25
zdef 45 levels
1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 85, 70, 60, 50, 40, 30, 20, 10, 7, 5, 3, 2, 1, 0.7, 0.3, 0.1, 0.03, 0.01
tdef 1464 linear 00Z01JAN2000 6hr

vars 67
u           45 99  Zonal Wind [m/s]
v           45 99  Meridional Wind [m/s]
pt          45 99  Potential Temperature [K]
t           45 99  Temperature Dagger derived from pt [K]
st          45 99  Mass Streamfunction [kg/s]
w           45 99  Vertical Velocity derived from st [m/s]
z           45 99  Geopotential Height [m]
epy         45 99  Meridional Component of EP Flux [kg/s^2]
depy        45 99  EP Flux Divergence dut to epy [m/s^2]
epz_form    45 99  Vertical Component of EP Flux (Form drag) [kg/s^2]
depz_form   45 99  EP Flux Divergence dut to epz_form [m/s^2]
epz_uv      45 99  Vertical Component of EP Flux (u'v') [kg/s^2]
depz_uv     45 99  EP Flux Divergence dut to epz_uv [m/s^2]
epz_ut      45 99  Vertical Component of EP Flux (u'T') [kg/s^2]
depz_ut     45 99  EP Flux Divergence dut to epz_ut [m/s^2]
epz         45 99  Vertical Component of EP Flux (Total) [kg/s^2]
depz        45 99  EP Flux Divergence dut to epz [m/s^2]
divf        45 99  EP Flux Divergence [m/s^2]
gy          45 99  Meridional Comp. of G (Meridional Momentum) Flux [kg/s^2]
dgy         45 99  G Flux Divergence due to gy [m/s^2]
gz          45 99  Vertical Component of G (Meridional Momentum) Flux [kg/s^2]
dgz         45 99  G Flux Divergence due to gz [m/s^2]
uux         45 99  Zonal mean (u'u') [m^2/s^2]
c_az_kz     45 99  C(Az->Kz) [J/(kg s)] or [W/kg]
c_kz_ae_u   45 99  C(Kz->Ae) by U (i.e. form drag) [J/(kg s)] or [W/kg]
c_kz_ae_v   45 99  C(Kz->Ae) by V [J/(kg s)] or [W/kg]
c_kz_ae     45 99  C(Kz->Ae) [J/(kg s)] or [W/kg]
c_ae_ke_u   45 99  C(Ae->Ke) by U [J/(kg s)] or [W/kg]
c_ae_ke_v   45 99  C(Ae->Ke) by V [J/(kg s)] or [W/kg]
c_ae_ke     45 99  C(Ae->Ke) [J/(kg s)] or [W/kg]
c_kz_ke_uy  45 99  C(Kz->Ke) by u * DFy [J/(kg s)] or [W/kg]
c_kz_ke_uz  45 99  C(Kz->Ke) by u * DFz^uw [J/(kg s)] or [W/kg]
c_kz_ke_vy  45 99  C(Kz->Ke) by v * DGy [J/(kg s)] or [W/kg]
c_kz_ke_vz  45 99  C(Kz->Ke) by v * DGz [J/(kg s)] or [W/kg]
c_kz_ke_tan 45 99  C(Kz->Ke) by tan [J/(kg s)] or [W/kg]
c_kz_ke     45 99  C(Kz->Ke) [J/(kg s)] or [W/kg]
c_kz_w      45 99  C(Kz->W) = C(Kz->Ke) + C(Kz->Ae) [J/(kg s)] or [W/kg]
q           45 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)]
ttswr       45 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by short wave radiation
ttlwr       45 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by long wave radiation
lrghr       45 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by large scale condensation
cnvhr       45 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by convection
vdfhr       45 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by vertical diffusion
qgz         45 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)]
ttswr_gz    45 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by short wave radiation
ttlwr_gz    45 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by long wave radiation
lrghr_gz    45 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by large scale condensation
cnvhr_gz    45 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by convection
vdfhr_gz    45 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by vertical diffusion
qe          45 99  Eddy Available Diabatic Heating [J/(kg s)]
ttswr_qe    45 99  Eddy Available Diabatic Heating [J/(kg s)] by short wave radiation
ttlwr_qe    45 99  Eddy Available Diabatic Heating [J/(kg s)] by long wave radiation
lrghr_qe    45 99  Eddy Available Diabatic Heating [J/(kg s)] by large scale condensation
cnvhr_qe    45 99  Eddy Available Diabatic Heating [J/(kg s)] by convection
vdfhr_qe    45 99  Eddy Available Diabatic Heating [J/(kg s)] by vertical diffusion
kz          45 99  Zonal Kinetic Energy [m^2/s^2] or [J/kg]
ke          45 99  Eddy Kinetic Energy [m^2/s^2] or [J/kg]
pz          45 99  Potential Energy (including groud state) [m^2/s^2] or [J/kg]
ae_total    45 99  Eddy Available Potential Energy (including surface effect)
dkzdt_vkz   45 99  Kz advection by v [m/s^3]
dkzdt_wkz   45 99  Kz advection by w [m/s^3]
dkedt_uy    45 99  Wave Energy Flux Div. d(u Fy)/dy [m/s^3]
dkedt_vy    45 99  Wave Energy Flux Div. d(v Gy)/dy [m/s^3]
dkedt_uz    45 99  Wave Energy Flux Div. d(u Fz^uw)/dz [m/s^3]
dkedt_vz    45 99  Wave Energy Flux Div. d(v Gz)/dz [m/s^3]
dkedt_vke   45 99  Ke advection by v [m/s^3]
dkedt_wke   45 99  Ke advection by w [m/s^3]
ENDVARS
