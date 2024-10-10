# MIM v-2024

Developed by Kosei Ohara  

MIM v-2024 is a modernized version of MIM.
This program is based on [MIM v0.35](https://github.com/mim-proj/mim)

## Namelist

- INPUT
    - INPUT\_UVT\_FILENAME  
      Optional  
      File name of the u (zonal wind [$`\mathrm{m\:s^{-1}}`$]), v (meridional wind [$`\mathrm{m\:s^{-1}}`$]), and t (temperature [$`\mathrm{K}`$]) on the pressure surfaces. Their order must be u &rarr; v &rarr; t.
    - INPUT\_U\_FILENAME  
      Optional  
      File name of the zonal wind [$`\mathrm{m\:s^{-1}}`$] on the pressure surfaces. This file is needed if INPUT\_UVT\_FILENAME is not specified. If the data is from JRA-3Q, anl_p125_ugrd.
    - INPUT\_V\_FILENAME  
      Optional  
      File name of the meridional wind [$`\mathrm{m\:s^{-1}}`$] on the pressure surfaces. This file is needed if INPUT\_UVT\_FILENAME is not specified. If the data is from JRA-3Q, anl_p125_vgrd.
    - INPUT\_T\_FILENAME  
      Optional  
      File name of the temperature [$`\mathrm{K}`$] on the pressure surfaces. This file is needed if INPUT\_UVT\_FILENAME is not specified. If the data is from JRA-3Q, anl_p125_tmp.
    - INPUT\_PS\_FILENAME  
      Optional  
      File name of the surface pressure [$`\mathrm{Pa}`$] / [$`\mathrm{hPa}`$]. If the data is from JRA-3Q, anl_surf125-PRES.
    - INPUT\_MSL\_FILENAME  
      Optional  
      File name of the mean sea level pressure [$`\mathrm{Pa}`$] / [$`\mathrm{hPa}`$]. This file is needed if INPUT\_PS\_FILENAME is not specified. If the data is from JRA-3Q, anl_surf125-PRMSL.
    - INPUT\_TS\_FILENAME  
      Optional  
      File name of the surface temperature [$`\mathrm{K}`$]. This file is needed if INPUT\_PS\_FILENAME is not specified. If the data is from JRA-3Q, anl_surf125-TMP.
    - INPUT\_Z\_FILENAME  
      Required  
      File name of the height [$`\mathrm{m}`$] or geopotential height [$`\mathrm{m^2\:s^{-2}}`$] on the pressure surfaces. If the data is from JRA-3Q, anl_p125_hgt.
    - INPUT\_OMEGA\_FILENAME  
      Optional  
      File name of the vertical velocity [$`\mathrm{Pa \: s^{-1}}`$]. If not specified, omega will be estimated from the continuity equation. If the data is from JRA-3Q, anl_p125_vvel.
    - INPUT\_TOPO\_FILENAME  
      Required  
      File name of the topography [$`\mathrm{m}`$] / [$`\mathrm{m^2 \: s^{-1}}`$]. If the data is from JRA-3Q, LL125_surf.
    - INPUT\_Q\_FILENAME  
      Optional  
      File name of the total diabatic heating [$`\mathrm{K \: s^{-1}}`$] / [$`\mathrm{K \: day^{-1}}`$]. If not specified, total diabatic heating is defied by the sum of diabatic heating by the short wave, long wave, large scale condensation, convective heating, and vertical diffusion. If one of them is not specified too, the diabatic heating will be estimated from the time derivative of the potential temperature.
    - INPUT\_SHORTWAVE\_FILENAME  
      Optional  
      File name of the diabatic heating by the short wave radiation [$`\mathrm{K \: s^{-1}}`$] / [$`\mathrm{K \: day^{-1}}`$]. If not specified, outputs related to the short wave radiation will be zero. If the data is from JRA-3Q, fcst_phyp125_ttswr.
    - INPUT\_LONGWAVE\_FILENAME  
      Optional  
      File name of the diabatic heating by the long wave radiation [$`\mathrm{K \: s^{-1}}`$] / [$`\mathrm{K \: day^{-1}}`$]. If not specified, outputs related to the long wave radiation will be zero. If the data is from JRA-3Q, fcst_phyp125_ttlwr.
    - INPUT\_LHR\_LARGE\_FILENAME  
      Optional  
      File name of the diabatic heating by the large scale condensation [$`\mathrm{K \: s^{-1}}`$] / [$`\mathrm{K \: day^{-1}}`$]. If not specified, outputs related to the large scale condensation will be zero. If the data is from JRA-3Q, fcst_phyp125_lrghr.
    - INPUT\_LHR\_CONV\_FILENAME  
      Optional  
      File name of the diabatic heating by the convective heating [$`\mathrm{K \: s^{-1}}`$] / [$`\mathrm{K \: day^{-1}}`$]. If not specified, outputs related to the convective heating will be zero. If the data is from JRA-3Q, fcst_phyp125_cnvhr.
    - INPUT\_DIFFUSION\_FILENAME  
      Optional  
      File name of the diabatic heating by the vertical diffusion [$`\mathrm{K \: s^{-1}}`$] / [$`\mathrm{K \: day^{-1}}`$]. If not specified, outputs related to the vertical diffusion will be zero. If the data is from JRA-3Q, fcst_phyp125_vdfhr.

- INPUT\_UNIT
    - INPUT\_UNIT\_Z  
      Optional (default : "m")  
      Unit of the data specified in INPUT\_Z\_FILENAME. "m" if the data unit is [$`\mathrm{m}`$] and "m^2/s^2" if [$`\mathrm{m^2 \: s^{-2}}`$]
    - INPUT\_UNIT\_PS  
      Optional (default : "hPa")  
      Unit of the data specified in INPUT\_PS\_FILENAME. "hPa" or "Pa".
    - INPUT\_UNIT\_MSL  
      Optional (default : "hPa")  
      Unit of the data specified in INPUT\_MSL\_FILENAME. "hPa" or "Pa".
    - INPUT\_UNIT\_TOPO  
      Optional (default : "m")  
      Unit of the data specified in INPUT\_TOPO\_FILENAME. "m" if the data unit is [$`\mathrm{m}`$] and "m^2/s^2" if [$`\mathrm{m^2 \: s^{-2}}`$].
    - INPUT\_UNIT\_Q  
      Optional (defalt : "K/s")  
      Unit of the data specified in INPUT\_Q\_FILENAME and the other datasets related to the diabatic heating. "K/s" if the data unit is [$`\mathrm{K \: s^{-1}}`$] and "K/day" if [$`\mathrm{K \: day^{-1}}`$].

- INPUT\_UNDEF  
    The value for the undefined grids. The default value is 9.999E+20. If INPUT\_UNDEF\_DEFAULT is specified by another value,  the undefined values in all files will be overwritten. However, if an undefined value has been individually specified for a particular dataset with a value other than 9.999E+20, that dataset will retain its unique undefined value.
    - INPUT\_UNDEF\_DEFAULT
    - INPUT\_UNDEF\_UVT
    - INPUT\_UNDEF\_U
    - INPUT\_UNDEF\_V
    - INPUT\_UNDEF\_T
    - INPUT\_UNDEF\_PS
    - INPUT\_UNDEF\_MSL
    - INPUT\_UNDEF\_TS
    - INPUT\_UNDEF\_Z
    - INPUT\_UNDEF\_OMEGA
    - INPUT\_UNDEF\_Q
    - INPUT\_UNDEF\_SHORTWAVE
    - INPUT\_UNDEF\_LONGWAVE
    - INPUT\_UNDEF\_LHR\_LARGE
    - INPUT\_UNDEF\_LHR\_CONV
    - INPUT\_UNDEF\_DIFFUSION

- INPUT\_ENDIAN  
    The endian for the input files. "native", "little_endian", and "big_endian" are valid. if "native", the endian is determined by the compiler option or the environment. If the endian for each file is not specified individually, it follows INPUT\_ENDIAN\_DEFAULT.
    - INPUT\_ENDIAN\_DEFAULT
    - INPUT\_ENDIAN\_UVT
    - INPUT\_ENDIAN\_U
    - INPUT\_ENDIAN\_V
    - INPUT\_ENDIAN\_T
    - INPUT\_ENDIAN\_PS
    - INPUT\_ENDIAN\_MSL
    - INPUT\_ENDIAN\_TS
    - INPUT\_ENDIAN\_Z
    - INPUT\_ENDIAN\_OMEGA
    - INPUT\_ENDIAN\_TOPO
    - INPUT\_ENDIAN\_Q
    - INPUT\_ENDIAN\_SHORTWAVE
    - INPUT\_ENDIAN\_LONGWAVE
    - INPUT\_ENDIAN\_LHR\_LARGE
    - INPUT\_ENDIAN\_LHR\_CONV
    - INPUT\_ENDIAN\_DIFFUSION

- INPUT\_XDEF
    - INPUT\_XDEF\_NUM  
      Required  
      Number of grids in x (longitudinal) direction.

- INPUT\_YDEF
    - INPUT\_YDEF\_TYPE  
      Required  
      "lat\_degree", "lat\_radian", and "linear" are valid. If specified by "lat\_degree", the latitude of each grid is determined by INPUT\_YDEF\_LEVEL and its unit is [degree]. If "lat\_radian", same as "lat\_degree", but its unit is [radian]. If "linear", the latitude is determined by INPUT\_YDEF\_NORTH, INPUT\_YDEF\_SOUTH, and INPUT\_YDEF\_NUM linearly.
    - INPUT\_YDEF\_NUM  
      Required  
      Number of grids in y (latitudinal) direction.
    - INPUT\_TDEF\_LEVEL  
      Optional  
      Latitudes of the grids. Both north &rarr; south and south &rarr; north orders are valid. Latitudes are determined only if INPUT\_YDEF\_TYPE is "lat\_degree" or "lat\_radian".
    - INPUT\_YDEF\_SOUTH  
      Optional  
      The latitude of the southern edge [degree]. This value is needed only if INPUT\_YDEF\_TYPE="linear".
    - INPUT\_YDEF\_NORTH  
      Optional  
      The latitude of the northern edge [degree]. This value is needed only if INPUT\_YDEF\_TYPE="linear".
    - INPUT\_YDEF\_YREV\_DEFAULT  
      Optional (default : .False.)  
      .True. and .False. are valid. If the input files are north &rarr; south (yrev), specify .True. If south &rarr; north, specify .False.
    - INPUT\_YDEF\_YREV\_TOPO  
      Optional (default : .False.)  
      Same as INPUT\_YDF\_YREV\_DEFAULT, but this option is for the topography file.

- INPUT\_ZDEF  
    - INPUT\_ZDEF\_NUM  
      Required  
      Number of levels.
    - INPUT\_ZDEF\_LEVEL  
      Required  
      Pressure of each surface [hPa]. Both upper &rarr; lower and lower &rarr; upper are valid.
    - INPUT\_ZDEF\_ZREV  
      Optional (default : .False.)  
      .True. and .False. are valid. If the input files are upper &rarr; lower (zrev), specify .True. If lower &rarr; upper, specify .False.
      
- INPUT\_TDEF
    - INPUT\_TDEF\_TYPE  
      Required  
      "tstep", "monthly", and "annual" are valid. If specified by "tstep", number of time steps is defined by INPUT\_TDEF\_TSTEP. If "monthly", number of time steps is defined by INPUT\_TDEF\_YEAR, INPUT\_TDEF\_MONTH, INPUT\_TDEF\_365DAY, and INPUT\_TDEF\_DAYNUM. "monthly" options is used to compute only for one month. If "annual", number of time steps is defined by INPUT\_TDEF\_YEAR, INPUT\_TDEF\_365DAY, and INPUT\_TDEF\_DAYNUM. "annual" options is used to compute only for one year.
    - INPUT\_TDEF\_DAYNUM  
      Requied  
      Number of time steps in each day.
    - INPUT\_TDEF\_TSTEP  
      Optional  
      Total number of time steps for the period of data analysis. Valid if INPUT\_TDEF\_TYPE="tstep".
    - INPUT\_TDEF\_YEAR  
      Optional  
      The year of the input data. Valid if INPUT\_TDEF\_TYPE is "monthly" or "annual".
    - INPUT\_TDEF\_MONTH  
      Optional  
      The month of the input data. Valid if INPUT\_TDEF\_TYPE="monthly".
    - INPUT\_TDEF\_365DAY  
      Optional  
      0 if leap year is assumed, 1 if not assumed. Valid if INPUT\_TDEF\_TYPE is "monthly" or "annual".

- WAVE
    - WAVE\_MAX\_NUMBER  
      Optional  
      Maxumum wavenumber for the spectral expansion of the form drag. Valid if OUTPUT\_WAVE\_FILENAME is specified.

- OUTPUT
    - OUTPUT\_ZONAL\_FILENAME  
      Optional  
      File name of the output zonal file. If specified, latitude-height crosssection zonal mean field will be output.
    - OUTPUT\_VINT\_FILENAME  
      Optional  
      File name of the output vint file. If specified, vertically integrated parameters will be output.
    - OUTPUT\_GMEAN\_FILENAME  
      Optional  
      File name of the output gmean file. If specified, global mean parameters will be output.
    - OUTPUT\_WAVE\_FILENAME  
      Optional  
      File name of the output wave file. If specified, the Form-Drag expanded in wavenumber will be output.
    - OUTPUT\_ERROR\_FILENAME  
      Required  
      File name of the output error file. Error of parameters will be output.
    - OUTPUT\_WARN\_FILENAME  
      Required  
      File name of the output warning file. Various warning will be output.

- OUTPUT\_ZDEF
    - OUTPUT\_ZDEF\_NUM  
      Required  
      Same as INPUT\_ZDEF\_NUM, but for the output. This value does not necessarily have to be equal to INPUT\_ZDEF\_NUM, however, specifying by the same value is strongly recommended to avoid crucial errors.
    - OUTPUT\_ZDEF\_LEVEL  
      Required  
      Same as INPUT\_ZDEF\_LEVEL, but for the output. This levels does not necessarily have to be equal to INPUT\_ZDEF\_LEVEL, however, specifying by the same value is strongly recommended to avoid crucial errors.



## Other Settings
Most settings should be applied in the namelists.
To optimize your requirement, edit the source codes.
- Numerical Precision  
To change the numerical precision, edit `src/params.f90`. `kp` is the kind parameter for real variables: `kp=4` for single precison, `kp=8` for double precision, and `kp=16` for quadruple precision. Similarly, `ckp` is the kind parameter for complex variables. `rkp` and `wkp` are kind parameters for input parameters and output parameters, respectively. For example, to read 8 byte real data and compute in double precision, set `rkp=8`, `kp=8`, and `ckp=8`. Execute `make re` to recomplie all source code.
- Output Parameters  
In the default source, the program outputs many parameters. To limit the outputs, edit `io_main.f90`. `write_zonal()`, `write_vint()`, `write_gmean()`, and `write_wave()` are subroutines for output.



## Output Parameters
Depending on the setting in the namelist, the program can generate 4 files. 4 are binary files, and the other is ASCII file. The endian of the binary files depends on your environment.

### ZONAL
Latitude-pressure crosssection data. The latitude of each grid is determined by the input files. The levels depends on OUTPUT\_ZDEF\_LEVEL. The data are yrev. 68 parameters are output.

- u  
Zonal mean zonal wind [$`\mathrm{m \: s^{-1}}`$].
- v  
Zonal mean meridional wind [$`\mathrm{m \: s^{-1}}`$].
- pt  
Potential temperature [$`\mathrm{K}`$].
- t  
Temperature in the zonal mean state [$`\mathrm{K}`$].
- st  
Mass streamfunction [$`\mathrm{kg \: s^{-1}}`$].
- w  
Zonal mean vetical wind [$`\mathrm{m \: s^{-1}}`$]. This parameter is computed from the mass streamfunction.
- z  
Zonal mean geopotential height [$`\mathrm{m}`$].
- epy  
Meridional component of the EP flux [$`\mathrm{kg \: s^{-2}}`$].
- depy  
Meridional divergence of the EP flux [$`\mathrm{m \: s^{-2}}`$].
- epz\_form  
Form Drag (One term of the vertical component of the EP flux) [$`\mathrm{kg \: s^{-2}}`$].
- depz\_form  
Vertical divergence of the Form Drag [$`\mathrm{m \: s^{-2}}`$].
- epz\_uv  
uv term of the vertical component of the EP flux [$`\mathrm{kg \: s^{-2}}`$].
- depz\_uv  
Vertical divergence of epz\_uv [$`\mathrm{m \: s^{-2}}`$].
- epz\_ut  
ut term of the vertical component of the EP flux [$`\mathrm{kg \: s^{-2}}`$].
- depz\_ut  
Vertical divergence of epz\_ut [$`\mathrm{m \: s^{-2}}`$].
- epz  
Vertical component of the EP flux (sum of epz\_form, epz\_uv, and epz\_ut).
- depz  
Vertical divergence of the EP flux.
- divf  
Divergence of the EP flux.
- gy  
Meridional component of the G flux [$`\mathrm{kg \: s^{-2}}`$].
- dgy  
Meridional divergence of the G flux [$`\mathrm{m \: s^{-2}}`$].
- gz  
Vertical component of the G flux [$`\mathrm{kg \: s^{-2}}`$].
- dgz  
Vertical divergence of the G flux [$`\mathrm{m \: s^{-2}}`$].
- uux  
Zonal mean of the square of the eddy zonal wind [$`\mathrm{m^2 \: s^{-2}}`$].
- c\_az\_kz  
Energy conversion from the zonal available potential energy to the zonal mean kinetic energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ae\_u  
Energy conversion from the zonal mean kinetic energy to the eddy available potential energy by the zonal wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ae\_v  
Energy conversion from the zonal mean kinetic energy to the eddy available potential energy by the meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ae  
Total energy conversion from the zonal mean kinetic energy to the eddy available potential energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_ae\_ke\_u  
Energy conversion from the eddy available potential energy to the eddy kinetic energy by the eddy zonal wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_ae\_ke\_v  
Energy conversion from the eddy available potential energy to the eddy kinetic energy by the eddy meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_ae\_ke  
Total energy conversion from the eddy available potential energy to the eddy kinetic energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_uy  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean zonal wind and the meridional divergence of the EP flux [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_uz  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean zonal wind and the depz\_uw [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_vy  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean meridional wind and the meridional divergence of the G flux [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_vz  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean meridional wind and the vertical divergence of the G flux [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_tan  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean meridional wind and the square of the eddy zonal wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke  
Total energy conversion from the zonal mean kinetic energy to the eddy kinetic energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_w  
Energy conversion from the zonal mean kinetic energy to the wave energy (sum of c\_kz\_ke and c\_kz\_ae) [$`\mathrm{m^2 \: s^{-3}}`$].
- q  
Zonal mean diabatic heating [$`\mathrm{m^2 \: s^{-3}}`$].
- shortwave  
Zonal mean diabatic heating by the short wave [$`\mathrm{m^2 \: s^{-3}}`$].
- longwave  
Zonal mean diabatic heating by the long wave [$`\mathrm{m^2 \: s^{-3}}`$].
- lhr\_large  
Zonal mean diabatic heating by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- lhr\_conv  
Zonal mean diabatic heating by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- diffusion  
Zonal mean diabatic heating by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz  
Generation rate to the zonal mean state [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_shortwave  
Generation rate to the zonal mean state by the short wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_longwave  
Generation rate to the zonal mean state by the long wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_lhr\_large  
Generation rate to the zonal mean state by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_lhr\_conv  
Generation rate to the zonal mean state by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_diffusion  
Generation rate to the zonal mean state by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- qe  
Eddy generation rate [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_shortwave  
Eddy generation rate by the short wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_longwave  
Eddy generation rate by the long wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_lhr\_large  
Eddy generation rate by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_lhr\_conv  
Eddy generation rate by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_diffusion  
Eddy generation rate by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- kz  
Zonal mean kinetic energy [$`\mathrm{m^2 \: s^{-2}}`$].
- ke  
Eddy kinetic energy [$`\mathrm{m^2 \: s^{-2}}`$].
- pz  
Zonal potential energy (NOT the available potential energy) [$`\mathrm{m^2 \: s^{-2}}`$].
- ae  
Eddy available potential energy [$`\mathrm{m^2 \: s^{-2}}`$].
- dkzdt\_vkz  
Advection of the zonal mean kinetic energy by the meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- dkzdt\_wkz  
Advection of the zonal mean kinetic energy by the vertical wind [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_uy  
Divergence of the wave energy flux d(u Fy)/dy [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_vy  
Divergence of the wave energy flux d(v Gy)/dy [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_uz  
Divergence of the wave energy flux d(u Fz\_uw)/dz [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_vz  
Divergence of the wave energy flux d(v Gz)/dz [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_vke  
Advection of the eddy kinetic energy by the meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_wke  
Advection of the eddy kinetic energy by the vertical wind [$`\mathrm{m^2 \: s^{-3}}`$].

### VINT
Virticaly integrated data (NOT the vertical mean). The latitude of each grid is determined by the input files. Grids under the ground are not used. The data are yrev. 50 parameters are output.

- kz  
Zonal mean kinetic energy [$`\mathrm{m^2 \: s^{-2}}`$].
- ke  
Eddy kinetic energy [$`\mathrm{m^2 \: s^{-2}}`$].
- az  
Zonal available potential energy [$`\mathrm{m^2 \: s^{-2}}`$].
- ae  
Eddy available potential energy [$`\mathrm{m^2 \: s^{-2}}`$].
- c\_az\_kz  
Energy conversion from the zonal available potential energy to the zonal mean kinetic energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ae\_u  
Energy conversion from the zonal mean kinetic energy to the eddy available potential energy by the zonal wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ae\_v  
Energy conversion from the zonal mean kinetic energy to the eddy available potential energy by the meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ae  
Total energy conversion from the zonal mean kinetic energy to the eddy available potential energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_ae\_ke\_u  
Energy conversion from the eddy available potential energy to the eddy kinetic energy by the eddy zonal wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_ae\_ke\_v  
Energy conversion from the eddy available potential energy to the eddy kinetic energy by the eddy meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_ae\_ke  
Total energy conversion from the eddy available potential energy to the eddy kinetic energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_uy  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean zonal wind and the meridional divergence of the EP flux [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_uz  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean zonal wind and the depz\_uw [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_vy  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean meridional wind and the meridional divergence of the G flux [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_vz  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean meridional wind and the vertical divergence of the G flux [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke\_tan  
Energy conversion from the zonal mean kinetic energy to the eddy kinetic energy by the zonal mean meridional wind and the square of the eddy zonal wind [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_ke  
Total energy conversion from the zonal mean kinetic energy to the eddy kinetic energy [$`\mathrm{m^2 \: s^{-3}}`$].
- c\_kz\_w  
Energy conversion from the zonal mean kinetic energy to the wave energy (sum of c\_kz\_ke and c\_kz\_ae) [$`\mathrm{m^2 \: s^{-3}}`$].
- q  
Zonal mean diabatic heating [$`\mathrm{m^2 \: s^{-3}}`$].
- shortwave  
Zonal mean diabatic heating by the short wave [$`\mathrm{m^2 \: s^{-3}}`$].
- longwave  
Zonal mean diabatic heating by the long wave [$`\mathrm{m^2 \: s^{-3}}`$].
- lhr\_large  
Zonal mean diabatic heating by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- lhr\_conv  
Zonal mean diabatic heating by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- diffusion  
Zonal mean diabatic heating by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz  
Generation rate to the zonal mean state [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_shortwave  
Generation rate to the zonal mean state by the short wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_longwave  
Generation rate to the zonal mean state by the long wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_lhr\_large  
Generation rate to the zonal mean state by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_lhr\_conv  
Generation rate to the zonal mean state by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- qgz\_diffusion  
Generation rate to the zonal mean state by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- qe  
Eddy generation rate [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_shortwave  
Eddy generation rate by the short wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_longwave  
Eddy generation rate by the long wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_lhr\_large  
Eddy generation rate by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_lhr\_conv  
Eddy generation rate by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- qe\_diffusion  
Eddy generation rate by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- qz  
Zonal generation rate [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_shortwave  
Zonal generation rate by the short wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_longwave  
Zonal generation rate by the long wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_lhr\_large  
Zonal generation rate by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_lhr\_conv  
Zonal generation rate by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_diffusion  
Zonal generation rate by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].
- dkzdt\_vkz  
Advection of the zonal mean kinetic energy by the meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- dkzdt\_wkz  
Advection of the zonal mean kinetic energy by the vertical wind [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_uy  
Divergence of the wave energy flux d(u Fy)/dy [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_vy  
Divergence of the wave energy flux d(v Gy)/dy [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_uz  
Divergence of the wave energy flux d(u Fz\_uw)/dz [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_vz  
Divergence of the wave energy flux d(v Gz)/dz [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_vke  
Advection of the eddy kinetic energy by the meridional wind [$`\mathrm{m^2 \: s^{-3}}`$].
- dkedt\_wke  
Advection of the eddy kinetic energy by the vertical wind [$`\mathrm{m^2 \: s^{-3}}`$].

### GMEAN
- az  
Zonal available potential energy [$`\mathrm{m^2 \: s^{-2}}`$].
- qz  
Zonal generation rate [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_shortwave  
Zonal generation rate by the short wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_longwave  
Zonal generation rate by the long wave radiation [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_lhr\_large  
Zonal generation rate by the large scale condensation [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_lhr\_conv  
Zonal generation rate by the convective heating [$`\mathrm{m^2 \: s^{-3}}`$].
- qz\_diffusion  
Zonal generation rate by the vertical diffusion [$`\mathrm{m^2 \: s^{-3}}`$].





