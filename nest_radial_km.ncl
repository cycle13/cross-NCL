;; -------------------------------------------------------------
;
; Script to calculate azimuthally-averaged radial velocity
; from Hurrican Rita WRF model run
;	Fix: time-step in timeloop, pressure levels
;
;; -------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


begin

; --- Read in file

path = "/nfs/a37/mm16lbb/data1/Lisa_thesis/models/"    ; where wrf output file lives
;a = addfile(path+"wrfout_d01_2005-09-20_00:00:00","r")
a = addfile(path+"wrfout_d02_2005-09-20_12:00:00","r")

;; ================================
;;   Variables
;; ================================

; --- What times and how many time steps are in the data set?

  time_0  = wrf_times_c(a->Times,0)
  time_1  = wrf_times_c(a->Times,1)

  ntimes = dimsizes(time_0)         ; number of times in the file

init_time = time_1(0)
init_time = (/ init_time - time_0(0) /)  ; incase of any time offset

; --- Lat/Lon in 1D and 2D

lat	= a->XLAT(0,:,:)   ; 2D lat/long
lon	= a->XLONG(0,:,:)  ; dims: [south_north | 232] x [west_east | 322] 

lat1	= lat(:,0)	   ; dims: [south_north | 232]
lon1	= lon(0,:)	   ; dims: [west_east | 322]

	;printVarSummary(lat)
	;printVarSummary(lon)
	;printVarSummary(lat1)
	;printVarSummary(lon1)

; --- 27 Actual pressure levels

  p_in	= wrf_user_getvar(a,"pressure",-1)	; all pressure levels
  p_level = p_in(0,:,0,0)
  level   = dimsizes(p_level)			; creating an array of levels
  p_level@units	= "hPa"
	;printVarSummary(p_level)

; --- Other variables

z	= wrf_user_getvar(a,"z",-1)		; geopotential height
u	= wrf_user_getvar(a,"ua",-1)
v	= wrf_user_getvar(a,"va",-1)		
;rh	= wrf_user_getvar(a,"rh",-1)		; relative humidity
	;printVarSummary(z)

printVarSummary(u)

; --- Interpolated variables

geo_plane    = wrf_user_intrp3d( z,p_in,"h",p_level,0.,False)
u_plane    = wrf_user_intrp3d( u,p_in,"h",p_level,0.,False)
v_plane    = wrf_user_intrp3d( v,p_in,"h",p_level,0.,False)
;rh_plane   = wrf_user_intrp3d( rh,p_in,"h",p_level,0.,False)


; --- Degrees to radians

pii	= 3.14159265
lat_rad	= lat1
lon_rad	= lon1
lat_rad	= (pii/180)*lat1
lon_rad	= (pii/180)*lon1	; 1D lat/lon
	;printVarSummary(lat_rad)
	;printVarSummary(lon_rad)
 
;; ================================ 
;;   Begin time loop
;; ================================

  do it = 24,ntimes-1,6		; every 6 hours
;  do it = 0,18,9          ; testing

   print("Working on time: " + cd_string(time_1(it),"") +" t+"+(tointeger(time_0(it))))

    title = cd_string(time_0(it),"")

format = "%d:%H"
time_str = cd_string(time_0(it),format)


;================================================================================
; Find cyclone centre on subset of global grid (use both height and vorticity)
;================================================================================

    ; Choose a level in the lower troposphere 
      ilev = 8	; 872 hPa

    ; First, concentrate on a subset of the global grid (not needed here?)

     ; u1   = u_plane(it,ilev,:,:)
     ; v1   = v_plane(it,ilev,:,:)
      geo1 = geo_plane(it,ilev,:,:)
     ; rh1  = rh_plane(it,ilev,:,:)

     ; u1!0 = "lat"
     ; u1!1 = "lon"
     ; u1&lat	= lat1
     ; u1&lon  = lon1
     ; copy_VarMeta(u1,v1)
	
	;printVarSummary(u1)
	;printVarSummary(v1)
	;printVarSummary(geo1)


    ; RELATIVE VORTICITY
	;vort_plane    = uv2vrF_Wrap(u1,v1)
     ; vort_plane    = uv2vr_cfd(u1,v1,lat1,lon1,0)	;different function
     ; vort_plane    = vort_plane * 10^5
     ; vort_plane@units = "10~S~-5~N~ s~S~-1~N~"
     ; vort_plane@description = "Relative vorticity"

     ; vort_max      = max(vort_plane)
     ; vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     ; vort_max_smth = max(vort_smth)

     ; dims          = dimsizes(vort_plane)
     ; vort1d        = ndtooned(vort_plane)
     ; inds          = ind_resolve(maxind(vort1d),dims)
     ; vort1d_sm     = ndtooned(vort_smth)
     ; inds_sm       = ind_resolve(maxind(vort1d_sm),dims)

     ; lat_max850 = 20.3587 + 0.02524*inds(0,0)  ; Distance between lat. points - 0.025
     ; lon_max850 = -92.00711 + 0.02603*inds(0,1) ; Distance between lon. points - 0.026
     ; print("Cyclone centre (vort): "+lat_max850+" degrees N, "+abs(lon_max850)+" degrees W")

	;printVarSummary(vort_plane)
	;printVarSummary(vort_smth)
	;printVarSummary(vort1d)

    ; GEOPOTENTIAL HEIGHT
      geo_min       = min(geo1)
      geo_smth      = smth9_Wrap(geo1, 0.5, 0.5, True)
      geo_min_smth  = min(geo_smth)
	
      dims_h        = dimsizes(geo1)
      geo1d         = ndtooned(geo1)
      inds_h        = ind_resolve(minind(geo1d),dims_h)
      geo1d_sm      = ndtooned(geo_smth)
      inds_h_sm     = ind_resolve(minind(geo1d_sm),dims_h)

      lat_min850 = lat1(inds_h(0,0))
      lon_min850 = lon1(inds_h(0,1))

      ;lat_min850 = 20.3587 + 0.02524*inds_h(0,0)  ; Distance between lat. points - 0.025
      ;lon_min850 = -92.00711 + 0.02603*inds_h(0,1) ; Distance between lon. points - 0.026
      ;print("Cyclone centre (hgt): "+lat_min850+" degrees N, "+abs(lon_min850)+" degrees W")

    ; Tidy up before looping over pressure levels
      ;delete(vort1d)
      ;delete(inds)
      ;delete(vort1d_sm)
      ;delete(vort_smth)
      ;delete(inds_sm)

      ;delete(u1)
      ;delete(v1)
;      delete(rh1)
      ;delete(vort_plane)
      delete([/geo1, inds_h/])

	;delete([/u1, v1, geo1, rh1/])

;=====================================================
; Set up coordinate arrays and create new variables
;=====================================================

    ; Loop over pressure levels ; 27 levels total
      p_start = 0   	; 27 levels total
      p_end   = level-1   ; initiated at begin
      p_int   = 1

      count_p = 0                   ; Counter for loop
      num_p   = (p_end - p_start)+1 ; Vertical array size

     radii  = new(25,float)   ; New arrays
     thetas = new(73,float)

     radii      = fspan(0.0,3,25)       ; Array of radial distances (3 deg from centre)
     thetas     = fspan(0.0,360.0,73)   ; Array of azimuth angles (0, 5, 10, etc)
     thetas_rad = (pii/180)*thetas      ; Degrees to radians
     rad_size   = dimsizes(radii)-1     ; Index for later calculation

     rad_km     = new(dimsizes(radii),float)

; Create new array of radii in meters vs degrees
 do k = 0,24
    klon = 0.125*k
    rad_km(k) = gc_latlon(lat_min850,lon_min850,lat_min850,lon_min850+klon,2,4)
    delete(klon)
 end do

   ; Define arrays for later use in radial velocity cross-section plots
     vrad_av0 = new((/num_p,dimsizes(radii)/),float) ; 0   - 90  deg (NE)
     vrad_av1 = new((/num_p,dimsizes(radii)/),float) ; 90  - 180 deg (NW)
     vrad_av2 = new((/num_p,dimsizes(radii)/),float) ; 180 - 270 deg (SW)
     vrad_av3 = new((/num_p,dimsizes(radii)/),float) ; 270 - 360 deg (SE)
     vrad_av4 = new((/num_p,dimsizes(radii)/),float) ; 0 - 360 (Whole cyclone)

    ; Add metadata for azimuthally averaged radial velocity
      vrad_av0!0     = "pres"
      vrad_av0!1     = "rad"
      vrad_av0&pres  = p_level(p_start:p_end)
      vrad_av0&rad   = radii
      vrad_av0@description = "Azimuthally averaged radial velocity (NE Quadrant)"
      vrad_av0@units       = "10~S~-5~N~ s~S~-1~N~"

      copy_VarMeta(vrad_av0,vrad_av1)
      copy_VarMeta(vrad_av0,vrad_av2)
      copy_VarMeta(vrad_av0,vrad_av3)
      copy_VarMeta(vrad_av0,vrad_av4)
      
      vrad_av1@description = "Azimuthally averaged radial velocity (NW Quadrant)"
      vrad_av2@description = "Azimuthally averaged radial velocity (SW Quadrant)"
      vrad_av3@description = "Azimuthally averaged radial velocity (SE Quadrant)"
      vrad_av4@description = "Azimuthally averaged radial velocity (Whole Cyclone)"

;	delete(ilev)	; ready for level loop?

;; === Begin level loop

   do ilev = p_start, p_end, p_int ; 850, 700, 500, 300, 200

     p = p_level(ilev)
     print("Working on: " + p + " hPa")

     u1   = u_plane(it,ilev,:,:)
     v1   = v_plane(it,ilev,:,:)
     	;vort_plane    = uv2vrF_Wrap(u_plane,v_plane)
     ;vort_plane    = uv2vr_cfd(u1,v1,lat1,lon1,0)
     ;vort_plane    = vort_plane * 10^5
     ;vort_plane@units = "10~S~-5~N~ s~S~-1~N~"
     ;vort_plane@description = "Relative vorticity"

     geo1 = geo_plane(it,ilev,:,:)
 ;    rh1  = rh_plane(it,ilev,:,:)
	;geo_plane     = ht(it,ilev,:,:)
	;rh_plane      = rh(it,ilev,:,:)
     ;vort_max      = max(vort_plane)
     ;vort_smth     = smth9_Wrap(vort_plane, 0.5, 0.5, True)
     ;vort_max_smth = max(vort_smth)

     ;dims          = dimsizes(vort_plane)
     ;vort1d        = ndtooned(vort_plane)
     ;inds          = ind_resolve(maxind(vort1d),dims)
     ;vort1d_sm     = ndtooned(vort_smth)
     ;inds_sm       = ind_resolve(maxind(vort1d_sm),dims)


    ; GEOPOTENTIAL HEIGHT
      geo_min       = min(geo1)
      geo_smth      = smth9_Wrap(geo1, 0.5, 0.5, True)
      geo_min_smth  = min(geo_smth)

      dims_h        = dimsizes(geo1)
      geo1d         = ndtooned(geo1)
      inds_h        = ind_resolve(minind(geo1d),dims_h)
      geo1d_sm      = ndtooned(geo_smth)
      inds_h_sm     = ind_resolve(minind(geo1d_sm),dims_h)

      lat_max = lat1(inds_h(0,0))
      lon_max = lon1(inds_h(0,1))


     ;delete(vort1d)
     ;delete(inds)
     ;delete(vort1d_sm)
     ;delete(vort_smth)
     ;delete(inds_sm)

   ; Account for cases in which another maximum is chosen in error
     lat_diff = abs(lat_max-lat_min850)
     lon_diff = abs(lon_max-lon_min850)

   ; Set to location of geopotential height min, rather than vorticity max
    ; if (lon_diff.gt.2.0.or.lat_diff.gt.2.0) then
        lat_max = lat_min850
        lon_max = lon_min850
    ; end if
     print("Cyclone centre: "+lat_max+" degrees N, "+abs(lon_max)+" degrees W")


;===========================================
; Switch to cylindrical coordinate system
;===========================================

   ; Create new arrays describing cylindrical coordinate system (2D for now)
     xpos = new((/dimsizes(radii),dimsizes(thetas)/),float)
     ypos = new((/dimsizes(radii),dimsizes(thetas)/),float)

   ; New arrays containing variables we want to plot
     u_int      = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Zonal velocity
     v_int      = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Meridional velocity

     v_rad   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Radial velocity
     v_tan   = new((/dimsizes(radii),dimsizes(thetas)/),float) ; Tangential velocity

   ; Add variable metadata
     copy_VarAtts(u1,u_int)
     u_int@description = "Zonal velocity"
     u_int@units       = "m s~S~-1~N~"

     copy_VarAtts(v1,v_int)
     v_int@description = "Meridional velocity"
     v_int@units       = "m s~S~-1~N~"

     v_rad@description = "Radial velocity"
     v_rad@units       = "m s~S~-1~N~"

     v_tan@description = "Tangential velocity"
     v_tan@units       = "m s~S~-1~N~"

     xpos!0   = "rad"
     xpos!1   = "azi"
     xpos&rad = radii
     xpos&azi = thetas
	
     ypos!0   = "rad"
     ypos!1   = "azi"
     ypos&rad = radii
     ypos&azi = thetas

     u_int!0   = "rad"
     u_int!1   = "azi"
     u_int&rad = radii
     u_int&azi = thetas

     v_tan!0 = "rad"
     v_tan!1 = "azi"
     v_tan&rad = radii
     v_tan&azi = thetas

     v_rad!0 = "rad"
     v_rad!1 = "azi"
     v_rad&rad = radii
     v_rad&azi = thetas

   ; Latitude/longitude (small grid)
     do irad = 0, dimsizes(radii)-1
       r = radii(irad)
       do iang = 0, dimsizes(thetas)-1
         theta = thetas_rad(iang)
           xpos(irad,iang) = lon_max + r*cos(theta)
           ypos(irad,iang) = lat_max + r*sin(theta)
       end do
     end do

  ; Loop over azimuth angles to fill entire circular array of data points
    do iang = 0, dimsizes(thetas_rad)-1
      u_int(:,iang)    = linint2_points(lon1,lat1,u1,False,\
                                         xpos(:,iang),ypos(:,iang),0)
      v_int(:,iang)    = linint2_points(lon1,lat1,v1,False,\
                                        xpos(:,iang),ypos(:,iang),0)
    end do

	;printVarSummary(u_int)

  ; Calculate radial and tangential windspeed at each point
    do iang = 0, dimsizes(thetas)-1
      theta = thetas_rad(iang)
        v_rad(:,iang) = u_int(:,iang)*cos(theta) + v_int(:,iang)*sin(theta)
       ; v_tan(:,iang) = -u_int(:,iang)*sin(theta) + v_int(:,iang)*cos(theta)
    end do

	;printVarSummary(v_rad)
	;print("min " + min(v_rad(:,:)))
	;print("max " + max(v_rad(:,:)))
	;printVarSummary(v_tan)

  ; Average the radial velocity over quadrants
    do irad = 0, dimsizes(radii)-1
      r = radii(irad)
      vrad_av0(count_p,irad) = avg(v_rad(irad,{0.0:89.5}))    ; NE Quadrant
      vrad_av1(count_p,irad) = avg(v_rad(irad,{90.0:179.5}))  ; NW Quadrant
      vrad_av2(count_p,irad) = avg(v_rad(irad,{180.0:269.5})) ; SW Quadrant
      vrad_av3(count_p,irad) = avg(v_rad(irad,{270.0:359.5})) ; SE Quadrant
      vrad_av4(count_p,irad) = avg(v_rad(irad,:))             ; Whole Cyclone
    end do

   	;printVarSummary(vrad_av0) 
	;printVarSummary(vrad_av1)
	;printVarSummary(vrad_av2)
	;printVarSummary(vrad_av3)
	;printVarSummary(vrad_av4)


   count_p = count_p + 1  ; Counter for loop

    delete(u_int)  ; Tidying up
    delete(v_int)
    delete(v_tan)
    delete(v_rad)
    delete(u1)
    delete(v1)
    ;delete(vort_plane)
    delete(geo1)

  end do                  ; End pressure level loop


vrad_km0        = vrad_av0
vrad_km1        = vrad_av1
vrad_km2        = vrad_av2
vrad_km3        = vrad_av3
vrad_km4        = vrad_av4

vrad_km0!0      = "pres"
vrad_km0!1      = "rad"
vrad_km0&pres   = p_level
vrad_km0&rad    = rad_km

copy_VarMeta(vrad_km0,vrad_km1)
copy_VarMeta(vrad_km0,vrad_km2)
copy_VarMeta(vrad_km0,vrad_km3)
copy_VarMeta(vrad_km0,vrad_km4)





;=================================
; Output file type and location
;=================================

plot_path = "/nfs/a37/mm16lbb/data1/Lisa_thesis/cross_sections/nest/plots/VRAD/"  ; where plots live
    output 	= (plot_path+"nest_vrad_km_Onegeo"+time_str)
   ; wks		= gsn_open_wks("png",output)
    wks         = gsn_open_wks("pdf",output)


;==========================
; Options for plotting
;==========================

     gsn_define_colormap(wks,"GMT_polar")

  ; Radial velocity
    opts_vel                              = True
    opts_vel@cnFillOn                     = True
    opts_vel@cnLineLabelInterval          = 5.0
    opts_vel@cnLineLabelFontHeightF       = 0.012
    opts_vel@cnLineLabelBackgroundColor   = "transparent"
    opts_vel@cnLineLabelPlacementMode     = "constant"
    opts_vel@cnLinesOn                    = False ; Contour lines off
    opts_vel@cnInfoLabelOn                = False
    opts_vel@cnLevelSelectionMode         = "ExplicitLevels"
    opts_vel@cnLevels                     = fspan(-20,20,21)

    opts_vel@gsnPaperOrientation         = "portrait"
    opts_vel@tiMainString                 = title
    opts_vel@tiMainFontHeightF            = 0.01
    opts_vel@gsnLeftString                = ""
    opts_vel@gsnRightString             = ""
    
    ;opts_vel@pmLabelBarOrthogonalPosF    = .10
    opts_vel@pmLabelBarHeightF            = 0.15
    opts_vel@pmLabelBarWidthF             = 0.6
    opts_vel@lbLabelFontHeightF           = 0.01
    opts_vel@lbPerimOn                    = False
    opts_vel@lbLabelFont                  = "Helvetica"
    opts_vel@lbTitleOn                    = True
    opts_vel@lbTitleString                = "m s~S~-1~N~"
    opts_vel@lbTitleFontHeightF           = 0.015
    ;opts_vel@lbTitleDirection            = "Across"
    opts_vel@lbLabelStride                = 10
    opts_vel@lbTitlePosition              = "Bottom"

  ; Additional plotting resources
    opts_vel@tiYAxisString                = "Pressure (hPa)"
    opts_vel@trYLog                       = False
   ; opts_vel@tiXAxisString                = "radius (degrees)"
   ; opts_vel@tiXAxisJust                  = "CenterRight"
   ; opts_vel@tiXAxisPosition              = "Right"
   ; opts_vel@tiXAxisFontHeightF                  = 0.01

    opts_vel@tmXTOn                       = "False" ; Turn off top x-axis TM
    opts_vel@tmYROn                       = "False" ; Turn off right y-axis TM
    opts_vel@gsnMaximize                  = True    ; Maximise plot size
    opts_vel@gsnAddCyclic                 = False
    opts_vel@tmXBMode        = "Manual"
    opts_vel@tmXBTickStartF  = 0
    opts_vel@tmXBTickEndF    = 310


;=============
; Plot data
;=============

  ; ---- Make titles for each plot

    opts_0                              = opts_vel
    opts_0@gsnCenterString              = "NE Quadrant"
    opts_0@gsnCenterStringFontHeightF   = opts_vel@tiMainFontHeightF

    opts_1                              = opts_vel
    opts_1@gsnCenterString              = "NW Quadrant"
    opts_1@gsnCenterStringFontHeightF   = opts_vel@tiMainFontHeightF 

    opts_2                              = opts_vel
    opts_2@gsnCenterString              = "SW Quadrant"
    opts_2@gsnCenterStringFontHeightF   = opts_vel@tiMainFontHeightF

    opts_3                              = opts_vel
    opts_3@gsnCenterString              = "SE Quadrant"
    opts_3@gsnCenterStringFontHeightF   = opts_vel@tiMainFontHeightF

    opts_4                              = opts_vel
    opts_4@gsnCenterString              = "Whole Hurricane"
    opts_4@gsnCenterStringFontHeightF   = opts_vel@tiMainFontHeightF

  ; Only plot data out to 20 radii from the cyclone center
    plot_vrad0 = gsn_csm_pres_hgt(wks,vrad_km0,opts_0)
    plot_vrad1 = gsn_csm_pres_hgt(wks,vrad_km1,opts_1)
    plot_vrad2 = gsn_csm_pres_hgt(wks,vrad_km2,opts_2)
    plot_vrad3 = gsn_csm_pres_hgt(wks,vrad_km3,opts_3)
    plot_vrad4 = gsn_csm_pres_hgt(wks,vrad_km4,opts_4)

        delete([/opts_0,opts_1,opts_2,opts_3,opts_4/])
         delete(opts_vel)
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Tidy up before next iteration
    
    delete([/ilev, time_str/])

end do    ; End time loop (do ifil = 0, numINPUT-1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end

