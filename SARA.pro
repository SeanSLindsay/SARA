function band_parameters,band_location,band_wave,band_reflect,band_err, band_continuum,poly_order,win,objname,output_dir,band2_depth
;This function reads in a defined band and it's linear continuum to calculate the
;band center, band depth, and band area.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Inputs:                                                                                ;;;
; band_location:  An array of the band indices                                          ;;;
; band_wave:      An array containing the wavelength values for the band                ;;;
; band_reflect:   An arrry containing the reflectance values for the band               ;;;
; band_err:       An arrry containing the reflectance 1-sigma errors for the band       ;;;
; band_continuum: An array containining the linear continuum for the band               ;;;
; poly_order:     The order of the polynomial fit for band center and area              ;;;
; win:            Integer to indicate which window to plot to [automatically set]       ;;;
; objname:        String holding object name                                            ;;;
; output_dir:     Directory for output [band_parameters creates JPGs of plots]          ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Outputs:                                                                               ;;;
; band_params:  A structure containing the band center, band depth, and band area       ;;;
;               with associated errors.                                                 ;;;
;  .bc - Band Center; .bc_err - Band Center Error                                       ;;;
;  .bd - Band Depth;  .bd_err - Band Depth Error                                        ;;;
;  .ba - Band Area;   .ba_err - Band Area Error                                         ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


band_size = n_elements(band_wave)
band_indices_sub=indgen(band_size)
band_smooth = smooth(band_reflect,10)

;Define structure to hold band center, band center error, band depth, band depth error, band area, and band area error
band_params = {bc:0.0, bc_err:0.0, bd:0.0, bd_err:0.0, ba:0.0, ba_err:0.0}

;Divide by the continuum
if win eq 1 then begin
  band_reflect_div = band_reflect/band_continuum
  band_continuum_div = band_continuum/band_continuum
  band_err_div = band_err/band_continuum
endif else begin
  band_reflect_div = band_reflect/band_continuum
  band_continuum_div = band_continuum/band_continuum
  band_err_div = band_err/band_continuum
endelse
  

;BAND CENTER
;Locate the band minimum
;The Band minimum will be defined by the minimum of the band when smoothed by a box 10 pixels in width
band_smooth_div = smooth(band_reflect_div,10)
band_minimum_loc_sub = where(band_smooth_div eq min(band_smooth_div))
band_minimum_loc = band_location[band_minimum_loc_sub]
band_minimum_wave = band_wave[band_minimum_loc_sub]
;print, band_minimum_loc_sub, band_minimum_loc, band_minimum_wave

;Calculate Band Depth
band_poly_coeff = poly_fit(band_wave,band_reflect_div,2,yfit=band_poly)
band_minimum_loc_sub = where(band_poly eq min(band_poly))
band_minimum_loc = band_location[band_minimum_loc_sub[0]]
if win eq 1 then band_depth_fraction = 1/2.  ;Change this for BAND I Depth Fraction
if win eq 2 then band_depth_fraction = band2_depth  ;Change this for BAND II Depth Fraction ;S-type => set to 1.; V-type => set to 0.5
band_depth_min_inverse = band_poly[band_minimum_loc_sub[0]]
band_depth_min = 1.0 - band_depth_min_inverse

if band_depth_fraction ne 1.0 then begin
  band_depth_boundary = band_depth_min_inverse[0] + band_depth_fraction*band_depth_min
  band_blue_boundary = band_minimum_loc_sub[0]
  band_red_boundary = band_minimum_loc_sub[0]
  
  ;if band_red_boundary eq n_elements(band_wave) then band_red_boundary = n_elements(band_wave)-1
  
  ;Radiate outward from minimum to find the first two points (+/-) that are outside the fractional band depth boundary
  while band_smooth_div[band_blue_boundary] le band_depth_boundary do band_blue_boundary = band_blue_boundary - 1
  whileflag=0
  while band_smooth_div[band_red_boundary] le band_depth_boundary and whileflag eq 0 do begin
    band_red_boundary = band_red_boundary + 1
    if band_red_boundary eq band_size then begin
      band_red_boundary = band_size - 1
      whileflag=1
    endif
  endwhile
endif else begin
  band_blue_boundary = 0
  band_red_boundary = band_size-1
endelse

;  band_red_boundary = band_red_boundary + 1
;  far_red_indices = where(band_wave gt 2.1)
;  far_red_reflect_div = band_reflect_div[far_red_indices]
;  far_red_max = where(far_red_reflect_div eq max(far_red_reflect_div))+min(far_red_indices)
;  if band_red_boundary eq far_red_max then whileflag=1
;endwhile
band_fit_bounds = band_indices_sub[band_blue_boundary:band_red_boundary]
band_fit_wave = band_wave[band_fit_bounds]
band_fit_reflect_div = band_reflect_div[band_fit_bounds]
band_fit_smooth_div = smooth(band_fit_reflect_div,10)
band_fit_err = band_err[band_fit_bounds]
band_fit_coeff = poly_fit(band_fit_wave,band_fit_reflect_div,poly_order,yfit=band_fit,status=statusnum,chisq=band_chisq)
band_center_location = where(band_fit eq min(band_fit),ncent)
if ncent gt 1 then begin
  print, 'More than one center found'
  band_center = mean(band_fit_wave[band_center_location])
endif else begin
band_center=band_fit_wave[band_center_location[0]]
endelse
;if win eq 2 then stop
band_center_err = band_center_montecarlo(band_fit_bounds, band_fit_wave, band_fit_reflect_div, band_fit_err, poly_order, 20000)
band_params.bc=band_center
band_params.bc_err = band_center_err
;BAND DEPTH 
;Method.  Take the range of band centers +/- the band_center_err, and compute the average noise.
;    Take this noise as the error in depth
band_depth_inverse = band_fit_smooth_div[band_center_location[0]]
band_depth = 1.0 - band_depth_inverse[0]
band_center_region = where(band_wave gt band_center-band_center_err and band_wave lt band_center+band_center_err, nbc_area)
band_depth_err = mean(band_err[band_center_region])
band_params.bd=band_depth
band_params.bd_err=band_depth_err
!p.font=1
window, win, title = 'Divided Band '+strtrim(string(win),2)
plot, band_wave, band_continuum_div,xr=[band_wave[0]-0.025,band_wave[band_size-1]+0.025],yr=[0.4,1.1], xs=1,ys=1, $
  xtitle = 'Wavelength (um)', ytitle='Continuum Divided Reflectance',$
  thick=2,charsize=3,xthick=3,ythick=3,charthick=3
oplot, band_wave, band_reflect_div,thick=2,psym=4,symsize=0.8
oplot, band_wave, band_smooth_div,thick=1, color=35
oplot, [band_minimum_wave,band_minimum_wave], [0.4,1.1], lines=2,thick=3
oplot, [band_center,band_center], [0.4,1.1],lines=2,color=50,thick=4
oploterr, band_wave, band_reflect_div,band_err,3
oplot, band_wave[band_fit_bounds], band_reflect_div[band_fit_bounds],color=250,thick=3
oplot, band_fit_wave,band_fit, color=90,thick=3
xyouts, band_center+0.01, 0.45, strtrim(string(band_center,format='(f6.3)'),2) + ' +/- ' + strtrim(string(band_center_err,format='(f6.3)'),2),charsize=3,charthick=2
xyouts, band_center+0.01, band_depth_inverse-0.1, 'Band Depth = ' + strtrim(string(band_depth,format='(f6.3)'),2) + ' +/- ' + strtrim(string(band_depth_err,format='(f6.3)'),2),charsize=3,charthick=2
xyouts, band_wave[0]+0.05, band_depth_inverse-0.15, 'Fit Chi-Square = ' +strtrim(string(band_chisq,format='(f6.3)'),2),charsize=3,charthick=2
xyouts, band_wave[0]+0.05,0.5,'Polynomial Order: '+strtrim(string(poly_order),2),charsize=3,charthick=2
;Save a .tiff image file of the band analyses
image = tvrd(true=1)
image_file = output_dir+'band_images/'+objname+'_band_'+strtrim(string(win),2)+'_po'+strtrim(string(poly_order),2)+'.jpg'
write_jpeg, image_file,tvrd(/true),quality=85,/true

;if win eq 2 then stop
;Band Area
;Band area is calculated using the trapazoidal rule on the 
;divided continuum minus the smooth (10 pixel box) of the divided reflectance
band_area=0
;HERE
;stop
area_function = band_continuum_div - band_reflect_div
for i_a=0,band_size-2,1 do band_area = band_area + 0.5*(area_function[i_a+1]+area_function[i_a])$
  *(band_wave[i_a+1] - band_wave[i_a])
;barea = tsum(band_wave,area_function)
;stop
;Band Area Error
;Method 1: Band Area Error is computed from a 20001 iteration Monte Carlo Simulation of false band data
band_area_err = band_area_montecarlo(band_wave,band_reflect_div, band_continuum_div, band_err,20000)
;
;;Method 2 - Moskovitz et al. (2010) Band area error is computed by determining the maximum and minimum areas based on the 1-sigma reflectance errors
;;Method:  Max. Area is computed by subtracting 1-sigma reflectance error from the divided band reflectances
;;         Min. Area is computed by adding 1-sigma reflectance errors to the divided band reflectances
;;   The 1-sigma error for the band area is then assumed to be half the difference between the maximum and minimum areas.
;band_area_max=0
;band_area_min=0
;area_function_max = band_continuum_div - smooth( (band_reflect_div-band_err),10)
;area_function_min = band_continuum_div - smooth( (band_reflect_div+band_err),10)
;for i_amax=0,band_size-2,1 do band_area_max = band_area_max + 0.5*(area_function_max[i_amax+1]+area_function_max[i_amax])$
;  *(band_wave[i_amax+1] - band_wave[i_amax])
;for i_amin=0,band_size-2,1 do band_area_min = band_area_min + 0.5*(area_function_min[i_amin+1]+area_function_min[i_amin])$
;  *(band_wave[i_amin+1] - band_wave[i_amin])
;band_area_err = (band_area_max - band_area_min)/2.
;print, 'Band Area: '+strtrim(string(band_area),2) + ' +/- ' + strtrim(string(band_area_err),2)
;print, band_area_min, band_area_max
band_params.ba=band_area
band_params.ba_err=band_area_err
return, band_params
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function band_center_montecarlo,fit_ind, fit_wave, fit_reflect_div, fit_err, poly_order, iterations
;This function performs a monte carlo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Inputs:                                                                                ;;;
; fit_ind:         An array containing the fit region indices                           ;;;
; fit_wave:        An array containing the wavelength values for the fit region         ;;;
; fit_reflect:     An arrry containing the reflectance values for the fit region        ;;;
; fit_err:         An array containining the 1-sigma errors for the fit region          ;;;
; poly_order:      The order for the polynomial fit                                     ;;;
; iterations:      The number of iterations the monte carlo simulation is done over     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Outputs:                                                                               ;;;
; band_center_err: The 1-sigma error associated with the monte carlo simulation         ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
centers=fltarr(iterations)
fit_pts = n_elements(fit_wave)

for i_mc=0, iterations-1,1 do begin
fit_reflect_div_new = fit_reflect_div + randomn(seed,fit_pts)*fit_err

fit_coeff = poly_fit(fit_wave,fit_reflect_div_new,poly_order,yfit=new_fit,status=sts)
fit_center_new = where(new_fit eq min(new_fit),ncent)
if ncent gt 1 then begin
;  print, 'More than one center found, taking the average'
  fit_center = mean(fit_wave[fit_center_new])
  centers[i_mc]=fit_center
endif else begin
fit_center=fit_wave[fit_center_new[0]]
centers[i_mc]=fit_center
endelse

endfor
center_stats=moment(centers,sdev=center_err)
center_diff = ( max(centers)-min(centers) )/2.
;print, center_err, center_diff
band_center_err = center_err
;band_center_err = max([center_err,center_diff])
;stop
return, band_center_err

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function band_area_montecarlo,band_wave,band_reflect,band_continuum, band_err,iterations
;This function performs a monte carlo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Inputs:                                                                                                    ;;;
; band_wave:           An array containing the wavelength values for the band region                        ;;;
; band_smooth_div:     An arrry containing the divided reflectance values for the band region               ;;;
; band_continuum_div:  An array containining the divided band continuum                                     ;;;
; band_err:            An array containining the 1-sigma errors for the band region                         ;;;
; iterations:          The number of iterations the monte carlo simulation is done over                     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Outputs:                                                                                                   ;;;
; band_area_err: The 1-sigma error associated with the monte carlo simulation                               ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
areas=fltarr(iterations)
band_pts = n_elements(band_wave)

for i_mc=0, iterations-1,1 do begin
;Generate False Band Data
band_reflect_new = band_reflect + randomn(seed,band_pts)*band_err
band_smooth_new = smooth(band_reflect_new,10)

areas[i_mc]=0
area_function = band_continuum - band_reflect_new
for i_a=0,band_pts-2,1 do areas[i_mc] = areas[i_mc] + 0.5*(area_function[i_a+1]+area_function[i_a]) $
  *(band_wave[i_a+1] - band_wave[i_a])

endfor
area_stats=moment(areas,sdev=area_err)
;print, center_err, center_diff
band_area_err = area_err
;band_center_err = max([center_err,center_diff])
;stop
return, band_area_err

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function band_slope_montecarlo, wave,cont1_ind,cont1_wave, cont1_reflect, cont1_err, cont2_ind,cont2_wave, cont2_reflect,cont2_err, iterations


slopes = fltarr(iterations)
n_cont1 = n_elements(cont1_wave)
n_cont2 = n_elements(cont2_wave)
for i_mc=0,iterations-1,1 do begin
  ;Generate False Continuum Region Data
  cont1_new_reflect = cont1_reflect + randomn(seed,n_cont1)*cont1_err
  cont2_new_reflect = cont2_reflect + randomn(seed,n_cont2)*cont2_err
  
  cont1_new_fit_coeff = poly_fit(cont1_wave,cont1_new_reflect,5,yfit=cont1_new_fit,status=sts_cont1)
  cont1_new_max_ind_sub = where(cont1_new_fit eq max(cont1_new_fit))
  cont1_new_max_ind = cont1_new_max_ind_sub[0]+cont1_ind[0]
  cont1_new_max_wave = cont1_wave[cont1_new_max_ind_sub]
  cont1_new_max_reflect = cont1_new_reflect[cont1_new_max_ind_sub]
  
  cont2_new_fit_coeff = poly_fit(cont2_wave,cont2_new_reflect,5,yfit=cont2_new_fit,status=sts_cont2)
  cont2_new_max_ind_sub = where(cont2_new_fit eq max(cont2_new_fit))
  cont2_new_max_ind = cont1_new_max_ind_sub[0]+cont2_ind[0]
  cont2_new_max_wave = cont2_wave[cont2_new_max_ind_sub] 
  cont2_new_max_reflect = cont2_new_reflect[cont2_new_max_ind_sub] 
  band_slope_p2p = (cont2_new_max_reflect[0] - cont1_new_max_reflect[0]) / (cont2_new_max_wave[0] - cont1_new_max_wave[0])
  ;Find the nearest point where  the slope of the continuum 2 polynomial fit is approximately the slope of the peak-to-peak continuum
  cont2_new_max_reflect_prime = deriv(cont2_wave, cont2_new_fit)
  cont2_new_max_reflect_dbl_prime = deriv(cont2_wave, cont2_new_max_reflect_prime)
  ;Note:  The second derivative tests to see if the potential tangent point is concave up or down.
  ; The tangent point we are looking for is the concave up point
  cont2_lt_slopes = where(cont2_new_max_reflect_prime lt band_slope_p2p and cont2_new_max_reflect_dbl_prime lt 0)
  tangent_point_sub = cont2_lt_slopes[0]
  tangent_point = cont2_lt_slopes[0]+min(cont2_new_max_ind)
  wave_tp = wave[tangent_point]
  reflect_tp = cont2_new_fit[tangent_point_sub]
  ;Find the new linear continuum for Band 1
  slopes[i_mc] = ( reflect_tp - cont1_new_max_reflect[0]) / ( wave_tp - cont1_new_max_wave[0] )
endfor

slope_stats = moment(slopes,sdev=band_slope_err)
return, band_slope_err

end
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro SARA,red_edge_median=red_edge_median, v_type=v_type, s_type=s_type ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This IDL routine is designed to work on visible + NIR (0.45 - 2.45um) spectral data.
;The routine computes the band centers, band areas, and band depths for the 1- and 2-um
;absorption bands due to Fe2+ in olivine+pyroxene and pyroxene, respectively.
; -Updated to work for data sets with NIR data only (12 Aug. 2013)
;
;The routine computes these values for 3rd, 4th, and 5th order polynomial fits to
;user defined percentage of band depth (I suggest full depth for noisy Band 2 data).
;For details see band_analysis_ssl_usersmanual.pdf
;
;Please reference Lindsay et al. 2013/2014, in prep. Icarus? if you make use of this routine
;Until the paper is released, reference:
;2013 Lindsay, S. S., J. P. Emery, F. Marchis, J. E. Enriquez, and M. Assafin. 
;“A spectroscopic and mineralogical study of multiple asteroid systems. 
; At The Division for Planetary Sciences (DPS) #45 Abstract# 112.04
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Inputs:  Three column spectral data with wavelengths, reflectances, and reflectance errors
;         Will locate and read in object from user prompt
;         Be sure to set directories below
;         3 Directories Needed:  Spectra location, Output, and output/band_images/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Output:  band_analysis:  Structure of Band Parameters 
;           (centers, areas, depths, Band 1 slope, and Band Area Ratio) in
;         band_analysis_avg:  Structure of the band parameter values averaged
;           over 3rd, 4th, and 5th polynomial fits
;   Note: band_analysis and band_analysis_avg are saved as binary IDL save files (.sav) in output directory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Author:  Sean S. Lindsay.  Earth and Planetary Sciences, University of Tennessee, Knoxville
;Date:  4 June 2013
;Last Modifiedy:  18 Feb. 2014 by SSL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

loadct,39,/silent
smoothwidth=10

;Generate directory paths
saradir='~/SARA/'

if file_search(saradir,/test_directory) eq '' then file_mkdir, saradir
cd, saradir
if file_search(saradir+'spectra/',/test_directory) eq '' then begin
  file_mkdir, saradir+'/spectra/'
  print, 'Please place spectral data in the following directory'
  print, 'before continuing'
  print, '~/SARA/spectra/'
  print, 'File format = ID#_OBJECTNAME'
  goto, done
endif
specdir=saradir+'spectra/'   ;Location of spectra with wavelengths, reflectances, and errors
if file_search(saradir+'/output/',/test_directory) eq '' then file_mkdir, 'output'
output_dir=saradir+'output/' ;Location for output: Files and band images in [output_path]/band_images/
objname=''
red_edge_po=5
red_edge_wave=2.44
default_band2_depth = 1.0
type_flag=0

;Check output directories for proper subdirectories and create them if needed
output_images = file_search(output_dir+'/band_images/',/test_directory)
output_defs = file_search(output_dir+'/band_definitions/',/test_directory)
output_master = file_search(output_dir+'/master_ouput/',/test_directory)
ast_temps = file_search(output_dir+'/ast_temps/',/test_directory)
if output_images eq '' then file_mkdir, output_dir+'band_images'
if output_defs eq '' then file_mkdir, output_dir+'band_definitions'
if output_master eq '' then file_mkdir, output_dir+'master_output'
if ast_temps eq '' then file_mkdir, output_dir+'ast_temps'
;stop
;Set keyword flag value
if keyword_set(s_type) eq 1 then type_flag = 2 
;band2_depth = 1.0 else band2_depth = default_band2_depth
if keyword_set(v_type) eq 1 then type_flag = 1 
case type_flag of
  0: begin
      band2_depth = default_band2_depth
      ast_type = 'U'
     end
  1: begin
      band2_depth = 0.5
      ast_type = 'V'
     end
  2: begin
      band2_depth = 1.0
      ast_type = 'S'
     end
endcase

;Input data
read,'Enter object name. [filename]: ',objname                    
spectrum=findfile(specdir+objname+'*.txt')
readcol,spectrum,wave1,reflect1,err1, FORMAT='F,F,F'
ndata = n_elements(wave1)
;Filter out potentially bad reflectance points
reflect_filter = where(reflect1 ge 0.0)
reflect = reflect1[reflect_filter]
wave = wave1[reflect_filter]
err = err1[reflect_filter]
reflect_smooth = smooth(reflect,smoothwidth)

;Create a structure to hold the entire band analysis for all polynomial fits
band_analysis_i = {band_analysis_i, poly_order:0.0,b1c:0.0, b1c_err:0.0, b1d:0.0, b1d_err:0.0, b1a:0.0, b1a_err:0.0,b1_slope:0.0, b1_slope_err:0.0,$
  b2c:0.0, b2c_err:0.0, b2d:0.0, b2d_err:0.0, b2a:0.0, b2a_err:0.0,bar:0.0,bar_err:0.0}
band_analysis = replicate(band_analysis_i,3)

;Loop over polynomial fits from 3rd to 5th
for poly_order=3,5,1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Stage 1: DEFINING THE CONTINUUM AND BAND BOUNDARIES                              ;;
;;Find the reflectance maxima near 0.7 and 1.4 um from polynomial fits to the      ;;
;;the wavelength regions 0.65 - 0.9 um [5th order] and 1.10 - 1.75 um [5th order]. ;;
;;The maxima are taken to be the maxima of the polynomial fits.                    ;;
;;This methodology follows that of Moskovitz et al. (2010, 2012).                  ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Locate the reflectance maxima in the data to be used to build object specific continuum fitting regions
if wave[0] lt 0.6 then begin
region1=where(wave gt 0.65 and wave lt 0.9)     ;General bounds for continuum maximum near 0.75 um
region1_wave = wave[region1]
region1_sm_reflect = reflect_smooth[region1]
region1_reflect_max_loc_sub = where(region1_sm_reflect eq max(region1_sm_reflect),nreg1)
region1_reflect_max_loc = region1_reflect_max_loc_sub[0] + min(region1)
region1_max_wave = wave[region1_reflect_max_loc[0]]
cont1_max_bounds = [region1_max_wave - 0.125,region1_max_wave + 0.125]

region2=where(wave gt 1.1 and wave lt 1.75)     ;General bounds for continuum maximum near 1.5 um
region2_wave = wave[region2]
region2_sm_reflect = reflect_smooth[region2]
region2_reflect_max_loc_sub = where(region2_sm_reflect eq max(region2_sm_reflect),nreg2)
region2_reflect_max_loc = region2_reflect_max_loc_sub[0] + min(region2)
region2_max_wave = wave[region2_reflect_max_loc[0]]
cont2_max_bounds = [region2_max_wave - 0.3,region2_max_wave + 0.3]
endif else begin  ;Else statement is for when there is only IR Data
  region1=where(wave gt wave[0] and wave lt 0.9)
  region1_wave = wave[region1]
  region1_sm_reflect = reflect_smooth[region1]
  region1_reflect_max_loc_sub = where(region1_sm_reflect eq max(region1_sm_reflect),nreg1)
  region1_reflect_max_loc = region1_reflect_max_loc_sub[0] + min(region1)
  region1_max_wave = wave[region1_reflect_max_loc[0]]
  cont1_max_bounds = [region1_max_wave - 0.125,region1_max_wave + 0.125]
  if objname eq '8116_Jeanperrin' then cont1_max_bounds = [region1_max_wave - 0.125,region1_max_wave + 0.25]
  region2=where(wave gt 1.1 and wave lt 1.75)
  region2_wave = wave[region2]
  region2_sm_reflect = reflect_smooth[region2]
  region2_reflect_max_loc_sub = where(region2_sm_reflect eq max(region2_sm_reflect),nreg2)
  region2_reflect_max_loc = region2_reflect_max_loc_sub[0] + min(region2)
  region2_max_wave = wave[region2_reflect_max_loc[0]]
  cont2_max_bounds = [region2_max_wave - 0.3,region2_max_wave + 0.3]
;  print, 'IR ONLY' & stop
endelse
;stop
;Continuum 1
;cont1_max_bounds = [0.65,0.9]
cont1_max_ind = where(wave gt cont1_max_bounds[0] and wave lt cont1_max_bounds[1])
 ;Some objects require user definition [Left in code as example]
 if objname eq '4607_Seilandfarm' then cont1_max_ind = where(wave gt 0.6 and wave lt 0.9) 
cont1_max_wave = wave[cont1_max_ind]
cont1_max_reflect = reflect[cont1_max_ind]
cont1_max_err = err[cont1_max_ind]
cont1_max_coeff = poly_fit(cont1_max_wave,cont1_max_reflect,5, yfit=cont1_max_fit)  ;Fit continuum peak with 5th order polynomial
cont1_max_loc_sub = where(cont1_max_fit eq max(cont1_max_fit))                      ;Index within band 1 region
cont1_fit_maximum_wave = cont1_max_wave[cont1_max_loc_sub]
cont1_fit_maximum_reflect = cont1_max_fit[cont1_max_loc_sub]
cont1_max_loc = where(cont1_max_fit eq max(cont1_max_fit)) + min(cont1_max_ind)     ;Index for entire spectrum

;stop
;Continuum 2
;cont2_max_bounds = [1.1,1.75]
cont2_max_ind = where(wave gt cont2_max_bounds[0] and wave lt cont2_max_bounds[1])
 ;Some objects require user definition [Left in code as example]
 if objname eq '939_Isberga' then cont2_max_ind = where(wave gt 1.275 and wave lt 1.7)
 if objname eq '7225_Huntress' then cont2_max_ind = where(wave gt 1.42 and wave lt 1.9)
 if objname eq '7225_sHuntress' then cont2_max_ind = where(wave gt 1.42 and wave lt 1.75)
 if objname eq '6708_Bobbievaile' then cont2_max_ind = where( wave gt 1.3 and wave lt 1.9)
 if objname eq '34706_2001op83' then cont2_max_ind = where( wave gt 1.2 and wave lt 1.6)
cont2_max_wave = wave[cont2_max_ind]
cont2_max_reflect = reflect[cont2_max_ind]
cont2_max_err = err[cont2_max_ind]
cont2_max_coeff = poly_fit(cont2_max_wave,cont2_max_reflect,5, yfit=cont2_max_fit)
;if objname eq '7225_Huntress' then cont2_max_coeff = poly_fit(cont2_max_wave,cont2_max_reflect,2, yfit=cont2_max_fit)
cont2_max_loc_sub = where(cont2_max_fit eq max(cont2_max_fit))                     ;Index within band 2 region
;if objname eq '7225_Huntress' then cont2_max_loc_sub = cont1_max_loc - min(cont2_max_ind[0])
cont2_fit_maximum_wave = cont2_max_wave[cont2_max_loc_sub]
cont2_fit_maximum_reflect = cont2_max_fit[cont2_max_loc_sub]
cont2_max_loc = where(cont2_max_fit eq max(cont2_max_fit)) + min(cont2_max_ind)    ;Index for entire spectrum
cont2_fit_function = cont2_max_coeff[0]+cont2_max_coeff[1]*wave+cont2_max_coeff[2]*wave^2;+cont2_max_coeff[3]*wave^3;+cont2_max_coeff[4]*wave^4
;stop
;;;;;This will most likely become a subroutine/function
;Find the linear continuum for Band 1 defined by:
;1) The maximum value for the 4th order polynomial fit to the blue edge Band 1 maximum
;2) The point tangent to the 3rd polynomial fit to the red edge Band 1 maximum
;Note:  A tangent point between the two intercepts of where the peak-to-peak continuum crosses
; the continuum 2 polynomial fit is gauranteed by the Mean Value Theorem 
band1_slope_p2p = (cont2_fit_maximum_reflect[0] - cont1_fit_maximum_reflect[0]) / (cont2_fit_maximum_wave[0] - cont1_fit_maximum_wave[0])
band1_intercept_p2p = cont1_fit_maximum_reflect[0] - cont1_fit_maximum_wave[0]*band1_slope_p2p
;Create continuum 1
continuum1_p2p = wave*band1_slope_p2p + band1_intercept_p2p

;Find the nearest point where  the slope of the continuum 2 polynomial fit is approximately the slope of the peak-to-peak continuum
cont2_max_reflect_prime = deriv(cont2_max_wave, cont2_max_fit)
cont2_max_reflect_dbl_prime = deriv(cont2_max_wave, cont2_max_reflect_prime)
;Note:  The second derivative tests to see if the potential tangent point is concave up or down.
; The tangent point we are looking for is the concave up point
cont2_lt_slopes = where(cont2_max_reflect_prime lt band1_slope_p2p and cont2_max_reflect_dbl_prime lt 0)
tangent_point_sub = cont2_lt_slopes[0]
tangent_point = cont2_lt_slopes[0]+min(cont2_max_ind)
wave_tp = wave[tangent_point]
reflect_tp = cont2_max_fit[tangent_point_sub]
;Find the new linear continuum for Band 1
band1_slope = ( reflect_tp - cont1_fit_maximum_reflect[0]) / ( wave_tp - cont1_fit_maximum_wave[0] )
band1_intercept = cont1_fit_maximum_reflect[0] - cont1_fit_maximum_wave[0]*band1_slope
continuum1 = wave*band1_slope + band1_intercept

;if objname eq '7225_Huntress' then begin
;  cont2_max_loc = tangent_point
;  cont2_fit_maximum_reflect = reflect[cont2_max_loc]
;endif

;Band 1 Slope Error
;Method: Monte-Carlo Simulation using false data generated for the two continuum pieces
band1_slope_err = band_slope_montecarlo(wave,cont1_max_ind,cont1_max_wave, cont1_max_reflect,cont1_max_err,cont2_max_ind,cont2_max_wave, cont2_max_reflect,cont2_max_err,20000)

;Define Bounds for Band 1
;Band definitions are the arguments for the band_parameters function.
b1_size = tangent_point - cont1_max_loc[0]+1
band1_indices=indgen(b1_size)+cont1_max_loc[0]
band1_wave = wave[band1_indices]
band1_reflect = reflect[band1_indices]
band1_err = err[band1_indices]
band1_continuum = continuum1[band1_indices]

band2_244_region = where(wave gt cont2_fit_maximum_wave[0] and wave le red_edge_wave,n244)
band2_244_wave = wave[band2_244_region]
band2_244_reflect = reflect[band2_244_region]
band2_polyfit_coeffs = poly_fit(band2_244_wave,band2_244_reflect,red_edge_po,yfit=band2_polyfit)
band2_red_edge_reflect = band2_polyfit[n244-1]
band2_red_wave = band2_244_wave
band2_red_reflect = band2_244_reflect
band2_red_edge_loc = band2_244_region[n244-1]
band2_red_edge_wave = band2_244_wave[n244-1]

;Definitions for Band 2
;Band definitions are the arguments for the band_parameters function.
b2_size = band2_red_edge_loc-cont2_max_loc[0]+1
band2_indices=indgen(b2_size)+cont2_max_loc[0]
band2_wave = wave[band2_indices]
band2_reflect = reflect[band2_indices]
band2_err = err[band2_indices]

band2_slope = ( band2_red_edge_reflect - cont2_fit_maximum_reflect[0] ) / (band2_red_edge_wave - band2_wave[0])
band2_intercept = cont2_fit_maximum_reflect[0] - cont2_fit_maximum_wave[0]*band2_slope
continuum2 = wave*band2_slope + band2_intercept
band2_continuum = continuum2[band2_indices]
;stop
!p.font=1
window,0
plot, wave, reflect,xr=[0.45,2.5],xs=1,yr=[0.4,max(reflect)+0.2],ys=1,charsize=3,thick=2,psym=4,symsize=0.5,$
  xtitle='Wavelength (um)', ytitle='Normalized Reflectance',xthick=4,ythick=5,charthick=2
oploterr, wave, reflect, err,3
oplot, wave[cont1_max_ind],cont1_max_fit,thick=3,color=150
oplot, wave[cont2_max_ind],cont2_max_fit,thick=3,color=150
;oplot, [wave[cont1_max_loc],wave[cont1_max_loc]],[0.4,max(reflect)+0.2], lines=2
;oplot, [wave[cont2_max_loc],wave[cont2_max_loc]],[0.4,max(reflect)+0.2], lines=2
;oplot, wave[cont1_max_loc:cont2_max_loc],continuum1_p2p[cont1_max_loc:cont2_max_loc],color=35,thick=2
oplot, wave[band1_indices],continuum1[band1_indices], thick=3,color=50
;oplot, band2_red_edge_wave,band2_red_edge_fit,color=250,thick=2
;oplot, band2_extrap_wave, band2_extrap_reflect,color=220,lines=2,thick=2
oplot, wave[band2_indices],continuum2[band2_indices],thick=3,color=50
oplot, band2_red_wave,band2_polyfit,color=90,thick=3
;oplot, band2_red_wave,band2_gaussfit,color=250,thick=3
xyouts, 0.7, max(reflect)-0.3, strtrim(objname,2),charsize=3,charthick=2
image = tvrd(true=1)
image_file = output_dir+'band_images/'+objname+'_fullrange'+'.jpg'
write_jpeg,image_file,tvrd(/true),quality=85,/true

;Calculate Band Area Centers and Band Areas
win=1 ;Window Index
;Band 1
band1_params = band_parameters(band1_indices, band1_wave, band1_reflect, band1_err, band1_continuum,poly_order,win,objname,output_dir,band2_depth)
win=2
;Band 2
band2_params = band_parameters(band2_indices, band2_wave, band2_reflect, band2_err, band2_continuum,poly_order,win,objname,output_dir,band2_depth)

 if objname eq '7225_Huntress' then band2_params.ba = band2_params.ba_err
 if objname eq '7225_sHuntress' then band2_params.ba = band2_params.ba_err
 if objname eq '6708_Bobbievaile' then band2_params.ba = band2_params.ba_err
 if objname eq '762_Pulcova' then band2_params.ba = band2_params.ba_err
 if objname eq '762_sPulcova' then band2_params.ba = band2_params.ba_err
;Collate results into single structure
band_analysis[poly_order-3].poly_order = poly_order
band_analysis[poly_order-3].b1c = band1_params.bc
band_analysis[poly_order-3].b1c_err = band1_params.bc_err 
band_analysis[poly_order-3].b1d = band1_params.bd
band_analysis[poly_order-3].b1d_err = band1_params.bd_err
band_analysis[poly_order-3].b1a = band1_params.ba
band_analysis[poly_order-3].b1a_err = band1_params.ba_err
band_analysis[poly_order-3].b1_slope = band1_slope
band_analysis[poly_order-3].b1_slope_err = band1_slope_err

band_analysis[poly_order-3].b2c = band2_params.bc
band_analysis[poly_order-3].b2c_err = band2_params.bc_err 
band_analysis[poly_order-3].b2d = band2_params.bd
band_analysis[poly_order-3].b2d_err = band2_params.bd_err
band_analysis[poly_order-3].b2a = band2_params.ba
band_analysis[poly_order-3].b2a_err = band2_params.ba_err


;Compute the Band Area Ratio and BAR Error
bar = band2_params.ba/band1_params.ba
bar_err = bar*sqrt( (band1_params.ba_err/band1_params.ba)^2 + (band2_params.ba_err/band2_params.ba)^2 )

band_analysis[poly_order-3].bar = bar
band_analysis[poly_order-3].bar_err = bar_err

print, 'x1 = '+strtrim(string(wave[cont1_max_loc]),2)
print, 'x2 = '+strtrim(string(wave_tp),2)
print, 'x3 = '+strtrim(string(band2_wave[0]),2)
print, 'x4 = '+strtrim(string(red_edge_wave),2)

;Print Band Analysis Parameters to Screen
print, '---------------'
print, 'Poly Order = '+strtrim(string(poly_order),2)
print, '---------------'
print, 'Band I Center: '+strtrim(string(band1_params.bc),2)+' +/- '+strtrim(string(band1_params.bc_err),2)
print, 'Band I Depth: '+strtrim(string(band1_params.bd),2)+' +/- '+strtrim(string(band1_params.bd_err),2)
print, 'Band I Area: '+strtrim(string(band1_params.ba),2)+' +/- '+strtrim(string(band1_params.ba_err),2)
print, 'Band I Slope: '+strtrim(string(band1_slope),2)+' +/- '+strtrim(string(band1_slope_err),2)
print, ''
print, 'Band II Center: '+strtrim(string(band2_params.bc),2)+' +/- '+strtrim(string(band2_params.bc_err),2)
print, 'Band II Depth: '+strtrim(string(band2_params.bd),2)+' +/- '+strtrim(string(band2_params.bd_err),2)
print, 'Band II Area: '+strtrim(string(band2_params.ba),2)+' +/- '+strtrim(string(band2_params.ba_err),2)
print, 'Band II Slope: '+strtrim(string(band2_slope),2)
print, ''
print, 'Band Area Ratio (BAR): '+strtrim(string(bar),2)+' +/- '+strtrim(string(bar_err),2)
print, '---------------'


;Write the Band Analysis parameters to a file
output= output_dir + objname + '_bandanalysis_' + strtrim(string(poly_order),2)+'.txt'
openw,lun11,output,/get_lun,width=300
printf,lun11,'Object','    Atype', '   BI Center  ','BI Cent. Err.', '   BI Depth  ', ' BI Depth Err.', ' BI Area ','  BI Area Err.','  BI Slope ', ' BI Slope Err.',$
                      ' BII Center  ',' BII Cent. Err.', '  BII Depth  ', ' BII Depth Err.', '  BII Area  ',' BII Area Err.', '  BAR     ', '    BAR Err.  ',$
                      ' Poly. Order'
printf,lun11,objname, '  '+ast_type, band1_params.bc,band1_params.bc_err,band1_params.bd,band1_params.bd_err,band1_params.ba,band1_params.ba_err, band1_slope, band1_slope_err,$
                      band2_params.bc,band2_params.bc_err,band2_params.bd,band2_params.bd_err,band2_params.ba,band2_params.ba_err, bar, bar_err, poly_order
close,lun11
free_lun,lun11
;if poly_order lt 5 then stop
endfor

point1=wave[cont1_max_loc[0]]
point2=wave_tp[0]
point3=band2_wave[0]
point4=red_edge_wave
output2 = output_dir+'/band_definitions/' + objname + '_band_definitions.txt'
openw,lun12,output2,/get_lun
printf,lun12,'Object',' Point (i) ', ' Point (ii) ', ' Point(iii) ', ' Point(iv) '
printf,lun12,objname, point1,point2,point3,point4
close,lun12
free_lun,lun12


;Create a new structure that contains the averaged band centers and band depths over the three polynomial fits.
;Errors are re-estimated as the standard deviation of the three computed centers & depths and are adopted if
; that error is larger than the average of the centers and depths errors.
band_analysis_avg = {b1c:0.0, b1c_err:0.0, b1d:0.0, b1d_err:0.0, b1a:0.0, b1a_err:0.0,b1_slope:0.0, b1_slope_err:0.0,b2c:0.0, b2c_err:0.0, b2d:0.0, b2d_err:0.0, b2a:0.0, b2a_err:0.0,bar:0.0,bar_err:0.0}

b1c_stats = moment(band_analysis.b1c, sdev=b1c_sdev)
b1c_err_stats = moment(band_analysis.b1c_err, sdev=b1c_err_sdev)
band_analysis_avg.b1c = b1c_stats[0]
band_analysis_avg.b1c_err = b1c_err_stats[0]
;if b1c_sdev ge b1c_err_stats[0] then band_analysis_avg.b1c_err = b1c_sdev else band_analysis_avg.b1c_err = b1c_err_stats[0]


b1d_stats = moment(band_analysis.b1d, sdev=b1d_sdev)
b1d_err_stats = moment(band_analysis.b1d_err, sdev=b1d_err_sdev)
band_analysis_avg.b1d = b1d_stats[0]
band_analysis_avg.b1d_err = b1d_err_stats[0]
;if b1d_sdev ge b1d_err_stats[0] then band_analysis_avg.b1d_err = b1d_sdev else band_analysis_avg.b1d_err = b1d_err_stats[0]

b2c_stats = moment(band_analysis.b2c, sdev=b2c_sdev)
b2c_err_stats = moment(band_analysis.b2c_err, sdev=b2c_err_sdev)
band_analysis_avg.b2c = b2c_stats[0]
band_analysis_avg.b2c_err = b2c_err_stats[0]
;if b2c_sdev ge b2c_err_stats[0] then band_analysis_avg.b2c_err = b2c_sdev else band_analysis_avg.b2c_err = b2c_err_stats[0]


b2d_stats = moment(band_analysis.b2d, sdev=b2d_sdev)
b2d_err_stats = moment(band_analysis.b2d_err, sdev=b2d_err_sdev)
band_analysis_avg.b2d = b2d_stats[0]
band_analysis_avg.b2d_err = b2d_err_stats[0]
;if b2d_sdev ge b2d_err_stats[0] then band_analysis_avg.b2d_err = b2d_sdev else band_analysis_avg.b2d_err = b2d_err_stats[0]

;Since the other band analyses values do not depend on the polynomial fit, simply fill the band_analysis_avg structure
;with the 0th realization of band_analysis structure
band_analysis_avg.b1a = band_analysis[0].b1a
band_analysis_avg.b1a_err = band_analysis[0].b1a_err
band_analysis_avg.b1_slope = band_analysis[0].b1_slope
band_analysis_avg.b1_slope_err = band_analysis[0].b1_slope_err
band_analysis_avg.b2a = band_analysis[0].b2a
band_analysis_avg.b2a_err = band_analysis[0].b2a_err
band_analysis_avg.bar = band_analysis[0].bar
band_analysis_avg.bar_err = band_analysis[0].bar_err

;stop
sav_file = output_dir+objname+'.sav'
sav_avg_file = output_dir+objname+'_avg.sav'
save, band_analysis, filename=sav_file
save, band_analysis_avg, filename=sav_avg_file
;stop
print, 'All Done!'
done:
end
