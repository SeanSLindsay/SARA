pro oplotbox, x1,x2,y1,y2,lin,thik,pcolor ;Plot an empty box 
;loadct,ct,/silent
oplot, [x1,x2],[y1,y1], linestyle=lin,thick=thik,color=pcolor
oplot, [x1,x2],[y2,y2], linestyle=lin,thick=thik,color=pcolor
oplot, [x1,x1],[y1,y2], linestyle=lin,thick=thik,color=pcolor
oplot, [x2,x2],[y1,y2], linestyle=lin,thick=thik,color=pcolor
end

pro plot_oc_boot
;Plot the OC Boot
loadct,17,/silent
oc_boot_x = [0.3,0.4,0.7,1.1,1.1,0.75,0.6,0.5]
oc_boot_y = [1.0,0.94,0.92,0.925,0.94,0.96,0.98,1.03]
oc_boot_n = n_elements(oc_boot_x)
for i=0,oc_boot_n-2, 1 do begin
  oplot, [oc_boot_x[i],oc_boot_x[i+1]],[oc_boot_y[i],oc_boot_y[i+1]],thick=3,color=240
endfor
oplot, [oc_boot_x[oc_boot_n-1],oc_boot_x[0]],[oc_boot_y[oc_boot_n-1],oc_boot_y[0]],thick=3,color=240
xyouts, 1.5,1.065,'Ordinary Chondrites',charsize=1.5,color=240
end

pro plot_ureilite_zone
;Plot Ureilite Zone
loadct,16,/silent
ur_x = [0.345,0.125,0.3,0.64]
ur_y = [0.975,0.96,0.91,0.9235]
for j=0,n_elements(ur_x)-2 do oplot, [ur_x[j],ur_x[j+1]],[ur_y[j],ur_y[j+1]],thick=3,color=250
xyouts,1.5,1.08,'Ureilites',charsize=1.5,color=250
end

pro plot_meteorite_boxes
loadct,17,/silent
oplotbox, 0.0,0.2,1.05,1.1,0,3,170      ;S(I) - Olivine Meteorites
oplotbox, 1.5,2.7,0.92,0.96,0,3,50     ;BA - Basaltic Achondrites
oplotbox, 0.65,1.5,0.91,0.925,0,3,120   ;Primitive Achondrites
xyouts, 1.5,1.125,'Olivine Meteorites',charsize=1.5,color=170
xyouts, 1.5,1.11, 'Basaltic Achondrites - HED', charsize=1.5, color=50
xyouts, 1.5,1.095,'Lodaranites/Acapulcoites', charsize=1.5,color=120
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro SARA_post3_analogs;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if file_search('~/SARA/plots/',/test_directory) eq '' then file_mkdir, '~/SARA/plots/'

cd,'~/SARA/plots/'

set_plot,'ps'
device, file = 'sara_meteorite_analog_zones.eps',/encapsulated,/color, xsize=9,ysize=9/1.5,font_size=6

mu_setup = "155B
mu = '!9' + string(mu_setup) + '!X'

output_dir = '~/Research/asteroids/binary_project/nir_band_analysis/output/v06/mineralogy/'
data_dir = '~/SARA/output/master_output/'
data_file = 'band_parameter_master_temp_corrections.txt'
if file_search(data_dir+data_file) eq '' then begin
  print, 'Post-processing product 2: band_parameter_master_temp_corrections.txt does not exist.'
  print, 'Please run sara_post2_mineralogy.pro before proceeding'
  goto,done
endif else print, 'Locating band_parameter_master_temp_corrections.txt'

UserSym, [ -1, 1, 1, -1, -1 ], [ 1, 1, -1, -1, 1 ], /Fill
readcol, data_dir+data_file, id,objname, atype,b1c, b1c_err,b1d,b1d_err,b1a,b1a_err,b1_slope,b1_slope_err,b2c,b2c_tcor,b2c_err,b2d,b2d_err,b2a,b2a_err,bar,bar_tcor,bar_err,$
  format='A,A,A'

n_ast = n_elements(objname)
symtype=intarr(n_ast)
atype = strtrim(atype,2)
bdmletter=strmid(atype,0,1)
for i=0,n_ast-1 do begin
  case atype[i] of
    'Sw': symtype[i]=8  ;Square
    'Sr': symtype[i]=8  ;Square
    'Sq': symtype[i]=8  ;Square
    'S' : symtype[i]=8  ;Square
    'Vw': symtype[i]=5  ;Triangle
    'V' : symtype[i]=5  ;Triangle
    'K' : symtype[i]=4  ;Diamond
 
 endcase
endfor


!p.font=1
!p.charsize=1
!p.thick=2
!p.charthick=1
atypes=['S','V']
asyms=[8,5]
;window, 0, title='BAR vs BI Center'
loadct,39,/silent
plot, [0,0],[0,0],xr=[-0.05,3.0],yr=[0.9,1.15],xs=1,ys=1, xtitle='Band Area Ratio',$
  ytitle='BI Center ('+mu+'m)',psym=symtype[0],symsize=1,charsize=2,$
  xthick=2,ythick=2

plot_meteorite_boxes
plot_ureilite_zone
plot_oc_boot

loadct,39,/silent
for i=0,n_ast-1 do oplot, [bar_tcor[i],bar_tcor[i]],[b1c[i],b1c[i]],psym=symtype[i],symsize=1
;for j=0,n_ast-1 do xyouts, bar_tcor[j]+0.05,b1c[j]-0.007,id[j];,charsize=2
oploterror,bar_tcor,b1c,bar_err,b1c_err,psym=3
for k=0,n_elements(atypes)-1 do oplot, [1.0,1.0],[1.13-0.015*k,1.13-0.015*k],psym=asyms[k],symsize=1,thick=1
for k=0,n_elements(atypes)-1 do xyouts, 1.075,1.1275-0.015*k, atypes[k]
;Manually place Huntress

;cgarrow, 1.3,0.96,1.3,0.93,/data,thick=4,clip=2

device,/close
;stop
print, 'All Done!'
done:
end
