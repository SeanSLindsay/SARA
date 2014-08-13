pro oplotbox, x1,x2,y1,y2,lin,thik,pcolor ;Plot an empty box 
;loadct,17,/silent
oplot, [x1,x2],[y1,y1], linestyle=lin,thick=thik,color=pcolor
oplot, [x1,x2],[y2,y2], linestyle=lin,thick=thik,color=pcolor
oplot, [x1,x1],[y1,y2], linestyle=lin,thick=thik,color=pcolor
oplot, [x2,x2],[y1,y2], linestyle=lin,thick=thik,color=pcolor
end

pro oplotellipse,x0,y0,a,b,theta,lin,thik,pcolor ; plot an empty ellipse
; x,y center coordinates
; a,b semi-major/minor axes
; theta clockwise rotation
loadct,0,/silent
ang=findgen(1400)*!dtor
thep=theta*!dtor
xp=a*cos(ang) & yp=b*sin(ang)
x=(xp*cos(thep))-(yp*sin(thep)) & y=((xp*sin(thep))+(yp*cos(thep)))/8.5
oplot,x+x0,y+y0,linestyle=lin,thick=thik,color=pcolor
end

pro plot_S4
;Plot Das Boot
loadct,17,/silent
das_boot_x = [0.25,0.41,0.62,1.16,1.15,0.72,0.52,0.4]
das_boot_y = [1.002,0.945,0.92,0.925,0.94,0.96,0.981,1.029]
das_boot_n = n_elements(das_boot_x)
for i=0,das_boot_n-2, 1 do begin
  oplot, [das_boot_x[i],das_boot_x[i+1]],[das_boot_y[i],das_boot_y[i+1]],$
         thick=2,color=240
endfor
oplot, [das_boot_x[das_boot_n-1],das_boot_x[0]],[das_boot_y[das_boot_n-1],$
                                                 das_boot_y[0]],thick=2,color=240
end

pro plot_S1BA
loadct,17,/silent
oplotbox, 0.0,0.2,1.05,1.1,0,2,170 ; S(I) - Olivine Meteorites
oplotbox, 1.5,2.7,0.92,0.96,0,2,50 ; BA - Basaltic Achondrites
end

pro plot_S2toS7
            ; x0,y0,a,b,theta,lin,thik,pcolor
oplotellipse,0.4,1.046,0.24,0.13,-30,0,2,0 ; S2
oplotellipse,0.26,0.965,0.32,0.08,-73.,0,2,0 ; S3 (lower)
oplotellipse,0.59,1.002,0.24,0.11,-74.,0,2,0 ; S3 (upper)
oplotellipse,1.268,0.957,0.379,0.134,-10.5,0,2,0 ; S5
oplotellipse,1.33,0.92,0.16,0.16,0.,0,2,0 ; S6
oplotellipse,2.1,0.91,0.6,0.07,0.,0,2,0 ; S7

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro SARA_post4_ssubtypes;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; DESCRIPTION:
;     re-creates the Gaffey (1993) S subtype plot (BI center/BAR) and
;     output automatically to an eps file
; CREATED BY: Eric MacClennan
; ADAPTED FOR SARA: Sean Lindsay

if file_search('~/SARA/plots/',/test_directory) eq '' then file_mkdir, '~/SARA/plots/'
cd,'~/SARA/plots/'

set_plot,'ps'
device, file = 'sara_meteorite_ssubtypes.eps',/encapsulated,/color, xsize=9,ysize=9/1.5,font_size=6

mu_setup = "155B
mu = '!9' + string(mu_setup) + '!X'

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
!p.thick=1
!p.charthick=1
thik=2
atypes=['S','V']
asyms=[8,5]
;window, 0, title='BAR vs BI Center'
loadct,39,/silent
plot, [0,0],[0,0],xr=[-0.05,3.0],yr=[0.9,1.15],xs=1,ys=1, xtitle='Band Area Ratio',$
  ytitle='BI Center ('+mu+'m)',psym=symtype[0],symsize=1,charsize=2,$
  xthick=2,ythick=2

;Add Meteorite Analog Zones
plot_S1BA
plot_S4
plot_S2toS7

; plot labels
xyouts, 0.3,1.11,'S(I)',charsize=2,color=0,charthick=4
xyouts, 0.5,1.07,'S(II)',charsize=2,color=0,charthick=4
xyouts, 0.7,1.04,'S(III)',charsize=2,color=0,charthick=4
xyouts, 0.85,0.99,'S(IV)',charsize=2,color=0,charthick=4
xyouts, 1.19,0.956,'S(V)',charsize=2,color=0,charthick=4
xyouts, 1.65,0.98,'S(VI)',charsize=2,color=0,charthick=4
xyouts, 2.25,0.97,'S(VII)',charsize=2,color=0,charthick=4
xyouts, 2.05,0.937, 'BA', charsize=2, color=0,charthick=4

; plot arrows
cgarrow,0.29,1.109,0.15,1.09,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S1
cgarrow,0.547,1.067,0.42,1.055,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S2
cgarrow,0.69,1.04,0.63,1.01,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S3
cgarrow,0.69,1.04,0.26,0.975,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S3
cgarrow,0.9,0.987,0.648,0.96,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S4
cgarrow,1.7,0.977,1.4,0.93,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S6
cgarrow,2.3,0.967,2.25,0.912,/solid,thick=2,/data,hsize=!D.X_SIZE/95.,hthick=0.5 ; S7

; plot ol/px mixing line
bar1=findgen(50)*.366/49.-0.02
oplot,bar1,0.903+0.1657*exp(-(bar1+0.02861)^2./(2*0.3695^2.)),thick=6.1
bar2=findgen(50)*.874/49.+.346 & oplot,bar2,0.9142+0.3249*exp(-3.776*bar2),thick=6.1
bar3=findgen(50)*0.48/49.+1.22 & oplot,bar3,0.9606-0.05269*bar3+0.01419*bar3^2.,thick=6.1
bar4=findgen(50)*1.3/49.+1.7 & oplot,bar4,0.9177-0.003324*bar4,thick=6.1

for i=0,n_ast-1 do oplot, [bar_tcor[i],bar_tcor[i]],[b1c[i],b1c[i]],psym=symtype[i],symsize=1
;for j=0,n_ast-1 do xyouts, bar_tcor[j]+0.05,b1c[j]-0.007,id[j];,charsize=2
oploterror,bar_tcor,b1c,bar_err,b1c_err,psym=3
for k=0,n_elements(atypes)-1 do oplot, [1.0,1.0],[1.13-0.015*k,1.13-0.015*k],psym=asyms[k],symsize=1,thick=1
for k=0,n_elements(atypes)-1 do xyouts, 1.075,1.1275-0.015*k, atypes[k]

device,/close
set_plot,'x'
print, 'All Done!'
done:
end
