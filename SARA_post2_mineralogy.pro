pro SARA_post2_mineralogy

cd, '~/SARA/output/master_output/
data_file = 'band_parameter_master_averages.txt'
ast_temps_path = '~/SARA/output/ast_temps/'
output_path = '~/SARA/output/mineralogy/'
if file_search(output_path,/test_directory) eq '' then file_mkdir, output_path

UserSym, [ -1, 1, 1, -1, -1 ], [ 1, 1, -1, -1, 1 ], /Fill
readcol, data_file, id,objname, atype,b1c, b1c_err,b1d,b1d_err,b1a,b1a_err,b1_slope,b1_slope_err,b2c,b2c_err,b2d,b2d_err,b2a,b2a_err,bar,bar_err,$
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

temp_file=''
read,'Please enter the filename for your asteroid temperatures.  [Enter 999 to skip temp. corrections] ',temp_file
file_string = file_search(ast_temps_path+temp_file+'*')
if file_string eq '' then temp_flag = temp_file else temp_flag = 'found'
case temp_flag of
  '999': begin
    b2c_tcor = b2c
    bar_tcor = bar
    end
  'found': begin
    readcol, file_string, temp
    n_temps = n_elements(temp)
    if n_temps ne n_ast then begin
      print, 'Number of Asteroid Temperatures and Number of Asteroids do not match.  Check Asteroid Temperature File'
      goto, done
      endif else begin
      bar_tcor = bar + 0.00075*temp - 0.23
      b2c_tcor = b2c + 0.06-0.0002*temp
    endelse
    end
    ELSE: begin
      print, 'File '+temp_file+' not found.  Please check filenames and paths.'
      goto, done
    end

endcase  

openw, lun3, 'band_parameter_master_temp_corrections.txt',/get_lun,width=300
printf, lun3,'ID   Object   Atype   BI_Center  BI_Cent_Err   BI_Depth   BI_Depth_Err BI_Area   BI_Area_Err  BI_Slope  BI_Slope_Err BII_Center  BIIC_Tcor BII_Cent_Err  BII_Depth   BII_Depth_Err  BII_Area   BII_Area_Err  BAR     BAR_Tcor      BAR_Err
for j=0,n_elements(id)-1,1 do begin
printf, lun3, strtrim((id[j]),2)+'  ',string(objname[j])+'  ',atype[j], b1c[j],b1c_err[j],b1d[j],b1d_err[j],b1a[j],b1a_err[j],b1_slope[j],b1_slope_err[j],b2c[j],b2c_tcor[j],b2c_err[j],b2d[j],b2d_err[j],b2a[j],b2a_err[j],bar[j],bar_tcor[j],bar_err[j]
endfor

;Create structures to hold mineralogical analysis
n_sample=n_elements(objname)
mineralogy_i={mineralogy_i, objname:'',Fs:0.0,Fs_gaffey:0.0,Wo:0.0,Wo_gaffey:0.0,En:0.0, En_gaffey:0.0, Fa:0.0,ol_ratio:0.0,ol_ratio_gaffey:0.0}
mineralogy=replicate(mineralogy_i,n_sample)
min_err_i={min_err_i, objname:'',Fs:0.0,Fs_gaffey:0.0,Wo:0.0,Wo_gaffey:0.0,En:0.0, En_gaffey:0.0, Fa:0.0,ol_ratio:0.0,ol_ratio_gaffey:0.0}
min_err=replicate(min_err_i,n_sample)
mineralogy[*].Fs = -999
mineralogy[*].Fs_gaffey = -999
mineralogy[*].Wo = -999
mineralogy[*].Wo_gaffey=-999
mineralogy[*].En = -999
mineralogy[*].En_gaffey = -999
mineralogy[*].Fa = -999
mineralogy[*].ol_ratio = -999
mineralogy[*].ol_ratio_gaffey = -999
for i_m=0,n_sample-1 do mineralogy[i_m].objname = objname[i_m]

;Dunn et al. (2010) - Applicable to S(IV) - Types
;siv_ind = where(objname eq '1089_Tama' or objname eq '2131_Mayall' or objname eq '3309_Brorfelde' or objname eq '4674_Pauling')
s_ind = where(bdmletter eq 'S' or bdmletter eq 'K')
ol_ratio = -0.242*bar_tcor[s_ind]+0.728  ;ol/(ol+px)
ol_ratio_err = 0.242*bar_err[s_ind]
mineralogy[s_ind].ol_ratio = ol_ratio
min_err[s_ind].ol_ratio = ol_ratio_err

Fa_ol = -1284.9*b1c[s_ind]^2+2656.5*b1c[s_ind]-1342.3
mineralogy[s_ind].Fa = Fa_ol
Fa_err = sqrt(b1c_err[s_ind]^2*(2*(-1284.9)*b1c[s_ind]+2656.5)^2)
min_err[s_ind].Fa=Fa_err

Fs_px = -879.1*b1c[s_ind]^2+1824.9*b1c[s_ind]-921.7
mineralogy[s_ind].Fs = Fs_px
Fs_err = sqrt(b1c_err[s_ind]^2*(2*(-879.1)*b1c[s_ind]+1824.9)^2)
min_err[s_ind].Fs=Fs_err
;stop

print, '----------------'
print, 'S-complex'
print, objname[s_ind]
print, ol_ratio
print, Fa_ol
print, Fs_px

;Burbine et al. (2007) - V-Type
vtype_ind = where(bdmletter eq 'V')
FsI_Burbine = 1023.4*b1c[vtype_ind]-913.82
FsII_Burbine = 205.86*b2c_tcor[vtype_ind]-364.3
WoI_Burbine = 396.13*b1c[vtype_ind]-360.55
WoII_Burbine = 79.905*b2c_tcor[vtype_ind]-148.3
EnI_Burbine = 100.0 - (FsI_Burbine+WoI_Burbine)
EnII_Burbine = 100.0 - (FsII_Burbine+WoII_Burbine)

mineralogy[vtype_ind].Fs = (FsI_Burbine+FsII_Burbine)/2
mineralogy[vtype_ind].Wo = (WoI_Burbine+WoII_Burbine)/2
mineralogy[vtype_ind].En = 100.0 - (mineralogy[vtype_ind].Fs + mineralogy[vtype_ind].Wo)
;stop
print, '----------------'
print, 'V-types'
print, objname[vtype_ind]
print, FsI_Burbine, FsII_Burbine, WoI_Burbine,WoII_Burbine,EnI_Burbine,EnII_Burbine
print, '---------------------------------'
print, mineralogy[vtype_ind].Fs, mineralogy[vtype_ind].Wo, mineralogy[vtype_ind].En
;;stop

;Gaffey et al. (2002) - Iterative Solution to all S or V-types
;First roughly determine if pyroxene or olivine dominate by checking
;if the Band I Center is > 1.0 (olivine-dominate) or < 1.0 (pyroxene-dominate)

;First pass uses the Equation 2c from Gaffey et al (2003) if V-type
;and uses Equation 2a if not V-type
aletter=strmid(atype,0,1)
pathway = strarr(n_sample)
for i=0,n_sample-1, 1 do begin
  converge='No'
  if aletter[i] eq 'V' then begin
    print, '---'
    
    ;Compute first-pass Wo Number and determine which Fs equation to appy
    wo = 418.9*b1c[i]-380.9
    if wo lt 11.0 then fs_flag = '3a'
    if wo ge 11.0 and wo lt 30.0 then fs_flag = '3b'
    if wo ge 30.0 and wo lt 45.0 then fs_flag = '3c'
    if wo ge 45 then fs_flag = '3d'
    pathway[i] = '2c '+ fs_flag+' '
    while converge eq 'No' do begin
      case fs_flag of
      
        '3a': begin
          fs = 268.2*b2c_tcor[i]-483.7 
          fs_err_nom = 5.0
          fs_err = 268.2*b2c_err[i]
         end
        '3b': begin
          fs = 57.5*b2c_tcor[i]-72.7
          fs_err_nom = 5.0
          fs_err = 57.5*b2c_err[i]
         end
        '3c': begin
          fs = -12.9*b2c_tcor[i]+45.9
          fs_err_nom = 4.0
          fs_err = 12.9*b2c_err[i]
        end
        '3d': begin
          fs = -118.0*b2c_tcor[i]+278.5
          fs_err_nom = 4.0
          fs_err = 118.0*b2c_err[i]
        end
     endcase
      
      ;Determine which Wo equation to apply
      if fs lt 10.0 then wo_flag = '2a'
      if fs ge 10.0 and fs lt 25.0 then wo_flag = '2b'
      if fs ge 25.0 and fs lt 50.0 then wo_flag = '2c'
      if fs gt 50.0 then begin
        print, 'Overestimated Fs: ', strtrim(string(fs,2))
        stop
      endif
      pathway[i]=pathway[i]+wo_flag+' '
      case wo_flag of
      
        '2a': begin
          wo_i = 347.9*b1c[i]-313.6
          wo_i_err_nom = 3.0
          wo_i_err = 347.9*b1c_err[i]
        end
        '2b': begin
          wo_i = 456.2*b1c[i]-416.9
          wo_i_err_nom = 3.0
          wo_i_err = 456.2*b1c_err[i]
        end
        '2c': begin
          wo_i = 418.9*b1c[i]-380.9
          wo_i_err_nom = 4.0
          wo_i_err = 418.9*b1c_err[i]
        end

      endcase
      ;Determine Fs equation and if convergence has occured
      if wo_i lt 11.0 then fs_i_flag = '3a'
      if wo_i ge 11.0 and wo lt 30.0 then fs_i_flag = '3b'
      if wo_i ge 30.0 and wo lt 45.0 then fs_i_flag = '3c'
      if wo_i ge 45 then fs_i_flag = '3d'
      pathway[i] = pathway[i]+fs_i_flag+' '
      if fs_i_flag eq fs_flag then begin
        print, pathway[i]
        print, wo, wo_i
;        stop
        mineralogy[i].Fs_gaffey = fs
        min_err[i].Fs_gaffey = fs_err
        mineralogy[i].Wo_gaffey = wo_i
        min_err[i].Wo_gaffey = wo_i_err
        converge = 'Yes'
        endif else fs_flag = fs_i_flag
        
    endwhile
 
  endif else begin
    
    print, '!!!'
    
    ;Compute first-pass Wo Number and determine which Fs equation to appy
    wo = 347.9*b1c[i]-313.6
    if wo lt 11.0 then fs_flag = '3a'
    if wo ge 11.0 and wo lt 30.0 then fs_flag = '3b'
    if wo ge 30.0 and wo lt 45.0 then fs_flag = '3c'
    if wo ge 45 then fs_flag = '3d'
    pathway[i] = '2a '+ fs_flag+' '
    while converge eq 'No' do begin
      case fs_flag of
      
        '3a': begin
          fs = 268.2*b2c_tcor[i]-483.7 
          fs_err = 5.0
         end
        '3b': begin
          fs = 57.5*b2c_tcor[i]-72.7
          fs_err = 5.0
         end
        '3c': begin
          fs = -12.9*b2c_tcor[i]+45.9
          fs_err = 4.0
        end
        '3d': begin
          fs = -118.0*b2c_tcor[i]+278.5
          fs_err = 4.0
        end
     endcase
      
      ;Determine which Wo equation to apply
      if fs lt 10.0 then wo_flag = '2a'
      if fs ge 10.0 and fs lt 25.0 then wo_flag = '2b'
      if fs ge 25.0 and fs lt 50.0 then wo_flag = '2c'
      if fs gt 50.0 then begin
        print, 'Overestimated Fs: ', strtrim(string(fs,2))
        fs = 50.
        ;converge_flag='Yes'
      endif
      pathway[i]=pathway[i]+wo_flag+' '
      case wo_flag of
      
        '2a': begin
          wo_i = 347.9*b1c[i]-313.6
          wo_i_err = 3.0
        end
        '2b': begin
          wo_i = 456.2*b1c[i]-416.9
          wo_i_err = 3.0
        end
        '2c': begin
          wo_i = 418.9*b1c[i]-380.9
          wo_i_err = 4.0
        end

      endcase
      ;Determine Fs equation and if convergence has occured
      if wo_i lt 11.0 then fs_i_flag = '3a'
      if wo_i ge 11.0 and wo lt 30.0 then fs_i_flag = '3b'
      if wo_i ge 30.0 and wo lt 45.0 then fs_i_flag = '3c'
      if wo_i ge 45 then fs_i_flag = '3d'
      pathway[i] = pathway[i]+fs_i_flag+' '
      if fs_i_flag eq fs_flag then begin
        print, pathway[i]
        print, objname[i]
        print, wo, wo_i
;        stop
        mineralogy[i].Fs_gaffey = fs
        mineralogy[i].Wo_gaffey = wo_i
        converge = 'Yes'
        endif else fs_flag = fs_i_flag
        
    endwhile
    
  endelse
  mineralogy[i].en_gaffey = 100.0 - (mineralogy[i].fs_gaffey + mineralogy[i].wo_gaffey)
  mineralogy[i].ol_ratio_gaffey = 1 - ( (0.417*bar_tcor[i]) + 0.052 )
endfor

;!p.font=-1
;window, 0, title='BAR vs BI Center'
;plot, [0,0],[0,0],xr=[0.0,3.0],yr=[0.9,1.15],xs=1,ys=1, xtitle='BAR',$
;  ytitle='BI Center (um)',charsize=2,psym=symtype[0],symsize=1.5
;for i=0,n_ast-1 do oplot, [bar[i],bar[i]],[b1c[i],b1c[i]],psym=symtype[i],symsize=1.5
;for j=0,n_ast-1 do xyouts, bar[j]+0.05,b1c[j]-0.007,objname[j],charsize=1.3
;oploterror,bar,b1c,bar_err,b1c_err,psym=3
;;Manually place Huntress
;;oplot, [0.0,0.0], [1.06,1.06], psym=6,symsize=1.5
;;xyouts, 0.05, 1.05, '7225_Huntress', charsize=1.3
;
;window, 1, title='BI Center vs BII Center'
;plot, [0,0],[0,0],xr=[1.86,2.08],yr=[0.91,0.96],xs=1,ys=1, xtitle='BI Center (um)',$
;  ytitle='BII Center (um)',charsize=2,psym=symtype[0],symsize=1.5
;for i=0,n_ast-1 do begin
;  if bdmletter[i] eq 'V' then begin
;    oplot, [b2c[i],b2c[i]],[b1c[i],b1c[i]],psym=symtype[i],symsize=1.5
;    xyouts, b2c[i]+0.003,b1c[i]-0.007,id[i],charsize=1.3
;    oploterror, b2c[i], b1c[i], b2c_err[i],b1c_err[i],psym=3
;  endif
;endfor
;
;window, 2, title='BAR vs BII Center'
;plot, [0,0],[0,0],xr=[0.7,3.0],yr=[1.86,2.08],xs=1,ys=1, xtitle='BAR',$
;  ytitle='BII Center (um)',charsize=2,psym=symtype[0],symsize=1.5
;for i=0,n_ast-1 do begin
;  if bdmletter[i] eq 'V' then begin
;    oplot, [bar[i],bar[i]],[b2c[i],b2c[i]],psym=symtype[i],symsize=1.5
;    xyouts, bar[i]-0.025,b2c[i]+0.003,id[i],charsize=1.3
;    oploterror, bar_err[i], b2c[i], bar_err[i],b2c_err[i],psym=3
;  endif
;endfor

;window, 3, title='Fs vs Wo - V-types'
;plot, [0,0],[0,0],xr=[8,60],yr=[0,17],xs=1,ys=1, xtitle='Fs (mol%)',$
;  ytitle='Wo (mol %)',charsize=2,psym=symtype[0],symsize=1.5
;for i=0,n_ast-1 do begin
;  if bdmletter[i] eq 'V' then begin
;    oplot, [mineralogy.Fs[i],mineralogy.Fs[i]],[mineralogy.Wo[i],mineralogy.Wo[i]],psym=symtype[i],symsize=1.5
;    xyouts, mineralogy.Fs[i]-1,mineralogy.Wo[i]+1,id[i],charsize=1.3
;;    oploterror, bar_err[i], b2c_tcor[i], bar_err[i],b2c_err[i],psym=3
;  endif
;endfor

;mineralogy_dir = data_dir+'/mineralogy/'
mineralogy_opx_file = 'mineralogy_opx.txt'
mineralogy_ol_file = 'mineralogy_ol.txt'
openw,lun1, output_path+mineralogy_opx_file,/get_lun,width=300
openw,lun2, output_path+mineralogy_ol_file,/get_lun, width=300
if temp_flag eq '999' then begin
  printf,lun1,'Temperature corrections to BAR and B2C have not been applied. Please see Sanchez et al. (2012)'
  printf,lun1,''
  printf,lun2,'Temperature corrections to BAR and B2C have not been applied. Please see Sanchez et al. (2012)'
  printf,lun2,''
endif
printf, lun1, 'Object  ', '  BDM_Taxon ', '     Fs   ', '    Wo   ', '  En   ', ' Fs(Gaffey) ', ' Wo(Gaffey) ', ' En(Gaffey) '
for ii=0, n_sample-1 do printf, lun1, objname[ii]+'  ', atype[ii]+'  ', mineralogy[ii].Fs, mineralogy[ii].Wo, mineralogy[ii].En,  mineralogy[ii].Fs_gaffey,mineralogy[ii].Wo_gaffey, mineralogy[ii].En_gaffey
printf, lun2, 'Object  ', '  BDM_Taxon ', '     Fa    ', '   Ol.Ratio  ', ' Ol.Ratio_G'
for jj=0, n_sample-1 do printf, lun2, objname[jj]+'  ', atype[jj]+'  ', mineralogy[jj].Fa, mineralogy[jj].ol_ratio, mineralogy[jj].ol_ratio_gaffey
close,/all
free_lun, lun1,lun2,lun3
print, 'All Done!'
done:
print, 'Exiting Program'
end
