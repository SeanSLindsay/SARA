pro SARA_post1_master

cd, '~/SARA/output/'
if file_search('~/SARA/output/master_output/',/test_directory) eq '' then begin
  print, 'ERROR: Missing master output directory
  print, 'Please run SARA before creating master output.'
  goto, done
endif
;results_dir = '~/Research/asteroids/binary_project/nir_band_analysis/output/v06_relab/H/2_40/'
filename = 'master_output/band_parameter_master.txt'
filename2 = 'master_output/band_parameter_master_averages.txt'
mas_result_files = findfile('*.txt',count=nfiles)

openw, lun1, filename,/get_lun,width=300
openw, lun2, filename2,/get_lun,width=300
printf, lun1,'ID   Object   Polyorder   Atype   BI_Center  BI_Cent_Err   BI_Depth   BI_Depth_Err BI_Area   BI_Area_Err  BI_Slope  BI_Slope_Err BII_Center   BII_Cent_Err  BII_Depth   BII_Depth_Err  BII_Area   BII_Area_Err  BAR         BAR_Err
printf, lun2,'ID   Object   Atype   BI_Center  BI_Cent_Err   BI_Depth   BI_Depth_Err BI_Area   BI_Area_Err  BI_Slope  BI_Slope_Err BII_Center   BII_Cent_Err  BII_Depth   BII_Depth_Err  BII_Area   BII_Area_Err  BAR         BAR_Err


b1c=(b1c_err=(b1d=(b1d_err=(b1a=(b1a_err=(b1s=(b1s_err=(b2c=(b2c_err=(b2d=(b2d_err=(b2a=(b2a_err=(bar=(bar_err=(poly_order=fltarr(nfiles)))))))))))))))))
id=lonarr(nfiles)
obj=(atype=strarr(nfiles))
b1c_avg=(b1c_err_avg=(b1d_avg=(b1d_err_avg=(b2c_avg=(b2c_err_avg=(b2d_avg=(b2d_err_avg=fltarr(nfiles/3))))))))
for i=0,nfiles-1,1 do begin
  readcol, mas_result_files[i], obj_idi,atypei, b1ci,b1c_erri,b1di,b1d_erri,b1ai,b1a_erri,b1si,b1s_erri,b2ci,b2c_erri,b2di,b2d_erri,b2ai,b2a_erri,bari,bar_erri,poly_orderi,format='A,A'
  object_id_split = strsplit(obj_idi,'_',/extract)
  id[i] = object_id_split[0]
  obj[i] = object_id_split[1]
  atype[i] = atypei
  b1c[i]=b1ci
  b1c_err[i]=b1c_erri
  b1d[i]=b1di
  b1d_err[i]=b1d_erri
  b1s[i]=b1si
  b1s_err[i]=b1s_erri
  b1a[i]=b1ai
  b1a_err[i]=b1a_erri
  b2c[i]=b2ci
  b2c_err[i]=b2c_erri
  b2d[i]=b2di
  b2d_err[i]=b2d_erri
  b2a[i]=b2ai
  b2a_err[i]=b2a_erri
  bar[i]=bari
  bar_err[i]=bar_erri
  poly_order[i]=poly_orderi
  if (i+1) mod 3 eq 0 then begin
    b1c_avg[((i+1)/3)-1] = mean(b1c[i-2:i])
    b1c_err_avg[((i+1)/3)-1] = mean(b1c_err[i-2:i])
    b1d_avg[((i+1)/3)-1] = mean(b1d[i-2:i])
    b1d_err_avg[((i+1)/3)-1] = mean(b1d_err[i-2:i])
    b2c_avg[((i+1)/3)-1] = mean(b2c[i-2:i])
    b2c_err_avg[((i+1)/3)-1] = mean(b2c_err[i-2:i])
    b2d_avg[((i+1)/3)-1] = mean(b2d[i-2:i])
    b2d_err_avg[((i+1)/3)-1] = mean(b2d_err[i-2:i])
  endif
endfor 
;stop
for j=0,nfiles-1,1 do begin
  printf, lun1, strtrim((id[j]),2)+'  ',string(obj[j])+'  ',atype[j], poly_order[j],b1c[j],b1c_err[j],b1d[j],b1d_err[j],b1a[j],b1a_err[j],b1s[j],b1s_err[j],b2c[j],b2c_err[j],b2d[j],b2d_err[j],b2a[j],b2a_err[j],bar[j],bar_err[j]
  if (j+1) mod 3 eq 0 then begin 
    printf, lun1,strtrim((id[j]),2)+'  ',string(obj[j])+'  ',atype[j],'        AVG  ',b1c_avg[((j+1)/3)-1],b1c_err_avg[((j+1)/3)-1],b1d_avg[((j+1)/3)-1],b1d_err_avg[((j+1)/3)-1],b1a[j],b1a_err[j],b1s[j],b1s_err[j],b2c_avg[((j+1)/3)-1],b2c_err_avg[((j+1)/3)-1],b2d_avg[((j+1)/3)-1],b2d_err_avg[((j+1)/3)-1],b2a[j],b2a_err[j],bar[j],bar_err[j]
    printf, lun2,strtrim((id[j]),2)+'  ',string(obj[j])+'  ',atype[j], b1c_avg[((j+1)/3)-1],b1c_err_avg[((j+1)/3)-1],b1d_avg[((j+1)/3)-1],b1d_err_avg[((j+1)/3)-1],b1a[j],b1a_err[j],b1s[j],b1s_err[j],b2c_avg[((j+1)/3)-1],b2c_err_avg[((j+1)/3)-1],b2d_avg[((j+1)/3)-1],b2d_err_avg[((j+1)/3)-1],b2a[j],b2a_err[j],bar[j],bar_err[j]
  endif
endfor
close, lun1,lun2
free_lun,lun1,lun2
print, 'All Done!'
done:
end
