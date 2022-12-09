pro resolve_beams

;before this:
;.r energy_ratio_imaging_new 

;First load in files, then make the nice three-method imaging figure (this also loads in needed functions)
foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files, new=2

;Use "new" option to use newer energy ratios method
;energy_ratio_imaging, files=event_files, erange=[19.,26.], video=0, nicefigure=1
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=0, nicefigure=1

files=event_files
erange=[17.,30.]


her1 = ratio_based_image(file=files[47], erange=erange, macrobins=1)
her2 = ratio_based_image(file=files[48], erange=erange, macrobins=1)
her3 = ratio_based_image(file=files[49], erange=erange, macrobins=1)
her4 = ratio_based_image(file=files[50], erange=erange, macrobins=1)

;_2files

cgPS_Open, 'resolve_beams_plot.ps', /nomatch
	
;cgLoadct, 8, NColors=bigmax, /reverse
cgLoadct, 1, /reverse

positions = cgLayout([2,2], OYMargin=[0,0], YGap=5, XGap=5, OXMargin=[5,0], aspect=0.85);, unit=0.5)
;couldn't get the layout the way I wanted, so adding a correction.
shift_pos = 0.175


x1=2100
x2=2279-60
y1=1440+60
y2=1619
;p[2] = p[2]-shift_pos

print, total(her2)		
her1 = her1[x1:x2, y1:y2]
her2 = her2[x1:x2, y1:y2]
her3 = her3[x1:x2, y1:y2]
her4 = her4[x1:x2, y1:y2]
print, total(her2)



;---------------------------------------
p = positions[*,0]

axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
				XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
				YTicks:(y2-y1)/60}
				



cgImage, her2*1000, position=p, /Axes, /NoErase, charsize=0.6, AXKEYWORDS=axis_format, stretch=1, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
			ytitle='Al-Side', xtitle='Pt-Side', title='Position 1'
			
print, 'MIN and MAX:', min(her2*1000), max(her2*1000)

num=0
for i=0, 4 do begin
	print, num
	cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor

num=0
for i=0, 4 do begin
	cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor

;---------------------------------------
p = positions[*,1]

axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
				XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
				YTicks:(y2-y1)/60}
				



cgImage, her3*1000, position=p, /Axes, /NoErase, charsize=0.6, AXKEYWORDS=axis_format, stretch=1, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
			ytitle='Al-Side', xtitle='Pt-Side', title='Position 2'
			
print, 'MIN and MAX:', min(her3*1000), max(her3*1000)

num=0
for i=0, 4 do begin
	print, num
	cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor

num=0
for i=0, 4 do begin
	cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor


;---------------------------------------
p = positions[*,2]

axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
				XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
				YTicks:(y2-y1)/60}
				



cgImage, her4*1000, position=p, /Axes, /NoErase, charsize=0.6, AXKEYWORDS=axis_format, stretch=1, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
			ytitle='Al-Side', xtitle='Pt-Side', title='Position 3'
			
print, 'MIN and MAX:', min(her4*1000), max(her4*1000)

num=0
for i=0, 4 do begin
	print, num
	cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor

num=0
for i=0, 4 do begin
	cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor

;---------------------------------------
p = positions[*,3]

cgLoadct, 3, /reverse

axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
				XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
				YTicks:(y2-y1)/60}
				



cgImage, (her4+her2)*1000, position=p, /Axes, /NoErase, charsize=0.6, AXKEYWORDS=axis_format,  stretch=1, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
			ytitle='Al-Side', xtitle='Pt-Side', title='Position 1+3 Overplot'
			
print, 'MIN and MAX:', min((her4+her2)*1000), max((her4+her2)*1000)

num=0
for i=0, 4 do begin
	print, num
	cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor

num=0
for i=0, 4 do begin
	cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
	num+=60
endfor



cgPS_Close
	
spawn, 'open resolve_beams_plot.ps'



stop

end