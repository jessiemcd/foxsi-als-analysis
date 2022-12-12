pro polarization_check_Lscan


;Goal: take different files from during the L scan, several hours or more apart in time, which each have 
;significant signal in the same strip. Then, plot the strip-spectra in that strip for all files 
;and see if there are any changes over time. 

;All the relevant (illuminated) strips from the overnight L-scan:
relevant_strips_p = indgen(10)+33
relevant_as_p = 2
relevant_strips_n = indgen(7)+25
relevant_as_n = 1

foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=2


print, ''
print, ''
print,'During N-Scan:'
print, ''
print, ''

files=event_files_n
for j=0, n_elements(files)-1 do begin
	restore, files[j]
	print, 'File: ', files[j]
	data = eventwise_spectra.EVENTDATA_SPEC
	;print, size(data)
	;print, j
	
	print, 'Counts in P Strips'
	for i=0, n_elements(relevant_strips_p)-1 do begin
		strip = relevant_strips_p[i]
		events = data[relevant_as_p,strip,*]
		events = events[where(events NE 0.)]
		print, 'Strip: ', strip, ' AISC: ', relevant_as_p, ' Number of Events: ', n_elements(events)
	endfor
	print,''
endfor

print, ''
print, ''
print,'During P-Scan:'
print, ''
print, ''

;FROM THE ABOVE, WE HAVE CONFIRMED THAT THROUGHOUT THE AL-SIDE SCAN, THE BEAM IS IN ASIC 2 STRIP 42.

files=event_files_p
for j=0, n_elements(files)-1 do begin
	restore, files[j]
	print, 'File: ', files[j]
	data = eventwise_spectra.EVENTDATA_SPEC
	;print, size(data)
	;print, j
	
	print, 'Counts in N Strips'
	for i=0, n_elements(relevant_strips_n)-1 do begin
		strip = relevant_strips_n[i]
		events = data[relevant_as_n,strip,*]
		events = events[where(events NE 0.)]
		print, 'Strip: ', strip, ' AISC: ', relevant_as_n, ' Number of Events: ', n_elements(events)
		;if j EQ 0 then print, strip, relevant_as_n, ' Number: ', n_elements(events)
		;if j EQ 90 then print, strip, relevant_as_n, ' Number: ', n_elements(events)
		;if j EQ 180 then print, strip, relevant_as_n, ' Number: ', n_elements(events)
	endfor
	
	print,''
	print, 'Counts in P Strips'
	for i=0, n_elements(relevant_strips_p)-1 do begin
		strip = relevant_strips_p[i]
		events = data[relevant_as_p,strip,*]
		events = events[where(events NE 0.)]
		print, 'Strip: ', strip, ' AISC: ', relevant_as_p, ' Number of Events: ', n_elements(events)
	endfor
	print, ''
	
	
	
endfor

print, ''
print, ''


;FROM THE ABOVE, WE HAVE CONFIRMED THAT DURING THE PT-SIDE SCAN, THE BEAM STARTS IN ASIC 1 STRIP 26
;AND GRADUALLY MOVES INTO ASIC 1 STRIP 25.

;ADDITIONALLY, THE P-SIDE SCAN STARTS IN STRIP 33 AND PROCEDES THROUGH THE STRIPS UNTIL IT REACHES 42.



;APPROACHING THE P-SIDE (DURING N SCAN, WHERE ILLUMINATED P-STRIP IS CONSISTENT) FIRST
file1 = event_files_n[0]
file2 = event_files_n[50]
file3 = event_files_n[100]
file4 = event_files_p[0]

print, file1
print, file2
print, file3

restore, file1
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_p,42,*]
;print, n_elements(events)
events1 = events[where(events NE 0.)]
;print, n_elements(events1)
hist1 = histogram(events1, binsize=0.25, min=0, max=35, locations=bins1)


restore, file2
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_p,42,*]
;print, n_elements(events)
events2 = events[where(events NE 0.)]
;print, n_elements(events2)
hist2 = histogram(events2, binsize=0.25, min=0, max=35, locations=bins2)

restore, file3
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_p,42,*]
;print, n_elements(events)
events3 = events[where(events NE 0.)]
;print, n_elements(events3)
hist3 = histogram(events3, binsize=0.25, min=0, max=35, locations=bins3)

restore, file4
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_p,33,*]
;print, n_elements(events)
events4 = events[where(events NE 0.)]
;print, n_elements(events3)
hist4 = histogram(events4, binsize=0.25, min=0, max=35, locations=bins4)


print, ''

;APPROACHING THE N-SIDE (DURING P SCAN) NEXT (WHERE THE ILLUMINATED N-STRIP IS (MORE) CONSISTENT)

file1_n = event_files_p[0]
file2_n = event_files_p[90]
file3_n = event_files_p[180]

restore, file1_n
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_n,25,*]
;print, n_elements(events)
events1 = events[where(events NE 0.)]
;print, n_elements(events1)
hist1_25 = histogram(events1, binsize=0.25, min=0, max=35, locations=bins1_25)

events = data[relevant_as_n,26,*]
;print, n_elements(events)
events1 = events[where(events NE 0.)]
;print, n_elements(events1)
hist1_26 = histogram(events1, binsize=0.25, min=0, max=35, locations=bins1_26)


restore, file2_n
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_n,25,*]
;print, n_elements(events)
events2 = events[where(events NE 0.)]
;print, n_elements(events2)
hist2_25 = histogram(events2, binsize=0.25, min=0, max=35, locations=bins2_25)

events = data[relevant_as_n,26,*]
;print, n_elements(events)
events2 = events[where(events NE 0.)]
;print, n_elements(events2)
hist2_26 = histogram(events2, binsize=0.25, min=0, max=35, locations=bins2_26)

restore, file3_n
data = eventwise_spectra.EVENTDATA_SPEC
events = data[relevant_as_n,25,*]
;print, n_elements(events)
events3 = events[where(events NE 0.)]
;print, n_elements(events3)
hist3_25 = histogram(events3, binsize=0.25, min=0, max=35, locations=bins3_25)

events = data[relevant_as_n,26,*]
;print, n_elements(events)
events3 = events[where(events NE 0.)]
;print, n_elements(events3)
hist3_26 = histogram(events3, binsize=0.25, min=0, max=35, locations=bins3_26)


;=========================================================================================
;=========================================================================================


popen, 'polarization_check.ps', $
		xsi=8, ysi=10
!Y.margin=4.
!X.margin=4.
ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
!p.multi=[0,1,4]
loadct, 2

labels = ['File 1 - Minute 0', 'File 50 - Minute ~140 (2 hr 20 min)', 'File 100 - Minute ~240 (4 hr)'];, $
			;'Minute 0 of Previous Branch, different strip (33)']
colors = [0,40, 100, 150]

;=========================================================================================
;=========================================================================================


plot, bins1, hist1, thick=th, /ylog, psym=10, yrange=[1, max(hist1)*1.5], ystyle=1, $
		title = 'Pt-side Spectra, all events in strip 42 (no event selection), three positions in Al-side scan'
oplot, bins2, hist2, color=40, thick=th, psym=10
oplot, bins3, hist3, color=100, thick=th, psym=10
;oplot, bins4, hist4, color=150, thick=th, psym=10

al_legend, labels , textcol=colors, box=1, /right, /top, /clear, charsi=0.6

;;previously commented section

labels = ['File 1 - Minute 0', 'File 90 - Minute ~180 (3 hr)', 'File 180 - Minute ~360 (6 hr)']

plot, bins1_25, hist1_25, thick=th, /ylog, psym=10, yrange=[1, max(hist1)*1.5], ystyle=1, $
	title = 'Al-side Spectra, all events in strip 25 (no event selection), three positions in Pt-side scan'
oplot, bins2_25, hist2_25, color=40, thick=th, psym=10
oplot, bins3_25, hist3_25, color=100, thick=th, psym=10

al_legend, labels , textcol=colors, box=1, /right, /top, /clear, charsi=0.6

plot, bins1_26, hist1_26, thick=th, /ylog, psym=10, yrange=[1, max(hist1)*1.5], ystyle=1, $
	title = 'Al-side Spectra, all events in strip 26 (no event selection), three positions in Pt-side scan'
oplot, bins2_26, hist2_26, color=40, thick=th, psym=10
oplot, bins3_26, hist3_26, color=100, thick=th, psym=10

al_legend, labels , textcol=colors, box=1, /right, /top, /clear, charsi=0.6

;;previously commented section end

;
;
;labels = ['First L-scan position 20190419 20:59', 'Efficiency scan midpoint (13) 20190418 18:58', $
;			'20190418 19:05', '20190418 19:11']
;colors = [0,40,100,150]
;
;plot, binsL, histL, thick=th, /ylog, psym=10, yrange=[1, max(hist1)*2], ystyle=1, $
;		title = 'Pt-side Spectra, all events in strip 33 (no event selection), positions in two scans'
;oplot, binsE, histE, color=40, thick=th, psym=10
;oplot, binsE5, histE5, color=100, thick=th, psym=10
;oplot, binsE7, histE7, color=150, thick=th, psym=10
;oplot, binsL, histL, thick=th, psym=10
;
;al_legend, labels , textcol=colors, box=1, /right, /top, /clear, charsi=0.6
;
;labels=['L-scan', 'Square scan']
;colors=[0,40]
;
;plot, binnpol, histnpol, thick=th, /ylog, psym=10, yrange=[1, max(histnpol)*3], ystyle=1, $
;		title = 'Al-side Spectra, all events in strip 25 (no event selection), positions in two scans'
;oplot, binnspol, histnspol, color=40, thick=th, psym=10
;
;al_legend, labels , textcol=colors, box=1, /right, /top, /clear, charsi=0.6


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
;DEVICE, /CLOSE
spawn, 'open polarization_check.ps'

;=========================================================================================
;=========================================================================================






stop
end

