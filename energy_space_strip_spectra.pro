pro energy_space_strip_spectra, singlestrip=singlestrip 

;GOALS:
;Make energy space spectra for every strip; fit peak locations

;badpar.pro and lclxtrem.pro should be in IDL path

default, singlestrip, 1
thold=2.

;Inputs
foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=event_files, new=2


;Pre-known relevant strips
relevant_strips_p = indgen(10)+33
relevant_as_p = 2
relevant_strips_n = indgen(7)+25
relevant_as_n = 1

;Defining strip arrays
for i=0, n_elements(relevant_strips_p)-1 do begin
	array_str = 'pstrip_'+strtrim(relevant_strips_p[i], 2)
	;print, array_str
	void = EXECUTE(array_str + '=[]')
endfor

for i=0, n_elements(relevant_strips_n)-1 do begin
	array_str = 'nstrip_'+strtrim(relevant_strips_n[i], 2)
	;print, array_str
	void = EXECUTE(array_str + '=[]')
endfor

doloop = 1

;cmns_p = []
;cmns_n = []

if doloop EQ 1 then begin


	for j=0, n_elements(event_files)-1 do begin
			restore, event_files[j]
			data = eventwise_spectra.EVENTDATA_SPEC
		
			;For each p-side strip, lets make an array of event energies (energy space) and
			;assign it to a strip-named variable.
			for i=0, n_elements(relevant_strips_p)-1 do begin
				strip = relevant_strips_p[i]
				events=[]
				
				;for each frame...
				for k=0, n_elements(data[0,0,*])-1 do begin
					vals = data[2,*,k]
					above = vals[where(vals GE thold, /NULL)]
					
					;If single strip is set, use spline to find a "threshold" channel value and
					;take only events which are the only above-threshold event in their frame
					if singlestrip EQ 1 then begin
						if n_elements(above) EQ 1 and above NE !NULL then begin
							;print, above
							mm = max(vals)
							if data[2,strip,k] EQ mm then events = [events, mm ]
						endif
					;if singlestrip is not set, add all values above 0 
					endif else begin
						events = [events, above] 
					endelse
				endfor
				;print, size(events)
				array_str = 'pstrip_'+strtrim(relevant_strips_p[i], 2)
				void = EXECUTE(array_str + '=['+array_str +', events]')
			endfor
			

		
			;repeat for n-side (see p-side documentation above)
			for i=0, n_elements(relevant_strips_n)-1 do begin
				strip = relevant_strips_n[i]
				events=[]
				for k=0, n_elements(data[0,0,*])-1 do begin
					vals = data[1,*,k]
					above = vals[where(vals GE thold, /NULL)]
					
					if singlestrip EQ 1 then begin
						if n_elements(above) EQ 1 and above NE !NULL then begin
							;print, above
							mm = max(vals)
							if data[1,strip,k] EQ mm then events = [events, mm ]
						endif
					;if singlestrip is not set, add all values above 0 
					endif else begin
						events = [events, above] 
					endelse
				endfor
				;events = data[*].data[1,strip]-cmn[*,1]
				array_str = 'nstrip_'+strtrim(relevant_strips_n[i], 2)
				void = EXECUTE(array_str + '=['+array_str +', events]')
			endfor	
	endfor 
	
endif

if doloop EQ 0 then print, 'Have not implemented a way to skip the loop'

;STOP

;=========================================================================================
;=========================================================================================

;Plot channel-space spectrum, and find peak locations (on each side, for each strip)

popen, 'strip_energy_spectrum.ps', $
	xsi=8, ysi=10
!Y.margin=4.
!X.margin=4.
ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
!p.multi=[0,2,4]
linecolors

;=========================================================================================
;=========================================================================================

	
binsize=0.1
coll=1
peak7s_p = []
peak21s_p = []
peak28s_p = []
;For each p-side strip...
for i=0, n_elements(relevant_strips_p)-1 do begin
	
	strip = relevant_strips_p[i]
	if coll EQ 5 then coll+=1
	array_str = 'pstrip_'+strtrim(strip, 2)
	void = EXECUTE('striparray = '+ array_str)
	;make and plot channel-space spectrum histogram
	hist = histogram(striparray, binsize=binsize, min=0, max=50, locations=sbins)
	plot, sbins, hist, thick=th, title=array_str, /ylog, psym=10, yrange=[1, max(hist)*1.5], $
			ystyle=1, ytitle='Counts', xtitle='Energy (keV)'
	oplot, sbins, hist, color=coll, thick=th, psym=10
	
	;via spline, estimate channels corresponding to 4, 7, 21 keV (two lines + commonly used threshold)
	;and plot for later comparison with found peaks. 
	val4 = 4
	val7 = 7
	val21 = 21
	val28 = 28
	oplot, [val7, val7], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
    oplot, [val21, val21], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
    oplot, [val28, val28], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [val4, val4], [0.1, max(hist)*2], linestyle=2, color=coll+2, thick=th
	
	;find maximum
	maxbin = sbins[where(hist EQ max(hist))]
	;oplot, [maxbin, maxbin], [0.1, max(hist)*2], linestyle=2, thick=2

	;find other maxima and mark them on plot
	maxes = lclxtrem(hist, 5, /maxima)
	maxes = maxes[where(sbins[maxes] LT 500)]
	maxes = maxes[where(sbins[maxes] GE maxbin[0])]
	;for m=0, n_elements(maxes)-1 do begin
	;	oplot, [sbins[maxes[m]], sbins[maxes[m]]], [0.1, max(hist)*2], linestyle=2
	;endfor

	iv=3

	;make data interval for fitting first peak + do gaussian fit
	fitbins = sbins[where(sbins GE maxbin[0]-iv)]
	fitbins2 = fitbins[where(fitbins LE maxbin[0]+iv)]
	fithist = hist[where(sbins GE maxbin[0]-iv)]
	fithist2 = fithist[where(fitbins LE maxbin[0]+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	;plot fit data interval, fit result, and fit center 
	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th

	;Add fit center to array of 7 keV peak channels 
	peak7s_p = [peak7s_p, coeff[1]]
	;Examine difference between peak-fit and spline values for 7 keV in channel space
	print, 'Sevens: ', val7, coeff[1], val7-coeff[1]

	;repeat the process above for 21 keV - note that the output of lclxtrem is sorted by peak height
	;which is why maxes[1] is the 21 keV peak in both p,n, despite different numbers of observed peaks
	fitbins = sbins[where(sbins GE sbins[maxes[1]]-iv)]
	fitbins2 = fitbins[where(fitbins LE sbins[maxes[1]]+iv)]
	fithist = hist[where(sbins GE sbins[maxes[1]]-iv)]
	fithist2 = fithist[where(fitbins LE sbins[maxes[1]]+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th

	peak21s_p = [peak21s_p, coeff[1]]
	print, 'Twenty-Ones: ', val21, coeff[1], val21-coeff[1]
	
	;Doing the 28 keV line (to find correct peak, sorting by channel bin)
	mbins = sbins[maxes]
	mbins = mbins[sort(mbins)]
	bin28 = mbins[3]
	if i EQ 0 then bin28 = mbins[6]
	if i EQ 7 or i EQ 1 or i EQ 6 then bin28 = mbins[5]
	if i EQ 4 or i EQ 5 or i EQ 9 then bin28 = mbins[4]
	bin28=28
	
	fitbins = sbins[where(sbins GE bin28-iv)]
	fitbins2 = fitbins[where(fitbins LE bin28+iv)]
	fithist = hist[where(sbins GE bin28-iv)]
	fithist2 = fithist[where(fitbins LE bin28+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	
	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th
	
	peak28s_p = [peak28s_p, coeff[1]]
	print, 'Twenty-Eights: ', val28, coeff[1], val28-coeff[1]
	
	al_legend, ['7: '+strtrim(peak7s_p[-1],2), '21: '+strtrim(peak21s_p[-1],2), '28: '+strtrim(peak28s_p[-1],2)], /top, /right
	
	

	coll+=1
endfor

print, ''
print, ''
print, ''

coll=1
peak7s_n = []
peak21s_n = []
peak28s_n = []
;For each n-side strip.... (see documentation in p-side loop. same deal here.)
for i=0, n_elements(relevant_strips_n)-1 do begin
	strip = relevant_strips_n[i]
	if coll EQ 5 then coll+=1
	array_str = 'nstrip_'+strtrim(strip, 2)
	;print, array_str
	void = EXECUTE('striparray = '+ array_str)
	hist = histogram(striparray, binsize=binsize, min=0, max=50, locations=sbins)
	plot, sbins, hist, thick=th, title=array_str, /ylog, psym=10, yrange=[1, max(hist)*1.5], $
			ystyle=1, ytitle='Counts', xtitle='Energy (keV)'
	oplot, sbins, hist, color=coll, thick=th, psym=10
	
	val4 = 4
	val7 = 7
	val21 = 21
	val28 = 28
	oplot, [val7, val7], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [val21, val21], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [val28, val28], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [val4, val4], [0.1, max(hist)*2], linestyle=2, color=coll+2, thick=th
	
	maxbin = sbins[where(hist EQ max(hist))]
	maxes = lclxtrem(hist, 5, /maxima)
	maxes = maxes[where(sbins[maxes] LT 500)]
	maxes = maxes[where(sbins[maxes] GE maxbin[0])]
	;for m=0, n_elements(maxes)-1 do begin
	;	oplot, [sbins[maxes[m]], sbins[maxes[m]]], [0.1, max(hist)*2], linestyle=2
	;endfor
	

	iv=3

	fitbins = sbins[where(sbins GE maxbin[0]-iv)]
	fitbins2 = fitbins[where(fitbins LE maxbin[0]+iv)]
	fithist = hist[where(sbins GE maxbin[0]-iv)]
	fithist2 = fithist[where(fitbins LE maxbin[0]+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th

	peak7s_n = [peak7s_n, coeff[1]]
	print, 'Sevens: ', val7, coeff[1], val7-coeff[1]

	fitbins = sbins[where(sbins GE sbins[maxes[1]]-iv)]
	fitbins2 = fitbins[where(fitbins LE sbins[maxes[1]]+iv)]
	fithist = hist[where(sbins GE sbins[maxes[1]]-iv)]
	fithist2 = fithist[where(fitbins LE sbins[maxes[1]]+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th

	peak21s_n = [peak21s_n, coeff[1]]
	print, 'Twenty-Ones: ', val21, coeff[1], val21-coeff[1]
	
	;Doing the 28 keV line (to find correct peak, sorting by channel bin)
	mbins = sbins[maxes]
	mbins = mbins[sort(mbins)]
	bin28 = 28
	;if i EQ 0 then bin28 = mbins[4]
	
	fitbins = sbins[where(sbins GE bin28-iv)]
	fitbins2 = fitbins[where(fitbins LE bin28+iv)]
	fithist = hist[where(sbins GE bin28-iv)]
	fithist2 = fithist[where(fitbins LE bin28+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	
	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th
	
	peak28s_n = [peak28s_n, coeff[1]]
	print, 'Twenty-Eights: ', val28, coeff[1], val28-coeff[1]
	
	al_legend, ['7: '+strtrim(peak7s_n[-1], 2), '21: '+strtrim(peak21s_n[-1],2), '28: '+strtrim(peak28s_n[-1],2)], /top, /right
	
	

	coll+=1
endfor

print, peak7s_n
print, peak7s_p


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
spawn, 'open strip_energy_spectrum.ps'

;=========================================================================================
;=========================================================================================


stop

end