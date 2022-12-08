pro boutique_cal, singlestrip=singlestrip, secondary=secondary, usepeaks=usepeaks, $
					measure_discrepancy=measure_discrepancy, thresholds=thresholds

;note: badpar.pro and lclxtrem.pro should be in IDL path or this will fail

;PURPOSE:
;Determine gain calibration for all strips crossed during 7 keV fine L-shaped ALS Scan
;Inputs: Original Peaks File, also the scan files themselves. Want to use both peaks values and scan peaks
;to do calibration for each strip. Also, want a version ONLY using ALS data.

;Method: Plot all peaks for each strip, for visual clarity (channel space). Fit gaussians to channel-space
;peaks to find locations. Then, fit linear (or quadratic?) functions to the peaks points to get a gain calibration.

;Output: .sav file with linear and quadratic gain function parameters

;set usepeaks=0 if you do not want to use the older peaksfile at all (this gives a better calibration for the scan, 
;due to differences in observed peak locations between the ALS data and what's recorded in the peaksfile).
default, usepeaks, 1

;set measure_discrepancy=1 if you want to use this procedure to measure the difference in channel space between
;the energies of the beam and the expected channel space energies from a linear fit to the peaksfile values 
;only. Note that this will just do that- no new gain calibration will be saved. 
;NOTE: designed to be done on one side of the detector at a time, edit below where files are loaded to switch
;between sides
default, measure_discrepancy, 0


;If singlestrip=1, use a threshold to include single strip events only in gain calibration
default, singlestrip, 1

;If secondary EQ 1, use a secondary threshold for finding double strip events (lower than the trigger threshold)
default, secondary, 1
if secondary EQ 1 and singlestrip EQ 0 then begin
	print, 'You cant use a secondary threshold without using a primary threshold!'
	stop
endif


;Inputs
peaksfile = '~/foxsi/routines/peaks_fe_3am_cdte_fec7_Al.sav'
foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=struct_files
restore, peaksfile

side=0

;Comment out the side you aren't using
if measure_discrepancy EQ 1 then begin
;	print, 'USING PT-SIDE BRANCH'
;	foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=struct_files
;	side=2
	print, 'USING AL-SIDE BRANCH'
	foxsi_cal_struct_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=struct_files
	side=1
endif

;Note: thresholding in channel space - for channel-space thresholds calculation, see thresholds.pro. Note that it 
;depends on an output result from this procedure. 
tholdn = 22. ;ADC, ~4 keV
thold2n = 11.	;ADC, ~2.5 keV
tholdp = 25. ;ADC, ~4 keV
thold2p = 14.	;ADC, ~2.5 keV

if keyword_set(thresholds) then begin
	tholdn = thresholds[0]
	thold2n = thresholds[1]
	tholdp = thresholds[2]
	thold2p = thresholds[3]
endif
	

;Pre-known relevant strips - which_strips.pro has code that determines strip coverage during ALS scans
relevant_strips_p = indgen(10)+33
relevant_as_p = 2
relevant_strips_n = indgen(7)+25
relevant_as_n = 1

;Defining event arrays for each strip
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


;making a channel-space array for one strip to export if desired for debugging
chan25 =[]

for j=0, n_elements(struct_files)-1 do begin
		;restore file (one beam position)
		restore, struct_files[j]
		n_evts = n_elements(data)
		
		;Making common mode array
		cmn = fltarr(n_evts, 4)
		for as=0, 3 do begin
			cmn[*,as] = data[*].cmn_median[as] + randomu(seed,n_evts*4+as) - 0.5
		endfor
	
		;For each p-side strip, let's make an array of event energies (channel space) and
		;assign it to a strip-named variable.
		for i=0, n_elements(relevant_strips_p)-1 do begin
			strip = relevant_strips_p[i]
			events=[]

			;for each frame...
			for k=0, n_evts-1 do begin
				vals = data[k].data[2,*]
				
				;If single strip is set, use spline to find a "threshold" channel value and
				;take only events which are the only above-threshold event in their frame
				if singlestrip EQ 1 then begin
					;old thresholding methods (spline, in particular, makes this take AGES)
					;thold = cmn[k,2] + (thresh-gain_coeffs[0])/gain_coeffs[1] 
					;thold = cmn[k,2] + spline(peaks[*,strip,2,1], peaks[*,strip,2,0], thresh)
					above = vals[where(vals GE cmn[k,2]+tholdp, /NULL)]
					
					;If there is only one strip above threshold and we are not using a secondary threshold:
					if n_elements(above) EQ 1 and above NE !NULL and secondary EQ 0 then begin
						mm = max(vals)
						if data[k].data[2,strip] EQ mm then events = [events, mm-cmn[k,2] ]
						;if data[k].data[2,strip] EQ mm then events = [events, mm]
					endif
					
					;If there is one strip above threshold and we ARE using a secondary threshold:
					if n_elements(above) EQ 1 and above NE !NULL and secondary EQ 1 then begin
						;thold2 = cmn[k,2] + (thresh2-gain_coeffs[0])/gain_coeffs[1]
						;thold2 = cmn[k,2] + spline(peaks[*,strip,2,1], peaks[*,strip,2,0], thresh2)
						above2 = vals[where(vals GE cmn[k,2]+thold2p, /NULL)]
						
						;If there is still only one event above the lower secondary threshold:
						if n_elements(above2) EQ 1 and above2 NE !NULL then begin
							mm = max(vals)
							if data[k].data[2,strip] EQ mm then events = [events, mm-cmn[k,2] ]
							
							;I think it is better to work in channel space with the common med already subtracted,
							;because if you don't do the subtraction here in the loop (where we are working on an event-by-event
							;basis) you then end up needing to use the mean/median/some other estimate of the common med
							;later, which seems to add some uncertainty. But if you want to work in raw ADC values, switch
							;to using the following (instead of the previous line) here and in all other instances.
							
							;if data[k].data[2,strip] EQ mm then events = [events, mm]
						endif
					endif
					
				;if singlestrip is not set, just take the maximum of each data frame on each side.
				endif else begin
					mm = max(vals)
					if data[k].data[2,strip] EQ mm then events = [events, mm-cmn[k,2] ] 
					;if data[k].data[2,strip] EQ mm then events = [events, mm ];-cmn[k,2]]
				endelse
			endfor
			array_str = 'pstrip_'+strtrim(relevant_strips_p[i], 2)
			void = EXECUTE(array_str + '=['+array_str +', events]')
			;print, strip, n_elements(events)
		endfor
			
	
		;repeat for n-side (see p-side documentation above for more details)
		for i=0, n_elements(relevant_strips_n)-1 do begin
			strip = relevant_strips_n[i]
			events=[]
			for k=0, n_evts-1 do begin
				vals = data[k].data[1,*]
				if singlestrip EQ 1 then begin
					;thold = cmn[k,1] + (thresh-gain_coeffs[0])/gain_coeffs[1]
					;thold = cmn[k,1] + spline(peaks[*,strip,1,1], peaks[*,strip,1,0], thresh)
					above = vals[where(vals GE  cmn[k,1]+tholdn, /NULL)]
					;If there is only one strip above threshold and we are not using a secondary threshold:
					if n_elements(above) EQ 1 and above NE !NULL and secondary EQ 0 then begin
						mm = max(vals)
						if data[k].data[1,strip] EQ mm then events = [events, mm-cmn[k,1] ]
						if data[k].data[1,strip] EQ mm and i EQ 0 then chan25 = [chan25, mm-cmn[k,1]]
						;if data[k].data[1,strip] EQ mm then events = [events, mm]
					endif
					;If there is one strip above threshold and we ARE using a secondary threshold:
					if n_elements(above) EQ 1 and above NE !NULL and secondary EQ 1 then begin
						;thold2 = cmn[k,1] + (thresh2_Al-gain_coeffs[0])/gain_coeffs[1]
						;thold2 = cmn[k,1] + spline(peaks[*,strip,1,1], peaks[*,strip,1,0], thresh2_Al)
						above2 = vals[where(vals GE cmn[k,1]+thold2n, /NULL)]
						;If there is still only one event above the lower secondary threshold:
						if n_elements(above2) EQ 1 and above2 NE !NULL then begin
							mm = max(vals)
							if data[k].data[1,strip] EQ mm then events = [events, mm-cmn[k,1] ]
							if data[k].data[1,strip] EQ mm and i EQ 0 then chan25 = [chan25, mm-cmn[k,1]]
							;if data[k].data[1,strip] EQ mm then events = [events, mm]
						endif
					endif
				endif else begin
					mm = max(vals)
					if data[k].data[1,strip] EQ mm then events = [events, mm-cmn[k,1] ] ;-cmn[k,1]]
					if data[k].data[1,strip] EQ mm and i EQ 0 then chan25 = [chan25, mm-cmn[k,1]]
					;if data[k].data[1,strip] EQ mm then events = [events, mm ]
				endelse
			endfor
			array_str = 'nstrip_'+strtrim(relevant_strips_n[i], 2)
			void = EXECUTE(array_str + '=['+array_str +', events]')
			;print, strip, n_elements(events)
		endfor	
		
		
		
		
	;stop
endfor 


;If desired, for exporting a one-strip channel-space array for debugging purposes
;save, chan25, filename='chan25.sav'

;STOP

;=========================================================================================
;=========================================================================================

;Plot channel-space spectrum, and find peak locations (on each side, for each strip)

popen, 'strip_channel_spectrum.ps', $
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

	
binsize=2
if side EQ 0 or side EQ 2 then begin
	coll=1
	peak7s_p = []
	peak21s_p = []
	peak28s_p = []
	;For each p-side strip...
	for i=0, n_elements(relevant_strips_p)-1 do begin
		strip = relevant_strips_p[i]
		print, 'Strip: ', strip
		if coll EQ 5 then coll+=1
		array_str = 'pstrip_'+strtrim(strip, 2)
		void = EXECUTE('striparray = '+ array_str)
		;make and plot channel-space spectrum histogram
		hist = histogram(striparray, binsize=binsize, max=600, locations=sbins)
		plot, sbins, hist, thick=th, title=array_str, /ylog, psym=10, yrange=[1, max(hist)*1.5], $
			ystyle=1, ytitle='Counts', xtitle='Channel Value'
		oplot, sbins, hist, color=coll, thick=th, psym=10
	
		;via spline, estimate channels corresponding to 4, 7, 21, 28 keV (lines + commonly used threshold)
		;and plot for later comparison with found peaks. 
		;if working in raw ADC space, need to add cmn:
	;	val4 = mean(cmn[*,2]) + spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 4.)
	;	val7 = mean(cmn[*,2]) + spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 7.)
	;	val21 = mean(cmn[*,2]) + spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 21.)
	;	val28 = mean(cmn[*,2]) + spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 28.)
		val4 = spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 4.)
		val7 = spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 7.)
		val21 = spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 21.)
		val28 = spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 28.)
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

		iv=15

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
		;print, 'Sevens: ', val7, coeff[1], val7-coeff[1]

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
		;print, 'Twenty-Ones: ', val21, coeff[1], val21-coeff[1]
	
		;Doing the 28 keV line (to find correct peak, sorting by channel bin)
		mbins = sbins[maxes]
		mbins = mbins[sort(mbins)]
		bin28 = 190
	
		fitbins = sbins[where(sbins GE bin28-iv)]
		fitbins2 = fitbins[where(fitbins LE bin28+iv)]
		fithist = hist[where(sbins GE bin28-iv)]
		fithist2 = fithist[where(fitbins LE bin28+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	
		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th
	
		peak28s_p = [peak28s_p, coeff[1]]
		;print, 'Twenty-Eights: ', val28, coeff[1], val28-coeff[1]
	
		al_legend, ['Total Single-Strip: '+strtrim(n_elements(striparray), 2), 'Chan 7: '+strtrim(peak7s_p[-1], 2), $
					'Chan 21: '+strtrim(peak21s_p[-1], 2), 'Chan 28: '+strtrim(peak28s_p[-1], 2)], color=[coll, 0,0,0], /top, /right
	

		coll+=1
	endfor
endif

if side EQ 0 or side EQ 1 then begin
	coll=1
	peak7s_n = []
	peak21s_n = []
	peak28s_n = []
	;For each n-side strip.... (see documentation in p-side loop. same deal here.)
	for i=0, n_elements(relevant_strips_n)-1 do begin
		strip = relevant_strips_n[i]
		print, 'Strip: ', strip
		if coll EQ 5 then coll+=1
		array_str = 'nstrip_'+strtrim(strip, 2)
		;print, array_str
		void = EXECUTE('striparray = '+ array_str)
		hist = histogram(striparray, binsize=binsize, max=600, locations=sbins)
		plot, sbins, hist, thick=th, title=array_str, /ylog, psym=10, yrange=[1, max(hist)*1.5], $
				ystyle=1, ytitle='Counts', xtitle='Channel Value'
		oplot, sbins, hist, color=coll, thick=th, psym=10
	;	
	;	val4 = mean(cmn[*,1]) + spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 4.)
	;	val7 = mean(cmn[*,1]) + spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 7.)
	;	val21 = mean(cmn[*,1]) + spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 21.)
	;	val28 = mean(cmn[*,1]) + spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 28.)
		val4 = spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 4.)
		val7 = spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 7.)
		val21 = spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 21.)
		val28 = spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 28.)
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
	

		iv=15

		fitbins = sbins[where(sbins GE maxbin[0]-iv)]
		fitbins2 = fitbins[where(fitbins LE maxbin[0]+iv)]
		fithist = hist[where(sbins GE maxbin[0]-iv)]
		fithist2 = fithist[where(fitbins LE maxbin[0]+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th

		peak7s_n = [peak7s_n, coeff[1]]
		;print, 'Sevens: ', val7, coeff[1], val7-coeff[1]

		fitbins = sbins[where(sbins GE sbins[maxes[1]]-iv)]
		fitbins2 = fitbins[where(fitbins LE sbins[maxes[1]]+iv)]
		fithist = hist[where(sbins GE sbins[maxes[1]]-iv)]
		fithist2 = fithist[where(fitbins LE sbins[maxes[1]]+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th

		peak21s_n = [peak21s_n, coeff[1]]
		;print, 'Twenty-Ones: ', val21, coeff[1], val21-coeff[1]
	
		;Doing the 28 keV line (to find correct peak, sorting by channel bin)
		mbins = sbins[maxes]
		mbins = mbins[sort(mbins)]
		bin28 = mbins[3]
		if i EQ 0 and mbins[3] LE 180 then bin28 = mbins[4]
		if i EQ 2 or i EQ 3 and n_elements(mbins) GE 5 then bin28 = mbins[4]
		if i EQ 6 then bin28 = val28
		bin28=210
	
		fitbins = sbins[where(sbins GE bin28-iv)]
		fitbins2 = fitbins[where(fitbins LE bin28+iv)]
		fithist = hist[where(sbins GE bin28-iv)]
		fithist2 = fithist[where(fitbins LE bin28+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	
		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th
	
		peak28s_n = [peak28s_n, coeff[1]]
		;print, 'Twenty-Eights: ', val28, coeff[1], val28-coeff[1]
	
		al_legend, ['Total Single-Strip: '+strtrim(n_elements(striparray), 2), 'Chan 7: '+strtrim(peak7s_n[-1], 2), $
					'Chan 21: '+strtrim(peak21s_n[-1], 2), 'Chan 28: '+strtrim(peak28s_n[-1], 2)], textcol=[coll, 0,0,0], /top, /right
	

		coll+=1
	endfor
endif

;print, peak7s_n
;print, peak7s_p


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
spawn, 'open strip_channel_spectrum.ps'

;=========================================================================================
;=========================================================================================


;=========================================================================================
;=========================================================================================

;Plotting peaks in channel-energy space and fitting the relationship

popen, 'strip_peaks.ps', $
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
symsize=1
psym=4

;=========================================================================================
;=========================================================================================


LINFITS_P = []
LINFITS_N = []
QUADFITS_P = []
QUADFITS_N = []

;If we are using this to look at line discrepancies, make arrays to hold values for each strip
if measure_discrepancy EQ 1 then begin
	disc7_n=[]
	disc21_n=[]
	disc28_n=[]
	disc7_p=[]
	disc21_p=[]
	disc28_p=[]
	
	edisc7_n=[]
	edisc21_n=[]
	edisc28_n=[]
	edisc7_p=[]
	edisc21_p=[]
	edisc28_p=[]
endif

if side EQ 0 or side EQ 2 then begin
	coll=1
	;For each p-side strip...
	for i=0, n_elements(relevant_strips_p)-1 do begin
		strip = relevant_strips_p[i]
		if coll EQ 5 then coll+=1
		array_str = 'pstrip_'+strtrim(strip, 2)
	

		if usepeaks EQ 1 then begin 
			;Add all original peaks, plus the new found peaks to new array
			all_chan = [peaks[*,strip,2,0], peak7s_p[i], peak21s_p[i], peak28s_p[i]]
			;all_chan = [peaks[*,strip,2,0], peak7s_p[i]-mean(cmn[*,2]), peak21s_p[i]-mean(cmn[*,2]), peak28s_p[i]-mean(cmn[*,2])]
			all_chan = all_chan[sort(all_chan)]
		
			;Same thing in energy space:
			all_e = [peaks[*,strip,2,1], 7., 21., 28.]
			all_e = all_e[sort(all_e)]
		
			;trimming to only consider the peaks in the vicinity of the range we are interested in
			all_e = all_e[2:7]
			all_chan = all_chan[2:7]
		endif 
		if usepeaks EQ 0 then begin
			
			;Noisepeak found from noise_counts.pro
			noisepeak_2=-1
			;Add only the new peaks
			all_chan = [noisepeak_2, peak7s_p[i], peak21s_p[i], peak28s_p[i]]
			all_e = [0., 7., 21., 28.]
		endif
		if measure_discrepancy EQ 1 then begin
			all_chan = peaks[2:5,strip,2,0]
			all_e = peaks[2:5,strip,2,1]
		endif
		

		;Plot energy-channel pairs we will use for fits
		plot, all_chan , all_e, title=array_str, psym=psym, symsize=symsize, thick=th, $
				xtitle='Channel-Space Value', ytitle='Energy-Space Value'
				;yrange=[0, 105], ystyle=1, xrange=[0, 800], xstyle=1, thick=th
				;yrange=[0, 30], thick=th, ystyle=1;, xrange=[0, 200], xstyle=1
		;Plot them again so they can be pretty colors
		oplot, all_chan , all_e, psym=psym, symsize=symsize, color=coll, thick=th
		;Overplot spline interpolation between points for comparison
		;oplot, spline(all_chan, all_e, findgen(400)), color=coll+1, thick=th
		;Overplot 7, 21 keV values found via spline method 
		;oplot, [spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 7.), spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 21.), $
		;		spline(peaks[*,strip,2,1], peaks[*,strip,2,0], 28.)], [7., 21., 28.], color=coll+2, psym=psym, thick=th 
	
		if measure_discrepancy EQ 1 then oplot, [peak7s_p[i], peak21s_p[i], peak28s_p[i]], [7., 21., 28.], psym=psym, symsize=symsize, color=coll+1, thick=th

		;Linear fit to values
		linfit_strip = LINFIT(all_chan, all_e)
		;print, 'Linfit parameters p:', linfit_strip
		linfit_arr = all_chan*linfit_strip[1]+linfit_strip[0]
		oplot, all_chan, linfit_arr, thick=th/2., linestyle=2, color=coll
	
		LINFITS_P = [[LINFITS_P], [linfit_strip]]
	
		if measure_discrepancy EQ 1 then begin
			;print, 'Expected channel value for 7 keV: ', (7. - linfit_strip[0])/linfit_strip[1]
			;print, 'ALS channel value for 7 keV: ', peak7s_p[i]
			;print, 'Difference: ', (7. - linfit_strip[0])/linfit_strip[1]-peak7s_p[i]
			disc = (7.-linfit_strip[0])/linfit_strip[1]-peak7s_p[i]
			disc7_p = [disc7_p, disc]
			dice = 7. - (peak7s_p[i]*linfit_strip[1]+linfit_strip[0])
			edisc7_p = [edisc7_p, dice]
			;print, 'Expected channel value for 21 keV: ', (21. - linfit_strip[0])/linfit_strip[1]
			;print, 'ALS channel value for 21 keV: ', peak21s_p[i]
			;print, 'Difference: ', (21. - linfit_strip[0])/linfit_strip[1]-peak21s_p[i]
			disc = (21. - linfit_strip[0])/linfit_strip[1]-peak21s_p[i]
			disc21_p = [disc21_p, disc]
			dice = 21. - (peak21s_p[i]*linfit_strip[1]+linfit_strip[0])
			edisc21_p = [edisc21_p, dice]
	;		print, 'Expected channel value for 28 keV: ', (28. - linfit_strip[0])/linfit_strip[1]
	;		print, 'ALS channel value for 28 keV: ', peak28s_p[i]
	;		print, 'Difference: ', (28. - linfit_strip[0])/linfit_strip[1]-peak28s_p[i]
			disc = (28. - linfit_strip[0])/linfit_strip[1]-peak28s_p[i]
			disc28_p = [disc28_p, disc]
			dice = 28. - (peak28s_p[i]*linfit_strip[1]+linfit_strip[0])
			edisc28_p = [edisc28_p, dice]
		endif
	
		;Quadratic fit to values
		result = POLY_FIT(all_chan, all_e, 2)
		;print, 'QUADfit parameters p:', result
		QUADfit_arr = all_chan^2*result[2]+all_chan*result[1]+result[0]
		;oplot, all_chan, QUADfit_arr, thick=th/2., linestyle=3, color=coll
	
		QUADFITS_P = [[QUADFITS_P], [reform(result)]]
	
		coll+=1
	
	endfor
endif


if side EQ 0 or side EQ 1 then begin
	coll=1
	;For each n-side strip.... (see p-side documentation for more detail)
	for i=0, n_elements(relevant_strips_n)-1 do begin
		strip = relevant_strips_n[i]
		if coll EQ 5 then coll+=1
		array_str = 'nstrip_'+strtrim(strip, 2)
	
		if usepeaks EQ 1 then begin 
			;Add all original peaks, plus the new found peaks to new array
			all_chan = [peaks[*,strip,1,0], peak7s_n[i], peak21s_n[i], peak28s_n[i]]
			;all_chan = [peaks[*,strip,1,0], peak7s_n[i]-mean(cmn[*,1]), peak21s_n[i]-mean(cmn[*,1]), peak28s_n[i]-mean(cmn[*,1])]
			all_chan = all_chan[sort(all_chan)]
		
			;Same thing in energy space:
			all_e = [peaks[*,strip,1,1], 7., 21., 28.]
			all_e = all_e[sort(all_e)]
		
			;trimming to only consider the peaks in the vicinity of the range we are interested in
			all_e = all_e[2:7]
			all_chan = all_chan[2:7]
		endif 
		if usepeaks EQ 0 then begin
			;Noisepeak found from noise_counts.pro
			noisepeak_1=-0.8
			;Add only the new peaks
			all_chan = [noisepeak_1, peak7s_n[i], peak21s_n[i], peak28s_n[i]]
			all_e = [0., 7., 21., 28.]
		endif
		if measure_discrepancy EQ 1 then begin
			all_chan = peaks[2:5,strip,1,0]
			all_e = peaks[2:5,strip,1,1]
		endif
	
		plot, all_chan , all_e, title=array_str, psym=psym, symsize=symsize, thick=th
				;yrange=[0, 105], ystyle=1, xrange=[0, 800], xstyle=1, thick=th
				;yrange=[0, 30], thick=th, ystyle=1;, xrange=[0, 200], xstyle=1
		oplot, all_chan , all_e, psym=psym, symsize=symsize, color=coll, thick=th
		;oplot, spline(all_chan, all_e, findgen(400)), color=coll+1, thick=th
		;oplot, [spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 7.), spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 21.) , $
		;		spline(peaks[*,strip,1,1], peaks[*,strip,1,0], 28.)], [7., 21., 28.], color=coll+2, psym=psym, thick=th 
	
		if measure_discrepancy EQ 1 then oplot, [peak7s_n[i], peak21s_n[i], peak28s_n[i]], [7., 21., 28.], psym=psym, symsize=symsize, color=coll+1, thick=th

	
		;-mean(cmn[*,1])
	
	;	plot, [peak7s_p[i],peak21s_p[i]] , [7., 21.], title=array_str, psym=psym, symsize=symsize, $
	;			yrange=[0, 105], ystyle=1, xrange=[0, 800], xstyle=1
	;	oplot, [peak7s_p[i],peak21s_p[i]] , [7., 21.], psym=psym, symsize=symsize, color=coll
	;	;Havent defined the common median to subtract yet. here, (line from test_spline) the mean value of all the
	;	;common median values for all the events in a file is used. Is this okay to represent ALL files?
	;	oplot, peaks[*,strip,2,0], peaks[*,strip,2,1], color=coll, psym=psym, symsize=symsize, thick=th



		;Linear fit to values
		linfit_strip = LINFIT(all_chan, all_e)
		;print, 'Linfit parameters n:', linfit_strip
		linfit_arr = all_chan*linfit_strip[1]+linfit_strip[0]
		oplot, all_chan, linfit_arr, thick=th/2., linestyle=2, color=coll
	
		LINFITS_N = [[LINFITS_N], [linfit_strip]]
	
		if measure_discrepancy EQ 1 then begin
	;		print, 'Expected channel value for 7 keV: ', (7. - linfit_strip[0])/linfit_strip[1]
	;		print, 'ALS channel value for 7 keV: ', peak7s_n[i]
	;		print, 'Difference: ', (7. - linfit_strip[0])/linfit_strip[1]-peak7s_n[i]
			disc = (7.-linfit_strip[0])/linfit_strip[1]-peak7s_n[i]
			disc7_n = [disc7_n, disc]
			dice = 7. - (peak7s_n[i]*linfit_strip[1]+linfit_strip[0])
			edisc7_n = [edisc7_n, dice]
	;		print, 'Expected channel value for 21 keV: ', (21. - linfit_strip[0])/linfit_strip[1]
	;		print, 'ALS channel value for 21 keV: ', peak21s_n[i]
	;		print, 'Difference: ', (21. - linfit_strip[0])/linfit_strip[1]-peak21s_n[i]
			disc = (21. - linfit_strip[0])/linfit_strip[1]-peak21s_n[i]
			disc21_n = [disc21_n, disc]
			dice = 21. - (peak21s_n[i]*linfit_strip[1]+linfit_strip[0])
			edisc21_n = [edisc21_n, dice]
	;		print, 'Expected channel value for 28 keV: ', (28. - linfit_strip[0])/linfit_strip[1]
	;		print, 'ALS channel value for 28 keV: ', peak28s_n[i]
	;		print, 'Difference: ', (28. - linfit_strip[0])/linfit_strip[1]-peak28s_n[i]
			disc = (28. - linfit_strip[0])/linfit_strip[1]-peak28s_n[i]
			disc28_n = [disc28_n, disc]
			dice = 28. - (peak28s_n[i]*linfit_strip[1]+linfit_strip[0])
			edisc28_n = [edisc28_n, dice]
		endif
	
	
		;Quadratic 
	
		result = POLY_FIT(all_chan, all_e, 2)
		;print, 'QUADfit parameters n:', result
		QUADfit_arr = all_chan^2*result[2]+all_chan*result[1]+result[0]
		oplot, all_chan, QUADfit_arr, thick=th/2., linestyle=3, color=coll
	
		QUADFITS_N = [[QUADFITS_N], [reform(result)]]
	
		coll+=1
	
	endfor
endif

;print, LINFITS_P
;print, size(LINFITS_P)
;print, QUADFITS_P
;print, size(QUADFITS_P)

;If we are using this to look at line discrepancies, make arrays to hold values for each strip
if measure_discrepancy EQ 1 then begin
	if side EQ 1 then begin
		print, 'N-Side:'
		print, disc7_n
		print, 'Mean discrepancy - nside 7 keV channel space: '  , mean(disc7_n)
		print, disc21_n
		print, 'Mean discrepancy - nside 21 keV channel space: ' , mean(disc21_n)
		print, disc28_n
		print, 'Mean discrepancy - nside 28 keV channel space: ' , mean(disc28_n)
		print, 'N-Side - keV:'
		print, edisc7_n
		print, 'Mean discrepancy - nside 7 keV energy space: ' , mean(edisc7_n)
		print, edisc21_n
		print, 'Mean discrepancy - nside 21 keV energy space: ' , mean(edisc21_n)
		print, edisc28_n
		print, 'Mean discrepancy - nside 28 keV energy space: ' , mean(edisc28_n)
	endif
	if side EQ 2 then begin
		print, 'P-Side:'
		print, disc7_p
		print, 'Mean discrepancy - pside 7 keV channel space: ' , mean(disc7_p)
		print, disc21_p
		print, 'Mean discrepancy - pside 21 keV channel space: ' , mean(disc21_p)
		print, disc28_p
		print, 'Mean discrepancy - pside 28 keV channel space: ' , mean(disc28_p)
		print, 'P-Side - keV:'
		print, edisc7_p
		print, 'Mean discrepancy - pside 7 keV energy space: ' ,mean(edisc7_p)
		print, edisc21_p
		print, 'Mean discrepancy - pside 21 keV energy space: ' ,mean(edisc21_p)
		print, edisc28_p
		print, 'Mean discrepancy - pside 28 keV energy space: ' ,mean(edisc28_p)
	endif
endif

if measure_discrepancy EQ 0 then save, LINFITS_P, LINFITS_N, QUADFITS_P, QUADFITS_N, filename='boutique_cal.sav'


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
;DEVICE, /CLOSE
spawn, 'open strip_peaks.ps'

;=========================================================================================
;=========================================================================================




stop

end