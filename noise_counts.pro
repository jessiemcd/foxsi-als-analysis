pro noise_counts

;Look at counts in many strips (including non-illuminated strips) to characterize noise pedestal

;This procedure is fairly hard-coded to this specific scan (ALS 2019 L-Scan) 
;(e.g. it does this only for ASICS 1 and 2, the illuminated ASICs during that scan, etc.)

;Beam spread also examined at the end (counts in all strips shown while beam is illuminating one)

;Starting channel-space thresholds, for reference. (t1: guess at 4 keV, t2: guess at extent of noise)
pcurrent_t1=25
pcurrent_t2=14
ncurrent_t1=22
ncurrent_t2=11

;find all the files
foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=struct_files_p
foxsi_cal_struct_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=struct_files_n

;to do with only a few non-illuminated strips
;relevant_strips_p = indgen(40)+20
;relevant_strips_n = indgen(45)+5

;to do with all strips
relevant_strips_p = indgen(64)
relevant_strips_n = indgen(64)

;To use a number of files
;noise_files_p = struct_files_n[0:75]
;noise_files_n = struct_files_p

;To use just one file per side
noise_files_p = struct_files_n[1]
noise_files_n = struct_files_p[1]

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



above_t2_p = []
above_t2_n = []

;For each file, put events into arrays for the strips in which they registered
for j=0, n_elements(noise_files_p)-1 do begin
		;restore file (one beam position)
		restore, noise_files_p[j]
		n_evts = n_elements(data)
		print, '# Events: ', n_evts
		
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
				val = data[k].data[2, strip] - cmn[k,2]
				events = [events, val]
			endfor
			array_str = 'pstrip_'+strtrim(relevant_strips_p[i], 2)
			void = EXECUTE(array_str + '=['+array_str +', events]')
			;print, strip, n_elements(events)
		endfor
endfor


for j=0, n_elements(noise_files_n)-1 do begin
		;restore file (one beam position)
		restore, noise_files_n[j]
		n_evts = n_elements(data)
		print, '# Events: ', n_evts
		
		;Making common mode array
		cmn = fltarr(n_evts, 4)
		for as=0, 3 do begin
			cmn[*,as] = data[*].cmn_median[as] + randomu(seed,n_evts*4+as) - 0.5
		endfor
	
		;For each n-side strip, let's make an array of event energies (channel space) and
		;assign it to a strip-named variable.
		for i=0, n_elements(relevant_strips_n)-1 do begin
			strip = relevant_strips_n[i]
			events=[]

			;for each frame...
			for k=0, n_evts-1 do begin
				val = data[k].data[1, strip] - cmn[k,1]
				events = [events, val]
			endfor
			array_str = 'nstrip_'+strtrim(relevant_strips_n[i], 2)
			void = EXECUTE(array_str + '=['+array_str +', events]')
			;print, strip, n_elements(events)
		endfor
endfor



;=========================================================================================
;=========================================================================================

;Plot channel-space spectrum, and find peak locations (on each side, for each strip)

popen, 'strip_channel_noise_spectrum.ps', $
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



above_t2_p = []
above_t2_n = []

above_t1_p = []
above_t1_n = []

p_fits = []
n_fits = []


	
binsize=2
coll=1
;For each p-side strip...
for i=0, n_elements(relevant_strips_p)-1 do begin
	strip = relevant_strips_p[i]
	print, 'Strip: ', strip
	if coll EQ 5 then coll+=1
	array_str = 'pstrip_'+strtrim(strip, 2)
	void = EXECUTE('striparray = '+ array_str)
	
	above = n_elements(where(striparray GT pcurrent_t2))
	above_t2_p = [above_t2_p, above]
	above = n_elements(where(striparray GT pcurrent_t1))
	above_t1_p = [above_t1_p, above]
	
	;make and plot channel-space spectrum histogram
	hist = histogram(striparray, binsize=binsize, max=600, min=-50, locations=sbins)
	plot, sbins, hist, thick=th, title=array_str, /ylog, psym=10, yrange=[1, max(hist)*1.5], ystyle=1,$
						xrange=[-50,150], xstyle=1
	oplot, sbins, hist, color=coll, thick=th, psym=10

	;plot 'guess' thresholds
	oplot, [pcurrent_t1, pcurrent_t1], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [pcurrent_t2, pcurrent_t2], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [8.5, 8.5], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	
	;find maximum
	maxbin = sbins[where(hist EQ max(hist))]
	iv=15

	;make data interval for fitting first peak + do gaussian fit
	fitbins = sbins[where(sbins GE maxbin[0]-iv)]
	fitbins2 = fitbins[where(fitbins LE maxbin[0]+iv)]
	fithist = hist[where(sbins GE maxbin[0]-iv)]
	fithist2 = fithist[where(fitbins LE maxbin[0]+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th
	
	;not including illuminated strips
	if strip LE 39 or strip GE 47 and strip NE 63 then p_fits = [[p_fits], [coeff, n_elements(striparray)]]

	coll+=1
endfor



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
	
	above = n_elements(where(striparray GT ncurrent_t2))
	above_t2_n = [above_t2_n, above]
	above = n_elements(where(striparray GT ncurrent_t1))
	above_t1_n = [above_t1_n, above]
	
	hist = histogram(striparray, binsize=binsize, max=600, min=-50, locations=sbins)
	plot, sbins, hist, thick=th, title=array_str, /ylog, psym=10, yrange=[1, max(hist)*1.5], ystyle=1,$
						xrange=[-50,150], xstyle=1
	oplot, sbins, hist, color=coll, thick=th, psym=10
;	
	oplot, [ncurrent_t1, ncurrent_t1], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [ncurrent_t2, ncurrent_t2], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	oplot, [9, 9], [0.1, max(hist)*2], linestyle=2, color=coll+1, thick=th
	
	;find maximum
	maxbin = sbins[where(hist EQ max(hist))]
	iv=15

	;make data interval for fitting first peak + do gaussian fit
	fitbins = sbins[where(sbins GE maxbin[0]-iv)]
	fitbins2 = fitbins[where(fitbins LE maxbin[0]+iv)]
	fithist = hist[where(sbins GE maxbin[0]-iv)]
	fithist2 = fithist[where(fitbins LE maxbin[0]+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(hist)*2], linestyle=2, col=100, thick=th
	
	if strip LE 7 or strip GE 34 and strip NE 63 then n_fits = [[n_fits], [coeff, n_elements(striparray)]]

	coll+=1
endfor


;print, peak7s_n
;print, peak7s_p


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
spawn, 'open strip_channel_noise_spectrum.ps'

;=========================================================================================
;=========================================================================================

sp = size(p_fits)

print, 'ASIC2 Mean Noise Center: ', mean(p_fits[1,*]), stdev(p_fits[1,*])
print, 'ASIC2 Mean Noise Sigma: ', mean(p_fits[2,*]), stdev(p_fits[2,*])
print, 'ASIC2 Mean Total Counts (non-illuminated): ', mean(p_fits[3,*]), stdev(p_fits[3,*])
print, 'ASIC2 1-Sigma Noise Threshold: ', mean(p_fits[1,*])+mean(p_fits[2,*]), ' Mean Counts Per File Above Estimate: ', 0.16*mean(p_fits[3,*])
print, 'ASIC2 2-Sigma Noise Threshold: ', mean(p_fits[1,*])+2*mean(p_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.025*mean(p_fits[3,*])
print, 'ASIC2 3-Sigma Noise Threshold: ', mean(p_fits[1,*])+3*mean(p_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.0015*mean(p_fits[3,*])
print, 'ASIC2 3.5-Sigma Noise Threshold: ', mean(p_fits[1,*])+3.5*mean(p_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.00023*mean(p_fits[3,*])
print, 'ASIC2 4-Sigma Noise Threshold: ', mean(p_fits[1,*])+4*mean(p_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.000032*mean(p_fits[3,*])

sn = size(n_fits)

print, ''
print, 'ASIC1 Noise Center: ', mean(n_fits[1,*]), stdev(n_fits[1,*])
print, 'ASIC1 Noise Sigma: ', mean(n_fits[2,*]), stdev(n_fits[2,*])
print, 'ASIC1 Mean Total Counts (non-illuminated): ', mean(n_fits[3,*]), stdev(n_fits[3,*])
print, 'ASIC1 1-Sigma Noise Threshold: ', mean(n_fits[1,*])+mean(n_fits[2,*]), ' Mean Counts Per File Above Estimate: ', 0.16*mean(n_fits[3,*])
print, 'ASIC1 2-Sigma Noise Threshold: ', mean(n_fits[1,*])+2*mean(n_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.025*mean(n_fits[3,*])
print, 'ASIC1 3-Sigma Noise Threshold: ', mean(n_fits[1,*])+3*mean(n_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.0015*mean(n_fits[3,*])
print, 'ASIC1 3.5-Sigma Noise Threshold: ', mean(n_fits[1,*])+3.5*mean(n_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.00023*mean(n_fits[3,*])
print, 'ASIC1 4-Sigma Noise Threshold: ', mean(n_fits[1,*])+4*mean(n_fits[2,*]), ' Mean Counts Per File Above Estimate: ',  0.000032*mean(n_fits[3,*])




;=========================================================================================
;=========================================================================================

;Plot channel-space spectrum, and find peak locations (on each side, for each strip)

popen, 'beam_spread.ps', $
	xsi=8, ysi=10
!Y.margin=4.
!X.margin=4.
ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
!p.multi=[0,1,4]
linecolors

;=========================================================================================
;=========================================================================================



plot, relevant_strips_p, above_t2_p, thick=th, psym=10, /ylog, xrange=[0.,64.], xstyle=1, $
		title='ASIC 2 Counts above threshold, log y scale', xtitle='Strip #', ytitle='Counts'
oplot,  relevant_strips_p, above_t1_p, thick=th, col=3, psym=10

al_legend, ['Primary Threshold', 'Secondary Threshold'], textcol=[3,0]

plot, relevant_strips_p, above_t2_p, thick=th, psym=10, xrange=[0.,64.], xstyle=1, $
		title='ASIC 2 Counts above threshold, linear y scale', xtitle='Strip #', ytitle='Counts'
oplot,  relevant_strips_p, above_t1_p, thick=th, col=3, psym=10

al_legend, ['Primary Threshold', 'Secondary Threshold'], textcol=[3,0]

plot, relevant_strips_n, above_t2_n, thick=th, psym=10, /ylog, xrange=[0.,64.], xstyle=1, $
		title='ASIC 1 Counts above threshold, log y scale', xtitle='Strip #', ytitle='Counts'
oplot,  relevant_strips_n, above_t1_n, thick=th, col=3, psym=10

al_legend, ['Primary Threshold', 'Secondary Threshold'], textcol=[3,0]

plot, relevant_strips_n, above_t2_n, thick=th, psym=10, xrange=[0.,64.], xstyle=1, $
		title='ASIC 1 Counts above threshold, linear y scale', xtitle='Strip #', ytitle='Counts'
oplot,  relevant_strips_n, above_t1_n, thick=th, col=3, psym=10

al_legend, ['Primary Threshold', 'Secondary Threshold'], textcol=[3,0]


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
spawn, 'open beam_spread.ps'

;=========================================================================================
;=========================================================================================












stop
end