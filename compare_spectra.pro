PRO compare_spectra, FILES=FILES, edit_correction=edit_correction

;TAKES IN 4 FILES ONLY, AND THEN MAKES PRETTY SPECTRA OF EACH IN ONE FIGURE. 
;FIRST TWO FILES ARE P-SIDE, SECOND TWO ARE NSIDE. 

;One center and one boundary in each L-scan branch, for examining calibration.
;deffiles = ['eventwisesp_20190419_224006.txt', 'eventwisesp_20190419_230246.txt', $
;			'eventwisesp_20190420_061854.txt', 'eventwisesp_20190420_064131.txt']

deffiles = ['eventwisesp_boutique_20190419_224006.txt', 'eventwisesp_boutique_20190419_230246.txt', $
			'eventwisesp_boutique_20190420_061854.txt', 'eventwisesp_boutique_20190420_063927.txt']

;n-side boundary slightly farther from actual boundary
;deffiles = ['eventwisesp_boutique_20190419_224006.txt', 'eventwisesp_boutique_20190419_230246.txt', $
;			'eventwisesp_boutique_20190420_061854.txt', 'eventwisesp_boutique_20190420_064131.txt']

default, files, deffiles

beam=7.


;Has functionality to add an energy-space correction to the spectra on the Al-side to have the CSA 
;peaks at the right energies (rather than the single-strip peaks used in calibration), under the 
;assumption that they should be the peaks actually registering the correct energy (due to sub-threshold
;charge loss in "single strip" events).
default, edit_correction, 1

if edit_correction EQ 1 then print, 'In mode for checking on the CSA-peak gain correction on Al-side'
if edit_correction EQ 0 then print, 'CSA-peak gain correction has been applied on Al-side!'
;=========================================================================================
;=========================================================================================


popen, 'pretty_spectrum.ps', $
		xsi=6, ysi=10
!Y.margin=5.
!X.margin=5.
ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
!p.multi=[0,2,3]
loadct, 2

;=========================================================================================
;=========================================================================================


lables = ['All Events', 'Single Strip', 'Double Strip', 'Double CSA']
colors = [0, 40, 80, 180]


for j=0, 1 do begin
	index = j
	restore, files[index]

	data = eventwise_spectra.EVENTDATA_SPEC
	pside_data = data[2:3, *, *]
	;All non-zero hits on p-side: to be 'total' spectrum
	nz_pside = pside_data[where(pside_data NE 0.)]


	SINGLE_FRAMES = []
	DOUBLE_FRAMES = []
	double_csa_energies = []

	pside_data = data[2:3, *, *]
	;All non-zero hits on p-side: to be 'total' spectrum
	nz = pside_data[where(pside_data NE 0.)]
	for i=0, n_elements(data[0,0,*])-1 do begin
		frame = data[*,*,i]
		;nside = where(frame[0:1,*] NE 0.)
		pside = where(frame[2:3,*] NE 0.)
		if n_elements(pside) EQ 1 and total(pside) NE -1. then begin
			SINGLE_FRAMES = [SINGLE_FRAMES, i]
			;print, pside
		endif
		if n_elements(pside) EQ 2 then begin
			as2 = where(frame[2,*] NE 0.)
			if n_elements(as2) EQ 2 then begin
				if as2[1]-as2[0] EQ 1 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies = [double_csa_energies, total(frame[2,*])]
				endif
			endif
			as3 = where(frame[3,*] NE 0.)
			if n_elements(as3) EQ 2 then begin
				if as3[1]-as3[0] EQ 1 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies = [double_csa_energies, total(frame[3,*])]
				endif
			endif		
		endif
	endfor
	singles_data = data[2:3,*, SINGLE_FRAMES]
	doubles_data = data[2:3,*, DOUBLE_FRAMES]
	
	
	nz_singles = singles_data[where(singles_data NE 0.)]
	nz_doubles = doubles_data[where(doubles_data NE 0.)]
	
	print, 'Min Doubles:', min(nz_doubles)
	
	binsize=0.5
	;Single event spectrum histogram
	s_hist = histogram(nz_singles, binsize=binsize, nbins=35*4, min=0, locations=sbins)
	;Double event spectrum histogram
	d_hist = histogram(nz_doubles, binsize=binsize, nbins=35*4, min=0, locations=dbins)
	;All event spectrum histogram
	p_hist = histogram(nz, binsize=binsize, nbins=35*4, min=0, locations=pbins)
	;CSA histogram from double pixel sums
	csa_hist = histogram(double_csa_energies, binsize=binsize, nbins=35*4, min=0, locations=cbins)
	
	;Getting the date and time from the file title
	filesplit = strsplit( files[index], '_.', /extract) 
	;date = strmid(filesplit[1],0,4)+'/'+strmid(filesplit[1],4,2)+'/'+strmid(filesplit[1],6,2)+' '+ $
     ;   strmid(filesplit[2],0,2)+':'+strmid(filesplit[2],2,2)+':'+strmid(filesplit[2],4,2)
    
    date = strmid(filesplit[2],0,4)+'/'+strmid(filesplit[2],4,2)+'/'+strmid(filesplit[2],6,2)+' '+ $
        strmid(filesplit[3],0,2)+':'+strmid(filesplit[3],2,2)+':'+strmid(filesplit[3],4,2)
	
	;Plotting all spectra with reference lines at expected beam location + harmonics
	plot, pbins, p_hist, xrange=[0.,35], yrange=[1., 1500.], thick=th, xthick=th, ythick=th, font=-1, $
			charthick=fth, charsize=charsize, title=date, xtitle='Energy (keV)', ytitle='Counts/bin', $
			/ylog, psym=10, xstyle=1
			;pos=positions[*,j], /ylog, psym=10, xstyle=1
	oplot, sbins, s_hist, color=40, thick=th, psym=10
	oplot, dbins, d_hist, color=80, thick=th, psym=10
	oplot, cbins, csa_hist, color=180, thick=th, psym=10
	oplot, [beam,beam,beam], [1.,100, 10000.], linestyle=2
	oplot, [beam*2,beam*2], [1,10000], linestyle=2
	oplot, [beam*3,beam*3], [1,10000], linestyle=2
	oplot, [beam*4,beam*4], [1,10000], linestyle=2
	al_legend, lables, textcol=colors, box=1, /right, /top, /clear, charsi=0.6
endfor

energiez = []
energies = [7.,21.,28.,7.,21.,28.]

for j=2, 3 do begin
	index = j
	restore, files[index]

	data = eventwise_spectra.EVENTDATA_SPEC
	pside_data = data[2:3, *, *]
	;All non-zero hits on p-side: to be 'total' spectrum
	nz_pside = pside_data[where(pside_data NE 0.)]


	SINGLE_FRAMES = []
	DOUBLE_FRAMES = []
	double_csa_energies1 = []
	double_csa_energies2 = []

	nside_data = data[0:1, *, *]
	;All non-zero hits on p-side: to be 'total' spectrum
	nz = nside_data[where(nside_data NE 0.)]
	for i=0, n_elements(data[0,0,*])-1 do begin
		frame = data[*,*,i]
		nside = where(frame[0:1,*] NE 0.)
		if n_elements(nside) EQ 1 and total(nside) NE -1. then begin
			SINGLE_FRAMES = [SINGLE_FRAMES, i]
			;print, nside
		endif
		if n_elements(nside) EQ 2 then begin
			as0 = where(frame[0,*] NE 0.)
			fframe = frame[0,*]
			if n_elements(as0) EQ 2 then begin
				if as0[1]-as0[0] EQ 1 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies1 = [double_csa_energies1, fframe[as0[0]]]
					double_csa_energies2 = [double_csa_energies2, fframe[as0[1]]]
				endif
			endif
			as1 = where(frame[1,*] NE 0.)
			fframe = frame[1,*]
			if n_elements(as1) EQ 2 then begin
				if as1[1]-as1[0] EQ 1 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies1 = [double_csa_energies1, fframe[as1[0]]]
					double_csa_energies2 = [double_csa_energies2, fframe[as1[1]]]
				endif
			endif		
		endif
	endfor
	singles_data = data[0:1,*, SINGLE_FRAMES]
	doubles_data = data[0:1,*, DOUBLE_FRAMES]
	
	double_csa_energies = double_csa_energies1+double_csa_energies2
	
	
	nz_singles = singles_data[where(singles_data NE 0.)]
	nz_doubles = doubles_data[where(doubles_data NE 0.)]
	
	print, 'Min Doubles:', min(nz_doubles)
	
	;CSA histogram from double pixel sums
	csa_hist = histogram(double_csa_energies, binsize=binsize, nbins=35*4, min=0, locations=cbins)
	
	;Getting the date and time from the file title
	filesplit = strsplit( files[index], '_.', /extract) 
	;date = strmid(filesplit[1],0,4)+'/'+strmid(filesplit[1],4,2)+'/'+strmid(filesplit[1],6,2)+' '+ $
     ;   strmid(filesplit[2],0,2)+':'+strmid(filesplit[2],2,2)+':'+strmid(filesplit[2],4,2)
    
    date = strmid(filesplit[2],0,4)+'/'+strmid(filesplit[2],4,2)+'/'+strmid(filesplit[2],6,2)+' '+ $
        strmid(filesplit[3],0,2)+':'+strmid(filesplit[3],2,2)+':'+strmid(filesplit[3],4,2)
        
    binsize=0.5
	;Single event spectrum histogram
	s_hist = histogram(nz_singles, binsize=binsize, nbins=35*4, min=0, locations=sbins)
	;Double event spectrum histogram
	d_hist = histogram(nz_doubles, binsize=binsize, nbins=35*4, min=0, locations=dbins)
	;All event spectrum histogram
	p_hist = histogram(nz, binsize=binsize, nbins=35*4, min=0, locations=pbins)
	;CSA histogram from double pixel sums
	csa_hist = histogram(double_csa_energies, binsize=binsize, nbins=35*4, min=0, locations=cbins)
	
	if edit_correction EQ 1 then begin
	
		;Plotting all spectra with reference lines at expected beam location + harmonics
		plot, pbins, p_hist, xrange=[0.,35], yrange=[1., 1500.], thick=th, xthick=th, ythick=th, font=-1, $
				charthick=fth, charsize=charsize, title=date+' fit CSA peak values', xtitle='Energy (keV)', ytitle='Counts/bin', $
				/ylog, psym=10, xstyle=1
				;pos=positions[*,j], /ylog, psym=10, xstyle=1
		oplot, sbins, s_hist, color=40, thick=th, psym=10
		oplot, dbins, d_hist, color=80, thick=th, psym=10
		oplot, cbins, csa_hist, color=180, thick=th, psym=10
	
	
		;If uncommented, the overplots in this loop will plot in the wrong window since the 
		;next spectrum has not been plotted yet.
		COLL=40
		maxbin = cbins[where(csa_hist EQ max(csa_hist))]
		maxes = lclxtrem(csa_hist, 5, /maxima)

		iv=3

		fitbins = cbins[where(cbins GE maxbin[0]-iv)]
		fitbins2 = fitbins[where(fitbins LE maxbin[0]+iv)]
		fithist = csa_hist[where(cbins GE maxbin[0]-iv)]
		fithist2 = fithist[where(fitbins LE maxbin[0]+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, 2000.], linestyle=2, col=100, thick=th


		print, 'Sevens: ', 7., coeff[1], 7.-coeff[1]
		energiez = [energiez, coeff[1]]

		fitbins = sbins[where(cbins GE cbins[maxes[1]]-iv)]
		fitbins2 = fitbins[where(fitbins LE cbins[maxes[1]]+iv)]
		fithist = csa_hist[where(cbins GE cbins[maxes[1]]-iv)]
		fithist2 = fithist[where(fitbins LE cbins[maxes[1]]+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, 2000.], linestyle=2, col=100, thick=th


		print, 'Twenty-Ones: ', 21, coeff[1], 21-coeff[1]
		energiez = [energiez, coeff[1]]
		
		bin28 = maxes[2]
		if j EQ 3 then bin28 = maxes[3]

		fitbins = cbins[where(cbins GE cbins[bin28]-iv)]
		fitbins2 = fitbins[where(fitbins LE cbins[bin28]+iv)]
		fithist = csa_hist[where(cbins GE cbins[bin28]-iv)]
		fithist2 = fithist[where(fitbins LE cbins[bin28]+iv)]
		fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

		oplot, fitbins2, fithist, thick=th, col=180,  psym=10
		oplot, fitbins2, fit, thick=th, col=180, linestyle=2
		oplot, [coeff[1],coeff[1]], [0.1, 2000.], linestyle=2, col=100, thick=th


		print, 'Twenty-Eights: ', 28, coeff[1], 28-coeff[1]
		energiez = [energiez, coeff[1]]
	endif
	if edit_correction EQ 0 then begin
		restore, 'CSA-based_correction.sav'
		
		;The logic here is as follows
		;
		;Correct CSA energy = Original_CSA_energy*M + B (M, B found from linear CSA peaks fit)
		;Original_CSA_energy = E1+E2
		;Correct CSA energy = (E1+E2)*M + B = E1*M + B/2 + E2*M + B/2
		
		;Define correction for ALL energies recorded as,
		;Corrected energy = E*M + B/2
		;This way we are not incorrectly "doubling the pedestal" by using CSA-line correction
		;on events which don't consist of two components.
		
		nz_singles = nz_singles*linfit_strip[1]+linfit_strip[0]/2
		nz_doubles = nz_doubles*linfit_strip[1]+linfit_strip[0]/2
		nz = nz*linfit_strip[1]+linfit_strip[0]/2
		double_csa_energies1 = double_csa_energies1*linfit_strip[1]+linfit_strip[0]/2
		double_csa_energies2 = double_csa_energies2*linfit_strip[1]+linfit_strip[0]/2
		double_csa_energies = double_csa_energies1+double_csa_energies2
	endif
	
	print, 'Min Doubles:', min(nz_doubles)
	
	
	binsize=0.5
	;Single event spectrum histogram
	s_hist = histogram(nz_singles, binsize=binsize, nbins=35*4, min=0, locations=sbins)
	;Double event spectrum histogram
	d_hist = histogram(nz_doubles, binsize=binsize, nbins=35*4, min=0, locations=dbins)
	;All event spectrum histogram
	p_hist = histogram(nz, binsize=binsize, nbins=35*4, min=0, locations=pbins)
	;CSA histogram from double pixel sums
	csa_hist = histogram(double_csa_energies, binsize=binsize, nbins=35*4, min=0, locations=cbins)
	
	
	;Plotting all spectra with reference lines at expected beam location + harmonics
	plot, pbins, p_hist, xrange=[0.,35], yrange=[1., 1500.], thick=th, xthick=th, ythick=th, font=-1, $
			charthick=fth, charsize=charsize, title=date, xtitle='Energy (keV)', ytitle='Counts/bin', $
			/ylog, psym=10, xstyle=1
			;pos=positions[*,j], /ylog, psym=10, xstyle=1
	oplot, sbins, s_hist, color=40, thick=th, psym=10
	oplot, dbins, d_hist, color=80, thick=th, psym=10
	oplot, cbins, csa_hist, color=180, thick=th, psym=10



	oplot, [beam,beam,beam], [1.,100, 10000.], linestyle=2
	oplot, [beam*2,beam*2], [1,10000], linestyle=2
	oplot, [beam*3,beam*3], [1,10000], linestyle=2
	oplot, [beam*4,beam*4], [1,10000], linestyle=2
	al_legend, lables, textcol=colors, box=1, /right, /top, /clear, charsi=0.6
endfor

if edit_correction EQ 1 then begin
	print, energies
	print, energiez
	plot, energiez, energies, thick=th, psym=4, xrange=[0,32], yrange=[0,32], $
			title='Fit CSA peaks vs. line energies - correction function', $
			xtitle='Fit Peak Energies', ytitle='Expected Lines'
	;Linear fit to values
	linfit_strip = LINFIT(energiez, energies)
	;print, 'Linfit parameters n:', linfit_strip
	linfit_arr = energiez*linfit_strip[1]+linfit_strip[0]
	oplot, energiez, linfit_arr, thick=th/2., linestyle=2, color=coll
	print, linfit_strip

	save, linfit_strip, filename='CSA-based_correction.sav'
endif

;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
;DEVICE, /CLOSE
spawn, 'open pretty_spectrum.ps'

;=========================================================================================
;=========================================================================================



STOP

END