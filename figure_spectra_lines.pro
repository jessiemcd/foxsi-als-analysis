PRO figure_spectra_lines

;plot four ALS files' spectra and fit their lines in energy space. Like a cleaner, more focused version of 
;energy_space_strip_spectra.pro. Hard coded (obv).

deffiles = ['eventwisesp_boutique_20190419_224006.txt', 'eventwisesp_boutique_20190419_230246.txt', $
			'eventwisesp_boutique_20190420_061854.txt', 'eventwisesp_boutique_20190420_064131.txt']


files = deffiles
number=n_elements(files)-1

beams=[7.,7.,7., 7.]
sides=['pside','pside', 'nside', 'nside']
fitranges1 = [[-3,3], [-3,3], [-2,3], [-3,3]]
fitranges3 = [[-3,2], [-3,+3], [-3,3], [-2,3]]

;=========================================================================================
;=========================================================================================

;For plotting the different spectra in a nice grid
positions = [[0.07,0.7,0.3,1],[0.4,0.7,0.63,1], $
			[0.07,0.3,0.3,0.6], [0.4,0.3,0.63,0.6]]


popen, 'line_spectrum.ps', $
		xsi=8, ysi=10
!Y.margin=4.
!X.margin=4.
ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
!p.multi=[0,1,40]
loadct, 2



;=========================================================================================
;=========================================================================================


for j=0, number do begin
	index = j
	restore, files[index]
	side=sides[index]
	beam=beams[index]
	
	SINGLE_FRAMES = []

	data = eventwise_spectra.EVENTDATA_SPEC
	if side EQ 'pside' then begin
		pside_data = data[2:3, *, *]
		;All non-zero hits on p-side: to be 'total' spectrum
		nz = pside_data[where(pside_data NE 0.)]
		
		strips=[]
		for i=0, n_elements(data[2,0,*])-1 do begin
			evt = data[2,*,i]
			pside = where(evt NE 0.)
			if n_elements(pside) EQ 1 and total(pside) NE -1. then SINGLE_FRAMES = [SINGLE_FRAMES, i]
			strips = [strips, where(evt EQ max(evt))]
		endfor
		
		strip = mode_val(strips)
		print, strip
	endif
	
	if side EQ 'nside' then begin
		nside_data = data[0:1, *, *]
		;All non-zero hits on p-side: to be 'total' spectrum
		nz = nside_data[where(nside_data NE 0.)]
		
		strips=[]
		for i=0, n_elements(data[1,0,*])-1 do begin
			evt = data[1,*,i]
			nside = where(evt NE 0.)
			if n_elements(nside) EQ 1 and total(nside) NE -1. then SINGLE_FRAMES = [SINGLE_FRAMES, i]
			strips = [strips, where(evt EQ max(evt))]
		endfor
		
		strip = mode_val(strips)
		print, strip
	endif
	
	binsize=0.5
	;All event spectrum histogram
	p_hist = histogram(nz, binsize=binsize, nbins=35*4, min=0, locations=pbins)
	
	if side EQ 'pside' then singles_data = data[2:3,*, SINGLE_FRAMES]
	if side EQ 'nside' then singles_data = data[0:1,*, SINGLE_FRAMES]
	nz_singles = singles_data[where(singles_data NE 0.)]
	
	;for single strip events specifically
	p_hist = histogram(nz_singles, binsize=binsize, nbins=35*4, min=0, locations=pbins)

	;Getting the date and time from the file title
	filesplit = strsplit( files[index], '_.', /extract) 
	date = strmid(filesplit[2],0,4)+'/'+strmid(filesplit[2],4,2)+'/'+strmid(filesplit[2],6,2)+' '+ $
        strmid(filesplit[3],0,2)+':'+strmid(filesplit[3],2,2)+':'+strmid(filesplit[3],4,2)
	
	;Plotting all spectra with reference lines at expected beam location + harmonics
	plot, pbins, p_hist, xrange=[0.,35], yrange=[1., 1500.], thick=th, xthick=th, ythick=th, font=-1, $
			charthick=fth, charsize=charsize, title=date, xtitle='Energy (keV)', ytitle='Counts/bin', $
			pos=positions[*,j], /ylog, psym=10, xstyle=1
	oplot, [beam,beam,beam], [1.,100, 10000.], linestyle=2
	oplot, [beam*2,beam*2], [1,10000], linestyle=2
	oplot, [beam*3,beam*3], [1,10000], linestyle=2
	oplot, [beam*4,beam*4], [1,10000], linestyle=2
	
	fitrange = fitranges1[*,j]
	print, fitrange
	
	fitbins = pbins[where(pbins GE beam+fitrange[0])]
	fitbins2 = fitbins[where(fitbins LE beam+fitrange[1])]

	fithist = p_hist[where(pbins GE beam+fitrange[0])]
	fithist2 = fithist[where(fitbins LE beam+fitrange[1])]
	
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	
	oplot, fitbins2, fithist, thick=th, col=80,  psym=10
	oplot, fitbins2, fit, thick=th, col=80, linestyle=2
	
	print, 'First Fit Range Coefficients: ', coeff
	print, 'FWHM: ', 2.355*coeff[2]
	
	labels = [side, 'Strip: '+strtrim(strip,2), 'Center: '+strtrim(coeff[1],2), 'Diff: '+strtrim((coeff[1]-beam), 2)]
	colors = [0, 0, 80, 80]
	
	fitrange = fitranges3[*,j]
	;print, fitrange
	
	fitbins = pbins[where(pbins GE beam*3+fitrange[0])]
	fitbins2 = fitbins[where(fitbins LE beam*3+fitrange[1])]

	fithist = p_hist[where(pbins GE beam*3+fitrange[0])]
	fithist2 = fithist[where(fitbins LE beam*3+fitrange[1])]
	
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)
	
	print, 'Second Fit Range Coefficients: ', coeff
	print, 'FWHM: ', 2.355*coeff[2]
	
	labels= [labels,'Center: '+strtrim(coeff[1],2), 'Diff: '+strtrim((coeff[1]-beam*3), 2)]
	colors=[colors,180, 180]
	
	;if j EQ 1 then labels[-1] = 'Diff: '+strtrim((coeff[1]-beam*2), 2)
	
	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	
	al_legend, labels, textcol=colors, /top,/right, charsi=0.7,/clear 
	
	print, ''
	
	
	
	

endfor

;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
;DEVICE, /CLOSE
spawn, 'open line_spectrum.ps'

;=========================================================================================
;=========================================================================================

;STOP
END