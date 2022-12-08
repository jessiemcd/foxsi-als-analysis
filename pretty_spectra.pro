PRO pretty_spectra, files=files, beam=beam, side=side

;GOAL: MAKE SPECTRUM PLOT INCLUDING
;-ALL PSIDE ENERGIES
;-SINGLE PSIDE ENERGIES
;-DOUBLE PSIDE ENERGIES (COMPONENTS)
;-> DOUBLE PSIDE ENERGIES (CSA)

;Only for making spectra based on p-side event selection so far (pretty_spectra_movie.pro does both sides)

default, files, 'eventwisesp_20190419_205914.txt'
default, side, 'pside'
default, beam, 7


number = n_elements(files)-1
if number GT 14 then begin
	number=14
	print, 'More than expected # of files, only the first few will plot!'
endif

;=========================================================================================
;=========================================================================================

;For plotting the different spectra in a nice grid
positions = [[0.07,0.85,0.3,1],[0.4,0.85,0.63,1],[0.73,0.85,0.96,1], $
			[0.07,0.65,0.3,0.8], [0.4,0.65,0.63,0.8], [0.73,0.65,0.96,0.8], $ 
			[0.07,0.45,0.3,0.6], [0.4,0.45,0.63,0.6], [0.73,0.45,0.96,0.6], $
			[0.07,0.25,0.3,0.4], [0.4,0.25,0.63,0.4], [0.73,0.25,0.96,0.4], $
			[0.07,0.05,0.3,0.2], [0.4,0.05,0.63,0.2], [0.73,0.05,0.96,0.2]]


popen, 'pretty_spectrum.ps', $
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

lables = ['All Events', 'Single Strip', 'Double Strip', 'Double CSA']
colors = [0, 40, 80, 180]



;=========================================================================================
;=========================================================================================


for j=0, number do begin
	index = j
	;index = j+60 ;center example, p-side Lbranch
	;index = j+49 ;boundary example, p-side Lbranch
	;index = j+71 ;boundary example, n-side Lbranch
	;index = j+82 ;center example, n-side branch
	;index=j+22
	restore, files[index]

	data = eventwise_spectra.EVENTDATA_SPEC
	pside_data = data[2:3, *, *]
	;All non-zero hits on p-side: to be 'total' spectrum
	nz_pside = pside_data[where(pside_data NE 0.)]


	SINGLE_FRAMES = []
	DOUBLE_FRAMES = []
	double_csa_energies = []

	if side EQ 'pside' then begin
		pside_data = data[2:3, *, *]
		;All non-zero hits on p-side: to be 'total' spectrum
		nz = pside_data[where(pside_data NE 0.)]
		for i=0, n_elements(data[0,0,*])-1 do begin
			frame = data[*,*,i]
			;nside = where(frame[0:1,*] NE 0.)
			pside = where(frame[2:3,*] NE 0.)
			if n_elements(pside) EQ 1 and total(pside) NE -1. then begin
				SINGLE_FRAMES = [SINGLE_FRAMES, i]
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
	endif
	
	if side EQ 'nside' then begin
		nside_data = data[0:1, *, *]
		;All non-zero hits on p-side: to be 'total' spectrum
		nz = nside_data[where(nside_data NE 0.)]
		for i=0, n_elements(data[0,0,*])-1 do begin
			frame = data[*,*,i]
			nside = where(frame[0:1,*] NE 0.)
			if n_elements(nside) EQ 1 and total(nside) NE -1. then begin
				SINGLE_FRAMES = [SINGLE_FRAMES, i]
			endif
			if n_elements(nside) EQ 2 then begin
				as0 = where(frame[0,*] NE 0.)
				if n_elements(as0) EQ 2 then begin
					if as0[1]-as0[0] EQ 1 then begin
						DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
						double_csa_energies = [double_csa_energies, total(frame[0,*])]
					endif
				endif
				as1 = where(frame[1,*] NE 0.)
				if n_elements(as1) EQ 2 then begin
					if as1[1]-as1[0] EQ 1 then begin
						DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
						double_csa_energies = [double_csa_energies, total(frame[1,*])]
					endif
				endif		
			endif
		endfor
		singles_data = data[0:1,*, SINGLE_FRAMES]
		doubles_data = data[0:1,*, DOUBLE_FRAMES]
	endif


	;print, size(SINGLE_FRAMES)
;	print, size(DOUBLE_FRAMES)
;
;	print, eventwise_spectra.singles
;	print, eventwise_spectra.doubles
;	print, ''

	

	nz_singles = singles_data[where(singles_data NE 0.)]
	nz_doubles = doubles_data[where(doubles_data NE 0.)]
	
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
    ;    strmid(filesplit[2],0,2)+':'+strmid(filesplit[2],2,2)+':'+strmid(filesplit[2],4,2)
    
    date = strmid(filesplit[2],0,4)+'/'+strmid(filesplit[2],4,2)+'/'+strmid(filesplit[2],6,2)+' '+ $
        strmid(filesplit[3],0,2)+':'+strmid(filesplit[3],2,2)+':'+strmid(filesplit[3],4,2)
	
	;Plotting all spectra with reference lines at expected beam location + harmonics
	plot, pbins, p_hist, xrange=[0.,35], yrange=[1., 1500.], thick=th, xthick=th, ythick=th, font=-1, $
			charthick=fth, charsize=charsize, title=date, xtitle='Energy (keV)', ytitle='Counts/bin', $
			pos=positions[*,j], /ylog, psym=10, xstyle=1
	oplot, sbins, s_hist, color=40, thick=th, psym=10
	oplot, dbins, d_hist, color=80, thick=th, psym=10
	oplot, cbins, csa_hist, color=180, thick=th, psym=10
	oplot, [beam,beam,beam], [1.,100, 10000.], linestyle=2
	oplot, [beam*2,beam*2], [1,10000], linestyle=2
	oplot, [beam*3,beam*3], [1,10000], linestyle=2
	oplot, [beam*4,beam*4], [1,10000], linestyle=2
	if j EQ 0 then al_legend, lables, textcol=colors, box=1, /right, /top, /clear, charsi=0.6

endfor

;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
;DEVICE, /CLOSE
spawn, 'open pretty_spectrum.ps'

;=========================================================================================
;=========================================================================================

;STOP
END