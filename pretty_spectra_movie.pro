PRO pretty_spectra_movie, files=files, nside_on=nside_on, beam=beam

;GOAL: MAKE SPECTRUM PLOT INCLUDING
;-ALL PSIDE ENERGIES
;-SINGLE PSIDE ENERGIES
;-DOUBLE PSIDE ENERGIES (COMPONENTS)
;-> DOUBLE PSIDE ENERGIES (COMPONENTS)

;ALSO HAVE OPTION TO SELECT N-SIDE AND DO THE SAME THING THERE

default, nside_on, 0
default, beam, 7

lables = ['All Events', 'Single Strip', 'Double Strip', 'Double CSA']
colors = [255, 40, 80, 180]


compile_opt idl2
data = dist(30)

;Name video file
if nside_on EQ 0 then video_file = 'beam_travel_pside.mp4'
if nside_on EQ 1 then video_file = 'beam_travel_nside.mp4'

;Code for making video taken from "Making Movies with IDL, Part 1"
;https://www.l3harrisgeospatial.com/Learn/Blogs/Blog-Details/ArtMID/10198/ArticleID/16793/Making-movies-with-IDL-part-I
video = idlffvideowrite(video_file)
framerate = 5
framedims = [640,512]
stream = video.addvideostream(framedims[0], framedims[1], framerate)


th=6
loadct, 2
set_plot, 'z', /copy
device, set_resolution=framedims, set_pixel_depth=24, decomposed=0
nframes = n_elements(files)
for j=0, nframes-1 do begin
	restore, files[j]

	data = eventwise_spectra.EVENTDATA_SPEC
	
	;Get all non-zero energies from pside, nside
	pside_data = data[2:3, *, *]
	nz_pside = pside_data[where(pside_data NE 0.)]
	
	nside_data = data[0:1, *, *]
	nz_nside = nside_data[where(nside_data NE 0.)]


	SINGLE_FRAMES = []
	DOUBLE_FRAMES = []
	double_csa_energies = []


	for i=0, n_elements(data[0,0,*])-1 do begin
		frame = data[*,*,i]
		;nside = where(frame[0:1,*] NE 0.)
		side = where(frame[2:3,*] NE 0.)
		if nside_on EQ 1 then side = where(frame[0:1,*] NE 0.)
		
		;Finding single strip events
		if n_elements(side) EQ 1 and total(side) NE -1. then begin
			SINGLE_FRAMES = [SINGLE_FRAMES, i]
		endif
		
		;Depending on selected side, finding double strip events (and making CSA array)
		if n_elements(side) EQ 2 then begin
			if nside_on EQ 0 then begin
				as2 = where(frame[2,*] NE 0.)
				if n_elements(as2) EQ 2 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies = [double_csa_energies, total(frame[2,*])]
				endif
				as3 = where(frame[3,*] NE 0.)
				if n_elements(as3) EQ 2 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies = [double_csa_energies, total(frame[3,*])]
				endif
			endif
			if nside_on EQ 1 then begin
				as0 = where(frame[0,*] NE 0.)
				if n_elements(as2) EQ 2 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies = [double_csa_energies, total(frame[0,*])]
				endif
				as1 = where(frame[1,*] NE 0.)
				if n_elements(as1) EQ 2 then begin
					DOUBLE_FRAMES = [DOUBLE_FRAMES, i]
					double_csa_energies = [double_csa_energies, total(frame[1,*])]
				endif
			endif		
		endif
	endfor


;	print, size(SINGLE_FRAMES)
;	print, size(DOUBLE_FRAMES)
;
;	print, eventwise_spectra.singles
;	print, eventwise_spectra.doubles
;	print, ''

	singles_data = data[2:3,*, SINGLE_FRAMES]
	doubles_data = data[2:3,*, DOUBLE_FRAMES]
	
	if nside_on EQ 1 then begin
		singles_data = data[0:1,*, SINGLE_FRAMES]
		doubles_data = data[0:1,*, DOUBLE_FRAMES]
	endif

	nz_singles = singles_data[where(singles_data NE 0.)]
	nz_doubles = doubles_data[where(doubles_data NE 0.)]
	
	;Making spectra (histograms)
	s_hist = histogram(nz_singles, binsize=0.5, nbins=35*4, min=0, locations=sbins)
	d_hist = histogram(nz_doubles, binsize=0.5, nbins=35*4, min=0, locations=dbins)
	main_hist = histogram(nz_pside, binsize=0.5, nbins=35*4, min=0, locations=mbins)
	if nside_on EQ 1 then main_hist = histogram(nz_nside, binsize=0.5, nbins=35*4, min=0, locations=mbins)
	csa_hist = histogram(double_csa_energies, binsize=0.5, nbins=35*4, min=0, locations=cbins)
	
	;Getting the date and time from the file title
	filesplit = strsplit( files[j], '_.', /extract) 
	date = strmid(filesplit[1],0,4)+'/'+strmid(filesplit[1],4,2)+'/'+strmid(filesplit[1],6,2)+' '+ $
        strmid(filesplit[2],0,2)+':'+strmid(filesplit[2],2,2)+':'+strmid(filesplit[2],4,2)
	
	;Plotting all spectra with reference lines at expected beam location + harmonics
	plot, mbins, main_hist, xrange=[0.,30.], yrange=[1., 1500.], thick=th, xthick=th, ythick=th, font=-1, $
			charthick=fth, charsize=charsize, title=date, xtitle='Energy (keV)', ytitle='Counts/bin', /ylog
			;pos=positions[*,j], /ylog
	oplot, sbins, s_hist, color=40, thick=th
	oplot, dbins, d_hist, color=80, thick=th
	oplot, cbins, csa_hist, color=180, thick=th
	oplot, [beam,beam,beam], [1.,100, 10000.], linestyle=2
	oplot, [beam*2,beam*2], [1,10000], linestyle=2
	oplot, [beam*3,beam*3], [1,10000], linestyle=2
	oplot, [beam*4,beam*4], [1,10000], linestyle=2
	al_legend, lables, textcol=colors, box=1, /right, /top, /clear, charsi=2, charth=3

	;Adding frame to movie
    timestamp = video.put(stream, tvrd(true=1))
endfor


device, /close
set_plot, strlowcase(!version.os_family) eq 'windows' ? 'win' : 'x'
video.cleanup
print, 'File "' + video_file + '" written to current directory.'
end


