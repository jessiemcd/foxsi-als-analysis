;"Charge Loss Region Indices" (N-side): 3,4,16,23,24,36,44,56,63,64,76,84,96
;"Charge Loss Region Indices" (P-side): 5,9,26,30,46,51,67,71,88,92,108,113,128,133,149,153,169,174

FUNCTION ratio_based_image, file=file, erange=erange, macrobins=macrobins

	;loadct, 1
	;reverse_ct

	default, macrobins, 0

	;We are currently assuming the beam is in ASICs 1 + 2
	asp = 2
	asn = 1
	
	;Energy threshold of detector
	thresh = 4.

	;For this method we will use 1 micron bins, AKA each strip has 60 bins
	p_crossing_width = 16.
	p_center = 60. - p_crossing_width
	n_crossing_width = 40.
	n_center = 60. - n_crossing_width
	
	;IN THIS NEW VERSION WE DO NOT ADD BACK (SEE energy_ratio_imaging.pro for this in use)
	;On n-side, single-strip events are really charge loss events this percent of the time:
	;(total below thresh)/(total below thresh + center)
	;n_losses = 12./32.
	
	;On p-side, single-strip events are really charge loss events this percent of the time:
	;p_losses = 5./50.
	
	
	restore, file
	data=eventwise_spectra.EVENTDATA_SPEC
	
	event_coords=[]

	
	singlecounts=0
	doublecounts=0
	
	initial_phits = []
	
	p_ers = []
	n_ers = []

	
	for i=0, n_elements(data[0,0,*])-1 do begin
		
		p_ran=[]
		n_ran=[]
	
		;For each event, find the MAXIMUM and TOTAL values on the p- and n-side 
		;ASICs selected
		mp=MAX(data[asp,*,i], p_strip)
		mn=MAX(data[asn,*,i], n_strip)
		pt=TOTAL(data[asp,*,i])
		nt=TOTAL(data[asn,*,i]) 
		
		;If the total on the n-side and the total on the p-side are both in the energy range
		if nt gt erange[0] and nt lt erange[1] and pt gt erange[0] and pt lt erange[1] then begin
		
			;Ultimately, we want to arrive at a set of coordinates, (p, n) defining where we are 
			;going to place this event, based on its properties and our charge-sharing motivated rules.
			
			;hit strips on the p-side
			phits = WHERE(data[asp,*,i] NE 0.)
			initial_phits = [initial_phits, n_elements(phits)]
			;If the event is single-strip on the p-side.
			if n_elements(phits) eq 1 then begin
				;Random position in center: strip bin + half crossing region width + random selection
				;within center region. 
				p_ran=p_strip*60.+p_crossing_width/2. + p_center*randomu(seed, 1)
				singlecounts+=1
				;For a certain percentage of events, assume charge loss. TURNED OFF IN THIS VERSION
;				if randomu(seed, 1) LE p_losses then begin
;					;Randomly assign either the left or right strip to the triggered one to have a new,
;					;sub-threshold energy.
;					if randomu(seed, 1) GE 0.5 then data[asp,p_strip+1,i] = randomu(seed, 1)*4.
;					if randomu(seed, 1) LT 0.5 then data[asp,p_strip-1,i] = randomu(seed, 1)*4.
;					singlecounts-=1
;				endif
			endif
			
			;print, 'P_Ran After Singles: ', p_ran
			
			;Re-defining, as some events may have become double strip events!
			phits = WHERE(data[asp,*,i] NE 0.)
			initial_phits = [initial_phits, n_elements(phits)]
			
			;If the event is double-strip on the p-side (and the second event occurs in an adjacent strip),
			;use the energy ratio to assign it to a specific region in the correct crossing region on one
			;side of the maximum-energy strip or the other. 
			if n_elements(phits) eq 2 then begin
				if data[asp,p_strip+1,i] NE 0 then begin
					ER = (data[asp,p_strip+1,i] - data[asp,p_strip,i])/(data[asp,p_strip+1,i] + data[asp,p_strip,i])
					p_ers = [p_ers, ER]
					p_ran=(p_strip+1)*60.+ER*p_crossing_width/2.
				endif
				if data[asp,p_strip-1,i] NE 0 then begin
					ER = (data[asp,p_strip,i] - data[asp,p_strip-1,i])/(data[asp,p_strip-1,i] + data[asp,p_strip,i])
					p_ers = [p_ers, ER]
					p_ran=(p_strip)*60.+ER*p_crossing_width/2.
				endif
				doublecounts+=1
			endif
			
			;print, 'P_Ran After Doubles: ', p_ran
			
		
			nhits = WHERE(data[asn,*,i] NE 0.)
			;If the event is single-strip on the n-side.
			if n_elements(nhits) eq 1 then begin
				;Random position in center: strip bin + half crossing region width + random selection
				;within center region. 
				n_ran=n_strip*60.+n_crossing_width/2. + n_center*randomu(seed, 1)
				singlecounts+=1
				;For a certain percentage of events, assume charge loss. TURNED OFF IN THIS VERSION
;				if randomu(seed, 1) LE n_losses then begin
;					;Randomly assign either the left or right strip to the triggered one to have a new,
;					;sub-threshold energy.
;					if randomu(seed, 1) GE 0.5 then data[asn,n_strip+1,i] = randomu(seed, 1)*4.
;					if randomu(seed, 1) LT 0.5 then data[asn,n_strip-1,i] = randomu(seed, 1)*4.
;				endif
			endif
			
			;print, 'N_Ran After Singles: ', n_ran
			
			;Re-defining, as some events may have become double strip events!
			nhits = WHERE(data[asn,*,i] NE 0.)
			
			;If the event is double-strip on the n-side (and the second event occurs in an adjacent strip),
			;use the energy ratio to assign it to a specific region in the correct crossing region on one
			;side of the maximum-energy strip or the other. 
			if n_elements(nhits) eq 2 then begin
				if data[asn,n_strip+1,i] NE 0 then begin
					ER = (data[asn,n_strip+1,i] - data[asn,n_strip,i])/(data[asn,n_strip+1,i] + data[asn,n_strip,i])
					n_ers = [n_ers, ER]
					n_ran=(n_strip+1)*60.+ER*n_crossing_width/2.
				endif
				if data[asn,n_strip-1,i] NE 0 then begin
					ER = (data[asn,n_strip,i] - data[asn,n_strip-1,i])/(data[asn,n_strip-1,i] + data[asn,n_strip,i])
					n_ers = [n_ers, ER]
					n_ran=(n_strip)*60.+ER*n_crossing_width/2.
				endif
				doublecounts+=1
			endif
			
			;print, 'N_Ran After Doubles: ', n_ran
			
			;print, 'Both at end: ', p_ran, n_ran
			;print, ''
			
			if p_ran ne !NULL and n_ran ne !NULL then event_coords=[[event_coords],[p_ran,n_ran]]
			
		endif	
	endfor
	
	;plot_hist, n_ers, bin=0.1;, xrange=[0.,4.]
	
	
	event_coords=transpose(event_coords)
	
	;print, event_coords[*,0]
	;print, event_coords[*,1]
	print, min(event_coords[*,0]), max(event_coords[*,0])
	print, min(event_coords[*,1]), max(event_coords[*,1])
	print, n_elements(data[0,0,*])
	

	min1=0.
    min2=0.
    max1=60.*64.
    max2=60.*64.
    bin1=1
    bin2=1
    
    print, 'Event coords size:', size(event_coords)
    
    ;stop

    h=hist2d(event_coords[*,0], event_coords[*,1], bin=[bin1,bin2],min=[min1,min2], max=[max1,max2])
    
    print, total(h)
    
    ;ranges for file 49
    ;plot_image, h, xrange=[2100., 2250.], $
    ;				yrange=[1500., 1650.]
    
    ;ranges for file 118
;    plot_image, h, xrange=[2300., 2450.], $
;    				yrange=[1450., 1600.]
;    
    				
;    lines= findgen(64)*60
;    for i=0, n_elements(lines)-1 do begin
;    	oplot, [lines[i],lines[i]], [0.,2000.], linestyle=2, thick=0.7
;    	oplot, [0.,3000.], [lines[i],lines[i]], linestyle=2, thick=0.7
;    endfor
    
    if macrobins EQ 1 then begin
		   print, "Total of histogram, pre averaging:", total(h)

		   ;artifically creating 'variable bin sizes' by averaging over bins in the non-double-strip-event region
		   ;averaging first in one direction, then the other. 
		   nh=float(h)
   
		   ;number of strips used in hist2D (max-min)
		   xx=max1-min1
		   yy=max2-min2
   
		   sh = size(nh)
   
		   half_p_crossing = round(p_crossing_width/2.)
		   ;print, half_p_crossing


;		   For every slice in p (column in plot)
		   for i=0, sh[1]-1 do begin
;			If there is at least one count in the column:
			if total(nh[*,i]) NE 0. then begin
				col = nh[*,i]
;				print, size(col)
;				print, n_elements(col)
;				for every strip in the column
				for j=0, n_elements(col)/60.-2 do begin
;					define one-strip array
					st = col[j*60:(j+1)*60-1]
;					if there are any counts in the array, continue
					if total(st) NE 0. then begin
;						find the middle region of the array
						mid = st[half_p_crossing+1:half_p_crossing+p_center]
;						Make a new array the same size containing the average in every position
						midarr = make_array(p_center)+total(mid)/p_center
;						print, midarr
;						Make the corresponding strip center in the full histogram that averaged array.
						nh[j*60+half_p_crossing+1:j*60+half_p_crossing+p_center,i] = midarr
					endif
				endfor
			endif
		   endfor


   
		   half_n_crossing = n_crossing_width/2
		   
		   ;For every slice in n (row in plot)
		   for i=0, sh[2]-1 do begin
			;If there is at least one count in the row:
			if total(nh[i,*]) NE 0. then begin
				row = nh[i,*]
				;print, size(row)
				;print, n_elements(row)
				;for every strip in the row
				for j=0, n_elements(row)/60.-2 do begin
					;define one-strip array
					st = row[j*60:(j+1)*60-1]
					;if there are any counts in the array, continue
					if total(st) NE 0. then begin
						;find the middle region of the array
						mid = st[half_n_crossing+1:half_n_crossing+n_center]
						;Make a new array the same size containing the average in every position
						midarr = make_array(n_center)+total(mid)/n_center
						;print, midarr
						;Make the corresponding strip center in the full histogram that averaged array.
						nh[i, j*60+half_n_crossing+1:j*60+half_n_crossing+n_center] = midarr
					endif
				endfor
			endif
		   endfor
		
		   
		   ;ranges for file 49
		   ;plot_image, nh, xrange=[2100., 2250.], $
;		   				yrange=[1500., 1650.]
;		   
		   ;ranges for file 118
		   ;plot_image, nh, xrange=[2300., 2450.], $
		   ;				yrange=[1450., 1600.]
		   				
;		   for i=0, n_elements(lines)-1 do begin
;		   	oplot, [lines[i],lines[i]], [0.,2000.], linestyle=2, thick=0.7
;		   	oplot, [0.,3000.], [lines[i],lines[i]], linestyle=2, thick=0.7
;		   endfor
		   
		   h=nh
	endif
    


	RETURN, h
	
END









FUNCTION zone_based_image, file=file, erange=erange
	
	;We are currently assuming the beam is in ASICs 1 + 2
	asp = 2
	asn = 1

	restore, file
	data=eventwise_spectra.EVENTDATA_SPEC

	event_coords=[]

	singlecounts=0
	doublecounts=0
	
	p_ran=[]
	n_ran=[]
	
	pwidth=20
	nwidth=42
	
	p_center = 60. - pwidth
	n_center = 60. - nwidth

	for i=0, n_elements(data[0,0,*])-1 do begin
		;For each event, find the MAXIMUM and TOTAL values on the p- and n-side ASICs selected
		mp=MAX(data[asp,*,i], p_strip)
		mn=MAX(data[asn,*,i], n_strip)
		pt=TOTAL(data[asp,*,i])
		nt=TOTAL(data[asn,*,i]) 
		
		;If the total on the n-side and the total on the p-side are both in the energy range
		if nt gt erange[0] and nt lt erange[1] and pt gt erange[0] and pt lt erange[1] then begin
		
			;Ultimately, we want to arrive at a set of coordinates, (p, n) defining where we are 
			;going to place this event, based on its properties and our charge-sharing motivated rules.
			
			phits = WHERE(data[asp,*,i] NE 0.)
			;If the event is single-strip on the p-side, randomly assign its location within the strip.
			if n_elements(phits) eq 1 then begin
				;p_ran=p_strip+randomu(seed,1)
				;p_ran=p_strip+0.1+(0.8*randomu(seed,1))
				p_ran=p_strip*60.+pwidth/2. + p_center*randomu(seed, 1)
				;singlecounts+=1
			endif
			
			;If the event is double-strip on the p-side (and the second event occurs in an adjacent strip),
			;randomly assign it within the appropriate double-strip event zone.
			if n_elements(phits) eq 2 then begin
				;if data[asp,p_strip+1,i] NE 0 then p_ran=p_strip+0.9+randomu(seed,1)/10.
				if data[asp,p_strip+1,i] NE 0 then p_ran=p_strip*60+pwidth/2+p_center+randomu(seed,1)*pwidth/2.
				;if data[asp,p_strip-1,i] NE 0 then p_ran=p_strip+randomu(seed,1)/10.
				if data[asp,p_strip-1,i] NE 0 then p_ran=p_strip*60+randomu(seed,1)*pwidth/2.
				;doublecounts+=1
			endif
			
			nhits = WHERE(data[asn,*,i] NE 0.)
			;If the event is single-strip on the n-side, randomly assign its location within the strip.
			if n_elements(nhits) eq 1 then begin
				;n_ran=n_strip+randomu(seed,1)
				;n_ran=n_strip+0.25+0.5*randomu(seed,1)
				n_ran=n_strip*60.+nwidth/2. + n_center*randomu(seed, 1)
				;singlecounts+=1
			endif
			
			;If the event is double-strip on the n-side (and the second event occurs in an adjacent strip),
			;randomly assign it within the appropriate double-strip event zone.
			if n_elements(nhits) eq 2 then begin
				if data[asn,n_strip+1,i] NE 0 then n_ran=n_strip*60+nwidth/2+n_center+randomu(seed,1)*nwidth/2.
				if data[asn,n_strip-1,i] NE 0 then n_ran=n_strip*60+randomu(seed,1)*nwidth/2.
				;doublecounts+=1
			endif
			
			
			if p_ran ne !NULL and n_ran ne !NULL then event_coords=[[event_coords],[p_ran,n_ran]]
			
		endif
	endfor
	
	event_coords=transpose(event_coords)

	min1=0.
    min2=0.
    max1=60.*64.
    max2=60.*64.
    bin1=1
    bin2=1
    
    print, 'Event coords size:', size(event_coords)
    
    ;stop

    h=hist2d(event_coords[*,0], event_coords[*,1], bin=[bin1,bin2],min=[min1,min2], max=[max1,max2])
    
    print, "Total of histogram, pre averaging:", total(h)

    ;artifically creating 'variable bin sizes' by averaging over bins in the non-double-strip-event region
    ;averaging first in one direction, then the other. 
    nh=float(h)

    ;number of strips used in hist2D (max-min)
    xx=max1-min1
    yy=max2-min2

	sh = size(nh)
   
	half_p_crossing = round(pwidth/2.)
	;print, half_p_crossing


	;For every slice in p (column in plot)
	for i=0, sh[1]-1 do begin
		;If there is at least one count in the column:
		if total(nh[*,i]) NE 0. then begin
			;isolate the column
			col = nh[*,i]
			;if total(col) GT 20 then print, col[35*60.:37*60.]
			;print, size(col)
			;print, n_elements(col)
			;for every strip in the column
			for j=0, n_elements(col)/60.-2 do begin
				;define one-strip array
				st = col[j*60:(j+1)*60-1]
				;if there are any counts in the one-strip array, continue
				if total(st) NE 0. then begin
					;print, st
					;find the middle region of the strip array
					mid = st[half_p_crossing+1:half_p_crossing+p_center]
					;print, n_elements(mid)
					;Make a new array the same size containing the average in every position
					midarr = make_array(p_center)+total(mid)/p_center
					;print, midarr
					;Make the corresponding strip center in the full histogram that averaged array.
					nh[j*60+half_p_crossing:j*60+half_p_crossing+p_center-1,i] = midarr
					;;;;print, n_elements(midarr)
				
					;find the first side region of the strip array, and repeat
					sid = st[0:half_p_crossing-1]
					;print, n_elements(sid)
					sidarr = make_array(half_p_crossing)+total(sid)/half_p_crossing
					nh[j*60:j*60+half_p_crossing-1,i] = sidarr
					;print, n_elements(sidarr)
					
					;find the second side region of the strip array, and repeat
					sid = st[half_p_crossing+p_center:-1]
					;print, n_elements(sid)
					sidarr = make_array(half_p_crossing)+total(sid)/half_p_crossing
					nh[j*60+half_p_crossing+p_center:(j+1)*60-1, i] = sidarr
					;print, n_elements(sidarr)
					
					;print, nh[j*60:(j+1)*60-1, i]
					
					;print, ''
				endif
			endfor
			col = nh[*,i]
			;if total(col) GT 20 then print, col[35*60.:37*60.]
			;print, ''
		endif
	endfor


   
   half_n_crossing = round(nwidth/2.)
   
	;For every slice in n (row in plot)
	for i=0, sh[2]-1 do begin
		;If there is at least one count in the row:
		if total(nh[i,*]) NE 0. then begin
			row = nh[i,*]
			;print, size(row)
			;print, n_elements(row)
			;for every strip in the row
			for j=0, n_elements(row)/60.-2 do begin
				;define one-strip array
				st = row[j*60:(j+1)*60-1]
				;print, size(st)
				;if there are any counts in the array, continue
				if total(st) NE 0. then begin
					;find the middle region of the array
					mid = st[half_n_crossing+1:half_n_crossing+n_center]
					;print, n_elements(mid)
					;Make a new array the same size containing the average in every position
					midarr = make_array(n_center)+total(mid)/n_center
					;print, midarr
					;Make the corresponding strip center in the full histogram that averaged array.
					nh[i, j*60+half_n_crossing:j*60+half_n_crossing+n_center-1] = midarr
				
					;find the first side region of the strip array, and repeat
					sid = st[0:half_n_crossing-1]
					;print, n_elements(sid)
					sidarr = make_array(half_n_crossing)+total(sid)/half_n_crossing
					nh[i, j*60:j*60+half_n_crossing-1] = sidarr
					
					;find the second side region of the strip array, and repeat
					sid = st[half_n_crossing+n_center:-1]
					;print, n_elements(sid)
					sidarr = make_array(half_n_crossing)+total(sid)/half_n_crossing
					nh[i, j*60+half_n_crossing+n_center:(j+1)*60-1] = sidarr
				
					;print, ''
				
				endif
			endfor
		endif
	endfor

	RETURN, nh
	

END


FUNCTION strip_based_image, file=file, erange=erange

	;GOAL: 
	;For each event in the file, extract sub-strip position based on whether it is a single or
	;double strip event
	;Plot 2D histogram image
	
	;We are currently assuming the beam is in ASICs 1 + 2
	asp = 2
	asn = 1

	restore, file
	data = eventwise_spectra.EVENTDATA_SPEC


	coordinates = []
	for i=0, n_elements(data[0,0,*])-1 do begin
			frame = data[*,*,i]
		
			;Get out quick if we are missing events on either side
			nside = where(frame[asn,*] NE 0., count)
			if count EQ 0 then continue
			if count GE 3 then continue
			pside = where(frame[asp,*] NE 0., count)
			if count EQ 0 then continue
			if count GE 3 then continue
			
			pt=TOTAL(frame[asp,*])
			nt=TOTAL(frame[asn,*]) 
			if nt gt erange[0] and nt lt erange[1] and pt gt erange[0] and pt lt erange[1] then begin
			

				;For p-side, n-side, find location of maximum energy hit, translate it to coordinate
				;In both ASIC, STRIP notation, as well as total strips notation (ie which of the 128
				;strips on the detector side).		
				pmax = where(frame[asp,*] EQ max(frame[asp,*]))
				pasic = asp
				pstrip = pmax
				if n_elements(pmax) GT 1 then pstrip = pmax[0]
				pcoord = (asp-2)*64 + pmax
				if n_elements(pmax) GT 1 then pcoord = (asp-2)*64 + pmax[0]
		
				nmax = where(frame[asn,*] EQ max(frame[asn,*]))
				nasic = asn
				nstrip = nmax
				if n_elements(nmax) GT 1 then nstrip = nmax[0]
				ncoord = asn*64 + nmax
				if n_elements(nmax) GT 1 then ncoord = asn*64 + nmax[0]
		
				coordinates = [[coordinates], [pcoord, ncoord]]
			endif
	endfor
	
	print, 'Coordinates size:', size(coordinates)
	;print, coordinates[*,0], coordinates[*,1]
	
	coordinates=transpose(coordinates)
	
	
	min1=0
    min2=0
    max1=128
    max2=128
    bin1=1
    bin2=1

    h=hist2d(coordinates[*,0], coordinates[*,1], bin=[bin1,bin2],min=[min1,min2], max=[max1,max2])
    
    print, "Total of histogram:", total(h)
    
    ;stop
	
	RETURN, float(h)
END

;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------------------------
	



PRO energy_ratio_imaging_new, files=files, filename=filename, xsize=xsize, ysize=ysize, fps=fps, $
								format=format, erange=erange, video=video, vidtype=vidtype, $
								nicefigure=nicefigure, branch=branch


;Format copied from vbins_histogram_imaging_new (which does strip-based and zone-based
;imaging). Goal here is to implement energy-ratio-based imaging as well. 

;Example call:
;foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
;energy_ratio_imaging_new, files=event_files_p, erange=[17.,24.], video=0, nicefigure=1


;For troubleshooting while updating imaging method
;newhist0 = zone_based_image(file=files[118], erange=erange)
;newhist0 = zone_based_image(file=files[180], erange=erange)
;nh0=newhist0

;newhist1 = zone_based_image(file=files[49], erange=erange)
;nh1=newhist1

;plot_image, nh0, xrange=[2400,2700], yrange=[1400,1600]
;plot_image, nh1, xrange=[2000,2300], yrange=[1500,1700]


									
loadct, 8
reverse_ct

xr=[35,42]
yr=[24,27]

default, filename, 'movie'
default, xsize, 1200
default, ysize, 400
default, fps, 10
default, format, 'mp4'
default, erange, [4.,30.]
default, video, 0
default, nicefigure, 1
default, branch, 'Pt'
default, vidtype, 'ratio'

filename = filename + '.' + format


if nicefigure EQ 1 then begin


		MAXES = []

		;newhist0 = zone_based_image(file=files[118], erange=erange)
		newhist0 = zone_based_image(file=files[180], erange=erange)
		nh0=newhist0
		;nh0 = BytScl(newhist0, TOP=255, MIN=0, MAX=256)
		
		MAXES = [MAXES, MAX(nh0)]
		
		newhist1 = zone_based_image(file=files[49], erange=erange)
		nh1=newhist1
		;nh1 = BytScl(newhist1, TOP=255, MIN=0, MAX=256)
		
		MAXES = [MAXES, MAX(nh1)]
		
		;newhist2 = strip_based_image(file=files[118], erange=erange)
		newhist2 = strip_based_image(file=files[180], erange=erange)
		nh2=newhist2
		;nh2 = BytScl(newhist2, TOP=255, MIN=0, MAX=256)
		
		MAXES = [MAXES, MAX(nh2/40.)]
		
		newhist3=strip_based_image(file=files[49], erange=erange)
		nh3=newhist3
		;nh3 = BytScl(newhist3, TOP=255, MIN=0, MAX=256)
		
		MAXES = [MAXES, MAX(nh3/40.)]
		
		bigmax = max(MAXES)
		
		print, 'bigmax:', bigmax
		
		;her1 = ratio_based_image(file=files[118], erange=erange, macrobins=1)
		her1 = ratio_based_image(file=files[180], erange=erange, macrobins=1)
		
		her2 = ratio_based_image(file=files[49], erange=erange, macrobins=1)
		;-----------------------------------------------------------------------------------------------------------
		;-----------------------------------------------------------------------------------------------------------
		
		cgPS_Open, 'multiannotate_plot.ps', /nomatch
	
		;cgLoadct, 8, NColors=bigmax, /reverse
		cgLoadct, 8, /reverse

		positions = cgLayout([2,3], OYMargin=[0,0], YGap=6, XGap=0, OXMargin=[0,0], aspect=0.5);, unit=0.5)
		;couldn't get the layout the way I wanted, so adding a correction.
		shift_pos = 0.175
		
		;-----------------------------------------------------------------------------------------------------------

		;ranges for larger, 5x4 array of intersections
;		x1=370
;		x2=419
;		y1=96
;		y2=111
		x1=380+10*3
		x2=409+10*3
		y1=96
		y2=107
		
		x1=2280+60*3
		x2=2459+60*3
		y1=1440
		y2=1619
		
		p = positions[*,2]
		p[0] = p[0]+shift_pos
		nh0 = nh0[x1:x2, y1:y2]
		print, size(nh0)
		axis_format = {XTickname:[indgen(2+(x2-x1)/10,start=x1/10)], YTickname:[indgen(2+(y2-y1)/4,start=y1/4)], $
						XTickV:[indgen(2+(x2-x1)/10)*10], YTickV:[indgen(2+(y2-y1)/4)*4], $
						XTicks:3, YTicks:3}
						
		axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
						XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
						YTicks:(y2-y1)/60}
						
	
	
		cgImage, nh0, position=p, /Axes, $
					ytitle='Al-Side', xtitle='Pt-Side', charsize=0.6, stretch=4, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
					AXKEYWORDS=axis_format, title='Sub-Strip Image, Beam at Center'
					
					
		num=0
		for i=0, 4 do begin
			cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
			num+=60
		endfor
		
		num=0
		for i=0, 4 do begin
			cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
			num+=60
		endfor
		
		;-----------------------------------------------------------------------------------------------------------
	
	    ;ranges for larger, 5x4 array of intersections
;		x1=340
;		x2=389
;		y1=96
;		y2=111
		x1=350
		x2=379
		y1=96
		y2=107
		
		x1=2100
		x2=2279
		y1=1440
		y2=1619
		
		p = positions[*,3]
		p[2] = p[2]-shift_pos
		nh1 = nh1[x1:x2, y1:y2]
		print, size(nh1)
		axis_format = {XTickname:[indgen(2+(x2-x1)/10,start=x1/10)], YTickname:[indgen(2+(y2-y1)/4,start=y1/4)], $
						XTickV:[indgen(2+(x2-x1)/10)*10], YTickV:[indgen(2+(y2-y1)/4)*4], $
						XTicks:3, YTicks:3}
						
		axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
						XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
						YTicks:(y2-y1)/60}
						
		
						
		cgImage, nh1, position=p, /Axes, /NoErase, charsize=0.6, stretch=4, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
					AXKEYWORDS=axis_format, ytitle='Al-Side', xtitle='Pt-Side', title='Sub-Strip Image, Beam at Boundaries'
				
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
		
		;-----------------------------------------------------------------------------------------------------------
		
		cgLoadct, 0, /reverse
		
		;-----------------------------------------------------------------------------------------------------------
		
	
;		x1=37
;		x2=41
;		y1=88
;		y2=91
		x1=38+3
		x2=40+3
		y1=88
		y2=90
		p = positions[*,0]
		p[0] = p[0]+shift_pos
		nh2 = nh2[x1:x2, y1:y2]
		print, size(nh2)
		axis_format = {XTickname:[indgen((x2-x1)+2,start=x1)], YTickname:[indgen((y2-y1)+2,start=y1-64)], $
						XTickV:[0,1,2,3], YTickV:[0,1,2,3], XTicks:3, YTicks:3}
	
		cgImage, nh2/40., position=p, /Axes, /NoErase, charsize=0.6, stretch=4, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
					AXKEYWORDS=axis_format, ytitle='Al-Side', xtitle='Pt-Side', title='Strip Image, Beam at Center'
					
		num=1
		for i=0, (y2-y1)+2 do begin
			;print, num
			cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
			num+=1
		endfor
		
		num=1
		for i=0, (y2-y1)+2 do begin
			;print, num
			cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
			num+=1
		endfor
		
		;-----------------------------------------------------------------------------------------------------------
					
	
;		x1=34
;		x2=38
;		y1=88
;		y2=91
		x1=35
		x2=37
		y1=88
		y2=90
		p = positions[*,1]
		p[2] = p[2]-shift_pos
		nh3_2 = nh3[x1:x2, y1:y2]
		print, size(nh3)
		axis_format = {XTickname:[indgen((x2-x1)+2,start=x1)], YTickname:[indgen((y2-y1)+2,start=y1-64)], $
						XTickV:[0,1,2,3], YTickV:[0,1,2,3], XTicks:3, YTicks:3}
	
		
		cgImage, nh3_2/40., position=p, /Axes, /NoErase, charsize=0.6, stretch=4, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
					AXKEYWORDS=axis_format, ytitle='Al-Side', xtitle='Pt-Side', title='Strip Image, Beam at Boundaries'
	
		num=1
		for i=0, (y2-y1)+2 do begin
			;print, num
			cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
			num+=1
		endfor
		
		num=1
		for i=0, (x2-x1)+2 do begin
			;print, num
			cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
			num+=1
		endfor
		
		;-----------------------------------------------------------------------------------------------------------
		
		;cgLoadct, 1, /reverse
		cgLoadct, 3, /reverse
		
		;-----------------------------------------------------------------------------------------------------------
		
		;-----------------------------------------------------------------------------------------------------------
					
	
		x1=2280+60*3
		x2=2459+60*3
		y1=1440
		y2=1619
		p = positions[*,4]
		p[0] = p[0]+shift_pos
		print, 'Before Trim Using: ', x1, x2, y1, y2
		print, total(her1)
		her1 = her1[x1:x2, y1:y2]
		print, 'After Trim', total(her1)
		axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
						XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
						YTicks:(y2-y1)/60}
		
		cgImage, her1*1000., position=p, /Axes, /NoErase, charsize=0.6, $
					AXKEYWORDS=axis_format, ytitle='Al-Side', xtitle='Pt-Side', title='Energy Ratio Image, Beam at Center'
	
		num=0
		for i=0, 4 do begin
			cgplot,[0,(x2-x1)*2],[num, num],thick=2,line=5,color=255, position=p, /overplot
			num+=60
		endfor
		
		num=0
		for i=0, 4 do begin
			cgplot,[num, num],[0, (y2-y1)*2],thick=2,line=5,color=255, position=p, /overplot
			num+=60
		endfor

		;-----------------------------------------------------------------------------------------------------------
	
		;-----------------------------------------------------------------------------------------------------------
					
	
		x1=2100
		x2=2279
		y1=1440
		y2=1619
		p = positions[*,5]
		p[2] = p[2]-shift_pos
		
		print, total(her2)		
		her2 = her2[x1:x2, y1:y2]
		print, total(her2)
		
		
		axis_format = {XTickname:[indgen((x2-x1)/60+2,start=x1/60)], YTickname:[indgen((y2-y1)/60+2,start=y1/60)], $
						XTickV:[indgen((x2-x1)/60+1)*60], YTickV:[indgen((y2-y1)/60+1)*60], XTicks:(x2-x1)/60, $
						YTicks:(y2-y1)/60}
						
		
	
		
		cgImage, her2*1000, position=p, /Axes, /NoErase, charsize=0.6, AXKEYWORDS=axis_format, $;, stretch=4;, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
					ytitle='Al-Side', xtitle='Pt-Side', title='Energy Ratio Image, Beam at Boundaries'
	
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
		
		;-----------------------------------------------------------------------------------------------------------
	
		
	
	
		;cGcolorbar, ncolors=bigmax, range=[0, bigmax], color=240, charsize=1, /XLOG
		cgPS_Close
	
		spawn, 'open multiannotate_plot.ps'

endif


if video EQ 1 and branch EQ 'Pt' then begin

	xr=[32,44]
	yr=[24,28]

	frames = n_elements(files)
	speed=1
	
	;-----------------------------------------------------------------------------------------------------------
	
	if vidtype EQ 'zone' then begin

		oVid=IDLffVideoWrite('Pt_zones_'+filename)
		vidStream = oVid.AddVideoStream(xsize, ysize, fps)

			window, xsize=xsize, ysize=ysize

		for i=0, n_elements(files)-1 do begin

			newhist3 = zone_based_image(file=files[i], erange=erange)

			plot_image, newhist3, xr=xr*60,yr=yr*60, origin=[0.5,0.5],$
			XTicklen=1.0, YTicklen=1.0,$
			XGridStyle=1, YGridStyle=1, ytickinterval=60, yminor=4, xtickinterval=60,POSITION=[0.1,0.15,0.94,0.94], $
			xtickname=indgen(13,start=xr[0]), $
			ytickname=indgen(9,start=yr[0]), ytitle='Al-Side', xtitle='Pt-Side'
			;colorbar

			time=oVid.Put(vidStream, tvrd(true=1) )

		endfor

		oVid = 0
		
	endif
			
	;-----------------------------------------------------------------------------------------------------------
;	

	if vidtype EQ 'strip' then begin
	
		loadct, 0
		reverse_ct


		oVid=IDLffVideoWrite('Pt_intersection_'+filename)
		vidStream = oVid.AddVideoStream(xsize, ysize, fps)

			window, xsize=xsize, ysize=ysize

		for i=0, n_elements(files)-1 do begin

			strip_image = strip_based_image(file=files[i], erange=erange)

			plot_image, strip_image, POSITION=[0.1,0.15,0.94,0.94], XTicklen=1.0, YTicklen=1.0, ytitle='Al-Side', xtitle='Pt-Side', $
			xr=xr,yr=yr+64, XGridStyle=1, YGridStyle=1, ytickinterval=1, xtickinterval=1, origin=[0.5,0.5], $
			ytickname=indgen(5,start=yr[0])
			;colorbar

			time=oVid.Put(vidStream, tvrd(true=1) )
			;write_png, '~/Desktop/images/image_'+strtrim(fix(i), 2)+'_'+strtrim(fix(erange[0]), 2)+'_to_'+strtrim(fix(erange[1]), 2)+'.png', tvrd(true=1)

		endfor

		oVid = 0
		
	endif

			
	;-----------------------------------------------------------------------------------------------------------

	if vidtype EQ 'ratio' then begin
	
		loadct, 3
		reverse_ct
	
		oVid=IDLffVideoWrite('Pt_ratio_'+filename)
		vidStream = oVid.AddVideoStream(xsize, ysize, fps)

			window, xsize=xsize, ysize=ysize

		for i=0, n_elements(files)-1 do begin
		
			ratio_image = ratio_based_image(file=files[i], erange=erange, macrobins=1)
		
			xr=[32,44]
			yr=[24,28]
	
			plot_image, ratio_image, POSITION=[0.1,0.15,0.94,0.94], XTicklen=1.0, YTicklen=1.0, $
			ytitle='Al-Side', xtitle='Pt-Side', origin=[0.5,0.5], xr=xr*60,yr=yr*60, XGridStyle=1, YGridStyle=1, $
			ytickinterval=60, xtickinterval=60, xtickname=[32,33,34,35,36,37,38,39,40,41,42,43,44],$
			ytickname=indgen(5,start=yr[0])
			;colorbar
		
			time=oVid.Put(vidStream, tvrd(true=1) )
			;write_png, '~/Desktop/images/image_'+strtrim(fix(i), 2)+'_'+strtrim(fix(erange[0]), 2)+'_to_'+strtrim(fix(erange[1]), 2)+'.png', tvrd(true=1)
	
	
		endfor

		oVid = 0	
	endif
	
	;-----------------------------------------------------------------------------------------------------------
	
	


endif

if video EQ 1 and branch EQ 'Al' then begin

	xr=[41,45]
	yr=[24,33]
	
	xsize=400
	ysize=900

	frames = n_elements(files)
	speed=1
	
	;-----------------------------------------------------------------------------------------------------------
	
	if vidtype EQ 'zone' then begin
	
		oVid=IDLffVideoWrite('Al_zones_'+filename)
		vidStream = oVid.AddVideoStream(xsize, ysize, fps)


			window, xsize=xsize, ysize=ysize

		for i=0, n_elements(files)-50 do begin

			newhist3 = zone_based_image(file=files[i], erange=erange)

			plot_image, newhist3, xr=xr*60, yr=yr*60, origin=[0.5,0.5],$
			XTicklen=1.0, YTicklen=1.0,$
			XGridStyle=1, YGridStyle=1, ytickinterval=60, yminor=4, xtickinterval=60,POSITION=[0.1,0.15,0.94,0.94], $
			xtickname=indgen(5,start=xr[0]), $
			ytickname=indgen(10,start=yr[0]), ytitle='Al-Side', xtitle='Pt-Side'
			;colorbar

			time=oVid.Put(vidStream, tvrd(true=1) )

		endfor

		oVid = 0
		
	endif
	
			
	;-----------------------------------------------------------------------------------------------------------

	if vidtype EQ 'strip' then begin
		loadct, 0
		reverse_ct


		oVid=IDLffVideoWrite('Al_intersection_'+filename)
		vidStream = oVid.AddVideoStream(xsize, ysize, fps)

			window, xsize=xsize, ysize=ysize

		for i=0, n_elements(files)-50 do begin

			strip_image = strip_based_image(file=files[i], erange=erange)

			plot_image, strip_image, POSITION=[0.1,0.15,0.94,0.94], XTicklen=1.0, YTicklen=1.0, ytitle='Al-Side', xtitle='Pt-Side', $
			xr=xr,yr=yr+64, XGridStyle=1, YGridStyle=1, ytickinterval=1, xtickinterval=1, origin=[0.5,0.5], $
			ytickname=indgen(10,start=yr[0])
			;colorbar

			time=oVid.Put(vidStream, tvrd(true=1) )
			;write_png, '~/Desktop/images/image_'+strtrim(fix(i), 2)+'_'+strtrim(fix(erange[0]), 2)+'_to_'+strtrim(fix(erange[1]), 2)+'.png', tvrd(true=1)

		endfor

		oVid = 0
	endif

	;-----------------------------------------------------------------------------------------------------------

	if vidtype EQ 'ratio' then begin
	
		loadct, 3
		reverse_ct
	
	
		oVid=IDLffVideoWrite('Al_ratio_'+filename)
		vidStream = oVid.AddVideoStream(xsize, ysize, fps)

			window, xsize=xsize, ysize=ysize

		for i=0, n_elements(files)-50 do begin
		
			ratio_image = ratio_based_image(file=files[i], erange=erange, macrobins=1)
		
			plot_image, ratio_image, POSITION=[0.1,0.15,0.94,0.94], XTicklen=1.0, YTicklen=1.0, $
			ytitle='Al-Side', xtitle='Pt-Side', origin=[0.5,0.5], xr=xr*60,yr=yr*60, XGridStyle=1, YGridStyle=1, $
			ytickinterval=60, xtickinterval=60, xtickname=[32,33,34,35,36,37,38,39,40,41,42,43,44],$
			ytickname=indgen(10,start=yr[0])
			;colorbar
		
			time=oVid.Put(vidStream, tvrd(true=1) )
			;write_png, '~/Desktop/images/image_'+strtrim(fix(i), 2)+'_'+strtrim(fix(erange[0]), 2)+'_to_'+strtrim(fix(erange[1]), 2)+'.png', tvrd(true=1)
	
		endfor

		oVid = 0
	endif	
	
	;-----------------------------------------------------------------------------------------------------------
	
	


endif

;stop


end
