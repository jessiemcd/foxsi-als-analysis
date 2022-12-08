function lookatlivetime

	;foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
	;foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=2

	foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=struct_files_p
	foxsi_cal_struct_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=struct_files_n

	files=struct_files_p
	p_livetimes=[]

	for j=0, n_elements(files)-1 do begin
		restore, files[j]
	
		livetime=0
		for i=0, n_elements(data[*].time[0])-1 do begin
			if min(data[i].time[*]) GE 0 then livetime+=min(data[i].time[*])
		endfor
	
		p_livetimes = [p_livetimes, livetime]
	
	endfor

	files=struct_files_n
	n_livetimes=[]

	for j=0, n_elements(files)-1 do begin
		restore, files[j]
	
		livetime=0
		for i=0, n_elements(data[*].time[0])-1 do begin
			if min(data[i].time[*]) GE 0 then livetime+=min(data[i].time[*])
		endfor
	
		n_livetimes = [n_livetimes, livetime]
	
	endfor

	plot, p_livetimes/1.0e6, yrange=[-10,10]
	oplot, n_livetimes/1.0e6, linestyle=2
	
	endval=100
	n_livetimes = [n_livetimes[0:37], make_array(21, val=1.), n_livetimes[39:endval], make_array(161-endval-1, val=1.)]
	
	print, size([[p_livetimes], [n_livetimes]])

	RETURN, [[p_livetimes], [n_livetimes]]
	
END

;==============================================================================================================================
;==============================================================================================================================
;==============================================================================================================================
;==============================================================================================================================
;==============================================================================================================================
;==============================================================================================================================
;==============================================================================================================================


PRO event_type_ratios_figure, files=files, side=side, energy_range=energy_range, $
								pickasic=pickasic, posfile=posfile, livetime_corr=livetime_corr

;GOAL: Find # of single, double, and triple+ events for each file/beam position
;-Plot ratios as a function of beam position
;-Do some sort of fitting to assign strip boundaries
;-Find widths of double-strip event regions

;for call:
;foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
;foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=2


;Additionally, plot single strip events in all scan-touched strips

default, side, 'pside'
default, energy_range, [0., 30.]
;nonsense ASIC, for if we don't want to pick a specific one on the side we are looking at
default, pickasic, 5
default, posfile, 'sample_pos_20190419_overnight_asic2.txt'
default, livetime_corr, 0

if livetime_corr EQ 1 then ff = lookatlivetime()

if side EQ 'nside' then posfile='sample_pos_20190419_overnight_asic1.txt'

;read in postions from separate file:
read_list_txt,file=posfile,tab=tab
if side EQ 'pside' then position=tab[*,1]
if side EQ 'nside' then position=tab[*,0]

;Converting positions to units of µm
position = position-float(min(position))
position *= 1000


if n_elements(files) EQ 0 then begin
	print, ''
	print, 'Hey, there are no files.'
	print, ''
endif

;These will be filled with the # of that type of event in each file (in selected energy range)
double_numbers = []
single_numbers = []
other_numbers = []
zero_numbers = []

hit_strips = []
single_counts = []

print, 'Number of Files:', n_elements(files)

for j=0, n_elements(files)-1 do begin
	restore, files[j]

	data = eventwise_spectra.EVENTDATA_SPEC

	doubles=0
	singles=0
	others=0
	zeroes=0
	
	single_hits = []
	
	if side EQ 'pside' then begin
		for i=0, n_elements(data[0,0,*])-1 do begin
			frame = data[*,*,i]
			
			;Take nonzero hits on pside (or, in specific selected ASIC)
			if pickasic EQ 5 then begin
				side_data = frame[2:3,*]
				pside = where(side_data NE 0.)
				tot = total(side_data[pside])
			endif
			if pickasic EQ 2 then begin
				data2 = frame[2,*]
				pside = where(data2 NE 0.)
				tot = total(data2[pside])
			endif
			if pickasic EQ 3 then begin
				data3 = frame[3,*]
				pside = where(data3 NE 0.)
				tot = total(data3[pside])
			endif
			
			totn=total(frame[1,*])
			if abs(tot-totn) GT 3 then begin
				;print, tot, totn
				continue
			endif
			
			;If there is ONLY ONE hit in the AISC or side chosen (single-strip event)
			;DOING BEFORE ENERGY SELECTION
			if n_elements(pside) EQ 1 and total(pside) NE -1. then begin
				single_hits = [single_hits, pside]
			endif
			
			
			;Energy selection condition. Note we are selecting in TOTAL ENERGY (SUM if 
			;more than one event in the side/asic selected)
			if tot LT energy_range[0] or tot GT energy_range[1] then continue
			
			;If there are NO hits in the ASIC or side chosen
			if total(pside) EQ -1 then zeroes+=1
			
			;If there is ONLY ONE hit in the AISC or side chosen (single-strip event)
			if n_elements(pside) EQ 1 and total(pside) NE -1. then begin
				singles+=1
				;single_hits = [single_hits, pside]
				continue
			endif
			
			;If there are TWO hits in the ASIC or side chosen - need to determine if they are
			;adjacent or nah
			if n_elements(pside) EQ 2 then begin
				if pickasic EQ 5 then begin
					as2 = where(frame[2,*] NE 0.)
					as3 = where(frame[3,*] NE 0.)
					if n_elements(as3) NE 2 and n_elements(as2) NE 2 then begin
						others+=1
						continue
					endif
					if n_elements(as2) EQ 2 then begin
						if as2[1]-as2[0] NE 1 then others+=1
						if as2[1]-as2[0] EQ 1 then doubles+=1
						continue
					endif
					if n_elements(as3) EQ 2 then begin
						if as3[1]-as3[0] NE 1 then others+=1
						if as3[1]-as3[0] EQ 1 then doubles+=1
					endif
				endif
				
				if pickasic EQ 2 then begin
					as2 = where(frame[2,*] NE 0.)
					if n_elements(as2) NE 2 then begin
						others+=1 
						continue
					endif
					if n_elements(as2) EQ 2 then begin
						if as2[1]-as2[0] NE 1 then others+=1
						if as2[1]-as2[0] EQ 1 then doubles+=1
					endif
				endif
				
				if pickasic EQ 3 then begin
					as3 = where(frame[3,*] NE 0.)
					if n_elements(as3) NE 2 then begin
						others+=1
						continue
					endif
					if n_elements(as3) EQ 2 then begin
						if as3[1]-as3[0] NE 1 then others+=1
						if as3[1]-as3[0] EQ 1 then doubles+=1
					endif
				endif
			endif
			if n_elements(pside) GT 2 then others+=1
		endfor
	endif
	
	if side EQ 'nside' then begin
		for i=0, n_elements(data[0,0,*])-1 do begin
			frame = data[*,*,i]
			
			;Take nonzero hits on nside (or, in specific selected ASIC)
			if pickasic EQ 5 then begin
				side_data = frame[0:1,*]
				nside = where(side_data NE 0.)
				tot = total(side_data[nside])
			endif
			if pickasic EQ 0 then begin
				data0 = frame[0,*]
				nside = where(data0 NE 0.)
				tot = total(data0[nside])
			endif
			if pickasic EQ 1 then begin
				data1 = frame[1,*]
				nside = where(data1 NE 0.)
				tot = total(data1[nside])
			endif
			
			totp=total(frame[2,*])
			if abs(tot-totp) GT 3 then begin
				;print, tot, totp
				continue
			endif
			
			;If there is ONLY ONE hit in the AISC or side chosen (single-strip event)
			;DOING BEFORE ENERGY SELECTION
			if n_elements(nside) EQ 1 and total(nside) NE -1. then begin
				single_hits = [single_hits, nside]
			endif
			
			;Energy selection condition. Note we are selecting in TOTAL ENERGY (SUM if 
			;more than one event in the side/asic selected)
			if tot LT energy_range[0] or tot GT energy_range[1] then continue
			
			;If there are NO hits in the ASIC or side chosen
			if total(nside) EQ -1 then zeroes+=1
			
			;If there is ONLY ONE hit in the AISC or side chosen (single-strip event)
			if n_elements(nside) EQ 1 and total(nside) NE -1. then begin
				singles+=1
				;single_hits = [single_hits, nside]
				continue
			endif
			
			
			;If there are TWO hits in the ASIC or side chosen - need to determine if they are
			;adjacent or nah
			if n_elements(nside) EQ 2 then begin
				if pickasic EQ 5 then begin
					as0 = where(frame[0,*] NE 0.)
					as1 = where(frame[1,*] NE 0.)
					if n_elements(as0) NE 2 and n_elements(as1) NE 2 then begin
						others+=1
						continue
					endif
					if n_elements(as0) EQ 2 then begin
						if as0[1]-as0[0] NE 1 then others+=1
						if as0[1]-as0[0] EQ 1 then doubles+=1
						continue
					endif
					if n_elements(as1) EQ 2 then begin
						if as1[1]-as1[0] NE 1 then others+=1
						if as1[1]-as1[0] EQ 1 then doubles+=1
					endif
				endif
				
				if pickasic EQ 0 then begin
					as0 = where(frame[0,*] NE 0.)
					if n_elements(as0) NE 2 then begin
						others+=1 
						continue
					endif
					if n_elements(as0) EQ 2 then begin
						if as0[1]-as0[0] NE 1 then others+=1
						if as0[1]-as0[0] EQ 1 then doubles+=1
					endif
				endif
				
				if pickasic EQ 1 then begin
					as1 = where(frame[1,*] NE 0.)
					if n_elements(as1) NE 2 then begin
						others+=1
						continue
					endif
					if n_elements(as1) EQ 2 then begin
						if as1[1]-as1[0] NE 1 then others+=1
						if as1[1]-as1[0] EQ 1 then doubles+=1
					endif
				endif
			endif
			if n_elements(nside) GT 2 then others+=1
		endfor
	endif
	
	;Defining "hit strip" for this file - strip with most single-pixel events
	hit_strip = mode_val(single_hits)
	hit_strips = [hit_strips, hit_strip]
	
	;Number of counts in that hit strip
	single_count = n_elements(where(single_hits EQ hit_strip))
	single_counts = [single_counts, single_count]
	
	
	double_numbers = [double_numbers, doubles]
	single_numbers = [single_numbers, singles]
	other_numbers = [other_numbers, others]
	zero_numbers = [zero_numbers, zeroes]
	
endfor

all_strips = hit_strips[sort(hit_strips)]
all_strips = all_strips[UNIQ(all_strips)]


if pickasic EQ 5 then begin
	strip_inds = array_indices(side_data, all_strips)
	strips = strip_inds[1,*]
	;print, strips
endif
if pickasic LE 3 then begin
	strips = all_strips
	;print, strips
endif


;The idea here is to make an array for every strip (strips) that has the number of single
;strip events in that strip at each location. Doing this AFTER energy selection. UNFINISHED 

strips_counts = fltarr(n_elements(strips), n_elements(files))
;print, size(strips_counts)

event_count=0.
add_count=0.

;Making strip-specific single strip event count arrays 
for j=0, n_elements(files)-1 do begin
	restore, files[j]
	data = eventwise_spectra.EVENTDATA_SPEC
	if side EQ 'pside' then begin
		for i=0, n_elements(data[0,0,*])-1 do begin
			event_count+=1
			frame = data[*,*,i]
			
			;Take nonzero hits on pside (or, in specific selected ASIC)
			if pickasic EQ 5 then begin
				side_data = frame[2:3,*]
				pside = where(side_data NE 0.)
				tot = total(side_data[pside])
			endif
			if pickasic EQ 2 then begin
				data2 = frame[2,*]
				pside = where(data2 NE 0.)
				tot = total(data2[pside])
			endif
			if pickasic EQ 3 then begin
				data3 = frame[3,*]
				pside = where(data3 NE 0.)
				tot = total(data3[pside])
			endif
			
			
			;print, pside
			
			
			;Energy selection condition. Note we are selecting in TOTAL ENERGY (SUM if 
			;more than one event in the side/asic selected)
			if tot LT energy_range[0] or tot GT energy_range[1] then continue
			
			;If there is ONLY ONE hit in the AISC or side chosen (single-strip event)
			if n_elements(pside) EQ 1 and total(pside) NE -1. then begin
				add_count+=1
				;print, 'Picked:', pside
				strip = pside
				if pickasic EQ 5 then strip = pside/2
				strips_counts[where(strips EQ strip[0]), j]+=1
			endif
		endfor
	endif
	
	if side EQ 'nside' then begin
		for i=0, n_elements(data[0,0,*])-1 do begin
			event_count+=1
			frame = data[*,*,i]
			
			;Take nonzero hits on pside (or, in specific selected ASIC)
			if pickasic EQ 5 then begin
				side_data = frame[0:1,*]
				nside = where(side_data NE 0.)
				tot = total(side_data[nside])
			endif
			if pickasic EQ 0 then begin
				data2 = frame[0,*]
				nside = where(data2 NE 0.)
				tot = total(data2[nside])
			endif
			if pickasic EQ 1 then begin
				data3 = frame[1,*]
				nside = where(data3 NE 0.)
				tot = total(data3[nside])
			endif
			
			
			;Energy selection condition. Note we are selecting in TOTAL ENERGY (SUM if 
			;more than one event in the side/asic selected)
			if tot LT energy_range[0] or tot GT energy_range[1] then continue
			
			;If there is ONLY ONE hit in the AISC or side chosen (single-strip event)
			if n_elements(nside) EQ 1 and total(nside) NE -1. then begin
				add_count+=1
				;print, 'Picked:', pside
				strip = nside
				if pickasic EQ 5 then strip = nside/2
				strips_counts[where(strips EQ strip[0]), j]+=1
			endif
		endfor
	endif
	
	

endfor

print, event_count, add_count

;STOP



total_numbers = double_numbers+single_numbers+other_numbers+zero_numbers

;Normalizing all the event type arrays by total events/file
dub_rat = double_numbers/float(total_numbers)
sing_rat = single_numbers/float(total_numbers)
ot_rat = other_numbers/float(total_numbers)
zero_rat = zero_numbers/float(total_numbers)

;Taking event type % over the whole scans biases based on how many edge/center regions are included 
;(Pt-side: too many centers, Al-side too many boundaries)
;Defining windows to reduce this bias:
if side EQ 'pside' then ew=[6:172]
if side EQ 'nside' then ew=[70:110]

print, 'Double Events Fraction: ', total(double_numbers[ew])/total(total_numbers[ew])
print, 'Single Events Fraction: ', total(single_numbers[ew])/total(total_numbers[ew])
print, 'Other Events Fraction: ', total(other_numbers[ew])/total(total_numbers[ew])


;If we are doing this for the second half of the overnight L-scan (with a beam outage and crossing some
;disabled strips), cut out those potions of the data so they don't mess things up later on. 
endval=100
if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
	single_counts = [single_counts[0:37], make_array(21, val=0.), single_counts[39:endval], make_array(161-endval-1, val=0.)]
	dub_rat= [dub_rat[0:37], make_array(21, val=0.), dub_rat[39:endval], make_array(161-endval-1, val=0.)]
	sing_rat= [sing_rat[0:37], make_array(21, val=1.), sing_rat[39:endval], make_array(161-endval-1, val=1.)]
	ot_rat= [ot_rat[0:37], make_array(21, val=0.), ot_rat[39:endval], make_array(161-endval-1, val=0.)]
	zero_rat= [zero_rat[0:37], make_array(21, val=0.), zero_rat[39:endval], make_array(161-endval-1, val=0.)]
endif



;FInd the maxima of the array
maxes = lclxtrem(dub_rat, 12, /maxima)

;For troubleshooting - can plot maxima locations over array to see which are being found.
;oplot, maxes, dub_rat[maxes]*100, psym=2, symsize=symsize*2, color=colors[2]

;Line for checking where in the array a certain index lies (e.g. whether a file is from the center
;or boundary of a strip)
;oplot, [position[82],position[82]],[-10,1000], color=2, thick=5 

ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
symsize=0.6
plot, position, dub_rat*100

;Defining gaussian fitting regions as 5 points in either direction around each maxima, fitting a 
;gaussian to said peak, overplotting the fit + included data points, and drawing a line at strip 
;boundary as estimated by the gaussian fit. Also, finds widths of gaussians to characterize the
;double-strip event region.
num = 10
if side EQ 'pside' then num=5
coll=1
widths=[]
bounds=[]
maxes = maxes[SORT(maxes)]
for i=0, n_elements(maxes)-1 do begin
	if coll EQ 5 then coll+=1
	if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' and i EQ 2 then continue
	peak_region = dub_rat[maxes[i]-num:maxes[i]+num]*100
	xx = position[maxes[i]-num:maxes[i]+num]
	oplot, xx, peak_region, psym=2, symsize=symsize*2, color=coll
	fit=GAUSSFIT(xx,peak_region,coeff,nterms=3)
    widths=[widths,coeff[2]]
    bounds=[bounds,coeff[1]]
    oplot, xx, fit, color=coll, thick=th
    oplot, [coeff[1], coeff[1]], [-10,1e5], thick=th/2, color=coll
	coll+=1
endfor

print, 'Widths of Double Strip Peaks: ', widths
print, 'Average Width of Double Strip Peaks: ', mean(widths)
print, 'Average Width FWHM: ', 2*SQRT(2*ALOG(2))*mean(widths)




;=========================================================================================
;=========================================================================================

;labels = ['Single Strip', 'Double Strip', 'Other Types', 'No Hits']
;colors = [2, 10, 7, 6]
labels = ['% Single-strip', '% Double-strip', '% Others']
colors = [2, 10, 7]

strip_labels=['Single-strip', 'counts in:']
for i=0, n_elements(strips)-1 do begin
	strip_labels=[strip_labels, 'Strip '+strtrim(strips[i],2)]
endfor


popen, 'event_ratios.ps', $
		xsi=8, ysi=10
!Y.margin=4.
!X.margin=4.
ch=1.1
th=4
lnth=4
fth=4
charsize=1.3
symsize=0.6
psym=4
!p.multi=[0,1,5]
linecolors
;=========================================================================================
;=========================================================================================

;Incredibly annoying process for printing the character "µ" on an IDL plot using .ps device, see:
;http://www.idlcoyote.com/ps_tips/greeksym.php
thisLetter = "155B
greekLetter = '!9' + String(thisLetter) + '!X'

;plot single strip counts in hit strip as a function of beam position
;plot, position, single_counts, thick=th, xtitle='Arbitrary Beam Position ('+greekletter+'m)', $
;		title='Single Strip Counts - '+side+' - All Energies', ytitle='Counts', charsi=charsize
;oplot, position, single_counts, thick=th, psym=psym, symsize=symsize

if posfile EQ 'sample_pos_20190419_overnight_asic2.txt' then plottitle='Pt-Side Branch of Scan - Selected Energy Range: '+strtrim(round(energy_range[0]), 2)+' to '+strtrim(round(energy_range[1]), 2)+' keV'
if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then plottitle='Al-Side Branch of Scan - Selected Energy Range: '+strtrim(round(energy_range[0]), 2)+' to '+strtrim(round(energy_range[1]), 2)+' keV'

;print, strips_counts

first_strip = strips_counts[0,*]

if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
	first_strip = [first_strip[0:37], make_array(21, val=0.), first_strip[39:endval], make_array(161-endval-1, val=0.)]
	if livetime_corr EQ 1 then first_strip = first_strip/(ff[*,1]/1.0e6)
endif 

if livetime_corr EQ 1 and posfile NE 'sample_pos_20190419_overnight_asic1.txt' then first_strip = first_strip/(ff[*,0]/1.0e6)

;instead, plotting single strip counts in ALL strips as function of beam postion.
if livetime_corr EQ 0 then begin
	plot, position, first_strip, thick=th,  xrange=[-620,0], xstyle=1,$ ;xtitle='Arbitrary Beam Position ('+greekletter+'m)', $
			title=plottitle, $
			ytitle='Counts', charsi=charsize, yrange=[0, max(strips_counts)+50], position=[0.05,0.8,0.9,0.95], xtickformat='(a1)'
endif

if livetime_corr EQ 1 then begin
	plot, position, first_strip, thick=th,  xrange=[-620,0], xstyle=1,$ ;xtitle='Arbitrary Beam Position ('+greekletter+'m)', $
			title=plottitle, $
			ytitle='Livetime-Corrected Counts/s', charsi=charsize, yrange=[0, max(first_strip)+50], position=[0.05,0.8,0.9,0.95], xtickformat='(a1)'
endif
;If we are doing this for the second half of the overnight L-scan (with a beam outage and crossing some
;disabled strips), highlight the no data regions in a calming blue. 
if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
	loadct,1
	hicol = 240
	oplot, [position[150], position[150]], [-10,1000], $
		  linestyle=0, thick=500, $
		  color=hicol
	oplot, [position[48], position[48]], [-10,1000], $
		  linestyle=0, thick=170, $
		  color=hicol
	linecolors	
endif
		
		
coll=1
strip_colors=[0,0]
for i=0, n_elements(strips)-1 do begin
	if coll EQ 5 then coll+=1
	strip_colors=[strip_colors, coll]
	strip_arr = strips_counts[i,*]
	if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
		strip_arr = [strip_arr[0:37], make_array(21, val=0.), strip_arr[39:endval], make_array(161-endval-1, val=0.)]
		if livetime_corr EQ 1 then strip_arr = strip_arr/(ff[*,1]/1.0e6)
	endif
	if livetime_corr EQ 1 and posfile NE 'sample_pos_20190419_overnight_asic1.txt' then begin
		strip_arr = strip_arr/(ff[*,0]/1.0e6)
	endif
	oplot, position, strip_arr, thick=th, color=coll
	oplot, position, strip_arr, thick=th, color=coll, psym=psym, symsize=symsize
	coll+=1
endfor


for i=0, n_elements(bounds)-1 do begin
	oplot, [bounds[i], bounds[i]], [-10, 1000], thick=th/2
endfor

if livetime_corr EQ 0 then begin
	if posfile EQ 'sample_pos_20190419_overnight_asic2.txt' then al_legend, strip_labels, textcol=strip_colors, box=1, /left, /top, /clear, charsi=0.65, pos=[-615,950]
	if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then al_legend, strip_labels[0:-1], textcol=strip_colors[0:-1], box=1, /left, /top, /clear, charsi=0.65, pos=[-615,750]
endif

if livetime_corr EQ 1 then begin
	if posfile EQ 'sample_pos_20190419_overnight_asic2.txt' then al_legend, strip_labels, textcol=strip_colors, box=1, /left, /top, /clear, charsi=0.65, pos=[-615,240]
	if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then al_legend, strip_labels[0:-1], textcol=strip_colors[0:-1], box=1, /left, /top, /clear, charsi=0.65, pos=[-615,240]
endif

;FInd the minima and maxima of the single strip counts array
mins = lclxtrem(single_counts, 6)
maxes = lclxtrem(single_counts, 7, /maxima)

;Select some of the minima (eliminating those that are close to maxima, as we don't expect
;those to actually be associated with strip boundaries)
selected_mins = []
for i=0, n_elements(mins)-1 do begin
	diffs = abs(maxes-mins[i])
	if min(diffs) GE 5 then selected_mins = [selected_mins, mins[i]]
endfor

;For troubleshooting - can plot maxima/minima locations over single counts array to see 
;which extrema are being found.
;oplot, position[mins], single_counts[mins], psym=2, symsize=symsize*2, color=colors[0]
;oplot, position[maxes], single_counts[maxes], psym=2, symsize=symsize*2, color=colors[2]
;oplot, position[selected_mins], single_counts[selected_mins], psym=2, symsize=symsize*2, color=colors[1]

;Defining gaussian fitting regions as the space between two adjacent minima, fitting a gaussian to said
;peak, overplotting the fit + included data points, and drawing a line at strip center as estimated by 
;the gaussian fit. 
coll=2
centers=[]
selected_mins = selected_mins[SORT(selected_mins)]
for i=0, n_elements(selected_mins)-2 do begin
	if coll EQ 5 then coll+=1
	if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' and i EQ 1 then continue
	peak_region = single_counts[selected_mins[i]:selected_mins[i+1]]
	xx = position[selected_mins[i]:selected_mins[i+1]]
	;oplot, xx, peak_region, psym=2, symsize=symsize*2, color=coll
	fit=GAUSSFIT(xx,peak_region,coeff,nterms=3)
	centers=[centers,coeff[1]]
    ;oplot, xx, fit, color=coll, thick=th
    oplot, [coeff[1], coeff[1]], [-10,1e5], thick=th/2, color=0, linestyle=2
    ;oplot, [selected_mins[i], selected_mins[i]], [-10, 1e5], thick=th/2, color=coll
	coll+=1
endfor



;Plot single, double, other strip event percentages as a function of beam position
plot, position, sing_rat*100, thick=th, yrange=[-10, 110], charsi=charsize, xrange=[-620,0], xstyle=1, $
		;title='Event Type Ratios - '+side+' - '+'Energy Range: '+strtrim(energy_range[0], 2)+' to '+strtrim(energy_range[1], 2)+' keV', $
		 ytitle='% of Events', ystyle=1, xtitle='Arbitrary Beam Position ('+greekletter+'m)', position=[0.05,0.65,0.9,0.8]

;If we are doing this for the second half of the overnight L-scan (with a beam outage and crossing some
;disabled strips), highlight the no data regions in a calming blue. 
if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
	loadct,1
	hicol = 240
	oplot, [position[150], position[150]], [-10,1000], $
		  linestyle=0, thick=500, $
		  color=hicol
	oplot, [position[48], position[48]], [-10,1000], $
		  linestyle=0, thick=170, $
		  color=hicol
	linecolors	
	
	;makes lines to indicate locations of beam spectra used in spectra figure
	;oplot, [position[71], position[71]], [-10,1000], linestyle=3, color=0, thick=4
	;oplot, [position[81], position[81]], [-10,1000], linestyle=3, color=0, thick=4
endif

if posfile NE 'sample_pos_20190419_overnight_asic1.txt' then begin
	;makes lines to indicate locations of beam spectra used in spectra figure
	;oplot, [position[60], position[60]], [-10,1000], linestyle=3, color=0, thick=4
	;oplot, [position[49], position[49]], [-10,1000], linestyle=3, color=0, thick=4
endif



oplot, position, sing_rat*100, thick=th, color=colors[0]
oplot, position, sing_rat*100, thick=th, color=colors[0], psym=psym, symsize=symsize
oplot, position, dub_rat*100, thick=th, color=colors[1]
oplot, position, dub_rat*100, thick=th, color=colors[1], psym=psym, symsize=symsize
oplot, position, ot_rat*100, thick=th, color=colors[2]
;oplot, position, zero_rat*100, thick=th, color=colors[3]
oplot, position, ot_rat*100, thick=th, color=colors[2], psym=psym, symsize=symsize
;oplot, position, zero_rat*100, thick=th, color=colors[3], psym=psym, symsize=symsize
oplot, [-1000,1000], [0,0], linestyle=2, thick=th
oplot, [-1000,1000], [100,100], linestyle=2, thick=th



cwidths=[]
;Overplot strip centers and boundaries as estimated by the two fitting processes above.
for i=0, n_elements(centers)-1 do begin
	oplot, [centers[i], centers[i]], [-10, 1000], thick=th/2, linestyle=2
	if i LT n_elements(centers)-1 then begin
		print, 'Center to center: ', centers[i]-centers[i+1]
		cwidths = [cwidths, abs(centers[i]-centers[i+1])]
	endif
endfor
bwidths=[]
for i=0, n_elements(bounds)-1 do begin
	oplot, [bounds[i], bounds[i]], [-10, 1000], thick=th/2
	if i LT n_elements(bounds)-1 then begin
		print, 'Bound to Bound: ', bounds[i]-bounds[i+1]
		bwidths = [bwidths, abs(bounds[i]-bounds[i+1])]
	endif
endfor

print, cwidths/60.
print, abs(cwidths/60.)-1.
if side EQ 'pside' then print, mean(abs(cwidths/60.)-1.)
if side EQ 'nside' then print, mean(abs(cwidths[1]/60.)-1.)

print, bwidths/60.
print, abs(bwidths/60.)-1.
if side EQ 'pside' then print, mean(abs(bwidths/60.)-1.)
if side EQ 'nside' then print, mean(abs([bwidths[0], bwidths[2], bwidths[3]]/60.)-1.)

al_legend, labels, textcol=colors, box=1, /left, /top, /clear, charsi=0.65, pos=[-615, 65]


;=========================================================================================
;=========================================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
;DEVICE, /CLOSE
spawn, 'open event_ratios.ps'

;=========================================================================================
;=========================================================================================





;STOP

END



