PRO simple_event_type_ratios, files=files, side=side, energy_range=energy_range, pickasic=pickasic, posfile=posfile

;GOAL: Edited version of event_type_ratios, with all fitting and boundary lines removed (helpful
;as a first pass when looking at new scan files or scans that may not have nice peaks)

default, side, 'pside'
default, energy_range, [0., 30.]
;nonsense ASIC, for if we don't want to pick a specific one on the side we are looking at
default, pickasic, 5
default, posfile, 'sample_pos_20190419_overnight_asic2.txt'


;if side EQ 'nside' then posfile='sample_pos_20190419_overnight_asic1.txt'

;read in postions from separate file:
read_list_txt,file=posfile,tab=tab
if side EQ 'pside' then position=tab[*,1]
if side EQ 'nside' then position=tab[*,0]

print, position

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
				pside_data = frame[2:3,*]
				pside = where(pside_data NE 0.)
				tot = total(pside_data[pside])
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
				nside_data = frame[0:1,*]
				nside = where(nside_data NE 0.)
				tot = total(nside_data[nside])
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

total_numbers = double_numbers+single_numbers+other_numbers+zero_numbers

;Normalizing all the event type arrays by total events/file
dub_rat = double_numbers/float(total_numbers)
sing_rat = single_numbers/float(total_numbers)
ot_rat = other_numbers/float(total_numbers)
zero_rat = zero_numbers/float(total_numbers)


;If we are doing this for the second half of the overnight L-scan (with a beam outage and crossing some
;disabled strips), cut out those potions of the data so they don't mess things up later on. 
endval=100
if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
	single_counts = [single_counts[0:37], make_array(21, val=0.), single_counts[39:endval], make_array(161-endval, val=0.)]
	dub_rat= [dub_rat[0:37], make_array(21, val=0.), dub_rat[39:endval], make_array(161-endval, val=0.)]
	sing_rat= [sing_rat[0:37], make_array(21, val=1.), sing_rat[39:endval], make_array(161-endval, val=1.)]
	ot_rat= [ot_rat[0:37], make_array(21, val=0.), ot_rat[39:endval], make_array(161-endval, val=0.)]
	zero_rat= [zero_rat[0:37], make_array(21, val=0.), zero_rat[39:endval], make_array(161-endval, val=0.)]
endif

;=========================================================================================
;=========================================================================================

lables = ['Single Strip', 'Double Strip', 'Other Types', 'No Hits']
colors = [2, 10, 7, 6]


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

;Incredibly annoying process for printing the character "µ" on an IDL plot, see:
;http://www.idlcoyote.com/ps_tips/greeksym.php
thisLetter = "155B
greekLetter = '!9' + String(thisLetter) + '!X'

;plot single strip counts in hit strip as a function of beam position
plot, position, single_counts, thick=th, xtitle='Arbitrary Beam Position ('+greekletter+'m)', $
		title='Single Strip Counts - All Energies', ytitle='Counts', charsi=charsize
oplot, position, single_counts, thick=th, psym=psym, symsize=symsize



;Plot double strip event percentage as a function of beam position
plot, position, dub_rat*100, thick=th, ytitle='% of Counts', charsi=charsize, $
		title='Percentage of Double Strip Events', yrange=[-10, 110], xtitle='Arbitrary Beam Position ('+greekletter+'m)'
oplot, position, dub_rat*100, thick=th, color=colors[1]
oplot, [-1000,1000], [0,0], linestyle=2, thick=th
oplot, [-1000,1000], [100,100], linestyle=2, thick=th


;Plot single, double, other strip event percentages as a function of beam position
plot, position, sing_rat*100, thick=th, yrange=[-10, 110], charsi=charsize, $
		title='Event Type Ratios - '+side+' - '+'Energy Range: '+strtrim(energy_range[0], 2)+' to '+strtrim(energy_range[1], 2)+' keV', $
		 ytitle='% of Events', ystyle=1, xtitle='Arbitrary Beam Position ('+greekletter+'m)'

;If we are doing this for the second half of the overnight L-scan (with a beam outage and crossing some
;disabled strips), highlight the no data regions in a calming blue. 
if posfile EQ 'sample_pos_20190419_overnight_asic1.txt' then begin
	loadct,1
	hicol = 240
	oplot, [position[150], position[150]], [-10,1000], $
		  linestyle=0, thick=570, $
		  color=hicol
	oplot, [position[48], position[48]], [-10,1000], $
		  linestyle=0, thick=200, $
		  color=hicol
	linecolors	
endif

oplot, position, sing_rat*100, thick=th, color=colors[0]
oplot, position, sing_rat*100, thick=th, color=colors[0], psym=psym, symsize=symsize
oplot, position, dub_rat*100, thick=th, color=colors[1]
oplot, position, dub_rat*100, thick=th, color=colors[1], psym=psym, symsize=symsize
oplot, position, ot_rat*100, thick=th, color=colors[2]
oplot, position, zero_rat*100, thick=th, color=colors[3]
oplot, position, ot_rat*100, thick=th, color=colors[2], psym=psym, symsize=symsize
oplot, position, zero_rat*100, thick=th, color=colors[3], psym=psym, symsize=symsize
oplot, [-1000,1000], [0,0], linestyle=2, thick=th
oplot, [-1000,1000], [100,100], linestyle=2, thick=th

al_legend, lables, textcol=colors, box=1, /left, /top, /clear, charsi=0.8


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



