PRO sharing_energetics

;This makes a bunch of .csv files which I then import to iPython notebooks because I could not 
;bring myself to wrestle with IDL plotting during the energy ratio analysis process. 

;energy range:
ER = [17.,24.]



;Import files
;foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p
;foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n

;event_files_p=file_search('quad_6new*')
;event_files_n=file_search('nquad_6new*')

;To use two-threshold calibrated data:
foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=2

;To use one-threshold calibrated data:
;foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=3
;foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=3
;=========================================================================================
;Looking at double-pixel events on the p-side, during the p-side branch of the scan
;=========================================================================================


p_rat_hists = fltarr(n_elements(event_files_p), 101)

all_energy_pairs = []
all_rats = []

print, 'Now, P Scan, p-side events:'

for j=0, n_elements(event_files_p)-1 do begin
		restore, event_files_p[j]
		
		array_str = 'file_no_'+strtrim(j, 2)
		;print, array_str
		void = EXECUTE(array_str + '=[]')
		
		p_energy_pairs = []
		data = eventwise_spectra.EVENTDATA_SPEC
		as = 2
		for i=0, n_elements(data[as,0,*])-1 do begin
			hits = where(data[as,*,i] NE 0.)
			
			;if n_elements(hits) EQ 2 then begin
			;if n_elements(hits) EQ 2 and total(data[as,*,i]) GE ER[0] and total(data[as,*,i]) LE ER[1] then begin
			if n_elements(hits) GE 2 and total(data[as,*,i]) GE ER[0] and total(data[as,*,i]) LE ER[1] then begin
				ss = sort(data[as,*,i])
				top2 = [ss[-1], ss[-2]]
				if abs(ss[-1]-ss[-2]) EQ 1 then p_energy_pairs = [[p_energy_pairs], [data[as,min(top2),i], data[as,max(top2),i]]]
				;if hits[1] EQ hits[0]+1 then p_energy_pairs = [[p_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
				;if hits[1] EQ hits[0]+1 then all_energy_pairs = [[all_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
			endif
		endfor
		
		;print, size(p_energy_pairs)
		
		;print, 'Number of Doubles: ', n_elements(p_energy_pairs[1,*])
		
		void = EXECUTE(array_str + '= p_energy_pairs')
		
		rats = (p_energy_pairs[0,*] - p_energy_pairs[1,*]) / (p_energy_pairs[0,*] + p_energy_pairs[1,*])
		all_rats = [all_rats, reform(rats)]
		
		;rat_hist = histogram(abs(rats), min=0, max=1, bin = 0.02)
		rat_hist = histogram(rats, min=-1, max=1, bin = 0.02)
		;print, size(rat_hist)
		
		p_rat_hists[j,*] = rat_hist
endfor

print, 'Rats min: ', min(all_rats), ' Rats max: ', max(all_rats)



;p_rat_hists_norm = fltarr(n_elements(event_files_p), 40)
;
;for j=0, n_elements(event_files_p)-1 do begin
;	mt = total(p_rat_hists[j,*])
;	p_rat_hists_norm[j,*] = 100*p_rat_hists[j,*]/mt
;
;endfor

;pep = transpose(all_energy_pairs)

;h=hist2d(pep[*,0], pep[*,1], bin=[0.5, 0.5],min=[0,0], max=[30,30])

;write_csv, 'all_double_energies_pscan.csv', all_energy_pairs

;write_csv, 'position_double_energies.csv', p_rat_hists
write_csv, 'highe_position_double_energies.csv', p_rat_hists

;=========================================================================================
;Looking at double-pixel events on the p-side, during the n-side branch of the scan
;=========================================================================================

;p_rat_hists_nside = fltarr(n_elements(event_files_p), 101)
;
;all_energy_pairs = []
;
;print, 'Now, N Scan, p-side events:'
;
;for j=0, n_elements(event_files_n)-1 do begin
;		restore, event_files_n[j]
;		
;		;stop
;		
;		counter=0
;
;		p_energy_pairs = []
;		data = eventwise_spectra.EVENTDATA_SPEC
;		as = 2
;		for i=0, n_elements(data[as,0,*])-1 do begin
;			hits = where(data[as,*,i] NE 0.)
;			;if n_elements(hits) EQ 2 then begin
;			if n_elements(hits) EQ 2 and total(data[as,*,i]) GE ER[0] and total(data[as,*,i]) LE ER[1] then begin
;				counter+=1
;				if hits[1] EQ hits[0]+1 then p_energy_pairs = [[p_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
;				if hits[1] EQ hits[0]+1 then all_energy_pairs = [[all_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
;			endif
;		endfor
;		
;		;print, counter
;		;print, size(p_energy_pairs)
;		
;		rats = (p_energy_pairs[0,*] - p_energy_pairs[1,*]) / (p_energy_pairs[0,*] + p_energy_pairs[1,*])
;		
;		rat_hist = histogram(rats, min=-1, max=1, bin = 0.02)
;		;print, size(rat_hist)
;		
;		p_rat_hists_nside[j,*] = rat_hist
;endfor
;
;;write_csv, 'all_double_energies_nscan_pside.csv', all_energy_pairs
;
;;write_csv, 'position_double_energies_nscan_pside.csv', p_rat_hists_nside
;write_csv, 'highe_position_double_energies_nscan_pside.csv', p_rat_hists_nside

;=========================================================================================
;Looking at double-pixel events on the n-side, during the n-side branch of the scan
;=========================================================================================


n_rat_hists = fltarr(n_elements(event_files_n), 101)

all_n_energy_pairs = []

all_rats = []

print, 'Now, N Scan, n-side events:'

for j=0, n_elements(event_files_n)-1 do begin
		restore, event_files_n[j]
		
		n_energy_pairs = []
		data = eventwise_spectra.EVENTDATA_SPEC
		as = 1
		for i=0, n_elements(data[as,0,*])-1 do begin
			hits = where(data[as,*,i] NE 0.)
			;if n_elements(hits) EQ 2 then begin
			;if n_elements(hits) EQ 2 and total(data[as,*,i]) GE ER[0] and total(data[as,*,i]) LE ER[1] then begin
			if n_elements(hits) GE 2 and total(data[as,*,i]) GE ER[0] and total(data[as,*,i]) LE ER[1] then begin
				ss = sort(data[as,*,i])
				top2 = [ss[-1], ss[-2]]
				if abs(ss[-1]-ss[-2]) EQ 1 then n_energy_pairs = [[n_energy_pairs], [data[as,min(top2),i], data[as,max(top2),i]]]
				;if hits[1] EQ hits[0]+1 then n_energy_pairs = [[n_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
				;if hits[1] EQ hits[0]+1 then all_n_energy_pairs = [[all_n_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
			endif
		endfor
		
		;print, size(n_energy_pairs)
		
		
		if N_ENERGY_PAIRS EQ !NULL then rat_hist = fltarr(1, 101)
		if N_ENERGY_PAIRS NE !NULL then begin
			rats = (n_energy_pairs[0,*] - n_energy_pairs[1,*]) / (n_energy_pairs[0,*] + n_energy_pairs[1,*])
			all_rats = [all_rats, reform(rats)]
		
			;rat_hist = histogram(abs(rats), min=0, max=1, bin = 0.02)
			rat_hist = histogram(rats, min=-1, max=1, bin = 0.02, locations=xbin)
			;print, 'xbin', xbin
			
			;stop
			;print, size(rat_hist)
		endif
		
		n_rat_hists[j,*] = rat_hist
endfor

print, 'Rats min: ', min(all_rats), ' Rats max: ', max(all_rats)


;write_csv, 'all_double_energies_nscan.csv', all_n_energy_pairs

;write_csv, 'position_double_energies_nside.csv', n_rat_hists
write_csv, 'highe_position_double_energies_nside.csv', n_rat_hists

;=========================================================================================
;Looking at double-pixel events on the n-side, during the p-side branch of the scan
;=========================================================================================

;n_rat_hists_pside = fltarr(n_elements(event_files_p), 101)
;
;all_energy_pairs = []
;
;print, 'Now, P Scan, n-side events:'
;
;for j=0, n_elements(event_files_p)-1 do begin
;		restore, event_files_p[j]
;
;		p_energy_pairs = []
;		data = eventwise_spectra.EVENTDATA_SPEC
;		as = 1
;		for i=0, n_elements(data[as,0,*])-1 do begin
;			hits = where(data[as,*,i] NE 0.)
;			;if n_elements(hits) EQ 2 then begin
;			if n_elements(hits) EQ 2 and total(data[as,*,i]) GE ER[0] and total(data[as,*,i]) LE ER[1] then begin
;				if hits[1] EQ hits[0]+1 then p_energy_pairs = [[p_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
;				if hits[1] EQ hits[0]+1 then all_energy_pairs = [[all_energy_pairs], [data[as,hits[0],i], data[as,hits[1],i]]]
;			endif
;		endfor
;		
;		;print, size(p_energy_pairs)
;		
;		rats = (p_energy_pairs[0,*] - p_energy_pairs[1,*]) / (p_energy_pairs[0,*] + p_energy_pairs[1,*])
;		
;		;rat_hist = histogram(abs(rats), min=0, max=1, bin = 0.02)
;		rat_hist = histogram(rats, min=-1, max=1, bin = 0.02)
;		;print, size(rat_hist)
;		
;		n_rat_hists_pside[j,*] = rat_hist
;endfor
;
;;write_csv, 'all_double_energies_pscan_nside.csv', all_energy_pairs
;
;;write_csv, 'position_double_energies_pscan_nside.csv', n_rat_hists_pside
;write_csv, 'highe_position_double_energies_pscan_nside.csv', n_rat_hists_pside
;

stop

;=========================================================================================
;=========================================================================================


;cgPS_Open, 'ratio_histogram.ps', /nomatch
;	
;cgLoadct, 8, /reverse
;
;positions = cgLayout([1,2]) ;, OYMargin=[10,12], YGap=8, XGap=2, OXMargin=[10,12], aspect=0.8)
;
;
;=========================================================================================
;=========================================================================================


;cgImage, h, charsize=0.6, /Axes

		
;cgImage, p_rat_hists_norm, charsize=0.6, position=positions[*,0], /Axes, $
;		ytitle='Ratio Abs((E1-E2)/(E1+E2))', xtitle='Beam Position #';, charsize=0.6;, stretch=4, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
;		;AXKEYWORDS=axis_format, title='Sub-Strip Image, Beam at Center'
;
;cgImage, p_rat_hists_norm[0:20, *], charsize=0.6, position=positions[*,1], /Axes, $
;		ytitle='Ratio Abs((E1-E2)/(E1+E2))', xtitle='Beam Position #';, charsize=0.6;, stretch=4, $ ;MinValue=0, MaxValue=bigmax, charsize=0.8, $
		;AXKEYWORDS=axis_format, title='Sub-Strip Image, Beam at Center'


;=========================================================================================
;=========================================================================================

;cgPS_Close
	
;spawn, 'open ratio_histogram.ps'

;=========================================================================================
;=========================================================================================











STOP 
END