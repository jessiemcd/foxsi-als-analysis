; Modified version of IDL procedure by Athiray
; Copyright (c) 2017, FOXSI Mission University of Minnesota.  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 08 Jul 2019 21:16
; Last Mod 	: 08 Jul 2019 23:52
; New Version : 03 Mar 2022 12:16 - Made by Jessie

;THIS NEW VERSION USES GAIN CALIBRATION FUNCTIONS FOR SPECIFIC STRIPS 
;(THOSE COVERED IN THE L-BRANCH 2019 ALS SCAN)
;STRIP FUNCTION PARAMETERS FOUND AND SAVED IN PROCEDURE BOUTIQUE_CAL.PRO

;Additionally, this new version allows for two different methods of input primary and secondary
;thresholds for event selection (e.g. determination of what counts as an event above the noise floor
;, and event multiplicity)


;-------------------  Details of the program --------------------------;
function eventdata_boutique_chanthresh, FILE, THRP = THRP, THRN = THRN , $
              BADCH=BADCH, NMAX = NMAX, SUBTRACT_COMMON = SUBTRACT_COMMON, $
              CMN_AVERAGE = CMN_AVERAGE, CMN_MEDIAN = CMN_MEDIAN, STOP = STOP, $
              CHANTHRESH=CHANTHRESH, THRESHOLDS=THRESHOLDS, FUNCFORM=FUNCFORM
              
  if not keyword_set(FUNCFORM) then FUNCFORM = 'linear'
              
  ;Pre-known relevant strips
	relevant_strips_p = indgen(10)+33
	relevant_as_p = 2
	relevant_strips_n = indgen(7)+25
	relevant_as_n = 1
	

  if not keyword_set(thrp) then thrp=4.0
  if not keyword_set(thrn) then thrn=5.0
  if not keyword_set(badch) then badch=intarr(4,64)
  
  if keyword_set(chanthresh) then begin
  		;to set thresholds like were used in Jessie's dissertation and SPIE paper
  		thrp1=25. ;ADC, ~4 keV
  		;thrp2= 14. ;ADC, ~2 keV
  		thrp2= 7.3 ;ADC, 3.5-sigma (from noise_counts)
  		;thrp2= -20 ;basically no threshold for adjacent events
  		thrn1=22. ;ADC, ~4 keV
  		;thrn2=11. ;ADC, ~2 keV
  		thrn2=7.9 ;ADC, 3.5-sigma (from noise_counts)
  		;thrn2= -20 ;basically no threshold for adjacent events
  endif
  
  if keyword_set(thresholds) then begin
    ;for user-input thresholds
	thrn1 = thresholds[0]
	thrn2 = thresholds[1]
	thrp1 = thresholds[2]
	thrp2 = thresholds[3]
  endif


  restore, file
  restore, 'boutique_cal.sav'
  ;Contents of this file : LINFITS_P, LINFITS_N, QUADFITS_P, QUADFITS_N
  ;linear and quadratic gain fit parameters for each side
  ;USE FUNCFORM command to pick which to use. 
  ;first strip in the relevant strips corresponds to the first set of coefficients in each file and so on
  
  n_evts = n_elements(data)
  print, n_evts, ' total events'

  cmn = fltarr(n_evts, 4)
  eventdata_spec = dblarr(4,64,n_evts) 

  for as=0, 3 do begin
      if keyword_set(subtract_common) then cmn[*,as] = data[*].common_mode[as] + randomu(seed,n_evts*4+as) - 0.5
      if keyword_set(cmn_average) then cmn[*,as] = data[*].cmn_average[as]
      if keyword_set(cmn_median) then cmn[*,as] = data[*].cmn_median[as] + randomu(seed,n_evts*4+as) - 0.5
  endfor


    if keyword_set(nmax) then n_evts = nmax
    ngood=long(0)

    hitch=0
    for evt = long(0), n_evts-1 do begin

        if (evt mod 1000) eq 0 then print, 'Event  ', evt, ' / ', n_evts
		if max( data[evt].data ) eq 0 then continue

        hitchnump=0
        hitchnumn=0

        if total(data[evt].packet_error) lt 10 then begin
            hitchnump=0
            hitchnumn=0
            as = 1 ; n-side
            ;for each of the L-scan strips on the Nside:
			for i = 0, n_elements(relevant_strips_n)-1 do begin
				ch = relevant_strips_n[i]
				if badch[as,ch] eq 0 then begin
					;Channel space energy in that strip in this frame:
					chane = data[evt].data[as,ch]-cmn[evt,as]
					;If said energy is above primary threshold:
					if chane GE thrn1 then begin
						if FUNCFORM EQ 'linear' then begin
							coeffs = LINFITS_N[*,i]
							edep = chane*coeffs[1]+coeffs[0]
						endif
						if FUNCFORM EQ 'quadratic' then begin
							coeffs = QUADFITS_N[*,i]
							edep = chane^2*coeffs[2]+chane*coeffs[1]+coeffs[0]
						endif
						hitchnumn+=1
						eventdata_spec[as,ch,evt]=edep
						;Channel-space energies in the two adjacent strips
						chane_pl = data[evt].data[as,ch+1]-cmn[evt,as]
						chane_m = data[evt].data[as,ch-1]-cmn[evt,as]
						;If either (or both - which is rare) are above the secondary threshold, 
						;then also add them to the events array
						if chane_pl GE thrn2 then begin
							eventdata_spec[as,ch+1,evt]=chane_pl^2*coeffs[2]+chane_pl*coeffs[1]+coeffs[0]
							;eventdata_spec[as,ch+1,evt]=chane_pl*coeffs[1]+coeffs[0]
							hitchnumn+=1
						endif
						if chane_m GE thrn2 then begin
							eventdata_spec[as,ch-1,evt]=chane_m^2*coeffs[2]+chane_m*coeffs[1]+coeffs[0]
							;eventdata_spec[as,ch-1,evt]=chane_m*coeffs[1]+coeffs[0]
							hitchnumn+=1
						endif
						;BREAK
					endif
				endif
			endfor
            as = 2 ; p-side - for more documentation of loop see n-side
			for i = 0, n_elements(relevant_strips_p)-1 do begin
				ch = relevant_strips_p[i]
				if badch[as,ch] eq 0 then begin
					chane = data[evt].data[as,ch]-cmn[evt,as]
					if chane GE thrp1 then begin
						if FUNCFORM EQ 'linear' then begin
							coeffs = LINFITS_P[*,i]
							edep = chane*coeffs[1]+coeffs[0]
						endif
						if FUNCFORM EQ 'quadratic' then begin
							coeffs = QUADFITS_P[*,i]
							edep = chane^2*coeffs[2]+chane*coeffs[1]+coeffs[0]
						endif
						hitchnump+=1
						eventdata_spec[as,ch,evt]=edep
						chane_pl = data[evt].data[as,ch+1]-cmn[evt,as]
						chane_m = data[evt].data[as,ch-1]-cmn[evt,as]
						if chane_pl GE thrp2 then begin
							if FUNCFORM EQ 'quadratic' then eventdata_spec[as,ch+1,evt]=chane_pl^2*coeffs[2]+chane_pl*coeffs[1]+coeffs[0]
							if FUNCFORM EQ 'linear' then eventdata_spec[as,ch+1,evt]=chane_pl*coeffs[1]+coeffs[0]
							hitchnump+=1
						endif
						if chane_m GE thrp2 then begin
							if FUNCFORM EQ 'quadratic' then eventdata_spec[as,ch-1,evt]=chane_m^2*coeffs[2]+chane_m*coeffs[1]+coeffs[0]
							if FUNCFORM EQ 'linear' then eventdata_spec[as,ch-1,evt]=chane_m*coeffs[1]+coeffs[0]
							hitchnump+=1
						endif
						;BREAK
					endif
				endif
			endfor
        
	    if hitchnump eq 1 and hitchnumn eq 1 then ngood += 1

        endif
    endfor

    print, 'good events: ', ngood, '/', n_evts
    

zeroles=intarr(4,1)
singles=intarr(4,1)
doubles=intarr(4,1)
triples=intarr(4,1)
allples=intarr(4,1)
for i =0,n_evts-1 do begin
     for as=0, 3 do begin
	 selection_index=where(EVENTDATA_SPEC[as,*,i] ne 0,count)
	 if count eq 0 then zeroles[as,0]+=1
	 if count eq 1 then singles[as,0]+=1
	 if count eq 2 then doubles[as,0]+=1
	 if count eq 3 then triples[as,0]+=1
	 if count gt 3 then allples[as,0]+=1
	 ;if count eq 3 and as EQ 1 then print, EVENTDATA_SPEC[as,relevant_strips_n,i]
	 ;if count eq 3 and as EQ 2 then print, EVENTDATA_SPEC[as,relevant_strips_p,i]
	 ;if count eq 3 then print, ''
	 if count ne 0 then begin
	     ;print,selection_index
	     ;print,doubles[as,*],triples[as,*]
	 endif
     endfor
endfor
print,singles,doubles,triples,allples,zeroles
eventwise_spectra = {eventdata_spec:dblarr(4,64,n_evts),	$	
		zeroles:intarr(4,1),	$	;
		singles:intarr(4,1),	$	;
		doubles:intarr(4,1), 	$	; 
		triples:intarr(4,1),    $
		allples:intarr(4,1),    $
	        ngood:intarr(1),	$
	        n_evts:intarr(1)}

eventwise_spectra.eventdata_spec=eventdata_spec
eventwise_spectra.zeroles=zeroles
eventwise_spectra.singles=singles
eventwise_spectra.doubles=doubles
eventwise_spectra.triples=triples
eventwise_spectra.allples=allples
eventwise_spectra.ngood=ngood
eventwise_spectra.n_evts=n_evts

return, eventwise_spectra 

END
