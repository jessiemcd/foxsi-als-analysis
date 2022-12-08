; Modified version of IDL procedure by Athiray
; Copyright (c) 2017, FOXSI Mission University of Minnesota.  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 08 Jul 2019 21:16
; Last Mod 	: 08 Jul 2019 23:52
; New Version : 03 Mar 2022 12:16

;THIS NEW VERSION USES GAIN CALIBRATION FUNCTIONS FOR SPECIFIC STRIPS (THOSE COVERED IN L-BRANCH 2019 ALS SCAN)
;STRIP FUNCTION PARAMETERS FOUND AND SAVED IN PROCEDURE BOUTIQUE_CAL.PRO


;-------------------  Details of the program --------------------------;
function eventdata_boutique, FILE, THRP = THRP, THRN = THRN , $
              BADCH=BADCH, NMAX = NMAX, SUBTRACT_COMMON = SUBTRACT_COMMON, $
              CMN_AVERAGE = CMN_AVERAGE, CMN_MEDIAN = CMN_MEDIAN, CHAN=CHAN, $
              FUNCFORM=FUNCFORM, STOP = STOP
              
  ;Pre-known relevant strips
	relevant_strips_p = indgen(10)+33
	relevant_as_p = 2
	relevant_strips_n = indgen(7)+25
	relevant_as_n = 1


  if not keyword_set(thrp) then thrp=4.0
  if not keyword_set(thrn) then thrn=5.0
  if not keyword_set(badch) then badch=intarr(4,64)
  
  if keyword_set(CHAN) then BEGIN
  		thrp1=25. ;ADC, ~4 keV
  		thrn1=22. ;ADC, ~4 keV
  endif
  

  restore,file
  restore, 'boutique_cal.sav'
  ;Contents of this file : LINFITS_P, LINFITS_N, QUADFITS_P, QUADFITS_N
  ;linear and quadratic gain fit parameters for each side
  ;Comment out lines in loops below to change which function is used.
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
			for i = 0, n_elements(relevant_strips_n)-1 do begin
				ch = relevant_strips_n[i]
				if badch[as,ch] eq 0 then begin
					chane = data[evt].data[as,ch]-cmn[evt,as]
					if FUNCFORM EQ 'linear' then begin
						coeffs = LINFITS_N[*,i]
						edep = chane*coeffs[1]+coeffs[0]
					endif
					if FUNCFORM EQ 'quadratic' then begin
						coeffs = QUADFITS_N[*,i]
						edep = chane^2*coeffs[2]+chane*coeffs[1]+coeffs[0]
					endif
					if CHAN EQ 1 and chane gt thrn1 then begin
						hitchnumn+=1
						hitchn=ch
						hitasicn=as
						eventdata_spec[as,ch,evt]=edep
					endif
					if CHAN EQ 0 and edep gt thrn then begin
						hitchnumn+=1
						hitchn=ch
						hitasicn=as
						eventdata_spec[as,ch,evt]=edep
					endif
				endif
			endfor
            as = 2 ; p-side
			for i = 0, n_elements(relevant_strips_p)-1 do begin
				ch = relevant_strips_p[i]
				if badch[as,ch] eq 0 then begin
					chane = data[evt].data[as,ch]-cmn[evt,as]
					if FUNCFORM EQ 'linear' then begin
						coeffs = LINFITS_P[*,i]
						edep = chane*coeffs[1]+coeffs[0]
					endif
					if FUNCFORM EQ 'quadratic' then begin
						coeffs = QUADFITS_P[*,i]
						edep = chane^2*coeffs[2]+chane*coeffs[1]+coeffs[0]
					endif
					if CHAN EQ 1 and chane gt thrp1 then begin
						hitchnumn+=1
						hitchn=ch
						hitasicn=as
						eventdata_spec[as,ch,evt]=edep
					endif
					if CHAN EQ 0 and edep gt thrp then begin
						hitchnumn+=1
						hitchn=ch
						hitasicn=as
						eventdata_spec[as,ch,evt]=edep
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
