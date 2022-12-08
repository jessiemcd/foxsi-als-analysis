; Modified version of  IDL procedure by Athiray
; Copyright (c) 2017, FOXSI Mission University of Minnesota.  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 08 Jul 2019 23:06
; Last Mod 	: 04 Sep 2019 09:16
; New Version : 03 Mar 2022 12:42 - Jessie

;THIS NEW VERSION USES EVENTDATA_BOUTIQUE FUNCTION, WHICH USES GAIN CALIBRATION FUNCTIONS 
;FOR SPECIFIC STRIPS (THOSE COVERED IN L-BRANCH 2019 ALS SCAN)
;STRIP FUNCTION PARAMETERS FOUND AND SAVED IN PROCEDURE BOUTIQUE_CAL.PRO, RESTORED IN 
;FUNCTION eventdata_boutique_chanthresh

;ALSO ALLOWS FOR INPUT CHOICE BETWEEN LINEAR AND QUADRATIC GAIN FUNCTIONS (FUNCFORM)
;ALSO ALLOWS FOR INPUT THRESHOLDS (THRESHOLDS)



;-------------------  Details of the program --------------------------;
PRO spectra_eventwise_boutique, files=files, chanthresh=chanthresh, $
								THRESHOLDS=THRESHOLDS, FUNCFORM=FUNCFORM

;One center and one boundary in each L-scan branch, for examining calibration.
deffiles = ['struct_data_20190419_230246.txt', 'struct_data_20190419_224006.txt', $
			'struct_data_20190420_061854.txt', 'struct_data_20190420_064131.txt']

default, files, deffiles
;Set chanthresh = 1 to use MULTIPLE CHANNEL SPACE THRESHOLDS
default, chanthresh, 0
default, FUNCFORM, 'linear'

;TO CHO0SE BETWEEN ONE CHANNEL SPACE THRESHOLD WHICH IS APPROXIMATELY 4 KEV AND AN ENERGY-SPACE THRESHOLD 
;OF 4 KEV, EDIT BELOW.

badch=intarr(4,64)
badch(0,0)=1
badch(0,1)=1
badch(0,2)=1
badch(0,3)=1
badch(0,4)=1
badch(0,5)=1

badch(0,8)=1
badch(0,9)=1
badch(0,10)=1
badch(0,11)=1

badch(0,18)=1
badch(0,19)=1

badch(0,21)=1
badch(0,22)=1

badch(0,26)=1
badch(0,27)=1

badch(0,30)=1
badch(0,31)=1
badch(0,32)=1
badch(0,33)=1

badch(0,35)=1
badch(0,36)=1

badch(0,42)=1
badch(0,43)=1
badch(0,44)=1
badch(0,45)=1

badch(1,32)=1
badch(1,33)=1

badch(1,63)=1

badch(2,63)=1

badch(3,0)=1
badch(3,1)=1

badch(3,17)=1
badch(3,18)=1

print, badch


pside_singles_fraction=dblarr(n_elements(files))
nside_singles_fraction=dblarr(n_elements(files))
pside_doubles_fraction=dblarr(n_elements(files))
nside_doubles_fraction=dblarr(n_elements(files))
pside_triples_fraction=dblarr(n_elements(files))
nside_triples_fraction=dblarr(n_elements(files))
for i =0,n_elements(files)-1 do begin
    f=files[i]
    ;FOR ENERGY SPACE 4 KEV THRESHOLD:
    ;if chanthresh EQ 0 then eventwise_spectra=eventdata_boutique(f, /cmn_med, thrn=4, badch=badch)
    ;FOR CHANNEL-SPACE "4 KEV" THRESHOLD:
    if chanthresh EQ 0 then eventwise_spectra=eventdata_boutique(f, /cmn_med, CHAN=1, badch=badch, $
    															FUNCFORM=FUNCFORM)
    ;FOR TWO-THRESHOLD CHANNEL SPACE METHOD (APPROX. 4, 2.5 KEV):
    if chanthresh EQ 1 then eventwise_spectra=eventdata_boutique_chanthresh(f, /cmn_med, /chanthresh, $
    															badch=badch, THRESHOLDS=THRESHOLDS, $
    															FUNCFORM=FUNCFORM)
    fname=strsplit(f,'_',/extract)
    f1=strsplit(fname[3],'.',/extract)
    if chanthresh EQ 1 then specfname = 'eventwisesp_boutique_'+fname[2]+'_'+f1[0]+'.txt'
    if chanthresh EQ 0 then specfname = 'eventwisesp_boutique1thresh_'+fname[2]+'_'+f1[0]+'.txt'
    ;specfname = 'new_eventwisesp_'+fname[2]+'_'+f1[0]+'.txt'
    save, eventwise_spectra, file = specfname
    zeroles=eventwise_spectra.zeroles
    singles=eventwise_spectra.singles
    doubles=eventwise_spectra.doubles
    triples=eventwise_spectra.triples
    allples=eventwise_spectra.allples
    n_evts=eventwise_spectra.n_evts
    ;;useable_evts=n_evts[0]-zeroles
    ;;pside_singles_fraction[i] = total(singles[2:3])/total(useable_evts[2:3])
    ;;nside_singles_fraction[i] = total(singles[0:1])/total(useable_evts[0:1])
    ;;pside_doubles_fraction[i] = total(doubles[2:3])/total(useable_evts[2:3])
    ;;nside_doubles_fraction[i] = total(doubles[0:1])/total(useable_evts[0:1])
    ;;pside_triples_fraction[i] = total(triples[2:3])/total(useable_evts[2:3])
    ;;nside_triples_fraction[i] = total(triples[0:1])/total(useable_evts[0:1])
endfor
stop

END

