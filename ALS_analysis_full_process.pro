;==============================================================================================================================
;What is this?
;==============================================================================================================================
;Documents process of performing analysis on L-Scan data from ALS:
;-Processing data into structure files
;-Examining noise + picking secondary thresholds for event selection
;-Gain calibration + processing data into energy space event files

;Many other parts of the analysis use the output of the above (energy-space event files):

;-Plotting nice spectra
;-Event ratios analysis
;-Energy ratios analysis


;Side investigations:
;-Do we see polarization? 

;General Comments: 

;The terms Al-side (n-side) and Pt-side (p-side) are often used interchangeably. The use of 
;the different terms in different scenarios does not necessarily mean that the code is meant 
;to be used with CdTe (Si) detectors. 


;==============================================================================================================================

;==============================================================================================================================
;Processing data into structure files
;==============================================================================================================================
;To start, note that an IDL session should be run in the directory containing the data files.
;Additionally, the file template_cal.sav should be in the data directory as well. 
;(https://github.com/foxsi/calsoft/blob/master/template_cal.sav)

;Find files from the fine L-Scan (entire duration)
;Note: there is a period of time between 04:30 and 05:12 on 4/20 where the beam went off- 
;remove these files from your files directory (they should be notably smaller than the rest).
foxsi_cal_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=files

;Make structure files:
;Note: this takes a while.
foxsi_cal_process_usb, files

;Note: once structure files are made, to identify them use foxsi_cal_struct_filename_timerange.pro
;like
;foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=struct_files
;==============================================================================================================================

;==============================================================================================================================
;Examining noise + picking secondary thresholds for event selection
;==============================================================================================================================

;Next: examine channel-space spectrum in each strip to characterize noise in non-illuminated
;strips + beam spread. 
;This procedure is fairly hard-coded to this specific scan (e.g. it does this only for
;ASICS 1 and 2, the illuminated ASICs during the L-scan, specific illuminated strips are hard-coded, etc.) 
noise_counts

;From the above, different sigma-values for the noise pedestal will be printed. 
;Moving forward, I am using the 3.5 sigma values, which (for both sides) corresponds 
;to a predicted ~1 noise count per file. 

thrp2= 7.3 ;ADC, 3.5-sigma (from noise_counts)
thrn2=7.9 ;ADC, 3.5-sigma (from noise_counts)

;Also adding estimates of 4 keV in ADC space:
thrn=22.
thrp=25.

thresholds = [thrn, thrn2, thrp, thrp2]
;==============================================================================================================================


;==============================================================================================================================
;Gain Calibration + making files with events in energy (keV) units
;==============================================================================================================================


;Gain calibration: there is a lot under the hood here, and a lot hard-coded to using this
;specific scan (It's called 'boutique' because it's highly specialized!). 
;However, there are several options to look through (e.g. using only ALS beam
;peaks vs. incorporating peaks from an older peaksfile as well). 
;
;Results will be two sets of images: 
;		channel-space spectrum (showing fits to line peak locations)
;		channel-space vs. energy-space plots (with fits, defining gain functions)
;
;Additionally, a .sav file (boutique_cal.sav) contains newly generated gain functions 
;for each illuminated strip. 

boutique_cal, singlestrip=1, secondary=1, usepeaks=0, thresholds=thresholds

;Now that we have our thresholds + saved strip-specific gain functions, it's time to use them to convert events
;from ADC space to energy space.

;find the files
foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=struct_files
;Note options for input thresholds, as well as functional form of gain to use. This will make a new set of files.
;It calls one of two functions (eventdata_boutique.pro, eventdata_boutique_chanthresh.pro) depending on whether you
;want to use multiple thresholds (chanthresh=1 version is for multiple thresholds). 
spectra_eventwise_boutique, files=struct_files, chanthresh=1, thresholds=thresholds, funcform='quadratic'

;To examine the results of the gain calibration in a way that is nice to compare to the output plots from boutique_cal.pro,
;run the following (makes spectrum over all files for each strip)
energy_space_strip_spectra

;Note that a final correction to the Al-side gain was made at the spectrum-plotting stage, to place CSA 
;(rather than single-strip) peaks at expected line locations.
;==============================================================================================================================


;==============================================================================================================================
;Plotting different formats of spectra
;==============================================================================================================================

;To find energy-space event files:
foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=2

;For a certain scan (input files in scan), generate spectra like those described in pretty_spectra.pro above, and then stitch
;them all together into a movie (.mp4) to help visualize spectral evolution over the course of the scan.
pretty_spectra_movie, files=event_files_p
pretty_spectra_movie, files=event_files_n, nside_on=1

;For a bunch of files, plot a bunch of individual beam position spectra in a contact-sheet format
pretty_spectra, files=event_files_p
pretty_spectra, files=event_files_n, side='nside'

;To make the figure from Jessie's dissertation (or the SPIE paper). The first time you run it, run it with 
;edit_correction=1 - this makes the additional gain correction (based on CSA peak locations) visible. 
;Then, run again with edit_correction=0 once you are satisfied that the correction is good (CSA peaks correctly
;found) - this makes the nice figure. 
compare_spectra, edit_correction=1
compare_spectra, edit_correction=0


;==============================================================================================================================

;==============================================================================================================================
;Event ratios analysis - what % of events are double-strip at each position?
;==============================================================================================================================

;To make the event ratios figure from Jessie's dissertation and/or ALS SPIE paper:

;Note: to make this figure, you need the position files corresponding to the scans, found in the FOXSI DRIVE
;Location: https://drive.google.com/drive/folders/1pbKWY9-Dwn16XFKUCHRApLQIOo8oZama?usp=sharing
;sample_pos_20190419_overnight_asic2.txt, sample_pos_20190419_overnight_asic1.txt

;To find energy-space event files:
foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files_p, new=2
foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files_n, new=2


;One side at a time - SPIE PAPER STYLE (for dissertation style, simply set livetime_corr=0):
event_type_ratios_figure, files=event_files_p, energy_range=[17.,24.], pickasic=2, livetime_corr=1
event_type_ratios_figure, files=event_files_n, energy_range=[17.,24.], pickasic=1, livetime_corr=1, side='nside'

;The above will also print sizes of double-strip event peaks (individually and on average)

;Also, a less-polished simple version which could be a good start to look at other scans. Note that the default
;position file is from one branch of the L-scan (Pt-side branch), so to use another scan make sure to specify the 
;position file (keyword: POSFILE)
simple_event_type_ratios, files=event_files_p, energy_range=[17.,24.], pickasic=2, side='pside' 
simple_event_type_ratios, files=event_files_n, energy_range=[17.,24.], pickasic=1, side='nside', posfile='sample_pos_20190419_overnight_asic1.txt'

;==============================================================================================================================

;==============================================================================================================================
;Energy ratios analysis - what is the distribution of ratios between the two energies in double-strip events 
;as a function of position?
;==============================================================================================================================

;This analysis was done partially in IDL and partially in python, because I started in IDL (as with everything else), but 
;became so frustrated with IDL's plotting options for representing 2D histograms with other functions overlayed that I 
;decided to export the array I was working with to .csv and do it in python instead. Apologies if this is an inconvenience to 
;anyone seeking to use this code for similar analysis in the future ¯\_(ツ)_/¯.

;Remember earlier, when we chose thresholds to use for analysis, including for making files full of events in energy-space? 
;By picking our secondary (lower) thresholds, we essentially bounded the possible energy ratio and eliminated the vast 
;majority of noise events (not actually caused by charge sharing). To look at energy ratios, we want a looser criteria – 
;we will lower the secondary thresholds, including more noise events but also making sure we aren't losing real double-strip
;events that have one very small energy.

;You might want to copy your struct files into a different directory to do this– otherwise you'll overwrite your existing 
;event files. That's not the end of the world, depending on whether having to re-run these procedures to re-make them with
;different thresholds is more or less annoying to you than having a fair amount of storage space taken up by all these files. 
;For example, all 342 energy-space event files for the full fine L-scan are collectively ~4 GB in size. 

;-20 ADC = basically no threshold for adjacent events
thrn=22.
thrp=25.
thresholds = [thrn, -20, thrp, -20]

;find the files
foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=struct_files
spectra_eventwise_boutique, files=struct_files, chanthresh=1, thresholds=thresholds, funcform='quadratic'


;Then, run:
sharing_energetics

;A good check to see if you are thresholding appropriately to your wishes is to note the minimum and maximum energy ratio
;values found by this procedure (printed). If you are trying to use essentially no secondary thresholds (as I did to make
;the figures in my dissertation + SPIE paper) the absolute values of the max/min should be > 0.99

;If you're using a different threshold, you can figure out what min/max energy ratio makes sense by doing it out:
;R_E_max = ((E_max - lower_threshold) - lower_threshold)/E_max
;R_E_min = (lower_threshold - (E_max - lower_threshold))/E_max
;Where E_max is the maximum incident photon energy expected from the beam (here, 28 keV)

;This makes two files:
;'highe_position_double_energies.csv'
;'highe_position_double_energies_nside.csv'
;
;Each containing an array which is a 1D histogram in energy-ratio at each beam position (so a single array with dimensions 
;101xn where 101 was the number of bins I chose to use in energy ratio space (ranging from -1 to 1), and n is the number of 
;beam positions on each side of the scan.)

;After this, I pop over to a jupyter notebook, import the CSV files, and use some handy functionalities of matplotlib to make
;figures displaying the energy ratio distribution arrays (and analyzing their behavior). A sample notebook is the below, which 
;includes a number of steps culminating with the generation of the figure from the SPIE paper (which is very similar to 
;the figure in Jessie's dissertation, but with color bars requested by a co-author in the revision process.)

;Example notebook: energy_ratio_analysis_and_figure_notebook.ipynb

;The notebook requires you to have installed numpy, matplotlib, pandas, and lmfit
;It works with python 3.10 (after minor changes to the matplotlib imshow() "norm" command), and it worked on python 3.6 as well 

;==============================================================================================================================

;==============================================================================================================================
;Imaging - using knowledge of charge sharing gained in other steps to make images of the ALS beam at various positions. 
;==============================================================================================================================

;Before starting imaging, go back to the old thresholds if you were just using the "no-threshold" thresholds to 
;look at energy ratios. (If not, you can skip this. Though it won't hurt.)

thrp2= 7.3 ;ADC, 3.5-sigma (from noise_counts)
thrn2=7.9 ;ADC, 3.5-sigma (from noise_counts)

;Also adding estimates of 4 keV in ADC space:
thrn=22.
thrp=25.

thresholds = [thrn, thrn2, thrp, thrp2]

;find the files
foxsi_cal_struct_filename_timerange, ['2019/04/19 20:59','2019/04/20 09:24'], file=struct_files
;Note options for input thresholds, as well as functional form of gain to use. This will make a new set of files.
;It calls one of two functions (eventdata_boutique.pro, eventdata_boutique_chanthresh.pro) depending on whether you
;want to use multiple thresholds (chanthresh=1 version is for multiple thresholds). 
spectra_eventwise_boutique, files=struct_files, chanthresh=1, thresholds=thresholds, funcform='quadratic'

;First load in files. Figure in SPIE paper/dissertation uses all files from one branch of the scan. Since
;the scan was slightly diagonal, this still allows us to explore boundary and center regions on both sides. 
foxsi_cal_eventdata_filename_timerange, ['2019/04/19 20:59','2019/04/20 03:11'], file=event_files, new=2

;To make figure:
;Note: pdf output will not appear correctly when opened in Mac Preview (looks weirdly blurry).
;If you have this problem, save the file and open it with Adobe Reader or something.
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=0, nicefigure=1

;To make energy ratio based video of Pt-side scan:
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=1, vidtype='ratio', nicefigure=0
;To make zone-based video of Pt-side scan:
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=1, vidtype='zone', nicefigure=0
;To make strip-based video of Pt-side scan:
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=1, vidtype='strip', nicefigure=0

;Now, let's make movies for the al-side branch of the scan.
foxsi_cal_eventdata_filename_timerange, ['2019/04/20 03:12','2019/04/20 09:24'], file=event_files, new=2

;To make energy ratio based video of Al-side scan:
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=1, vidtype='ratio', nicefigure=0, branch='Al'
;To make zone-based video of Pt-side scan:
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=1, vidtype='zone', nicefigure=0, branch='Al'
;To make strip-based video of Pt-side scan:
energy_ratio_imaging_new, files=event_files, erange=[17.,30.], video=1, vidtype='strip', nicefigure=0, branch='Al'

;Note: in the Al-side branch movies, you can visually see the "skip" when the beam encounters strips that
;are not enabled (oops!)

;Another note: energy_ratio_imaging_new.pro includes three functions which each make a single image from a
;FOXSI data file, using the three different imaging methods. These functions can be used separately from the
;procedure as a whole to make images at any beam position if desired. 

;==============================================================================================================================





