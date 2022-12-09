# foxsi-als-analysis
Code written for analysis of data from FOXSI detector tests at the Advanced Light Source, mostly in IDL. There are a few helper routines sourced from various non-FOXSI IDL libraries included here in addition to the FOXSI-specific code (documentation from the creators remains unmodified in those files).
This analysis is detailed in an SPIE paper https://doi.org/10.1117/12.2629443, as well as in Jessie Duncan's dissertation https://hdl.handle.net/11299/241752. 

The ALS data that it is designed to work with should be downloaded, and everything is designed to run from the data directory. You can get the FOXSI ALS data from the google drive, or from the UMN DRUM Repository. 

DRUM Link: https://conservancy.umn.edu/handle/11299/228019 

Google Drive Link (for FOXSI team members - need access): https://drive.google.com/drive/folders/1joI_Zi_4FK68YSnat7qNZb1SO_KhMBso?usp=sharing 

Note that if you want to use "real" (arbitrary) beam positions for plotting how any properties of the data vary along a scan, these are availible in the FOXSI drive in .txt files with name formats like "sample_pos_201904*.txt". The same information (x,y positions for every data file) contained in these files is also in the summary CSV files in the DRUM repository, so if using DRUM to get your data you'll need to make your own .txt files containing all the beam positions. 

Place this FOXSI calibration template file in the directory with your data: https://github.com/foxsi/calsoft/blob/master/template_cal.sav

You will need to have the FOXSI CALSOFT and FOXSI SCIENCE codes in your IDL path (you can find these at https://github.com/foxsi/calsoft and https://github.com/foxsi/foxsi-science). 

Finally, to make the figure comparing different imaging methods (Figure 5 in https://doi.org/10.1117/12.2629443), you will need the Coyote IDL Library. This is not necessary for making FOXSI images in general, it was just used to make a nice multi-pannel figure. If you don't have coyote-idl yet, you can find it here: http://www.idlcoyote.com/code_tips/installcoyote.php

Instructions for using all other procedures/etc. are contained in ALS_analysis_full_process.pro. It is set up so you can read the notes and then run things a few lines at a time, depending on which parts of the analysis process you're interested in replicating. Everything is set up to be used specifically with the fine, overnight, L-shaped scan that occured from '2019/04/19 20:59' to '2019/04/20 09:24'. 

