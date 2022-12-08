# foxsi-als-analysis
Code written for analysis of data from FOXSI detector tests at the Advanced Light Source, mostly in IDL. There are a few helper routines sourced from various non-FOXSI IDL libraries included here in addition to the FOXSI-specific code (documentation from the creators remains unmodified in those files).

The ALS data that it is designed to work with should be downloaded, and everything is designed to run from the data directory. You can get the FOXSI ALS data from the google drive, or from the UMN DRUM Repository. 

DRUM Link: https://conservancy.umn.edu/handle/11299/228019 

Google Drive Link (for FOXSI team members - need access): https://drive.google.com/drive/folders/1joI_Zi_4FK68YSnat7qNZb1SO_KhMBso?usp=sharing 

Note that if you want to use "real" (arbitrary) beam positions, these are availible in the FOXSI drive, file formats like sample_pos_20190419_overnight_asic2.txt. The same information contained in these files is also in the summary CSV files in the DRUM repository, so if using DRUM to get your data you'll need to make your own .txt files containing all the positions. 

You will also need this FOXSI calibration template file in the directory with your data: https://github.com/foxsi/calsoft/blob/master/template_cal.sav

Finally, you will need to have the FOXSI CALSOFT and FOXSI SCIENCE codes in your IDL path (you can find these at https://github.com/foxsi/calsoft and https://github.com/foxsi/foxsi-science). 

Instructions for using all other procedures/etc. are contained in ALS_analysis_full_process.pro. It is set up so you can read the notes and then run things a few lines at a time, depending on which parts of the analysis process you're interested in replicating. Everything is set up to be used specifically with the fine, overnight, L-shaped scan that occured from '2019/04/19 20:59' to '2019/04/20 09:24'. 

