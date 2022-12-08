PRO foxsi_cal_process_usb, files
  
  ;+
  ;  PROJECT
  ;    Foxsi Rocket Detector Calibration
  ;
  ;  DESCRIPTION
  ;    This procedure will process several usb files to foxsi structure files
  ;    in the current folder
  ;
  ;  INPUTS
  ;    files: an array of usb files of the form data_yyyymmdd_hhmmss.txt
  ;
  ;  EXAMPLE
  ;    foxsi_cal_filename_timerange, '2019/04/16 '+['18:00','22:00'], file=files
  ;    foxsi_cal_process_usb, files
  ;
  ;  HISTORY
  ;    2019/04/18 SMusset (UMN) - initial release
  ;-
  
  FOR k=0, n_elements(files)-1 DO BEGIN
    f=files[k]
    data = read_data_struct_cal(f)
    str = 'struct_'+f
    save, data, file=str
  END
  
END