PRO foxsi_cal_struct_filename_timerange, timerange, dir=dir, filenames=filenames
  
  ;+
  ;  PROJECT
  ;    Foxsi Rocket Detector Calibration
  ;    
  ;  DESCRIPTION
  ;		Modified version of foxsi_cal_filename_timerange, this version for struct_ files.
  ;		Made by Jessie December 2021
  ;
  ;		Old Documentation:
  ;
  ;    This procedure will find all the data files in a given time range and can return the list of
  ;    filenames in an array of strings
  ;    This works only for USB files right now, since it is looking for filenames
  ;    of this kind: data_yyyymmdd_*.txt
  ;    
  ;  INPUTS
  ;    timerange: a 2-element array containing a time range
  ;    
  ;  KEYWORDS
  ;    dir (input), string: directory where to look for the files. Default is current directory
  ;    filenames (output), strarr: array of filenames found
  ;    
  ;  CALLS
  ;    make_time_stamp
  ;    
  ;  EXAMPLE
  ;    foxsi_cal_filename_timerange, '2019/04/16 '+['18:00','22:00'], file=files
  ;    
  ;  HISTORY
  ;    2019/04/18 SMusset (UMN) - initial release
  ;-
 
  cd, current=current
  DEFAULT, dir, current
  cd, dir
  
  date = make_time_stamp(timerange[0],/date)
  datefin = make_time_stamp(timerange[-1],/date)
  print, date, datefin
  ;IF date EQ datefin THEN files = file_search('data_'+date+'_*.txt') ELSE files = file_search('data_*_*.txt')
  IF date EQ datefin THEN files = file_search('struct_data_'+date+'_*.txt') ELSE files = file_search('struct_data_*_*.txt')
  dates = strarr(n_elements(files))
  
  
  IF n_elements(files) EQ 1 THEN BEGIN
    IF files EQ '' THEN BEGIN
      print, 'no file found'
      filenames=0
    ENDIF
  ENDIF ELSE BEGIN
  
    FOR k=0, n_elements(files)-1 DO BEGIN
      filesplit = strsplit( files[k], '_.', /extract)  
      ;print, filesplit          
      dates[k] = strmid(filesplit[2],0,4)+'/'+strmid(filesplit[2],4,2)+'/'+strmid(filesplit[2],6,2)+' '+ $
        strmid(filesplit[3],0,2)+':'+strmid(filesplit[3],2,2)+':'+strmid(filesplit[3],4,2)
      ;print, dates[k]
      ;stop
    ENDFOR
  
    select = where(anytim(dates) GE anytim(timerange[0]) AND anytim(dates) LE anytim(timerange[1]))
    print, n_elements(select), ' files found in time range'
    filenames = files[select]
  ENDELSE
  
END