FUNCTION make_time_stamp, tim, date=date, time=time
; this function make time stamp in the form '20110215_011248'
; 

time_array = anytim(tim, /ext)
time_string = string(time_array)
time_str = time_string
FOR k=0, n_elements(time_string)-1 do TIME_STR = STRSPLIT(TIME_STRING, /EXTRACT)
; correct for number under 10
FOR k=0, 5 DO IF time_array[k] lt 10 THEN TIME_STR[k] = '0'+TIME_STR[k]

TIME_STAMPP = time_str[6]+time_str[5]+time_str[4]+'_'+time_str[0]+time_str[1]+time_str[2]
IF keyword_set(date) THEN TIME_STAMPP=time_str[6]+time_str[5]+time_str[4]
IF keyword_set(time) THEN TIME_STAMPP=time_str[0]+time_str[1]+time_str[2]

return, TIME_STAMPP
end