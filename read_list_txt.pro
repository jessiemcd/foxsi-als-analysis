PRO read_list_txt, file=file, tab=tab

default, file , 'C:\Users\SMusset\Documents\GitHub\Jets\flare_list_oct18.txt'
;data = READ_ASCII(file, delimiter=[' ',':','/'])
;listjets = data.field1

openr, lun, file, /get_lun
array = ''
line = ''
liste = list()
WHILE NOT EOF(lun) DO BEGIN
  readf, lun, line
  array = [array, line]
  sptline = strsplit(line,/extract)
  liste.add, sptline
ENDWHILE

tab = liste.toarray()

FREE_LUN, lun

;dates = listjets[0,*]
;tb = listjets[1,*]
;te = listjets[2,*]
END
