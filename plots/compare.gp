set datafile  separator comma
set term qt persist
patients = "HG001 HG002 HG003 HG004 HG005 HG006 HG007"
kits = "agilent idt"
# Convert patient to integere
patients_ids(s) = sum [i=1:words(patients)] (word(patients,i) eq s ? i : 0)
kits_ids(s) = sum [i=1:words(kits)] (word(kits,i) eq s ? i : 0)
# plot 'analysis/all-giab.csv' using 16:17:(patients_ids(stringcolumn("patient"))) w p pt 7 lc variable
plot 'analysis/all-giab.csv' using 16:17:(kits_ids(stringcolumn("kit"))) w p pt 7 lc variable
