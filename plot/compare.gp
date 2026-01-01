set datafile  separator comma
patients = "HG001 HG002 HG003 HG004 HG005 HG006 HG007"
# Convert patient to integere
code(s) = sum [i=1:words(patients)] \
          (word(patients,i) eq s ? i : 0)
plot 'analysis/all-giab.csv' using 16:17:(code(stringcolumn("patient"))) w p pt 7 lc variable
