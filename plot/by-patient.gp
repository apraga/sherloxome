set terminal qt persist
set key outside
set datafile  separator comma
plot for [j=8:14] 'by-patient-giab.csv' u 7:j title columnheader(j) pt 7
