set terminal qt persist
set key outside
set datafile  separator comma
plot for [j=8:10] 'by-kit-giab.csv' u 7:j title columnheader(j) pt 7
