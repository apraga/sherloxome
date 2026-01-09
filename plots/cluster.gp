set terminal qt persist size 800,1000
set datafile  separator comma
set multiplot layout 4,2

array files[4] = ["by-patient-giab.csv", "by-kit-giab.csv", "by-sequencer-giab.csv", "by-depth-giab.csv"]
array jends[4] = [14, 10, 9, 10]
array types[2] = ["INDEL", "SNP"]

do for [i=1:4] {
    do for [t=1:2] {
        if (i == 1) {  set title types[t] } else { set title }
        if (t==2) { set key right top } else { set key off }
        if (t==1) {
            set xrange [0.7:1]
            set yrange [0.7:1]
        } else {
            set xrange [0.95:1]
            set yrange [0.95:1]
        }
        if ( i > 1) { set xtics format ""}

        plot for [j=8:jends[i]] files[i] u (stringcolumn(5) eq types[t] ? $7 : NaN):j  title columnheader(j) pt 7
    }
}

#set key off
#set title "INDEL"
#plot for [j=8:14] 'by-patient-giab.csv' u (stringcolumn(5) eq "INDEL" ? $7 : NaN):j  title columnheader(j) pt 7
#set key outside
#set title "SNP"
#plot for [j=8:14] 'by-patient-giab.csv' u (stringcolumn(5) eq "SNP" ? $7 : NaN):j  title columnheader(j) pt 7
#
##plot for [j=8:9] 'by-sequencer-giab.csv' u 7:j title columnheader(j) pt 7
##plot for [j=8:10] 'by-depth-giab.csv' u 7:j title columnheader(j) pt 7
unset multiplot
