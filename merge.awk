# Must be called with awk -f merge.awk analysis/*.summary.csv
NR==1 { print "patient,sequencer,kit,depth,caller,"$0}
!/METRIC/ {
  # Fileformat PATIENT-SEQUENCER-KIT-DEPTH.CALLER.summary.csv
  # Get metadata from filename
  split(FILENAME, meta, "-", seps)
  # Get variant caller
  split(meta[4], meta2, ".", seps)
  # Remove leading directory
  sub("analysis/", "",meta[1])
  # patient,sequencer,kit,depth,caller
  print meta[1]","meta[2]","meta[3]","meta2[1]","meta2[2]","$0
}
