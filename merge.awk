# Must be called with find analysis -name \*.summary.csv -execdir awk -f ../merge.awk {}
NR==1 { print "patient,sequencer,kit,depth,caller,"$0}
NR>1
{
  # Fileformat PATIENT-SEQUENCER-KIT-DEPTH.CALLER.summary.csv
  # Get metadata from filename
  split(FILENAME, meta, "-", seps)
  # Get variant caller
  split(meta[4], meta2, ".", seps)
  # Remove leading dots
  sub("./", "",meta[1])
  # patient,sequencer,kit,depth,caller
  print meta[1]","meta[2]","meta[3]","meta2[1]","meta2[2]","$0
}
