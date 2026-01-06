all: by-patient-giab.csv by-kit-giab.csv by-sequencer-giab.csv by-depth-giab.csv

by-patient-giab.csv: ../analysis/all-giab.csv
	nu -c "open $? | where FILTER == PASS | polars into-df \
	| polars pivot -o [patient] \
	-i  [sequencer kit depth caller Type Filter "METRIC.Recall"] \
	--values ["METRIC.Precision"] | polars save $@"

by-kit-giab.csv: ../analysis/all-giab.csv
	nu -c "open $? | where FILTER == PASS | polars into-df \
	| polars pivot -o [kit] \
	-i  [sequencer patient depth caller Type Filter "METRIC.Recall"] \
	--values ["METRIC.Precision"] | polars save $@"

by-sequencer-giab.csv: ../analysis/all-giab.csv
	nu -c "open $? | where FILTER == PASS | polars into-df \
	| polars pivot -o [sequencer] \
	-i  [kit patient depth caller Type Filter "METRIC.Recall"] \
	--values ["METRIC.Precision"] | polars save $@"

by-depth-giab.csv: ../analysis/all-giab.csv
		nu -c "open $? | where FILTER == PASS | polars into-df \
		| polars pivot -o [depth] \
		-i  [kit patient sequencer caller Type Filter "METRIC.Recall"] \
		--values ["METRIC.Precision"] | polars save $@"
