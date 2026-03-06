run-agilent:
	nextflow run nf-core/sarek --input samplesheet.csv  -r 3.5.1 --outdir giab -c conf/slurm.config \
	--tools haplotypecaller --skip_tools haplotypecaller_filter	 \
	--wes --intervals ./capture/Agilent_SureSelect_All_Exons_v7_hg38_Regions.bed -bg

run-idt:
	nextflow run nf-core/sarek --input samplesheet.csv  -r 3.5.1 --outdir giab -c conf/slurm.config \
	--tools haplotypecaller --skip_tools haplotypecaller_filter	 \
	--wes --intervals ./capture/xgen-exome-hyb-panel-v2-targets-hg38.bed -bg

run-truseq:
	nextflow run nf-core/sarek --input samplesheet.csv  -r 3.5.1 --outdir giab -c conf/slurm.config \
	--tools haplotypecaller --skip_tools haplotypecaller_filter	 \
	--wes --intervals ./capture/truseq-dna-exome-targeted-regions-manifest-v1-2-lifted-grch38.bed -bg

compare:
	bash compare.sh

merge:
	awk -f merge.awk analysis/*.summary.csv > analysis/all-giab.csv
