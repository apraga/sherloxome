HG001 = HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG001_GRCh38_1_22_v4.2.1_benchmark.bed
HG002 = HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG003 = HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG004 = HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG005 = HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG005_GRCh38_1_22_v4.2.1_benchmark.bed
HG006 = HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG006_GRCh38_1_22_v4.2.1_benchmark.bed
HG007 = HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG007_GRCh38_1_22_v4.2.1_benchmark.bed

GIAB = ${HG001} ${HG002} ${HG003} ${HG004} ${HG005} ${HG006} ${HG007}

# Store into suddirectory ref
REF = $(GIAB:%=ref/%)
all: truth

truth: $(REF)

# Custom rule to download each patient
ref/HG001_%:
	wget -P ref https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_$*
ref/HG002_%:
	wget -P ref https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_$*
ref/HG003_%:
	wget -P ref https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_$*
ref/HG004_%:
	wget -P ref https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_$*
ref/HG005_%:
	wget -P ref https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_$*
ref/HG006_%:
	wget -P ref https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG006_NA24694_father/latest/GRCh38/HG006_$*
ref/HG007_%:
	wget -P ref https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG007_NA24695_mother/latest/GRCh38/HG007_$*

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
