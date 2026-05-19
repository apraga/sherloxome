`sherloxome` is a small command-line utility to help you benchmark your genomics pipeline. It is aimed at exome or targeted sequencing.
A setup phase download FASTQ for reference patient and generate in silico data according to [a configuration file](03-configuration.md).
A samplesheet is generated, aimed at `nf-core/sarek` pipeline, but any pipeline can be used.
Then all VCF file in a folder are analyzed regarding to their reference and a neat summary and plots are shown. This step is quite flexible, providing file follows [our filenaming scheme](./06-filenaming.md).


