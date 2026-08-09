A simple filenaming convention is used for VCF and FASTQ :

    SAMPLE_SEQUENCER_CAPTURE_DEPTH{_SILICO}.*

Value are free text for maximum compatibility, with the exception of DEPTH that must be an integer followed by the 'x' character. Here is an example

    HG002_novaseq_agilent_100x_simuscop.haplotypecaller.vcf.gz

Where SAMPLE is HG002, SEQUENCER is novaseq, DEPTH is 100. It's in silico data as SILICO is simuscop. The rest of the filename is discarded.
A set of values is defined in according to Baid 2020, with
- Agilent, IDT, Truseq as the capture kit
- HG002... HG007 as patient (GIAB data),
- 50x, 75x and 100x as depth.

Note that not all combinations are available (see [baid2020.rs](https://github.com/apraga/sherloxome/blob/ci/src/baid2020.rs)).
