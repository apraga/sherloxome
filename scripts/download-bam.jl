function url(sequencer::String, capture::String, coverage::String, patient::String)
    root =  "https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam"
    dir = "$sequencer/wes_$capture/$coverage"
    "$root/$dir/$patient.$sequencer.wes-$capture.$coverage.dedup.grch38.bam"
end

function download(sequencer::String, capture::String, coverage::String)
    outdir = "../baid2020/grch38/bam/$sequencer/wes_$capture/$coverage"
    mkpath(outdir)
    patients = map(x -> "HG00$x", 1:7)
    urls = map(x -> url(sequencer, capture, coverage, x), patients)
    for url in urls
        run(Cmd(`sdf get $url`, dir=outdir))
    end
end

seqs = ["hiseq4000", "novaseq"]
capts = ["agilent", "idt", "truseq"]
covs = ["50x"]
_ = [download(seq, capt, cov) for (seq, capt, cov) in Iterators.product(seqs, capts, covs)]
