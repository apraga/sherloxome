# Install dependencies

## Rust toolchain

Install Rust via [rustup](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then compile sherloxome:

```bash
cargo build --release
```

The binary is placed at `target/release/sherloxome`. Add it to your `PATH` or call it directly.

## External tools

The following tools must be available in your `PATH`. Only `hap.py` and `rtg` are required for the `analyze` step; `bwa`, `samtools`, and `muteditor` are only needed when generating in silico controls.

| Tool | Minimum version | Purpose |
|------|----------------|---------|
| `samtools` | 1.15 | BAM/FASTQ processing |
| `bwa` | 0.7 | Read alignment for in silico BAM editing |
| `muteditor` | any | Variant injection (in silico controls) |
| `hap.py` | 0.3.15 | Variant call comparison against truth sets |
| `rtg` (RTG Tools) | 3.12 | vcfeval engine used internally by hap.py |

### Installing samtools and bwa

Most package managers carry recent versions:

```bash
# conda
conda install -c bioconda samtools bwa

# apt (Debian/Ubuntu)
sudo apt install samtools bwa
```

### Installing hap.py and RTG Tools

```bash
conda install -c bioconda hap.py rtg-tools
```

### Installing muteditor (varben)

Follow the instructions in the [muteditor repository](https://github.com/bioinformatics-centre/ART).

### Checking your environment

When generating in silico controls, sherloxome verifies at startup that `bwa`, `samtools`, and `muteditor` are present and prints their paths at `debug` log level:

```bash
RUST_LOG=debug sherloxome setup -c config.toml
```
