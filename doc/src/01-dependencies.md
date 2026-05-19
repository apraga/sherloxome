# Install dependencies

`sherloxome` is a command-line tool that can be installed with `cargo`, the Rust package manager.  Install Rust via [rustup](https://rustup.rs/):
Then compile sherloxome:

```bash
cargo build --release
```

Running it with
```bash
cargo run --release -- setup
```

Or install it with `cargo install` and run `sherloxome` on its own afterwards.

## External tools

We provide an easy setup with nix. `nix profile add .` will install and add to the PATH the following

| Tool | Minimum version | Purpose |
|------|----------------|---------|
| `samtools` | 1.15 | BAM/FASTQ processing |
| `bwa` | 0.7 | Read alignment for in silico BAM editing |
| `muteditor` | any | Variant injection (in silico controls) |
| `hap.py` | 0.3.15 | Variant call comparison against truth sets |
| `rtg` (RTG Tools) | 3.12 | vcfeval engine used internally by hap.py |

Otherwise, install them manually and make them available in your `PATH`. Only `hap.py` and `rtg` are required for the `analyze` step; `bwa`, `samtools`, and `muteditor` are only needed when generating in silico controls.

When generating in silico controls, sherloxome verifies at startup that `bwa`, `samtools`, and `muteditor` are present and prints their paths at `debug` log level:
```bash
RUST_LOG=debug sherloxome setup -c config.toml
```
