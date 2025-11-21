# Virasca

Virasca is a CLI tool for viral genome assembly and identification, using Snakemake for workflow management.

## Installation

It is recommended to use `micromamba` (or `conda`/`mamba`) to manage dependencies.

1. **Clone virasca repository and submodules**:
   ```bash
   git clone https://github.com/WallauBioinfo/virasca.git
   cd virasca
   ```

2. **Create the environment**:
   ```bash
   micromamba env create -f environment.yml
   ```

3. **Activate the environment**:
   ```bash
   micromamba activate virasca
   ```

## Usage

### Configure Database

This command downloads and configures the necessary databases (NCBI Virus genomes and taxonomy).

```bash
virasca configure-database
```

Options:
- `--cores`: Number of cores to use (default: 1).
- `--dry-run`: Print the workflow without executing.

### Run Analysis

This command runs the viral assembly and identification pipeline.

```bash
virasca run --input <scaffolds.fasta> --reads-r1 <R1.fastq.gz> --reads-r2 <R2.fastq.gz>
```

Arguments:
- `--input`: Input FASTA file with contigs/scaffolds.
- `--reads-r1`: Path to R1 reads (FASTQ).
- `--reads-r2`: Path to R2 reads (FASTQ).

Options:
- `--cores`: Number of cores to use (default: 1).
- `--dry-run`: Print the workflow without executing.

## Testing

To verify the installation and workflow structure (dry-run):

```bash
# Ensure you are in the virasca environment
micromamba activate virasca

# Run the test script
python3 tests/test_cli.py
```
