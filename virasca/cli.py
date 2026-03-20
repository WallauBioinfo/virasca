import click
import snakemake
import os
from pathlib import Path

def get_snakefile(file_name):
    return Path(__file__).parent / "workflow" / file_name

def run_snakemake(target, config=None, cores=1, dryrun=False):
    snakefile = get_snakefile("Snakefile")
    
    # Snakemake 7.x API
    success = snakemake.snakemake(
        snakefile=snakefile,
        targets=[target] if target else None,
        config=config,
        cores=cores,
        dryrun=dryrun,
        printshellcmds=True
    )
    return success

@click.group()
def cli():
    pass

@cli.command()
@click.option("--taxon-id", default=10239, help="NCBI taxon ID for virus database download (default: 10239 for all viruses)")
@click.option("--cores", default=4, help="Number of cores to use (default: 4)")
@click.option("--dry-run", is_flag=True, help="Dry run without executing")
def configure_database(taxon_id, cores, dry_run):
    """Configure the database for virasca."""
    click.echo("Configuring database...")
    
    config = {
        "datasets": {
            "taxon_id": taxon_id
        },
        "threads": cores
    }
    
    run_snakemake("configure_database_all", config=config, cores=cores, dryrun=dry_run)

@cli.command()
@click.option("--input", required=True, type=click.Path(exists=True), help="Path to contigs FASTA file")
@click.option("--reads-r1", required=True, type=click.Path(exists=True), help="Path to R1 FASTQ file")
@click.option("--reads-r2", required=True, type=click.Path(exists=True), help="Path to R2 FASTQ file")
@click.option("--database-seq", required=True, type=click.Path(exists=True), help="Path to blast database FASTA")
@click.option("--database-metadata", required=True, type=click.Path(exists=True), help="Path to blast database metadata file")
@click.option("--output", required=True, type=click.Path(), help="Output directory for analysis results")
@click.option("--use-virseqimprover", is_flag=True, help="Run Virseqimprover step")
@click.option("--tax-level", default="species", type=click.Choice(["phylum", "class", "order", "family", "genus", "species"]), 
              help="Taxonomic level for reference selection (default: species)")
@click.option("--blast-task", default="blastn", help="BLAST task to use (default: blastn)", type=click.Choice(["blastn", "megablast", "dc-megablast", "blastn-short", "rmblastn"]))
@click.option("--blast-word-size", default=11, help="BLAST word size (default: 11)")
@click.option("--fastp-threads", default=4, help="Number of threads for fastp (default: 4)")
@click.option("--bwa-threads", default=4, help="Number of threads for BWA (default: 4)")
@click.option("--min-len", default=50, help="Minimum read length for fastp (default: 50)")
@click.option("--trim-len", default=0, help="Trim length for fastp (default: 0)")
@click.option("--mapping-quality", default=20, help="Minimum mapping quality for iVar consensus (default: 20)")
@click.option("--cores", default=4, help="Number of cores for Snakemake (default: 4)")
@click.option("--ragtag-threads", default=1, help="Number of threads for RagTag (default: 1)")
@click.option("--ragtag-mm2-preset", default="asm5", type=click.Choice(["asm5", "asm10", "asm20"]), help="Minimap2 preset for RagTag (default: asm5)")
@click.option("--ragtag-min-unique-len", default=1000, help="Minimum unique alignment length (default: 1000)")
@click.option("--ragtag-min-mapq", default=10, help="Minimum MAPQ for alignments (default: 10)")
@click.option("--ragtag-infer-gaps", is_flag=True, help="Infer gap sizes (-r flag)")
@click.option("--ragtag-remove-small", is_flag=True,help="Remove unique alignments shorter than --ragtag-min-unique-len")
@click.option("--dry-run", is_flag=True, help="Dry run without executing")
def run(input, reads_r1, reads_r2, database_seq, database_metadata, output,  blast_task, blast_word_size, use_virseqimprover, 
        tax_level, fastp_threads, bwa_threads, min_len, trim_len, mapping_quality, cores, ragtag_threads, ragtag_mm2_preset,
        ragtag_min_unique_len, ragtag_min_mapq, ragtag_infer_gaps, ragtag_remove_small, dry_run):
    """Run the virasca analysis pipeline."""
    click.echo("Running analysis...")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.abspath(output)
    os.makedirs(output_dir, exist_ok=True)
    
    config = {
        "input_fasta": os.path.abspath(input),
        "reads_R1": os.path.abspath(reads_r1),
        "reads_R2": os.path.abspath(reads_r2),
        "database_seq": os.path.abspath(database_seq),
        "database_metadata": os.path.abspath(database_metadata),
        "output_dir": output_dir,
        "blast_task": blast_task,
        "blast_word_size": blast_word_size,
        "use_virseqimprover": use_virseqimprover,
        "tax_level": tax_level,
        "threads": cores,
        "params": {
            "fastp_threads": fastp_threads,
            "bwa_threads": bwa_threads,
            "minLen": min_len,
            "trimLen": trim_len,
            "mapping_quality": mapping_quality,
            "ragtag_threads": ragtag_threads,
            "ragtag_mm2_preset": ragtag_mm2_preset,
            "ragtag_min_unique_len": ragtag_min_unique_len,
            "ragtag_min_mapq": ragtag_min_mapq,
            "ragtag_infer_gaps": ragtag_infer_gaps,
            "ragtag_remove_small": ragtag_remove_small
        }
    }
    
    run_snakemake("run_all", config=config, cores=cores, dryrun=dry_run)

if __name__ == "__main__":
    cli()
