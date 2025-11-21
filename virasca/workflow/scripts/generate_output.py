import pandas as pd
import argparse

def main(blast_classified, ref_mapping, metadata, ragtag_done, output_file):
    """
    Generate final output TSV combining BLAST results, taxonomy, and reference information.
    
    Args:
        blast_classified: Path to classified BLAST results TSV
        ref_mapping: Path to reference-to-contigs mapping file
        metadata: Path to database metadata TSV with taxonomy
        ragtag_done: Path to ragtag completion marker (not used, just for dependency)
        output_file: Path to output TSV file
    """
    # Read BLAST classified results
    blast_df = pd.read_csv(blast_classified, sep="\t")
    
    # Read metadata
    metadata_df = pd.read_csv(metadata, sep="\t")
    
    # Read reference mapping (reference<tab>contigs)
    ref_map = {}
    with open(ref_mapping, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    ref_acc, contigs = parts
                    contig_list = contigs.split(',')
                    for contig in contig_list:
                        ref_map[contig] = ref_acc
    
    # Merge BLAST results with metadata to get taxonomy
    merged_df = blast_df.merge(
        metadata_df[["accession", "species_name", "genus_name", "family_name", "order_name", "class_name", "phylum_name"]], 
        left_on="sseqid", 
        right_on="accession",
        how="left"
    )
    
    # Add reference used in ragtag
    merged_df["reference_used"] = merged_df["qseqid"].map(ref_map)
    
    # Select and reorder columns for final output
    output_columns = [
        "qseqid",           # Input sequence name
        "sseqid",           # BLAST subject (hit)
        "status",           # Classification status
        "pident",           # Percent identity
        "qcovs",            # Query coverage
        "bitscore",         # Bitscore
        "species_name",     # Taxonomy
        "genus_name",
        "family_name",
        "order_name",
        "class_name",
        "phylum_name",
        "reference_used"    # Reference used in ragtag
    ]
    
    # Create final output
    final_df = merged_df[output_columns]
    
    # Write to output
    final_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final output summary")
    parser.add_argument("blast_classified", help="Classified BLAST results TSV")
    parser.add_argument("ref_mapping", help="Reference-to-contigs mapping file")
    parser.add_argument("metadata", help="Database metadata TSV")
    parser.add_argument("ragtag_done", help="Ragtag completion marker")
    parser.add_argument("output_file", help="Output TSV file")
    
    args = parser.parse_args()
    
    main(args.blast_classified, args.ref_mapping, args.metadata, args.ragtag_done, args.output_file)
