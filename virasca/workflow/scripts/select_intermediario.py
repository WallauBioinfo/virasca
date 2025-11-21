import pandas as pd
import argparse

def main(blast_classified, metadata_file, output_file, tax_level="species"):
    """
    Select references for "intermediario" genomes.
    
    For each taxonomic group with intermediario status:
    1. Select the reference (sseqid) with highest bitscore
    2. Output: reference_accession (one per line)
    
    Args:
        blast_classified: Path to classified BLAST results TSV
        metadata_file: Path to database metadata TSV with taxonomy
        output_file: Path to output file for selected references
        tax_level: Taxonomic level to group by (default: species)
    """
    # Valid taxonomic levels
    valid_levels = ["phylum", "class", "order", "family", "genus", "species"]
    if tax_level not in valid_levels:
        raise ValueError(f"Invalid tax_level: {tax_level}. Must be one of {valid_levels}")
    
    # Column name in metadata for this taxonomic level
    tax_col = f"{tax_level}_name"
    
    # Read BLAST classified results
    blast_df = pd.read_csv(blast_classified, sep="\t")
    
    # Read metadata
    metadata_df = pd.read_csv(metadata_file, sep="\t")
    
    # Merge BLAST results with metadata to get taxonomy
    merged_df = blast_df.merge(
        metadata_df[["accession", tax_col]], 
        left_on="sseqid", 
        right_on="accession",
        how="left"
    )
    
    # Filter for intermediario hits only
    intermediario_hits = merged_df[merged_df["status"] == "intermediario"]
    
    if intermediario_hits.empty:
        # No intermediario hits, write empty file
        with open(output_file, "w") as f:
            f.write("")
        return
    
    # For each taxonomic group, select reference with highest bitscore
    results = []
    
    for tax_group in intermediario_hits[tax_col].unique():
        # Skip NA or empty taxonomic groups
        if pd.isna(tax_group) or tax_group == "na":
            continue
            
        tax_group_hits = intermediario_hits[intermediario_hits[tax_col] == tax_group]
        
        # Select the reference (sseqid) with highest bitscore in this group
        best_ref_row = tax_group_hits.loc[tax_group_hits["bitscore"].idxmax()]
        best_ref = best_ref_row["sseqid"]
        
        results.append(best_ref)
    
    if not results:
        # No valid references found
        with open(output_file, "w") as f:
            f.write("")
        return
    
    # Remove duplicates and write to output
    unique_refs = list(set(results))
    
    with open(output_file, "w") as f:
        f.write("\n".join(unique_refs))
        f.write("\n")  # Ensure trailing newline for shell scripts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select references for intermediario genomes")
    parser.add_argument("blast_classified", help="Classified BLAST results TSV")
    parser.add_argument("metadata_file", help="Database metadata TSV with taxonomy")
    parser.add_argument("output_file", help="Output file for selected references")
    parser.add_argument("--tax-level", default="species", 
                        choices=["phylum", "class", "order", "family", "genus", "species"],
                        help="Taxonomic level to group by (default: species)")
    
    args = parser.parse_args()
    
    main(args.blast_classified, args.metadata_file, args.output_file, args.tax_level)
