import pandas as pd
import argparse

def main(blast_file, metadata_file, output_file, tax_level="species"):
    """
    Select the best reference for each query based on taxonomic grouping.
    
    For each taxonomic group:
    1. Select the reference with highest bitscore
    2. Collect all contigs/scaffolds that map to this taxonomic group
    3. Output: reference_accession <tab> contig1,contig2,contig3
    
    Args:
        blast_file: Path to classified BLAST results TSV
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
    
    # Read BLAST results
    blast_df = pd.read_csv(blast_file, sep="\t")
    
    # Read metadata
    metadata_df = pd.read_csv(metadata_file, sep="\t")
    
    # Merge BLAST results with metadata to get taxonomy information
    # sseqid in BLAST corresponds to accession in metadata
    merged_df = blast_df.merge(
        metadata_df[["accession", tax_col]], 
        left_on="sseqid", 
        right_on="accession",
        how="left"
    )
    
    # Filter for good hits (incompleto as per user request)
    good_hits = merged_df[merged_df["status"].isin(["incompleto"])]
    
    if good_hits.empty:
        # No good hits, write empty file
        with open(output_file, "w") as f:
            f.write("")
        return
    
    # For each taxonomic group, select reference with highest bitscore and collect contigs
    results = []
    
    for tax_group in good_hits[tax_col].unique():
        # Skip NA or empty taxonomic groups
        if pd.isna(tax_group) or tax_group == "na":
            continue
            
        tax_group_hits = good_hits[good_hits[tax_col] == tax_group]
        
        # Select the reference (sseqid) with highest bitscore in this group
        best_ref_row = tax_group_hits.loc[tax_group_hits["bitscore"].idxmax()]
        best_ref = best_ref_row["sseqid"]
        
        # Collect all contigs (qseqid) that belong to this taxonomic group
        contigs = tax_group_hits["qseqid"].unique().tolist()
        
        if len(contigs) >= 2:
            contigs_str = ",".join(contigs)
            results.append(f"{best_ref}\t{contigs_str}")
    
    if not results:
        # No valid references found
        with open(output_file, "w") as f:
            f.write("")
        return
    
    # Write to output
    with open(output_file, "w") as f:
        f.write("\n".join(results))
        f.write("\n")  # Ensure trailing newline for shell scripts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select best reference genome based on taxonomy")
    parser.add_argument("--blast", required=True, help="Classified BLAST results TSV")
    parser.add_argument("--metadata", required=True, help="Database metadata TSV with taxonomy")
    parser.add_argument("--output", required=True, help="Output file for selected references")
    parser.add_argument("--tax-level", default="species", 
                        choices=["phylum", "class", "order", "family", "genus", "species"],
                        help="Taxonomic level to group by (default: species)")
    
    args = parser.parse_args()
    
    main(args.blast, args.metadata, args.output, args.tax_level)
