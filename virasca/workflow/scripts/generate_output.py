import pandas as pd
import argparse
import os

def main(blast_classified, ref_mapping, intermediario_mapping, metadata, output_dir, output_tsv, output_log):
    """
    Generate final output TSV combining BLAST results, taxonomy, and reference information,
    and produce a detailed tracking log (output_log) describing the fate of each sequence.
    """
    # Read BLAST classified results
    blast_df = pd.read_csv(blast_classified, sep="\t")
    
    # Read metadata
    metadata_df = pd.read_csv(metadata, sep="\t")
    
    # Process reference_mapping for ragtag (incompleto)
    ragtag_ref_map = {}
    ragtag_line_map = {} 
    ragtag_contig_counts = {}
    
    if os.path.exists(ref_mapping) and os.path.getsize(ref_mapping) > 0:
        line_num = 1
        with open(ref_mapping, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        ref_acc, contigs = parts
                        contig_list = contigs.split(',')
                        ragtag_line_map[ref_acc] = line_num
                        ragtag_contig_counts[ref_acc] = len(contig_list)
                        for contig in contig_list:
                            ragtag_ref_map[contig] = ref_acc
                line_num += 1

    # Process intermediario_mapping for read_assembly
    inter_line_map = {}
    if os.path.exists(intermediario_mapping) and os.path.getsize(intermediario_mapping) > 0:
        line_num = 1
        with open(intermediario_mapping, 'r') as f:
            for line in f:
                ref_acc = line.strip()
                if ref_acc:
                    inter_line_map[ref_acc] = line_num
                line_num += 1

    # Merge BLAST results with metadata to get taxonomy
    merged_df = blast_df.merge(
        metadata_df[["accession", "species_name", "genus_name", "family_name", "order_name", "class_name", "phylum_name"]], 
        left_on="sseqid", 
        right_on="accession",
        how="left"
    )
    
    # Add reference used in ragtag to TSV
    merged_df["reference_used"] = merged_df["qseqid"].map(ragtag_ref_map)
    
    output_columns = [
        "qseqid", "sseqid", "status", "pident", "qcovs", "bitscore", 
        "species_name", "genus_name", "family_name", "order_name", 
        "class_name", "phylum_name", "reference_used"
    ]
    
    final_df = merged_df[output_columns]
    final_df.to_csv(output_tsv, sep="\t", index=False)
    
    # Prepare best reference map for intermediario
    inter_best_refs = {}
    inter_hits = final_df[final_df["status"] == "intermediario"]
    if not inter_hits.empty:
        for species in inter_hits["species_name"].unique():
            if pd.isna(species) or species == 'na': 
                continue
            group = inter_hits[inter_hits["species_name"] == species]
            best_ref = group.loc[group["bitscore"].idxmax()]["sseqid"]
            inter_best_refs[species] = best_ref

    # Prepare best reference map for incompleto singletons
    incomp_singleton_refs = {}
    incomp_hits = final_df[final_df["status"] == "incompleto"]
    if not incomp_hits.empty:
        for species in incomp_hits["species_name"].unique():
            if pd.isna(species) or species == 'na': 
                continue
            group = incomp_hits[incomp_hits["species_name"] == species]
            if len(group["qseqid"].unique()) < 2:
                best_ref = group.loc[group["bitscore"].idxmax()]["sseqid"]
                incomp_singleton_refs[species] = best_ref
            
    # Generate Output Log tracking what happened to each sequence
    fate_groups = {}
    consensus_table = {}
    
    for _, row in final_df.iterrows():
        qseqid = row["qseqid"]
        status = row["status"]
        sseqid = row["sseqid"]
        
        consensus_path = None
        
        if status == "completo":
            fate_msg = f"{{seqs}} considered COMPLETE against {sseqid}. No further action was needed."
            
        elif status == "intermediario":
            species = row["species_name"]
            best_ref = inter_best_refs.get(species)
            if best_ref and best_ref in inter_line_map:
                line_n = inter_line_map[best_ref]
                fate_msg = f"{{seqs}} considered INTERMEDIATE. The respective taxonomic group triggered Reference-Guided Assembly using reference '{best_ref}', generating outputs in '{output_dir}/read_assembly_{line_n}_{best_ref}'."
                consensus_path = f"{output_dir}/read_assembly_{line_n}_{best_ref}/consensus.fa"
            else:
                fate_msg = f"{{seqs}} considered INTERMEDIATE, but lacked sufficient taxonomic resolution for Reference-Guided Assembly."
                
        elif status == "incompleto":
            ref_used = row["reference_used"]
            if pd.notna(ref_used) and ref_used in ragtag_line_map:
                line_n = ragtag_line_map[ref_used]
                count = ragtag_contig_counts[ref_used]
                fate_msg = f"{{seqs}} considered INCOMPLETE. Submitted to RagTag against reference '{ref_used}' (with {count} total contigs), composing a scaffold output in '{output_dir}/ragtag_{line_n}_{ref_used}'."
                consensus_path = f"{output_dir}/ragtag_{line_n}_{ref_used}/ragtag_output/ragtag.scaffold.fasta"
            else:
                species = row["species_name"]
                best_ref = incomp_singleton_refs.get(species)
                if best_ref and best_ref in inter_line_map:
                    line_n = inter_line_map[best_ref]
                    fate_msg = f"{{seqs}} considered INCOMPLETE against reference '{sseqid}'. However, since it was the ONLY contig assigned to its taxonomic group, it was redirected to Reference-Guided Assembly using reference '{best_ref}', generating outputs in '{output_dir}/read_assembly_{line_n}_{best_ref}'."
                    consensus_path = f"{output_dir}/read_assembly_{line_n}_{best_ref}/consensus.fa"
                else:
                    fate_msg = f"{{seqs}} considered INCOMPLETE against '{sseqid}', but no valid redirection mapping was found."
        else:
            fate_msg = f"{{seqs}} status '{status}'. No specific refinement steps taken."
            
        if fate_msg not in fate_groups:
            fate_groups[fate_msg] = []
        fate_groups[fate_msg].append(qseqid)
        
        if consensus_path:
            if consensus_path not in consensus_table:
                consensus_table[consensus_path] = []
            consensus_table[consensus_path].append(qseqid)
            
    log_lines = [
        "# Virasca Workflow Sequence Tracking Detail Log",
        "",
        "## Sequences Processed",
        ""
    ]
    
    for fate_msg, seqs in fate_groups.items():
        if len(seqs) == 1:
            seq_str = f"Sequence `{seqs[0]}` was"
        elif len(seqs) == 2:
            seq_str = f"Sequences `{seqs[0]}` and `{seqs[1]}` were"
        else:
            quoted_seqs = [f"`{s}`" for s in seqs]
            seq_str = f"Sequences {', '.join(quoted_seqs[:-1])} and {quoted_seqs[-1]} were"
            
        log_lines.append("- " + fate_msg.format(seqs=seq_str))
        
    if consensus_table:
        log_lines.extend([
            "",
            "## Consensus Mapping",
            "",
            "| Consensus Path | Sequences in Consensus |",
            "|---|---|"
        ])
        
        for c_path, seqs in consensus_table.items():
            quoted_seqs = [f"`{s}`" for s in seqs]
            seq_str = ", ".join(quoted_seqs)
            log_lines.append(f"| `{c_path}` | {seq_str} |")
            
    with open(output_log, 'w') as f:
        f.write("\n".join(log_lines) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final output summary and detailed trace log")
    parser.add_argument("--blast", required=True, help="Classified BLAST results TSV")
    parser.add_argument("--ragtag-mapping", required=True, help="Reference-to-contigs mapping file")
    parser.add_argument("--inter-mapping", required=True, help="Intermediario References mapping file")
    parser.add_argument("--metadata", required=True, help="Database metadata TSV")
    parser.add_argument("--output-dir", required=True, help="Output Directory")
    parser.add_argument("--output-tsv", required=True, help="Output TSV file")
    parser.add_argument("--output-log", required=True, help="Detailed tracking log file")
    
    args = parser.parse_args()
    main(args.blast, args.ragtag_mapping, args.inter_mapping, args.metadata, args.output_dir, args.output_tsv, args.output_log)
