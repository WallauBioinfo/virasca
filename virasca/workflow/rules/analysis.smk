from virasca import PKG_PATH
OUTPUT_DIR = config.get("output_dir", ".").strip()



rule run_all:
    input:
        f"{OUTPUT_DIR}/output.tsv",
        f"{OUTPUT_DIR}/read_assembly_done.txt"

rule virseqimprove:
    input:
        scaffold = config.get("input_fasta", ""),
        r1 = config.get("reads_R1", ""),
        r2 = config.get("reads_R2", "")
    output:
        f"{OUTPUT_DIR}/scaffold.fasta"
    params:
        out_dir = f"{OUTPUT_DIR}/result_virseq"
    shell:
        """
        mkdir -p {params.out_dir}
        
        python3 {PKG_PATH}/submodules/virseqimprover/virseqimprover/Virseqimprover.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -scaffold {input.scaffold} \
            -o {params.out_dir}

        if [ -f {params.out_dir}/pilon_out.fasta ]; then
            cp {params.out_dir}/pilon_out.fasta {output}
        else
            cp {params.out_dir}/scaffold.fasta {output}
        fi
        
        grep '^>' {output} \
            | sed 's/^>//; s/ .*//; s/:.*//; s/_pilon$//' \
            > {params.out_dir}/treated_ids.txt

        
        awk -v ids={params.out_dir}/treated_ids.txt '
            BEGIN {{
                while ((getline < ids) > 0) seen[$1]=1
            }}
            /^>/ {{
                id=$0
                sub(/^>/,"",id)
                sub(/ .*/,"",id)
                sub(/:.*$/,"",id)
                sub(/_pilon$/,"",id)
                keep = !(id in seen)
            }}
            keep {{ print }}
        ' {input.scaffold} > {params.out_dir}/not_treated.fasta

        cat {output} {params.out_dir}/not_treated.fasta \
            > {params.out_dir}/final_tmp.fasta

        mv {params.out_dir}/final_tmp.fasta {output}
        """

def get_blastn_input(wildcards):
    if config.get("use_virseqimprover", False):
        return f"{OUTPUT_DIR}/scaffold.fasta"
    else:
        return config.get("input_fasta")

rule blastn:
    input:
        query = get_blastn_input,
        db = config.get("database_seq")
    output:
        f"{OUTPUT_DIR}/blastn_results.tsv"
    threads: config.get("threads")
    shell:
        """
        blastn -db {input.db} \
            -query {input.query} \
            -out {output} \
            -task blastn \
            -evalue 0.001 \
            -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp length" \
            -max_hsps 1 \
            -max_target_seqs 1 \
            -num_threads {threads}
        """


rule classify_blast:
    input:
        f"{OUTPUT_DIR}/blastn_results.tsv"
    output:
        f"{OUTPUT_DIR}/blastn_classified.tsv"
    shell:
        "python3 {PKG_PATH}/workflow/scripts/classify.py {input} {output}"

rule select_reference:
    input:
        blast = f"{OUTPUT_DIR}/blastn_classified.tsv",
        metadata = config.get("database_metadata")
    output:
        f"{OUTPUT_DIR}/selected_reference.acc"
    params:
        tax_level = config.get("tax_level", "species")
    shell:
        "python3 {PKG_PATH}/workflow/scripts/select_reference.py {input.blast} {input.metadata} {output} --tax-level {params.tax_level}"

rule select_intermediario:
    input:
        blast = f"{OUTPUT_DIR}/blastn_classified.tsv",
        metadata = config.get("database_metadata", "")
    output:
        f"{OUTPUT_DIR}/selected_intermediario.acc"
    params:
        tax_level = config.get("tax_level", "species")
    shell:
        "python3 {PKG_PATH}/workflow/scripts/select_intermediario.py {input.blast} {input.metadata} {output} --tax-level {params.tax_level}"

rule read_assembly:
    input:
        refs = f"{OUTPUT_DIR}/selected_intermediario.acc",
        database_seq = config.get("database_seq", ""),
        reads_r1 = config.get("reads_R1", ""),
        reads_r2 = config.get("reads_R2", "")
    output:
        f"{OUTPUT_DIR}/read_assembly_done.txt"
    params:
        output_dir = OUTPUT_DIR,
        fastp_threads = config.get("params", {}).get("fastp_threads", 4),
        bwa_threads = config.get("params", {}).get("bwa_threads", 4),
        minLen = config.get("params", {}).get("minLen", 50),
        trimLen = config.get("params", {}).get("trimLen", 0),
        mapping_quality = config.get("params", {}).get("mapping_quality", 20)
    shell:
        """
        # Check if refs file is empty
        if [ ! -s {input.refs} ]; then
            echo "No intermediario references found" > {output}
            exit 0
        fi
        
        # Process each reference
        line_num=0
        while read -r ref_acc || [ -n "$ref_acc" ]; do
            # Skip empty lines
            [ -z "$ref_acc" ] && continue
            
            line_num=$((line_num + 1))
            
            # Create output directory for this reference
            ref_dir="{params.output_dir}/read_assembly_${{line_num}}_${{ref_acc}}"
            mkdir -p "$ref_dir"
            
            # Extract reference genome from database
            ref_file="${{ref_dir}}/reference.fasta"
            echo "$ref_acc" > "${{ref_dir}}/ref_acc.txt"
            seqkit grep -f "${{ref_dir}}/ref_acc.txt" {input.database_seq} > "$ref_file"
            
            # Run fastp
            fastp -i {input.reads_r1} -I {input.reads_r2} \
                --detect_adapter_for_pe \
                --thread {params.fastp_threads} \
                -o "${{ref_dir}}/clean_R1.fq.gz" \
                -O "${{ref_dir}}/clean_R2.fq.gz" \
                -j "${{ref_dir}}/fastp.json" \
                -l {params.minLen} -f {params.trimLen} -t {params.trimLen} \
                -F {params.trimLen} -T {params.trimLen} \
                --cut_front --cut_tail --qualified_quality_phred 20
            
            # Build BWA index
            bwa index "$ref_file"
            
            # Map reads with BWA
            bwa mem \
                -t {params.bwa_threads} \
                "$ref_file" \
                "${{ref_dir}}/clean_R1.fq.gz" \
                "${{ref_dir}}/clean_R2.fq.gz" | \
                samtools sort -o "${{ref_dir}}/mapped.bam" -
            
            # Index BAM
            samtools index "${{ref_dir}}/mapped.bam"
            
            # Generate consensus with iVar
            samtools mpileup \
                -aa \
                -d 50000 \
                --reference "$ref_file" \
                -A \
                -Q 0 \
                "${{ref_dir}}/mapped.bam" | \
                ivar consensus \
                    -p "${{ref_dir}}/consensus" \
                    -q {params.mapping_quality} \
                    -t 0 \
                    -m 10 \
                    -n N \
                    -c 0.51
            
        done < {input.refs}
        
        # Mark as done
        echo "Read assembly completed for all intermediario references" > {output}
        """
def get_ragtag_input_fasta(wildcards):
    if config.get("use_virseqimprover", False):
        return f"{OUTPUT_DIR}/scaffold.fasta"
    else:
        return config.get("input_fasta", "")

rule ragtag:
    input:
        refs_mapping = f"{OUTPUT_DIR}/selected_reference.acc",
        database_seq = config.get("database_seq", ""),
        input_fasta = get_ragtag_input_fasta
    output:
        f"{OUTPUT_DIR}/ragtag_done.txt"
    params:
        output_dir = OUTPUT_DIR, 
        ragtag_threads = config["params"].get("ragtag_threads", 1),
        ragtag_mm2_preset = config["params"].get("ragtag_mm2_preset", "asm5"),
        ragtag_min_unique_len = config["params"].get("ragtag_min_unique_len", 1000),
        ragtag_min_mapq = config["params"].get("ragtag_min_mapq", 10),

        infer_gaps_flag = "-r" if config["params"].get("ragtag_infer_gaps", False) else "",
        remove_small_flag = "--remove-small" if config["params"].get("ragtag_remove_small", False) else ""
    shell:
        """
        # Check if refs_mapping is empty
        if [ ! -s {input.refs_mapping} ]; then
            echo "No references found" > {output}
            exit 0
        fi
        
        # Process each line in the reference mapping file
        # Note: Using while read with || [ -n "$ref_acc" ] to handle files without trailing newline
        line_num=0
        while IFS=$'\t' read -r ref_acc contigs || [ -n "$ref_acc" ]; do
            # Skip empty lines
            [ -z "$ref_acc" ] && continue
            
            line_num=$((line_num + 1))
            
            # Create output directory for this reference
            ref_dir="{params.output_dir}/ragtag_${{line_num}}_${{ref_acc}}"
            mkdir -p "$ref_dir"
            
            # Extract reference genome from database
            ref_file="${{ref_dir}}/reference.fasta"
            echo "$ref_acc" > "${{ref_dir}}/ref_acc.txt"
            seqkit grep -f "${{ref_dir}}/ref_acc.txt" {input.database_seq} > "$ref_file"
            
            # Convert comma-separated contigs to one per line
            contig_list="${{ref_dir}}/contig_list.txt"
            echo "$contigs" | tr ',' '\\n' > "$contig_list"
            
            # Extract contigs from input fasta
            contigs_file="${{ref_dir}}/contigs.fasta"
            seqkit grep -f "$contig_list" {input.input_fasta} > "$contigs_file"
            
            # Run ragtag for this reference
            cd "$ref_dir"
            ragtag.py scaffold "$ref_file" "$contigs_file" -t {params.ragtag_threads} --mm2-params "-x {params.ragtag_mm2_preset}" -f {params.ragtag_min_unique_len} -q {params.ragtag_min_mapq} {params.infer_gaps_flag} {params.remove_small_flag} 
            cd -
            
        done < {input.refs_mapping}
        
        # Mark as done
        echo "Ragtag completed for all references" > {output}
        """

rule generate_output:
    input:
        blast_classified = f"{OUTPUT_DIR}/blastn_classified.tsv",
        ref_mapping = f"{OUTPUT_DIR}/selected_reference.acc",
        metadata = config.get("database_metadata", ""),
        ragtag_done = f"{OUTPUT_DIR}/ragtag_done.txt"
    output:
        f"{OUTPUT_DIR}/output.tsv"
    shell:
        "python3 {PKG_PATH}/workflow/scripts/generate_output.py {input.blast_classified} {input.ref_mapping} {input.metadata} {input.ragtag_done} {output}"
