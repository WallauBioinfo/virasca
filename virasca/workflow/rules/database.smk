rule configure_database_all:
    input:
        "database/metadata.tsv",
        "database/sequences.fasta",
        expand("{home}/.taxonkit/names.dmp", home=os.environ["HOME"])

rule download_datasets:
    output:
        "database/ncbi_dataset.zip"
    params:
        taxon_id = config["datasets"]["taxon_id"]
    shell:
        """
        datasets download virus genome taxon {params.taxon_id} \
            --refseq \
            --include genome \
            --fast-zip-validation \
            --filename {output}
        """

rule unzip_datasets:
    input:
        "database/ncbi_dataset.zip"
    output:
        "database/ncbi_dataset/data/data_report.jsonl",
        "database/ncbi_dataset/data/genomic.fna"
    shell:
        "unzip -o {input} -d database"

rule process_metadata:
    input:
        "database/ncbi_dataset/data/data_report.jsonl"
    output:
        "database/tmp_metadata.tsv"
    shell:
        """
        echo -e "accession\\tsegment\\tvirus_name\\tvirus_tax_id\\tphylum_name\\tphylum_tax_id\\tclass_name\\tclass_tax_id\\torder_name\\torder_tax_id\\tfamily_name\\tfamily_tax_id\\tgenus_name\\tgenus_tax_id\\tspecies_name\\tspecies_tax_id" > {output}
        dataformat tsv virus-genome \\
            --inputfile {input} \\
            --fields accession,segment,virus-name,virus-tax-id | \\
            grep -v "Accession" >> {output}
        """

rule setup_taxonomy:
    output:
        expand("{home}/.taxonkit/{file}", home=os.environ["HOME"], file=["names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp"])
    params:
        taxonkit_dir = expand("{home}/.taxonkit", home=os.environ["HOME"])
    shell:
        """
        mkdir -p tmp_taxdump
        cd tmp_taxdump
        wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -zxvf taxdump.tar.gz
        mkdir -p {params.taxonkit_dir}
        cp names.dmp nodes.dmp delnodes.dmp merged.dmp {params.taxonkit_dir}
        cd ..
        rm -rf tmp_taxdump
        """

rule map_taxonomy:
    input:
        "database/tmp_metadata.tsv",
        expand("{home}/.taxonkit/names.dmp", home=os.environ["HOME"])
    output:
        "database/taxid_mapping.tsv"
    shell:
        """
        echo -e "virus_tax_id\\tphylum_name\\tphylum_tax_id\\tclass_name\\tclass_tax_id\\torder_name\\torder_tax_id\\tfamily_name\\tfamily_tax_id\\tgenus_name\\tgenus_tax_id\\tspecies_name\\tspecies_tax_id" > {output}
        cut -f 4 {input[0]} | grep -v "virus_tax_id" | taxonkit lineage \
            | taxonkit reformat2 -t -f "{{phylum}};{{class}};{{order}};{{family}};{{genus}};{{species}}" \
            | cut -f 1,3,4 \
            | sed 's/;/\t/g' \
            | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$8,$3,$9,$4,$10,$5,$11,$6,$12,$7,$13}}' >> {output}
        """

rule join_metadata:
    input:
        meta="database/tmp_metadata.tsv",
        tax="database/taxid_mapping.tsv"
    output:
        "database/metadata.tsv"
    run:
        shell("""
        awk -F'\\t' '
        BEGIN {{ FS = OFS = "\\t" }}
        FNR==NR {{
            map[$1]=$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13
            next
        }}
        FNR==1 {{
            print $0
            next
        }}
        {{
            key=$4
            if(key in map){{
                print $1"\\t"$2"\\t"$3"\\t"$4"\\t"map[key]
            }} else {{
                print $1"\\t"$2"\\t"$3"\\t"$4"\\tna\\tna\\tna\\tna\\tna\\tna\\tna\\tna\\tna\\tna\\tna\\tna"
            }}
        }}
        ' {input.tax} {input.meta} > {output}
        """)

rule make_blast_db:
    input:
        sequence = "database/ncbi_dataset/data/genomic.fna",
        metadata = "database/metadata.tsv"
    output:
        "database/sequences.fasta"
    shell:
        """
        cp {input.sequence} {output}
        makeblastdb -in {output} -dbtype nucl
        rm -rf database/ncbi_dataset database/md5sum.txt database/ncbi_dataset.zip database/README.md database/taxid_mapping.tsv database/tmp_metadata.tsv
        """