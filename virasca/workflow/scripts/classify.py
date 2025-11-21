import pandas as pd
import sys

def get_scov(row):
    slen = row["slen"]
    length = row["length"]
    return (length / slen) * 100

def classify(row, complet_threshold=98, intermediate_threshold=90, incomplete_threshold=60):
    scov = row["scov"]
    pident = row["pident"]

    if scov >= complet_threshold:
        return "completo"
    elif pident >= intermediate_threshold:
        return "intermediario"
    elif pident >= incomplete_threshold:
        return "incompleto"
    else:
        return "desconhecido"

def main(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", names=["qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "pident", "qcovs", "qcovhsp", "length"])
    df["scov"] = df.apply(get_scov, axis=1)
    df["status"] = df.apply(classify, axis=1)

    df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python classify.py <input_tsv> <output_tsv>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])
