import pandas as pd

def parse_blast_output(blast_file: str) -> pd.DataFrame:
    """
    Parse BLAST tabular output (outfmt 6) into a DataFrame.
    """
    cols = [
        "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"
    ]
    df = pd.read_csv(blast_file, sep="\t", names=cols)
    return df