from pygenomeviz import GenomeViz
import mysql.connector
import pandas as pd
from typing import List

CONTIG_ID = [4339, 2082, 1991]
DB_CONFIG = {
    "host": "localhost",
    "user": "gmg",
    "database": "kes2021"
}
AMR_COLOR = "#33ffff"

MAKEBLASTDB_PARAMS = {
    "dbtype": "prot",
    "num_threads": 1
}

BLASTP_PARAMS = {
    "outfmt": 3,
    "evalue": 1e-9,
    "word_size": 11,
    "min_identity_percentage": 60,
    "max_target_seqs": 1000000,
    "num_threads": 1
}

def get_contig_data(contig_id: int, conn, aro_index_df: pd.DataFrame) -> pd.DataFrame:
    """
    Fetch all proteins for a given contig_ID, their associated AMR information, and merge it with the ARO index.
    Args:
        contig_id (int): The ID of the contig to fetch data for.
        conn: A MySQL database connection object.
        aro_index_df (pd.DataFrame): A DataFrame containing the ARO index.
    Returns:
        pd.DataFrame: A DataFrame containing the merged data for the given contig_ID.
    """
    proteins_df = fetch_proteins(contig_id, conn)
    amr_df = fetch_amr_for_proteins(proteins_df['id'].tolist(), conn)
    if not amr_df.empty:
        amr_df['ARO'] = amr_df['ARO'].astype(int)
        df = proteins_df.merge(amr_df, how='left', on='id')
        df = df.merge(aro_index_df, how='left', left_on='ARO', right_on='ARO Accession')
    else:
        df =proteins_df.copy()
    return df

def load_aro_index(file_path: str) -> pd.DataFrame:
    """
    Load the ARO index from a TSV file and return it as a pandas DataFrame.
    Args:
        file_path (str): The path to the ARO index TSV file.
    Returns:
        pd.DataFrame: A DataFrame containing the ARO index.
    """
    df = pd.read_csv(file_path, sep='\t')
    df['ARO Accession'] = df['ARO Accession'].str.replace("ARO:", "").astype(int)
    return df



def main():
    conn = mysql.connector.connect(**DB_CONFIG)
    aro_index_df = load_aro_index("card-data/aro_index.tsv")
    for contig_id in CONTIG_ID:
        contig_df = get_contig_data(contig_id, conn, aro_index_df)
        all_contigs_df = pd.concat([all_contigs_df, contig_df], ignore_index=True) if 'all_contigs_df' in locals() else contig_df
    conn.close()
    gv = GenomeViz()
    track_dict = {}
    for contig_id in CONTIG_ID:
        contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
        genome_size = int(contig_df[['start', 'end']].max().max())
        name = f"Contig {contig_id}"
        track = gv.add_feature_track(name, genome_size)
        track_dict[contig_id] = track
    for contig_id, track in track_dict.items():
        contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
        for _, row in contig_df.iterrows():
            is_amr = pd.notna(row.get('ARO'))
            if row['start'] > row['end']:
                track.add_feature(
                    start = row['start'],
                    end = genome_size,
                    strand = 1 if row['strand'] == '+' else -1,
                    label = row['CARD Short Name'] if is_amr else "",
                    plotstyle = "bigarrow",
                    text_kws = {"size": 5},
                    facecolor = AMR_COLOR if is_amr else "orange",
                )
                track.add_feature(
                    start = 0,
                    end = row['end'],
                    strand = 1 if row['strand'] == '+' else -1,
                    label = row['CARD Short Name'] if is_amr else "",
                    plotstyle = "bigarrow",
                    text_kws = {"size": 5},
                    facecolor = AMR_COLOR if is_amr else "orange",
                )
            else:
                track.add_feature(
                    start = row['start'],
                    end = row['end'],
                    strand = 1 if row['strand'] == '+' else -1,
                    label = row['CARD Short Name'] if is_amr else "",
                    plotstyle = "bigarrow",
                    text_kws = {"size": 5},
                    facecolor = AMR_COLOR if is_amr else "orange",
                )
    gv.savefig(f"results/contig_{CONTIG_ID}_plot.png", dpi=300)

if __name__ == "__main__":
    main()