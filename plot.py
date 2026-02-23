from pygenomeviz import GenomeViz
import pandas as pd
from typing import List, Dict, Tuple

AMR_COLOR = "#33ffff"
REP_COLOR = "#d81159"
DEF_COLOR = "#ff8800"

def plot_contigs(all_contigs_df: pd.DataFrame, contig_ids: List[int], blast_results: Dict[Tuple[int, int], pd.DataFrame], min_identity: float, output_file: str):
    """
    Plot the contigs using pygenomeviz.
    Args:
        all_contigs_df (pd.DataFrame): A DataFrame containing the merged data for all contigs.
        contig_ids (List[int]): A list of contig IDs to plot.
        blast_results (Dict[Tuple[int, int], pd.DataFrame]): A dictionary containing BLAST results for each contig pair.
        min_identity (float): Minimum percent identity to consider for plotting BLAST links.
        output_file (str): The path to the output file.
    Returns:
        None
    """
    gv = GenomeViz(track_align_type = "center")
    track_dict = {}
    # ---------- Create tracks for each contig ----------
    for contig_id in contig_ids:
        contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
        genome_size = int(contig_df[['start', 'end']].max().max())
        name = f"Contig {contig_id}"
        track = gv.add_feature_track(name, genome_size)
        track_dict[contig_id] = track

    # ---------- Add features to tracks ----------
    for contig_id, track in track_dict.items():
        contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
        for _, row in contig_df.iterrows():
            is_amr = pd.notna(row.get('ARO'))
            is_replicon = pd.notna(row.get('replicon'))
            if is_replicon:
                color = REP_COLOR
            elif is_amr:
                color = AMR_COLOR
            else:
                color = DEF_COLOR
            strand = 1 if row['strand'] == '+' else -1
            if row['start'] > row['end']:  # wrap-around for circular genome
                track.add_feature(
                    start=row['start'], end=int(contig_df[['start', 'end']].max().max()),
                    strand=strand,
                    plotstyle="bigarrow",
                    label=row['replicon'] if is_replicon else (row['ARO Name'] if is_amr else None),
                    text_kws={"size": 5},
                    facecolor=color 
                    )
                track.add_feature(
                    start=0, end=row['end'],
                    strand=strand,
                    plotstyle="bigarrow",
                    label=row['replicon'] if is_replicon else (row['ARO Name'] if is_amr else None),
                    text_kws={"size": 5},
                    facecolor=color
                    )
            else:
                track.add_feature(
                    start=row['start'], end=row['end'],
                    strand=strand,
                    plotstyle="bigarrow",
                    label=row['replicon'] if is_replicon else (row['ARO Name'] if is_amr else None),
                    text_kws={"size": 5},
                    facecolor=color
                    )

    # ---------- Add links for BLAST hits ----------
    if blast_results:
        for (query, subject), blast_df in blast_results.items():
            filtered_blast_df = blast_df[blast_df['pident'] >= min_identity]  # Filter for high identity hits
            for _, hit in filtered_blast_df.iterrows():
                try:
                    q_row = all_contigs_df[
                        (all_contigs_df['contig_ID'] == query) & 
                        (all_contigs_df['id'] == hit['qseqid'])
                    ].iloc[0]
                    s_row = all_contigs_df[
                        (all_contigs_df['contig_ID'] == subject) & 
                        (all_contigs_df['id'] == hit['sseqid'])
                    ].iloc[0]
                except IndexError:
                    continue
                gv.add_link(
                    target1=(f"Contig {query}", q_row['start'], q_row['end']),
                    target2=(f"Contig {subject}", s_row['start'], s_row['end']),
                    color="grey",
                    v=hit['pident'],
                    vmin=filtered_blast_df['pident'].min(),
                    vmax=100,
                    curve=True
                )
                gv.set_colorbar(
                    colors=["grey", "grey"],
                    vmin=filtered_blast_df['pident'].min(),
                    vmax=100,
                )
    gv.savefig(output_file)