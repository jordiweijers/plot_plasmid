from pygenomeviz import GenomeViz
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from typing import List, Dict, Tuple

AMR_COLOR = "#33ffff"
REP_COLOR = "#5e2bff"
CONJUGATION_COLOR = "#8b0000"
PLASANN_COLOR = "#efa0bd"
DEF_COLOR = "#ff8800"
WRAP_COLOR = "#ff0000"

def plot_contigs(
        all_contigs_df: pd.DataFrame,
        contig_ids: List[int],
        blast_results: Dict[Tuple[int, int], pd.DataFrame],
        plasann_replicons: Dict[int, pd.DataFrame],
        min_identity: float,
        output_file: str,
        contig_categories: Dict[int, str],
    ) -> Tuple[GenomeViz, Dict[int, object]]:
    """
    Plot the contigs using pygenomeviz.
    Args:
        all_contigs_df (pd.DataFrame): A DataFrame containing the merged data for all contigs.
        contig_ids (List[int]): A list of contig IDs to plot.
        blast_results (Dict[Tuple[int, int], pd.DataFrame]): A dictionary containing BLAST results for each contig pair.
        plasann_replicons (Dict[int, pd.DataFrame]): A dictionary containing PlasAnn replicon data for each contig.
        min_identity (float): Minimum percent identity to consider for plotting BLAST links.
        output_file (str): The path to the output file.
        contig_categories (Dict[int, str]): A dictionary mapping contig_id to a category label for track-label coloring. Contigs absent from this dictionary get the default label color.
    Returns:
        Tuple[GenomeViz, Dict[int, FeatureTrack]]: the GenomeViz figure object and a dict mapping
            contig_id to its FeatureTrack -- callers that want to add more to the figure (e.g. an
            extra subtrack) can add it to a track here and call `gv.plotfig()` again themselves to
            re-render including the addition.
    """
    gv = GenomeViz(track_align_type = "center")
    gv.set_scale_bar(ymargin=0.5)
    contig_categories = contig_categories or {}
    # ---------- Create a color mapping for the categories ----------
    unique_categories = sorted(set(contig_categories.values()))
    cmap = plt.get_cmap("tab10")
    category_colors = {category: mcolors.to_hex(cmap(i % cmap.N)) for i, category in enumerate(unique_categories)}
    # ---------- Create tracks for each contig ----------
    track_dict = {}
    for contig_id in contig_ids:
        contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
        cds_max = int(contig_df[['start', 'end']].max().max())
        rep_df = plasann_replicons.get(contig_id)
        if rep_df is not None and not rep_df.empty:
            rep_max = int(rep_df[['start', 'end']].max().max())
            genome_size = max(cds_max, rep_max)
        else:
            genome_size = cds_max
        name = f"Contig {contig_id}"
        category = contig_categories.get(contig_id)
        label_kws = {"color": category_colors[category]} if category is not None else None
        track = gv.add_feature_track(name, genome_size, label_kws=label_kws)
        track_dict[contig_id] = track


    # ---------- Add features to tracks ----------
    for contig_id, track in track_dict.items():
        contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
        for _, row in contig_df.iterrows():
            is_amr = pd.notna(row.get('ARO'))
            is_replicon = pd.notna(row.get('replicon'))
            is_conjugation = pd.notna(row.get('conjugation'))
            if is_replicon:
                color = REP_COLOR
            elif is_conjugation:
                color = CONJUGATION_COLOR
            elif is_amr:
                color = AMR_COLOR
            else:
                color = DEF_COLOR
            alpha = 1.0
            edgecolor = None
            linewidth = 0.0
            if row['type'] == 'pseudo':
                alpha = 0.2
                edgecolor = color
                linewidth = 0.5
            strand = 1 if row['strand'] == '+' else -1
            is_wrap = row['start'] > row['end'] # wrap around for circular contigs
            label = f"{str(row['replicon']).split(' ')[0]}_{int(row['mcl_id'])}" if is_replicon and pd.notna(row['replicon']) and pd.notna(row['mcl_id']) else None
            if is_wrap: 
                track.add_feature(
                    start=row['start'], end=int(contig_df[['start', 'end']].max().max()),
                    strand=strand,
                    plotstyle="bigarrow",
                    label=label,
                    text_kws={"size": 5},
                    facecolor=color,
                    edgecolor=WRAP_COLOR,
                    alpha=alpha,
                    lw=0.5,
                    )
                track.add_feature(
                    start=0, end=row['end'],
                    strand=strand,
                    plotstyle="bigarrow",
                    label=label,
                    text_kws={"size": 5},
                    facecolor=color,
                    edgecolor=WRAP_COLOR,
                    alpha=alpha,
                    lw=0.5,
                    )
            else:
                track.add_feature(
                    start=row['start'], end=row['end'],
                    strand=strand,
                    plotstyle="bigarrow",
                    label=label,
                    text_kws={"size": 5},
                    facecolor=color,
                    edgecolor=edgecolor,
                    alpha=alpha,
                    lw=linewidth,
                    )
    
    # ---------- Add links for BLAST hits ----------
    if blast_results:
        for (query, subject), blast_df in blast_results.items():
            filtered_blast_df = blast_df[blast_df['pident'] >= min_identity]  # Filter for high identity hits
            if filtered_blast_df.empty:
                continue
            vmin = filtered_blast_df['pident'].min()
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
                q_wrap = q_row['start'] > q_row['end']
                s_wrap = s_row['start'] > s_row['end']
                q_size = int(all_contigs_df[all_contigs_df['contig_ID'] == query][['start', 'end']].max().max())
                s_size = int(all_contigs_df[all_contigs_df['contig_ID'] == subject][['start', 'end']].max().max())
                def wrap_segments(row, size):
                    if row['start'] <= row['end']:
                        yield (row['start'], row['end'])
                    else:
                        yield (row['start'], size)
                        yield (0, row['end'])
                for q_seg_start, q_seg_end in wrap_segments(q_row, q_size):
                    for s_seg_start, s_seg_end in wrap_segments(s_row, s_size):
                        gv.add_link(
                            target1=(f"Contig {query}", q_seg_start, q_seg_end),
                            target2=(f"Contig {subject}", s_seg_start, s_seg_end),
                            color="grey",
                            v=hit['pident'],
                            vmin=vmin,
                            vmax=100,
                            curve=True
                        )
            gv.set_colorbar(
                colors=["grey", "grey"],
                vmin=vmin,
                vmax=100,
            )
    fig = gv.plotfig()
    if plasann_replicons:
        for contig_id, track in track_dict.items():
            rep_df = plasann_replicons.get(contig_id)
            if rep_df is None or rep_df.empty:
                continue
            contig_df = all_contigs_df[all_contigs_df['contig_ID'] == contig_id]
            genome_size = int(contig_df[['start', 'end']].max().max())
            for _, row in rep_df.iterrows():
                start = int(row['start'])
                end = int(row['end'])
                label = row['gene_name']
                def plot_region(start, end):
                    ts = track.transform_coord(start)
                    te = track.transform_coord(end)
                    x, y = (ts, te, te, ts), (-1, -1, 1, 1)
                    track.ax.fill(x, y, fc = PLASANN_COLOR, edgecolor='none', zorder=-1)
                    text_x, text_y = (ts + te) / 2, - 1.5
                    track.ax.text(text_x, text_y, s=label, ha='center', va='bottom', size=5, color="black")
                if start <= end:
                    plot_region(start, end)
                else:
                    plot_region(start, genome_size)
                    plot_region(0, end)
    if category_colors:
        handles = [mpatches.Patch(color=color, label=category) for category, color in sorted(category_colors.items())]
        fig.legend(handles=handles, loc="upper right", title="Category", fontsize=10, title_fontsize=12, handlelength=1.5, handleheight=1.5)
    fig.savefig(output_file)
    return gv, track_dict