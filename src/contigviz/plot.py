from typing import Dict, Iterator, List, Tuple, Optional
from pathlib import Path
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
from pygenomeviz import GenomeViz

from contigviz.models import Contig, FamilyID, Feature, PartialContig, PseudoID

DEF_COLOR = "#ff8800"
WRAP_COLOR = "#ff0000"
LINK_COLOR = "grey"


def pygenomeviz_range(start: int, end: int, length: int) -> Iterator[Tuple[int, int]]:
    """
    Turn the 1-based, inclusive coordinates of a feature into the 0-based, half-open ones of pygenomeviz, and split a
    feature that wraps around the end of a circular contig into the parts before and after the origin.
    Args:
        start (int): The start of the feature.
        end (int): The end of the feature, smaller than start if the feature wraps around.
        length (int): The length of the contig.
    Yields:
        Tuple[int, int]: The 0-based start and end of each part: one part, or two if the feature wraps around.
    """
    if start <= end:
        yield start - 1, end
    else:
        yield start - 1, length
        yield 0, end


def find_links(contig1: Contig, contig2: Contig) -> List[Tuple[Feature, Feature]]:
    """
    Find the pairs of features in the same protein family on two contigs. A family with several copies links every copy
    on the first contig to every copy on the second.
    Args:
        contig1 (Contig): The first contig.
        contig2 (Contig): The second contig.
    Returns:
        List[Tuple[Feature, Feature]]: Each pair of a feature on contig1 and a feature on contig2 in the same protein family.
    """
    family_features: Dict[FamilyID, List[Feature]] = {}
    for feature in contig2.features:
        if feature.protein_family is not None:
            family_features.setdefault(feature.protein_family, []).append(feature)
    return [
        (feature1, feature2)
        for feature1 in contig1.features if feature1.protein_family is not None
        for feature2 in family_features.get(feature1.protein_family, [])
    ]


def add_legend(fig: Figure, contigs: List[Contig]) -> None:
    """
    Add a legend of the categories of the contigs and their features to a figure, in the order they first appear from
    the top contig down. Categories with the same label but another color get their own entry.
    Args:
        fig (Figure): The figure to add the legend to.
        contigs (List[Contig]): The contigs in the figure, in order from top to bottom.
    """
    # dict.fromkeys drops duplicates like a set, but keeps the order in which the categories first appear.
    used = list(dict.fromkeys(
        category
        for contig in contigs
        for category in (contig.category, *(feature.category for feature in contig.features))
        if category is not None
    ))
    if used:
        handles = [mpatches.Patch(color=category.color, label=category.label) for category in used]
        fig.legend(handles=handles, loc="upper right", title="Category", fontsize=10, title_fontsize=12, handlelength=1.5, handleheight=1.5)


def plot_contigs(contigs: List[Contig], output_path: Optional[Path] = None, default_color: str = DEF_COLOR) -> Figure:
    """
    Plot the contigs from top to bottom with pygenomeviz, coloring features and contig labels by category and linking
    the features of neighboring contigs that are in the same protein family.
    Args:
        contigs (List[Contig]): The contigs to plot, in order from top to bottom.
        output_path (Optional[Path]): The path to save the figure to. If None, the figure is not saved.
        default_color (str): The color of features without a category.
    Returns:
        Figure: The matplotlib Figure object containing the plotted contigs.
    """
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar(ymargin=0.5)

    # ---------- Create tracks for each contig ----------
    track_names = {}
    track_dict = {}
    for contig in contigs:
        if isinstance(contig, PartialContig):
            fallback = f"Contig {contig.id.source_id} ({contig.start}-{contig.end})"
        else:
            fallback = f"Contig {contig.id.source_id}"
        track_names[contig.id] = contig.label if contig.label is not None else fallback
        label_kws = {"color": contig.category.color} if contig.category is not None else None
        track_dict[contig.id] = gv.add_feature_track(track_names[contig.id], (contig.start - 1, contig.end), label_kws=label_kws)

    # ---------- Add features to tracks ----------
    for contig in contigs:
        track = track_dict[contig.id]
        for feature in contig.features:
            color = feature.category.color if feature.category is not None else default_color
            style = {"facecolor": color, "alpha": 1.0, "edgecolor": None, "lw": 0.0}
            if isinstance(feature.id, PseudoID):
                style.update(alpha=0.2, edgecolor=color, lw=0.5)
            if feature.start > feature.end:  # wraps around the end of a circular contig
                style.update(edgecolor=WRAP_COLOR, lw=0.5)
            for start, end in pygenomeviz_range(feature.start, feature.end, contig.length):
                track.add_feature(start, end, 1 if feature.strand == "+" else -1, plotstyle="bigarrow", **style)

    # ---------- Add links between protein families of neighboring contigs ----------
    for contig1, contig2 in zip(contigs, contigs[1:]):
        for feature1, feature2 in find_links(contig1, contig2):
            for start1, end1 in pygenomeviz_range(feature1.start, feature1.end, contig1.length):
                for start2, end2 in pygenomeviz_range(feature2.start, feature2.end, contig2.length):
                    gv.add_link(
                        (track_names[contig1.id], start1, end1),
                        (track_names[contig2.id], start2, end2),
                        color=LINK_COLOR,
                        curve=True,
                    )
    # ---------- Plot the figure and add the legend ----------
    fig = gv.plotfig()
    add_legend(fig, contigs)

    # ---------- Save the figure ----------
    if output_path is not None:
        fig.savefig(output_path)
    return fig
