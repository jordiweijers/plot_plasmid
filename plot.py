from typing import Dict, Iterator, List, Tuple

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from pygenomeviz import GenomeViz
from pygenomeviz.track import FeatureTrack

from plot_plasmid.models import Category, Contig, Feature

DEF_COLOR = "#ff8800"
WRAP_COLOR = "#ff0000"
LINK_COLOR = "grey"


def wrap_segments(start: int, end: int, length: int) -> Iterator[Tuple[int, int]]:
    """
    Split a feature that wraps around the end of a circular contig into the parts before and after the origin.
    Args:
        start (int): The start of the feature.
        end (int): The end of the feature, smaller than start if the feature wraps around.
        length (int): The length of the contig.
    Yields:
        Tuple[int, int]: The start and end of each part: one part, or two if the feature wraps around.
    """
    if start <= end:
        yield start, end
    else:
        yield start, length
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
    family_features: Dict[str, List[Feature]] = {}
    for feature in contig2.features:
        if feature.protein_family is not None:
            family_features.setdefault(feature.protein_family, []).append(feature)
    return [
        (feature1, feature2)
        for feature1 in contig1.features if feature1.protein_family is not None
        for feature2 in family_features.get(feature1.protein_family, [])
    ]


def plot_contigs(
        contigs: List[Contig],
        categories: List[Category],
        output_file: str,
    ) -> Tuple[GenomeViz, Dict[str, FeatureTrack]]:
    """
    Plot the contigs from top to bottom with pygenomeviz, coloring features and contig labels by category and linking
    the features of neighboring contigs that are in the same protein family.
    Args:
        contigs (List[Contig]): The contigs to plot, in order from top to bottom.
        categories (List[Category]): The categories, whose colors are used for the features and contig labels in them.
            Categories with the same name share the color of the first one. Contig categories that are not in this
            list, such as those from a contig file, get a color from the tab10 colormap.
        output_file (str): The path to the output file.
    Returns:
        Tuple[GenomeViz, Dict[str, FeatureTrack]]: the GenomeViz figure object and a dict mapping
            contig ID to its FeatureTrack -- callers that want to add more to the figure (e.g. an
            extra subtrack) can add it to a track here and call `gv.plotfig()` again themselves to
            re-render including the addition.
    """
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar(ymargin=0.5)

    # ---------- Create a color mapping for the categories ----------
    category_colors: Dict[str, str] = {}
    for category in categories:
        category_colors.setdefault(category.name, category.color)
    other_categories = sorted({contig.category for contig in contigs if contig.category is not None} - set(category_colors))
    cmap = plt.get_cmap("tab10")
    for i, name in enumerate(other_categories):
        category_colors[name] = mcolors.to_hex(cmap(i % cmap.N))

    # ---------- Create tracks for each contig ----------
    track_names = {contig.id: f"Contig {contig.id}" for contig in contigs}
    track_dict = {}
    for contig in contigs:
        label_kws = {"color": category_colors[contig.category]} if contig.category is not None else None
        track_dict[contig.id] = gv.add_feature_track(track_names[contig.id], contig.length, label_kws=label_kws)

    # ---------- Add features to tracks ----------
    for contig in contigs:
        track = track_dict[contig.id]
        for feature in contig.features:
            color = category_colors.get(feature.category, DEF_COLOR)
            style = {"facecolor": color, "alpha": 1.0, "edgecolor": None, "lw": 0.0}
            if feature.type == "pseudo":
                style.update(alpha=0.2, edgecolor=color, lw=0.5)
            if feature.start > feature.end:  # wraps around the end of a circular contig
                style.update(edgecolor=WRAP_COLOR, lw=0.5)
            for start, end in wrap_segments(feature.start, feature.end, contig.length):
                track.add_feature(start, end, 1 if feature.strand == "+" else -1, plotstyle="bigarrow", **style)

    # ---------- Add links between protein families of neighboring contigs ----------
    for contig1, contig2 in zip(contigs, contigs[1:]):
        for feature1, feature2 in find_links(contig1, contig2):
            for start1, end1 in wrap_segments(feature1.start, feature1.end, contig1.length):
                for start2, end2 in wrap_segments(feature2.start, feature2.end, contig2.length):
                    gv.add_link(
                        (track_names[contig1.id], start1, end1),
                        (track_names[contig2.id], start2, end2),
                        color=LINK_COLOR,
                        curve=True,
                    )

    fig = gv.plotfig()
    used = {feature.category for contig in contigs for feature in contig.features} | {contig.category for contig in contigs}
    handles = [mpatches.Patch(color=color, label=name) for name, color in category_colors.items() if name in used]
    if handles:
        fig.legend(handles=handles, loc="upper right", title="Category", fontsize=10, title_fontsize=12, handlelength=1.5, handleheight=1.5)
    fig.savefig(output_file)
    return gv, track_dict
