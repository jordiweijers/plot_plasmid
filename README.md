# contigviz

> [!IMPORTANT]
> This tool is made for members of the Dagan group. It reads its data from the group's MySQL databases,
> which are installed locally on the HPC (`cauclusterg`), so it only works there. The code is public for
> reference, but it will not run elsewhere without your own databases.

Plot contigs from the group's databases with [pyGenomeViz](https://github.com/moshi4/pyGenomeViz). Each contig is
drawn as a track of its proteins and pseudogenes, features of the same MCL protein family on neighboring contigs are
connected by links, and you can color features and contigs by category, for example AMR genes or replication proteins.

<!-- TODO: example image -->

There are two ways to use it:

- **`plot-contigs`**: a command that plots contigs with up to one category. Enough for most plots.
- **The `contigviz` Python API**: the building blocks behind the command, for plots with several categories, cut-out regions or anything else you want to do in your own script.

## Installation

Requires Python 3.11 or newer and an account on the HPC. Install it in a virtual environment:

```bash
python3 -m venv ~/venvs/contigviz
source ~/venvs/contigviz/bin/activate
pip install git+https://github.com/jordiweijers/contigviz.git
plot-contigs -h
```

While the environment is active, `plot-contigs` is on your `PATH` and `import contigviz` works in Python. In a new session, activate it again with `source ~/venvs/contigviz/bin/activate`.

## Quick start

Plot three contigs from kes2021, a hybrid plasmid between its two parents:

```bash
plot-contigs -d kes2021 -c 1991 2082 4338 -o hybrid.svg
```

To color the AMR genes, download the example query [`examples/amr_proteins.sql`](examples/amr_proteins.sql) and give it as a category:

```bash
curl -O https://raw.githubusercontent.com/jordiweijers/contigviz/main/examples/amr_proteins.sql
plot-contigs -d kes2021 -c 1991 2082 4338 -o hybrid_amr.svg \
    --category amr_proteins.sql --id-type protein --label AMR --color "#33ffff" --cache amr_proteins.txt
```

Some SQL queries take a long time. With `--cache amr_proteins.txt`, the IDs are saved in that file, and the next plot reads them from there instead of rerunning the query.

## `plot-contigs`

```bash
plot-contigs -d DATABASE (-c ID [ID ...] | -f FILE) -o OUTPUT [options]
```

| Option | Description |
|---|---|
| `-d`, `--database` | Database to load the contigs from: `kes2021`, `ESKAPEE25`, `pantoea2025_v2` or `acinetobacter_2023` (required). |
| `-c`, `--contig_ids` | Contig IDs to plot, separated by spaces. |
| `-f`, `--contig_file` | File with one contig ID per line; anything after a tab is ignored. Give either `-c` or `-f`. |
| `-o`, `--output_file` | Output file (required). The extension sets the format: `.svg`, `.png`, `.pdf`, ... A log file `<name>.log` is written next to it. |
| `-s`, `--order-by-similarity` | Reorder the contigs so that contigs sharing many protein families are next to each other (hierarchical clustering on Jaccard similarity). Without it, contigs are plotted in the given order. |

### One category

A category colors the features or contigs whose IDs are in a file:

| Option | Description |
|---|---|
| `--category` | A `.sql` file with a query that returns one column of IDs from the whole database, or any other file with one ID per line (blank lines and lines starting with `#` are skipped). |
| `--id-type` | What the IDs are (required with `--category`): `protein` or `pseudo` (feature IDs), `family` (MCL family IDs: colors every protein and pseudogene in the family) or `contig` (colors the contig's name). |
| `--label` | Name of the category in the legend. Default: the file name without extension. |
| `--color` | Color of the category, any matplotlib color such as `red` or `"#ff0000"`. Default: `#d62728`. |
| `--cache` | Only for a `.sql` category: the file to save the IDs of the query in. |
| `--refresh-cache` | Rerun the query and overwrite the IDs in `--cache`, even if they are up to date. |

Without `--cache`, the query runs on every plot. That is fine for fast queries, but queries on protein families can take more than ten minutes on the larger databases. With `--cache`, the IDs are saved together with the database, the query and the time; the query only runs again when one of them no longer matches, or when the tables it reads have changed since. A cache file is itself a list of IDs, so you can also give it to `--category` directly.

For more than one category, use the API.

### Reading the plot

- Each track is one contig, labeled `Contig <ID>`.
- Arrows are proteins, pointing in the direction of their strand. Pseudogenes are drawn transparent.
- Features in a category get the category's color; other features are orange. A contig in a category gets its name
  in the category's color.
- A feature that crosses the end of a circular contig is drawn at both ends, with a red outline.
- Grey links connect features of the same MCL protein family on neighboring tracks. A family with several copies
  links every copy to every copy.

## The `contigviz` API

With the API you can combine as many categories as you like. Below, the hybrid plasmid from the quick start is built up step by step; all steps together are the example script [`examples/plot_categories.py`](examples/plot_categories.py). To follow along, download it with its two data files into one folder:

```bash
for file in plot_categories.py replication_families.txt amr_proteins.sql; do
    curl -O https://raw.githubusercontent.com/jordiweijers/contigviz/main/examples/$file
done
```

### 1. Load and plot contigs

```python
from contigviz import ContigID, Database, fetch_contig, plot_contigs
from contigviz.schemas import kes2021

with Database(kes2021) as database:
    contigs = [fetch_contig(database, ContigID(contig_id)) for contig_id in ("1991", "2082", "4338")]

plot_contigs(contigs, "hybrid_plasmid.svg")
```

`Database(kes2021)` connects to the kes2021 database, and the `with` block closes the connection afterwards. `fetch_contig` loads a contig with its proteins, pseudogenes and their protein families. `plot_contigs` plots the contigs from top to bottom in the given order, so the hybrid in the middle is linked to both parents.

After each of the next steps, `plot_contigs(contigs, "hybrid_plasmid.svg")` shows the result.

### 2. Assign categories: a rule from IDs in your code to color contigs

```python
from contigviz import Category, CategoryRule, assign_category

parent = CategoryRule(Category("Parent", "#1f77b4"), {ContigID("1991"), ContigID("4338")})
hybrid = CategoryRule(Category("Hybrid", "#d62728"), {ContigID("2082")})
assign_category(contigs, parent)
assign_category(contigs, hybrid)
```

A `Category` is a label and a color for the legend. A `CategoryRule` gives its category to every contig or feature with one of its IDs, and `assign_category` applies it. A `ContigID` colors the label of a contig.

### 3. Assign categories: a rule from a file of IDs to color proteins based on their protein family ID

```python
from pathlib import Path
from contigviz import FamilyID

replication = CategoryRule.from_file(Category("Replication", "#5e2bff"), FamilyID, Path("replication_families.txt"))
assign_category(contigs, replication)
```

[`replication_families.txt`](examples/replication_families.txt) has one ID per line. A `FamilyID` colors every protein and pseudogene in that MCL protein family.

### 4. Assign categories: a rule from an SQL query to color proteins based on their protein ID

```python
from contigviz import ProteinID

with Database(kes2021) as database:
    amr = CategoryRule.from_query(
        database, Category("AMR", "#33ffff"), ProteinID, Path("amr_proteins.sql").read_text(), cache=Path("amr_proteins.txt"),
    )
assign_category(contigs, amr)
```

[`amr_proteins.sql`](examples/amr_proteins.sql) runs on the whole database and returns one column of IDs, here protein IDs. With `cache`, the IDs are saved in that file and the query only runs again when the database changes, like `--cache` of `plot-contigs`.

A rule overwrites the category that an earlier rule gave to the same contig or feature, so assign the rule you most want to see last.

### Building blocks

Everything below can be imported from `contigviz`; the docstrings describe all arguments.

| | |
|---|---|
| `contigviz.schemas` | The schemas of the group's databases: `kes2021`, `ESKAPEE25`, `pantoea2025_v2` and `acinetobacter_2023`. |
| `Database(schema)` | A connection to one database. Use it in a `with` block, so it is closed afterwards. It only reads, so it cannot change the database. |
| `fetch_contig(database, contig_id)` | Loads a contig with its proteins, pseudogenes and their protein families. |
| `ContigID`, `ProteinID`, `PseudoID`, `FamilyID` | The kinds of ID, e.g. `ContigID("1991")`. |
| `Category(label, color)` | A label and color for the legend. |
| `CategoryRule(category, ids)` | Gives `category` to the contigs and features with one of the `ids`. Also made with `CategoryRule.from_file(category, id_class, path)` and `CategoryRule.from_query(database, category, id_class, query, cache=None)`. |
| `assign_category(contigs, rule)` | Gives the rule's category to the matching contigs and features. Call it once per rule; a later rule overwrites the category of an earlier one, with a warning. |
| `contig.slice(start, end)` | Cuts a region out of a contig, keeping its coordinates and features. Assign categories before slicing. |
| `order_by_clustering(contigs)` | Returns the contigs ordered by protein family similarity, like `-s`. |
| `plot_contigs(contigs, output_path=None)` | Plots the contigs from top to bottom, saves the plot if you give a path, and returns the matplotlib `Figure`, so you can change it further. |

### Slice contigs

To plot only the first 20 kb of a contig:

```python
with Database(kes2021) as database:
    contig = fetch_contig(database, ContigID("4338"))
plot_contigs([contig.slice(1, 20_000)], "4338_start.svg")
```

Coordinates are 1-based and inclusive, like in GenBank. 

### Logging

contigviz reports what it does through Python's `logging` module. To see its messages in your own script, add these lines (they are also at the top of the example script):

```python
import logging

logging.basicConfig(format="%(message)s")
logging.getLogger("contigviz").setLevel(logging.INFO)
```

## Databases

The schemas in [`src/contigviz/schemas.py`](src/contigviz/schemas.py) describe where the contigs and features are in each database, so the same code can read databases with different table and column names. They connect as the MySQL user `gmg` on `localhost`.

To add a database, copy one of the schemas in that file and change the names; `plot-contigs -d` picks it up by itself. Each feature table needs the columns `start`, `end` and `strand` (`+` or `-`), and `is_pseudo` is an SQL expression that tells whether a feature is a pseudogene, for example `"TRUE"` for a table of only pseudogenes.

## Project layout

```
src/contigviz/      the package
  cli/              the commands, e.g. plot_contigs.py for plot-contigs
  schemas.py        the database schemas
  database.py       Database: the read-only connection
  loader.py         fetch_contig
  category.py       CategoryRule and assign_category
  cache.py          saves and checks the IDs of category queries
  models.py         Contig, Feature, Category and the ID classes
  order.py          order_by_clustering
  plot.py           plot_contigs
examples/           an example script, query and ID file
```

To work on contigviz itself, clone the repository and install it with `pip install -e .`, so your changes take
effect without reinstalling. Planned work is in [`TODO.md`](TODO.md).
