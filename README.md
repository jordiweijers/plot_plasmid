# plot_plasmid

> [!IMPORTANT]
> This tool is made for members of the Dagan group. It reads its data from the group's MySQL databases,
> which are installed locally on the HPC (`cauclusterg`), so it only works there and with access to
> those databases. The code is public for reference, but it will not run elsewhere without your own
> databases and an adapted `config.yaml`.

Plot plasmid contigs from a MySQL database with [pyGenomeViz](https://github.com/moshi4/pyGenomeViz).
Each contig is drawn as a track of its proteins and pseudogenes. Features are colored by category
(for example AMR, replication or conjugation), and the features of neighboring contigs that belong
to the same MCL protein family are connected by links.

Supported out of the box: the `kes2021`, `ESKAPEE25` and `pantoea2025_v2` databases. Other databases
can be added in `config.yaml` without changing the code.

## Installation

Requires Python 3.9 or newer, and an account on the HPC with access to the group's MySQL databases.

```bash
git clone https://github.com/jordiweijers/plot_plasmid.git
cd plot_plasmid
python3 -m venv venv
source venv/bin/activate
pip install -e .
plot-plasmid -h
```

While the venv is active, the `plot-plasmid` command is on your `PATH`. In a new session, activate it
again with `source /path/to/plot_plasmid/venv/bin/activate`, and leave it with `deactivate`.

Install with `-e`: the command then reads `config.yaml` and `queries/` from your clone, and saves its
category lists in `cache/` there. 

The first plot of a database runs its category queries and saves the results in `cache/`. For the
protein family categories this takes a while: about 12 minutes per query on kes2021, and longer on
ESKAPEE25. Later plots read the saved lists and take seconds.

## Usage

```bash
# Plot three contigs from kes2021, in the given order
plot-plasmid -d kes2021 -c 1991 2028 4338 -o plots/kes2021_example.svg

# Plot the contigs listed in a file, ordered by protein family similarity, as PNG
plot-plasmid -d ESKAPEE25 -f my_contigs.tsv -s -o my_contigs.png
```

| Option | Description |
|---|---|
| `-d`, `--database` | Database to load the contigs from, as named in `config.yaml` (required). |
| `-c`, `--contig_ids` | Contig IDs to plot, separated by spaces. |
| `-f`, `--contig_file` | File with one contig ID per line (see below). Give either `-c` or `-f`. |
| `-o`, `--output_file` | Output file (required). The extension sets the format: `.svg`, `.png`, `.pdf`, ...  |
| `-s`, `--order-by-similarity` | Reorder the contigs so that contigs sharing many protein families are next to each other (hierarchical clustering on Jaccard similarity). Without it, contigs are plotted in the given order. |
| `--config` | Use another config file instead of the `config.yaml` of the installation. |
| `--refresh-cache` | Rerun the category queries and overwrite their saved lists, even if they are up to date. |

A log file with the same name as the plot (`<name>.log`) is written next to it.

### Contig file

One contig ID per line. Optionally add a tab and a label: the contig's name in the plot then gets a
color per label, which is shown in the legend.

```
1991	parent
2028	hybrid
4338	parent
```

### Reading the plot

- Each track is one contig, labeled `Contig <ID>`.
- Arrows are proteins, pointing in the direction of their strand. Pseudogenes are drawn transparent.
- Features in a category get the category's color; other features are orange.
- A feature that crosses the end of a circular contig is drawn at both ends, with a red outline.
- Grey links connect features of the same MCL protein family on neighboring tracks. A family with
  several copies links every copy to every copy.


## Configuration

Everything database-specific lives in `config.yaml`: where the contigs and features are in each
database, and the categories.

### Databases

Each database describes where its contigs and features are, so the same code can read databases with
different table and column names.

```yaml
databases:
  kes2021:                         # the MySQL database; also the name to give with -d
    contigs:
      table: contigs               # table with one row per contig
      id_column: contig_ID         # contig ID
      length_column: Size          # contig length in bp
    features:                      # one entry per table of features
      - table: proteins            # table with one row per feature
        id_column: id              # feature ID
        contig_column: contig_ID   # contig the feature is on
        is_pseudo: "FALSE"         # SQL expression: is the feature a pseudogene?
        mcl_table: mcl30           # table that assigns features to MCL families
        mcl_id_column: id          # feature ID in mcl_table
        mcl_family_column: clust   # MCL family ID in mcl_table
      - table: pseudogenes
        id_column: pseudo_id
        contig_column: contig_ID
        is_pseudo: "TRUE"
        mcl_table: pmcl
        mcl_id_column: pseudo_id
        mcl_family_column: clust
```

- **features** lists every table with features. kes2021 and ESKAPEE25 keep proteins and pseudogenes in
  separate tables, so they have two entries; pantoea2025_v2 has both in one table, `features_cds`.
- Each feature table also needs the columns `start`, `end` and `strand` (`+` or `-`). These names are fixed.
- **is_pseudo** is an SQL expression that MySQL evaluates for each feature, in quotes. Use `"TRUE"` or
  `"FALSE"` for a table with only pseudogenes or only proteins. For a table with both, use an
  expression on its columns, for example pantoea2025_v2:
  `"FIND_IN_SET('pseudo', REPLACE(attributes, ';', ',')) > 0"`.

To add a database: copy an existing entry, change the database name and the table and column names, and
add categories with `database:` set to the new name.

### Categories

A category colors the proteins, pseudogenes, protein families or contigs whose IDs it contains.

```yaml
categories:
  - name: Replication                                  # label in the legend
    color: '#5e2bff'                                   # any matplotlib color, in quotes
    type: family                                       # protein, pseudo, family or contig
    database: kes2021                                  # the database it applies to
    query: queries/kes2021/replication_families.sql    # SQL that returns one column of IDs
    path: cache/kes2021/replication_families.txt       # where the IDs are saved
```

- **type** says what the IDs are: `protein` or `pseudo` (feature IDs), `family` (MCL family IDs: colors
  all proteins and pseudogenes in the family) or `contig` (colors the contig's name).
- Give a **query**, a **path**, or both:
  - *query only*: the query runs on every plot;
  - *path only*: the IDs are read from the file, for example a list you made yourself
    (one ID per line; blank lines and lines starting with `#` are skipped);
  - *both*: the query runs once and its IDs are saved at `path`; later plots read the file.
- Queries run on the whole database and must return exactly one column.
- Paths are relative to `config.yaml`.
- If a feature is in several categories, the first one in the list wins.
- One category applies to one database. For the same category in another database, add another entry;
  categories with the same name share their legend entry and the color of the first one.

The included categories:

| Category | Databases | Rule |
|---|---|---|
| AMR | kes2021, ESKAPEE25 | Proteins in the `amr` table. |
| Replication | all three | MCL families whose most common protein name contains "replication initiat". |
| Conjugation | all three | MCL families whose most common protein symbol is `tra` plus one capital letter (`traA`, `traB`, ...). |

A category with both a query and a path saves its IDs in a list under `cache/`.

## Project layout

```
bin/            the Python package (imported as plot_plasmid)
  main.py       command line and run_plot()
  config.py     reads and checks config.yaml
  database.py   loads categories, contigs and features from MySQL
  cache.py      reads, writes and checks the saved category lists
  order.py      orders contigs by protein family similarity
  plot.py       draws the plot
  models.py     data classes: Contig, Feature, Category
config.yaml     databases and categories
queries/        SQL files of the categories, per database
cache/          saved category lists (created on first use, not in git)
```

`run_plot()` in `main.py` does the same as the command, for use from another Python script or a
notebook.
