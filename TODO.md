# TODO

## Scripts

- [ ] Move `scripts/plot_contigs.py` into the package as `contigviz.cli` and remove `scripts/` (shared options such as `-d`, the cache and logging in `cli/common.py`), with commands in `[project.scripts]` of `pyproject.toml`. `queries/` stays in the repo, outside the package.
- [ ] Generic plot command: one optional category with `--category FILE` (a `.sql` query, run with the cache, or a `.txt` with one ID per line), `--id-type` (protein, pseudo, family or contig; required with `--category`), `--label` (default: the file name) and `--color`. More categories means writing a script on the API.
- [ ] Synteny workflow: find the loci where a given set of protein families (an operon) lie next to each other, cut them out with a few flanking genes and plot them. Each operon family gets its own color; `--category` adds one more.

## Package

- [ ] `loader.py`: replace the `ContigLoader` class with a function `fetch_contig(database, contig_id)`.
- [ ] `plot.py`: a feature that wraps around the origin gets an arrowhead on both parts; it should only get one, at the end it points to.
- [ ] `plot.py`: give a clear error when two contigs get the same track name (pygenomeviz raises a ValueError).
- [ ] Load feature and contig labels from the database.
- [ ] Support slices that cross the origin of a circular contig (`start > end` on `Contig.slice`).
- [ ] Tests under `tests/` with pytest.

## Repository

- [ ] Update `README.md` (still describes `plot_plasmid`, `bin/` and `config.yaml`).
- [ ] Rename the GitHub repository, then update the Repository URL in `pyproject.toml` and the git remote.
