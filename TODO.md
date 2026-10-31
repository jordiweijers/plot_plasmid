# TODO

## Scripts

- [ ] Synteny workflow: find the loci where a given set of protein families (an operon) lie next to each other, cut them out with a few flanking genes and plot them. Each operon family gets its own color; `--category` adds one more.

## Package

- [ ] `plot.py`: a feature that wraps around the origin gets an arrowhead on both parts; it should only get one, at the end it points to.
- [ ] `plot.py`: give a clear error when two contigs get the same track name (pygenomeviz raises a ValueError).
- [ ] Load feature and contig labels from the database.
- [ ] Support slices that cross the origin of a circular contig (`start > end` on `Contig.slice`).
- [ ] Tests under `tests/` with pytest.

## Repository

- [ ] Update `README.md` (still describes `plot_plasmid`, `bin/` and `config.yaml`).
- [ ] Rename the GitHub repository, then update the Repository URL in `pyproject.toml` and the git remote.
