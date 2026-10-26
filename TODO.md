# TODO

## Scripts

- [ ] `scripts/utils.py`: let `setup_logging` configure the root (or `contigviz`) logger instead of `"main_logger"`, so the messages of the package reach its handlers.
- [ ] `scripts/main.py`: rewrite on the new API (argparse, logging, `Database(schema)`, `ContigLoader`, `CategoryRule`s per database, `assign_categories` → `order_by_clustering` → `plot_contigs`).
- [ ] Decide whether `queries/` moves to `scripts/queries/`.

## Package

- [ ] `__init__.py`: decide whether to export the models (e.g. `from contigviz import Contig`).
- [ ] `loader.py`: mention the length-0 contig error under Raises in the `fetch_contig` docstring.
- [ ] `plot.py`: check `default_color` with `is_color_like`, like `Category` does.
- [ ] `plot.py`: a feature that wraps around the origin gets an arrowhead on both parts; it should only get one, at the end it points to.
- [ ] `plot.py`: give a clear error when two contigs get the same track name (pygenomeviz raises a ValueError).
- [ ] Load feature and contig labels from the database.
- [ ] Support slices that cross the origin of a circular contig (`start > end` on `Contig.slice`).
- [ ] Tests under `tests/` with pytest.

## Repository

- [ ] Update `README.md` (still describes `plot_plasmid`, `bin/` and `config.yaml`).
- [ ] Rename the GitHub repository, then update the Repository URL in `pyproject.toml` and the git remote.
