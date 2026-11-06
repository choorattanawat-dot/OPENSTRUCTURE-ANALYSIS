# Contributing

OPENSTRUCTURE-ANALYSIS is a structural-engineering library. A change is useful when a reviewing engineer can follow the arithmetic and the tests still match a known result.

## Setup

```bash
python -m pip install -e ".[dev]"
ruff check src tests
pytest
```

The analysis tests need OpenSeesPy. Design tests do not call the solver.

## Layout

- `src/openstructure/model.py` — nodes, members, loads, combinations. No OpenSees calls.
- `src/openstructure/analysis.py` — linear-static OpenSees engine and section-force diagrams.
- `src/openstructure/design/` — code plugins and RC calculations. No file or plot I/O.
- `src/openstructure/io.py`, `pipeline.py`, `plotting.py`, `cli.py` — files, reports, figures, command line.

Keep new analysis types (pushover, fiber sections) behind the same model and combination objects. Do not fold them into the YAML loader.

## Units

Internal analysis uses N and mm. Public results use kN, m, kN·m, and MPa, and the names say so. `1 kN/m` is the same number as `1 N/mm`.

## Design codes

Do not invent a clause number or a factor. If a value is taken from ACI 318-19, cite the section. If a วสท. 1008 value is only known from a textbook rounding or from the ACI edition that standard was based on, mark it `ASSUMPTION` in `docs/CODE_ASSUMPTIONS.md` and in the code plugin's `assumptions()` list. The calculation report prints that list.

## Benchmarks

A new analysis or design check needs a test whose expected number is written out as arithmetic, with the source named in the test docstring and in `docs/BENCHMARKS.md`. Matching the library against a second copy of the same function is not a benchmark.

## Notebooks

`examples/two_bay_two_storey.ipynb` should keep running after a change to the public API or the example file.
