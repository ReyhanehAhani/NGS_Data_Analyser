# Pillar 2 — Reproducibility (Papermill)

Parameterized runs for **batch reports** and audit-friendly reruns.

## Layout

| File | Purpose |
|------|---------|
| `params/example.yaml` | Default parameters; copy to `params/local.yaml` (gitignored optional). |
| `notebooks/ngs_run_report.ipynb` | Tagged `parameters` cell + run checklist + optional stats over `input_root`. |
| `run_papermill.sh` | One command to execute the notebook into `reports/`. |

## Quick start

```bash
pip install -r papermill/requirements-papermill.txt
bash papermill/run_papermill.sh
```

Artifacts land in `papermill/reports/` (see `.gitkeep`).

## Tie-in to existing Python tools

Your repo centers on `filter_all.py`, `filter_big_data.py`, etc. This notebook does **not** replace them—it records **which run_id / paths / notes** accompanied a batch, and previews row counts if CSVs exist under `input_root`.

## Resume

- *Wrapped genomic filtering batches in **Papermill**-parameterized notebooks for reproducible reporting.*
