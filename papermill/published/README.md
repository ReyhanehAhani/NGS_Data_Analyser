# Published Papermill runs (actually executed)

These files were produced by `papermill` + `jupyter nbconvert` so it is obvious the parameterized pipeline runs end-to-end.

| File | Description |
|------|-------------|
| `ngs_run_report_executed.ipynb` | Executed notebook (parameters from `params/example.yaml`). |
| `ngs_run_report_executed.html` | Same run as static HTML for quick viewing in a browser. |

## Reproduce

```bash
pip install -r papermill/requirements-papermill.txt
bash papermill/run_papermill.sh papermill/params/example.yaml
python3 -m jupyter nbconvert --to html papermill/reports/ngs_run_report_*.ipynb --output-dir papermill/published
```

(Or rerun `papermill` targeting `papermill/published/` as in the committed history.)
