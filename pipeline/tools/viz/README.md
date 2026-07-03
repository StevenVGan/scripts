# scripts/pipeline/tools/viz/

Shared figure helpers, promoted from the lab's multiome analyses (~30 style
copies, a handful for the profile plotter). Existing per-project copies are
stale-by-design (like the pipeline-step copies) — new work imports these.

- **`_figure_style.py`** — `apply_publication_style()` sets the lab's canonical
  matplotlib rcParams: the **big-canvas** style (font 16/20/17, sized for
  14-22in figures, `pdf.fonttype=42`, sans-serif). A smaller variant
  (font 13/11, dpi 200) exists in one project's analysis; reconcile there if you
  switch the house style.
- **`_profile_plot.py`** — `plot_profile()`, a deepTools computeMatrix profile
  plotter (legend in the rightmost subplot, honors `_figure_style` rcParams, one
  subplot per region group). Pairs with `tools/heatmap.sh`.

Use from an analysis `script/`:

```python
import sys, matplotlib
matplotlib.use("Agg")
sys.path.insert(0, "/mnt/home/digan/work/scripts/pipeline/tools/viz")   # adjust to this node's WORK root
from _figure_style import apply_publication_style; apply_publication_style()
from _profile_plot import plot_profile
```

(or copy the two files into the analysis `script/` for a self-contained dir.)
