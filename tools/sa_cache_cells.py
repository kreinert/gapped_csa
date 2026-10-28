"""List the (dataset, shape, blob) cells of a run_suite.py SA cache.

Each sidecar in <sa-cache>/meta/ names one blob. Gapped shapes have one sidecar
each. Integer shapes (gcsa's regular mode) share one sidecar per dataset,
<dataset>__plain.json with "shape": "plain", because the plain SA is the same
for every k -- so a plain sidecar expands to one cell per k in `plain_ks`.
"""
import glob, json, os, sys

PLAIN_SA_KEY = "plain"


def default_plain_ks():
    """The integer shapes in bench/datasets.py SHAPES (empty if none)."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(root, "bench"))
    try:
        import datasets
        return [s for s in datasets.SHAPES if str(s).isdigit()]
    except Exception:
        return []


def sidecar_cells(sa_cache, plain_ks, only=None, warn=print):
    cells, skipped = [], []
    for p in sorted(glob.glob(os.path.join(sa_cache, "meta", "*.json"))):
        try:
            m = json.load(open(p))
        except (OSError, ValueError):
            continue
        ds, sh = m.get("dataset"), m.get("shape")
        blob = os.path.join(sa_cache, m.get("blob", ""))
        if not ds or not sh or not os.path.isfile(blob):
            continue
        if only and ds not in only:
            continue
        if sh == PLAIN_SA_KEY or sh.isdigit():   # isdigit: sidecars written before sharing
            if not plain_ks:
                skipped.append(ds)
            for k in plain_ks:
                cells.append((ds, str(k), blob))
        else:
            cells.append((ds, sh, blob))
    if skipped:
        warn(f"[note] plain-SA sidecars for {len(skipped)} dataset(s) skipped "
             f"({', '.join(skipped[:3])}{'...' if len(skipped) > 3 else ''}): "
             f"pass --plain-k K ... to say which k to run on them")
    return sorted(set(cells))
