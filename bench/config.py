"""Where things live on *your* machine -- kept out of everything tracked by
git, on purpose: different people keep the benchmark data in different
places, and no real filesystem path should end up in a commit or a GitHub
diff.

Three settings, resolved in this order (first one found wins):

  1. CLI flags -- run_suite.py / fetch_data.py accept --data-dir / --bin-dir
     / --tmp-dir directly, for a one-off override.
  2. Environment variables -- GCSA_BENCH_DATA_DIR / GCSA_BENCH_BIN_DIR /
     GCSA_BENCH_TMP_DIR, for CI or a shell profile.
  3. bench/config.local.py -- a gitignored file (copy config.local.py.example
     to create it) defining DATA_DIR / BIN_DIR / TMP_DIR as plain strings,
     for a setting you want to persist across sessions on one machine
     without exporting env vars every time.
  4. The built-in defaults below, which assume the common layout: bench/ and
     src/ both directly under the gapped_csa checkout, with a data folder as
     a sibling of that checkout. Works out of the box for that layout;
     everyone else uses #2 or #3.

Usage (from run_suite.py / fetch_data.py):
    import config
    data_dir = config.resolve("DATA_DIR", args.data_dir)
"""
import importlib.util
import os
from pathlib import Path

_HERE = Path(__file__).resolve().parent          # .../gapped_csa/bench
_REPO_ROOT = _HERE.parent                          # .../gapped_csa

_DEFAULTS = {
    "DATA_DIR": _REPO_ROOT.parent / "data",
    "BIN_DIR": _REPO_ROOT,
    "TMP_DIR": Path("/tmp/gcsa_bench"),
}

_SETTINGS = ("DATA_DIR", "BIN_DIR", "TMP_DIR")


def _load_local_overrides() -> dict:
    local_path = _HERE / "config.local.py"
    if not local_path.exists():
        return {}
    spec = importlib.util.spec_from_file_location("gcsa_bench_config_local", local_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return {k: getattr(mod, k) for k in _SETTINGS if hasattr(mod, k)}


_LOCAL = _load_local_overrides()


def resolve(name: str, cli_value=None) -> Path:
    """Resolve one of DATA_DIR/BIN_DIR/TMP_DIR through the precedence chain
    described above. `cli_value` is whatever argparse got for the
    corresponding --data-dir/--bin-dir/--tmp-dir flag (None if not passed)."""
    if name not in _SETTINGS:
        raise ValueError(f"unknown setting {name!r}, expected one of {_SETTINGS}")
    if cli_value is not None:
        return Path(cli_value)
    env_key = f"GCSA_BENCH_{name}"
    if env_key in os.environ:
        return Path(os.environ[env_key])
    if name in _LOCAL:
        return Path(_LOCAL[name])
    return _DEFAULTS[name]


if __name__ == "__main__":
    # `python3 config.py` -- print what each setting currently resolves to,
    # and which source it came from, without running anything else.
    for name in _SETTINGS:
        env_key = f"GCSA_BENCH_{name}"
        if env_key in os.environ:
            source = f"env {env_key}"
        elif name in _LOCAL:
            source = "config.local.py"
        else:
            source = "built-in default"
        print(f"{name:10} = {resolve(name)}   ({source})")
