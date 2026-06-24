"""Single source of truth: load config.yaml. CLI flags override these defaults."""
from __future__ import annotations
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config.yaml"


def load_config(path=CONFIG_PATH) -> dict:
    return yaml.safe_load(Path(path).read_text())


CONFIG = load_config()
