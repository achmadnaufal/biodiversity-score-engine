"""Shared pytest configuration.

Prepends the repo root to ``sys.path`` so tests can ``import src.*``
whether pytest is invoked from the repo root or from ``tests/``.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
