"""Make the repository root importable under pytest and unittest."""
from pathlib import Path
import sys
ROOT = str(Path(__file__).parents[1])
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
