"""Locate writable portable-app files independently of the current directory."""
from pathlib import Path
import sys

def application_root():
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[1]

VERSION = "0.1.0"
