import sys
import os

# Add the project root to sys.path so tests can import app-level modules
# (parsers, and examples.inference_demo) without installing the package.
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
