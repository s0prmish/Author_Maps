import os
from pathlib import Path

home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project", "group5")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
DATA_DIR = os.path.join(PROJECT_DIR, "data")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)


