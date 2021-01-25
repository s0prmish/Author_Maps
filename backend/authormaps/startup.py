import os
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Default Directory Paths
# home_dir = os.path.expanduser("~")  # Also acceptable
home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project", "group5")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
DATA_DIR = os.path.join(PROJECT_DIR, "data")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Logging Configuration
LOG_FILE_PATH = os.path.join(LOG_DIR, "group5_log.log")
logging.basicConfig(filename=LOG_FILE_PATH,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
