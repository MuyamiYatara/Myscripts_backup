

import os

from pathlib import Path
import numpy as np
import yaml

class ConfigManager:
    def __init__(self, config_path: str):
        self.config = self._load_config(config_path)

    def _load_config(self, file_path: str) -> dict:
        """Load configuration from a YAML file."""
        with open(file_path, 'r') as file:
            return yaml.safe_load(file)

    def update_config(self, new_config: dict):
        """Update the current configuration with the provided dictionary."""
        for key,value in new_config.items():
            if key not in self.config:
                self.config[key] = {}

            self.config[key]['value'] = value
            
    def get_config(self) -> dict:
        """Retrieve the current configuration."""
        return self.config
    
# Important directory paths
FILE = Path(__file__).resolve()
PKG_DIR = str(FILE.parents[1])  # pyprocar
ROOT = str(FILE.parents[2])  # PyProcar
LOG_DIR=os.path.join(ROOT,'logs')
DATA_DIR=os.path.join(ROOT,'data')
CONFIG_FILE=os.path.join(PKG_DIR,'cfg','package.yml')

# Load config from yaml file
with open(CONFIG_FILE, 'r') as f:
    CONFIG = yaml.safe_load(f)