#TODO string co robi ten skrypt
import argparse
import yaml
from pathlib import Path
from pubmed_analysis import *

# Parsowanie
parser = argparse.ArgumentParser()
parser.add_argument("--year", help="Rok do analizy")
parser.add_argument("--config", help="Ścieżka do pliku config.yaml")

args = parser.parse_args()
year = int(args.year)

config_file = Path(args.config)
with open(config_file) as f:
    config = yaml.safe_load(f)

user_entrez_email = config["user_entrez_email"]
user_entrez_api_key = config["user_entrez_api_key"]
queries_list = config["queries"]

# Analiza


