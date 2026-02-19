"""
Pobiera metadane z archiwum GIOS archive and saves it
as a local .xlsx file.

Wywoływany przez Snakemake rule `download_metadata`.
Dostaje od snakemake object:
    config      -zawartość config.yaml
    output[0]   -ścieżka do wynikowego pliku .xlsx
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import requests

config = snakemake.config

base_url = config["pm25_data_link"]
metadata_id = config["metadata_id"]
output_path = snakemake.output[0]

os.makedirs(os.path.dirname(output_path), exist_ok=True)

print(f"[download_metadata] Pobieramy metadane")

url = f"{base_url}{metadata_id}"
response = requests.get(url)
response.raise_for_status()

with open(output_path, "wb") as out:
    out.write(response.content)

print(f"[download_metadata] Metadane zapisane do {output_path}\n\n\n")