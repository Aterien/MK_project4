"""
Pobiera surowe dane stężenia PM2.5 dla danego roku z archiwum GIOS i zapisuje go lokalnie do pliku .xlsx.

Wywoływany przez Snakemake rule `download_year`.
Dostaje od snakemake object:
    wildcards.year      -rok
    config              -zawartość config.yaml
    output[0]           -ścieżka do wynikowego pliku .xlsx
"""

import sys
import os

import requests
import zipfile
import io

year = int(snakemake.wildcards.year)
config = snakemake.config

base_url = config["pm25_data_link"]
archive_id = config["years"][year]["archive_id"]
filename = config["years"][year]["filename"]
output_path = snakemake.output[0]

print(f"[download_year] Pobieramy dane za rok {year} plik: {filename})")

response = requests.get(f"{base_url}{archive_id}")
response.raise_for_status()

with zipfile.ZipFile(io.BytesIO(response.content)) as z:
    with z.open(filename) as f:
        content = f.read()

with open(output_path, "wb") as out:
    out.write(content)

print(f"[download_year] Zapisany do {output_path}")