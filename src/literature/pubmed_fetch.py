#TODO string co robi ten skrypt

# Skrypt pobiera dane z PubMed dla podanego roku i zapisuje wyniki do katalogu {output}/{year}/
# Może być uruchamiany samodzielnie lub przez Snakemake (rule pubmed_year)

import argparse
import yaml
import os
from pathlib import Path
from pubmed_analysis import *

parser = argparse.ArgumentParser()
parser.add_argument("--year", required=True, help="Rok do analizy")
parser.add_argument("--config", required=True, help="Ścieżka do pliku config.yaml")
parser.add_argument("--output", required=True, help="Katalog wyjściowy (np. results/literature)")

args = parser.parse_args()
year = int(args.year)

config_file = Path(args.config)
with open(config_file) as f:
    config = yaml.safe_load(f)

# Katalog wyjściowy dla danego roku
output_dir = Path(args.output) / str(year)
os.makedirs(output_dir, exist_ok=True)
print(f"[pubmed_year][{year}] Katalog wyjściowy: {output_dir}")

user_entrez_email = config["user_entrez_email"]
user_entrez_api_key = config["user_entrez_api_key"]
queries_list = config["pubmed_queries"]
retmax = int(config["retmax"])
max_top_journals_sample_size = config["sample_size"]

publications_per_year = {}

# --- Analiza zapytań ---
for query in queries_list:
    result_df = search_for_papers(
        query=query,
        year=year,
        entrez_email=user_entrez_email,
        entrez_api_key=user_entrez_api_key,
        retmax=retmax
    )
    query_filename = query.replace(" ", "_")
    out_path = output_dir / f"{query_filename}_search_result.csv"
    result_df.to_csv(out_path, index=False)
    publications_per_year[query] = result_df.shape[0]
    print(f"[pubmed_year][{year}] Zapisano wyniki zapytania '{query}' -> {out_path}")

# Summary
publications_per_year_df = pd.DataFrame(publications_per_year, index=[year])
summary_path = output_dir / "summary_by_year.csv"
publications_per_year_df.to_csv(summary_path, index=True)
print(f"[pubmed_year][{year}] Zapisano summary -> {summary_path}")

# Top journals
top_journals_df = top_journals(
    year=year,
    entrez_email=user_entrez_email,
    entrez_api_key=user_entrez_api_key,
    top=10,
    sample_size=max_top_journals_sample_size
)
top_journals_path = output_dir / "top_journals.csv"
top_journals_df.to_csv(top_journals_path, index=False)
print(f"[pubmed_year][{year}] Zapisano top journals -> {top_journals_path}")

# --- Wizualizacja ---
figure_path = output_dir / "figure_top_journals_per_year.png"
figure_top_journals_per_year(top_journals_df, str(figure_path))
print(f"[pubmed_year][{year}] Zapisano wykres -> {figure_path}")

print(f"[pubmed_year][{year}] Zakończono pomyślnie.")