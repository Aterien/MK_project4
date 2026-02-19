#TODO string co robi ten skrypt

import argparse
import yaml
from pathlib import Path
from pubmed_analysis import *

# W zadaniu było wymagane podanie roku i configu jako argumenty skryptu, a nie przez Snakemake
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
queries_list = config["pubmed_queries"]
retmax = int(config["retmax"])
max_top_journals_sample_size = config["sample_size"]

publications_per_year = {}

# Analiza
for query in queries_list:
    result_df = search_for_papers(query = query,
                                  year = year,
                                  entrez_email= user_entrez_email,
                                  entrez_api_key = user_entrez_api_key,
                                  retmax = retmax)
    result_df.to_csv(f"../../results/literature/{year}/{query.replace(" ","_")}_search_result.csv", index=False)

    # Dodajemy liczbę wyników zapytania
    publications_per_year[query] = result_df.shape[0]
    print(f"[pubmed_year][{year}] Wyniki zapytania {query} zapisano do results/literature/{year}/{query.replace(" ","_")}_search_result.csv")

# Tworzymy tabelę z liczbą publikacji dla wszystkich zapytań za rok
publications_per_year_df = pd.DataFrame(publications_per_year, index = [year])
publications_per_year_df.to_csv(f"../../results/literature/{year}/summary_by_year.csv", index=True)

top_journals_df, total_found = top_journals(year = year,
                               entrez_email = user_entrez_email,
                               entrez_api_key = user_entrez_api_key,
                               top = 10,
                               sample_size=max_top_journals_sample_size)
top_journals_df.to_csv(f"../../results/literature/{year}/top_journals.csv", index=True)

# Wizualizacja (jedyny wykres, który ma sens przy analizie danych za pojedyńczy rok)
figure_top_journals_per_year(top_journals_df,f"../../results/literature/{year}/figure_top_journals_per_year.png")