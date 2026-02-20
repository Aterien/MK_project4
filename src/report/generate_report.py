"""
generate_report.py

Skrypt generujący raport Markdown podsumowujący wyniki analizy PM2.5
oraz przeglądu literatury PubMed dla wszystkich lat.
Uruchamiany przez regułę Snakemake generate_report.
"""

import os
import sys
import pandas as pd
from pathlib import Path

years   = list(snakemake.config["years"].keys())
queries = snakemake.config["pubmed_queries"]
cities  = snakemake.config["cities"]
years_str = "_".join(str(y) for y in years)

output_md = snakemake.output[0]
os.makedirs(os.path.dirname(output_md), exist_ok=True)

# Ścieżki do plików - pobieramy z snakemake.input
city_figures         = snakemake.input.city_figures          # list, jeden na rok
heatmap              = snakemake.input.heatmap
voivodeship_plot     = snakemake.input.voivodeship_plot
query_summary_plot   = snakemake.input.query_summary_plot
query_csvs           = snakemake.input.query_csvs            # list: wszystkie {query}_search_result.csv
top_journal_figures  = snakemake.input.top_journal_figures   # list, jeden na rok


def rel(path):
    """Zwraca ścieżkę relatywną względem katalogu wynikowego raportu."""
    return os.path.relpath(path, os.path.dirname(output_md))


lines = []

# ============================================================
# Nagłówek
# ============================================================
lines += [
    f"# Raport zadanie 4",
    f"",
    f"Analiza obejmuje lata: **{', '.join(str(y) for y in years)}**",
    f"",
]

# ============================================================
# Sekcja PM2.5
# ============================================================
lines += [
    f"## Analiza stężenia PM2.5",
    f"",
    f"### Średnie miesięczne stężenie dla wybranych miast ({', '.join(cities)})",
    f"",
]

for year, fig_path in zip(years, city_figures):
    lines += [
        f"#### Rok {year}",
        f"",
        f"![Porównanie miast {year}]({rel(fig_path)})",
        f"",
    ]

lines += [
    f"### Heatmapa stężeń dla wszystkich wspólnych miast",
    f"",
    f"![Heatmapa średnich miesięcznych]({rel(heatmap)})",
    f"",
    f"### Liczba dni z przekroczeniem normy po województwach",
    f"",
    f"![Przekroczenia normy po województwach]({rel(voivodeship_plot)})",
    f"",
]

# ============================================================
# Sekcja literatura
# ============================================================
lines += [
    f"## Przegląd literatury",
    f"",
    f"### Liczba wyników wyszukiwań dla poszczególnych zapytań",
    f"",
    f"Analizowane zapytania:",
    f"",
]

for q in queries:
    lines.append(f"- `{q}`")
lines.append("")

lines += [
    f"![Liczba wyników wyszukiwań]({rel(query_summary_plot)})",
    f"",
    f"#### Pierwszy wynik dla każdego zapytania",
    f"",
]

# Dla każdego zapytania i roku - pierwszy wiersz z CSV
for query in queries:
    query_slug = query.replace(" ", "_")
    lines += [
        f"##### `{query}`",
        f"",
    ]
    for year in years:
        # Szukamy pasującego pliku w liście query_csvs
        matching = [f for f in query_csvs if f"{year}/{query_slug}_search_result" in f]
        if not matching:
            lines.append(f"*Rok {year}: brak danych*")
            lines.append("")
            continue
        try:
            df = pd.read_csv(matching[0])
            if df.empty:
                lines.append(f"*Rok {year}: brak wyników*")
                lines.append("")
                continue
            row = df.iloc[0]
            lines += [
                f"**Rok {year}**",
                f"",
                f"- **Tytuł:** {row.get('title', 'brak')}",
                f"- **Autorzy:** {row.get('authors', 'brak')}",
                f"- **Czasopismo:** {row.get('journal', 'brak')}",
                f"- **PMID:** {row.get('PMID', 'brak')}",
                f"",
            ]
        except Exception as e:
            lines.append(f"*Rok {year}: błąd odczytu — {e}*")
            lines.append("")

lines += [
    f"### Czasopisma z największą liczbą publikacji",
    f"",
]

for year, fig_path in zip(years, top_journal_figures):
    lines += [
        f"#### Rok {year}",
        f"",
        f"![Top czasopisma {year}]({rel(fig_path)})",
        f"",
    ]

# ============================================================
# Zapis pliku
# ============================================================
with open(output_md, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(f"[generate_report] Raport zapisany -> {output_md}")
