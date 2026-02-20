"""
report_assembly.py

Skrypt do generowania wykresów porównawczych dla wielu lat.
Uruchamiany przez regułę Snakemake report_task4.
"""

import sys
import os
import matplotlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from report.multi_year_plots import (
    concat_monthly_means_by_city,
    concat_exceedance_days_by_voivodeship
)
from pm25.wizualizacja import wykres_heatmap_srednie, plot_exceedence_by_voivodeship

years = list(snakemake.config["years"].keys())
pm25_daily_limit = float(snakemake.config["pm25_daily_limit"])
output_heatmap = snakemake.output[0]
output_voivodeship = snakemake.output[1]

# --- 1. Heatmap średnich miesięcznych po miastach ---
print(f"[report_task4] Wczytywanie monthly_means_by_city dla lat {years}...")
monthly_means_by_city = concat_monthly_means_by_city(years=years)

print(f"[report_task4] Rysowanie heatmapy...")
wykres_heatmap_srednie(
    srednie_po_miejscach=monthly_means_by_city,
    lata=years,
    output_file_name=output_heatmap,
    show=False
)
print(f"[report_task4] Zapisano heatmapę -> {output_heatmap}")

# --- 2. Wykres przekroczeń normy po województwach ---
print(f"[report_task4] Wczytywanie exceedance_days_by_voivodeship dla lat {years}...")
exceedance_by_voivodeship = concat_exceedance_days_by_voivodeship(years=years)

print(f"[report_task4] Rysowanie wykresu przekroczeń po województwach...")
plot_exceedence_by_voivodeship(
    df=exceedance_by_voivodeship,
    daily_norm=pm25_daily_limit,
    output_file_name=output_voivodeship,
    show=False
)
print(f"[report_task4] Zapisano wykres -> {output_voivodeship}")

print(f"[report_task4] Raport zakończony pomyślnie.")