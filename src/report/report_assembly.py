"""
report_assembly.py

Skrypt do generowania wykresów porównawczych dla wielu lat.
Uruchamiany przez regułę Snakemake report_task4.
"""

import sys
import os
import matplotlib
matplotlib.use("Agg")  # przed importem pyplot - tryb bez GUI

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from report.multi_year_plots import (
    concat_monthly_means_by_city,
    concat_exceedance_days_by_voivodeship,
    all_years_query_summary,
    year_summary_by_queries_plot
)
from pm25.wizualizacja import wykres_heatmap_srednie, plot_exceedence_by_voivodeship

years = list(snakemake.config["years"].keys())
pm25_daily_limit = float(snakemake.config["pm25_daily_limit"])

output_heatmap       = snakemake.output[0]
output_voivodeship   = snakemake.output[1]
output_query_summary = snakemake.output[2]

os.makedirs(os.path.dirname(output_heatmap), exist_ok=True)

# --- 1. Heatmapa średnich miesięcznych po miastach ---
print(f"[report] Wczytywanie monthly_means_by_city dla lat {years}...")
monthly_means_by_city = concat_monthly_means_by_city(
    years=years,
    input_files=snakemake.input.monthly_means_by_city
)
wykres_heatmap_srednie(
    srednie_po_miejscach=monthly_means_by_city,
    lata=years,
    output_file_name=output_heatmap,
    show=False
)
print(f"[report] Zapisano heatmapę -> {output_heatmap}")

# --- 2. Wykres przekroczeń normy po województwach ---
print(f"[report] Wczytywanie exceedance_days_by_voivodeship dla lat {years}...")
exceedance_by_voivodeship = concat_exceedance_days_by_voivodeship(
    years=years,
    input_files=snakemake.input.exceedance_by_voivodeship
)
plot_exceedence_by_voivodeship(
    df=exceedance_by_voivodeship,
    daily_norm=pm25_daily_limit,
    output_file_name=output_voivodeship,
    show=False
)
print(f"[report] Zapisano wykres przekroczeń -> {output_voivodeship}")

# --- 3. Wykres liczby publikacji PubMed po zapytaniach i latach ---
print(f"[report] Wczytywanie summary_by_year dla lat {years}...")
query_summary_df = all_years_query_summary(
    years=years,
    input_files=snakemake.input.query_summary
)
year_summary_by_queries_plot(
    all_years_query_summary_df=query_summary_df,
    output_file_name=output_query_summary,
    show=False
)
print(f"[report] Zapisano wykres zapytań -> {output_query_summary}")

print(f"[report] Wykresy zakończone pomyślnie.")