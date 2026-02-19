import sys
import os
import pandas as pd
from wczytaj_wyczysc import *
from analiza import *
from wizualizacja import *

year = int(snakemake.wildcards.year)
year_raw_data = snakemake.input[0]

config = snakemake.config
metadata_file = snakemake.input[1]
pm25_daily_limit = float(config["pm25_daily_limit"])
cities = config["cities"]

# Wczytywanie
print(f"[{year}] Wczytywanie danych surowych z {year_raw_data}...")
raw_df = pd.read_excel(year_raw_data, header=None, decimal=",")
print(f"[{year}] Wczytywanie metadanych z {metadata_file}...")
metadata = pd.read_excel(metadata_file, engine='openpyxl')

# Czyszczenie
print(f"[{year}] Czyszczenie i ujednolicanie danych...")
data = df_gotowy(raw_df_dict={year: raw_df}, metadane=metadata)
print(f"[{year}] Wymiar danych: {data.shape}")

# Analiza
# Średnie miesięczne
monthly_means = srednie_miesieczne(data)
monthly_means.to_csv(snakemake.output[0],index=True)
print(f"[{year}] Średnie miesięczne po stacjach zapisane do: {snakemake.output[0]}")

# Średnie miesięczne pogrupowane po miastach
monthly_means_by_city = srednie_po_miastach(monthly_means)
monthly_means_by_city.to_csv(snakemake.output[1],index=True)
print(f"[{year}] Średnie miesięczne po miastach zapisane do: {snakemake.output[1]}")

# Liczba dni z przekroczeniem po stacjach
exceedance_days_by_station = dni_przekroczenia_normy(df_pomiary=data,norma_dobowa = pm25_daily_limit, years=[year])
exceedance_days_by_station.to_csv(snakemake.output[2],index=True)
print(f"[{year}] Liczba dni z przekroczeniem po stacjach zapisana do: {snakemake.output[2]}")

# Stacje z największym i najmniejszym przekroczeniem
top3_max_min_stations, days_over_norm = wybierz_stacje_max_min(ile_dni_wiecej_normy = exceedance_days_by_station, rok = year)
days_over_norm.to_csv(snakemake.output[3],index=True)
print(f"[{year}] Stacje z największym i najmniejszym przekroczeniem zapisane do: {snakemake.output[3]}")

# Liczba dni z przekroczeniem po województwach
exceedance_days_by_voivodeship = overnorm_by_voivodeship(df = data, metadata=metadata, daily_norm=pm25_daily_limit, years=[year])
exceedance_days_by_voivodeship.to_csv(snakemake.output[4],index=True)
print(f"[{year}] Liczba dni z przekroczeniem po województwach zapisana do: {snakemake.output[4]}")

# Wizualizacja
# wykres_porownanie_miast(srednie_miast=monthly_means_by_city,lata = [year],miasta=cities)

print(f"[{year}] Analiza zakończona.\n\n\n")
