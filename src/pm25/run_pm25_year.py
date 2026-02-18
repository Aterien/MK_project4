import sys
import os
import pandas as pd

from wczytaj_wyczysc import *
from analiza import *
from wizualizacja import *

# Umożliwia import z src/pm25/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

year = int(snakemake.wildcards.year)
year_raw_data = snakemake.input[0]

config = snakemake.config
metadata_file = snakemake.input[1]
pm25_daily_limit = float(config["pm25_daily_limit"])
cities = config["cities"]

# Wczytywanie
raw_df = pd.read_excel(year_raw_data, header=None, decimal=",")
metadata = pd.read_excel(metadata_file,engine='openpyxl')

# Czyszczenie
data = df_gotowy(raw_df_dict={year:raw_df},metadane=metadata)

# Analiza
# Średnie miesięczne
monthly_means = srednie_miesieczne(data)
monthly_means.to_csv(snakemake.output[0],index=True)

# Średnie miesięczne pogrupowane po miastach
monthly_means_by_city = srednie_po_miastach(monthly_means)
monthly_means_by_city.to_csv(snakemake.output[1],index=True)

# Liczba dni z przekroczeniem po stacjach
exceedance_days_by_station = dni_przekroczenia_normy(df_pomiary=data,norma_dobowa = pm25_daily_limit, years=[year])
exceedance_days_by_station.to_csv(snakemake.output[2],index=True)

# Stacje z największym i najmniejszym przekroczeniem
top3_max_min_stations, days_over_norm = wybierz_stacje_max_min(ile_dni_wiecej_normy = exceedance_days_by_station, rok = year)
days_over_norm.to_csv(snakemake.output[3],index=True)

# Liczba dni z przekroczeniem po województwach
exceedance_days_by_voivodeship = overnorm_by_voivodeship(df = data, metadata=metadata, daily_norm=pm25_daily_limit, years=[year])
exceedance_days_by_voivodeship.to_csv(snakemake.output[4],index=True)

# Wizualizacja
# wykres_porownanie_miast(srednie_miast=monthly_means_by_city,lata = [year],miasta=cities)



