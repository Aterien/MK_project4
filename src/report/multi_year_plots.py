"""
multi_year_plots.py

Funkcje do wczytywania i łączenia wyników analizy PM2.5 z wielu lat,
oraz rysowania wykresów porównawczych.
"""

import sys
import os
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../pm25"))
from wczytaj_wyczysc import wspolne_stacje


def concat_monthly_means_by_city(years: list[int], results_dir: str = "results/pm25") -> pd.DataFrame:
    """
    Wczytuje monthly_means_by_city.csv dla każdego roku, znajduje miasta
    wspólne dla wszystkich lat i zwraca połączony DataFrame.

    Indeks wierszy jest dwupoziomowy: (Rok, Miesiąc).
    Kolumny to miasta obecne we wszystkich latach jednocześnie.

    Parameters
    ----------
    years : list of int
        Lista lat do wczytania.
    results_dir : str
        Ścieżka bazowa do wyników (domyślnie "results/pm25").

    Returns
    -------
    pd.DataFrame
        Połączony DataFrame z przecięciem miast, indeks (Rok, Miesiąc).
    """
    dfs = {}
    for year in years:
        path = os.path.join(results_dir, str(year), "monthly_means_by_city.csv")
        df = pd.read_csv(path, index_col=[0, 1])
        df.index.names = ["Rok", "Miesiąc"]
        dfs[year] = df

    # Przecięcie kolumn (miast) obecnych we wszystkich latach
    common_cities = set(dfs[years[0]].columns)
    for year in years[1:]:
        common_cities &= set(dfs[year].columns)
    common_cities = sorted(common_cities)

    print(f"[concat_monthly_means_by_city] Wspólne miasta ({len(common_cities)}): {common_cities}")

    combined = pd.concat([dfs[year][common_cities] for year in years])
    combined.index.names = ["Rok", "Miesiąc"]
    return combined

def concat_exceedance_days_by_voivodeship(years: list[int], results_dir: str = "results/pm25") -> pd.DataFrame:
    """
    Wczytuje exceedance_days_by_voivodeship.csv dla każdego roku,
    znajduje województwa wspólne dla wszystkich lat
    i zwraca połączony DataFrame.

    Indeks wierszy to lata.
    Kolumny to województwa obecne we wszystkich latach jednocześnie.

    Parameters
    ----------
    years : list of int
        Lista lat do wczytania.
    results_dir : str
        Ścieżka bazowa do wyników (domyślnie "results/pm25").

    Returns
    -------
    pd.DataFrame
        Połączony DataFrame z przecięciem województw.
    """
    dfs = []
    for year in years:
        path = os.path.join(results_dir, str(year), "exceedance_days_by_voivodeship.csv")
        df = pd.read_csv(path, index_col=0)
        df.index = [int(year)]
        df.index.name = "Rok"
        dfs.append(df)

    # Przecięcie kolumn (województw) obecnych we wszystkich latach
    common_voiv = set(dfs[0].columns)
    for df in dfs[1:]:
        common_voiv &= set(df.columns)
    common_voiv = sorted(common_voiv)


    combined = pd.concat([df[common_voiv] for df in dfs])
    return combined
