"""
multi_year_plots.py

Funkcje do wczytywania i łączenia wyników analizy PM2.5 i PubMed z wielu lat,
oraz rysowania wykresów porównawczych.
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../pm25"))
from wczytaj_wyczysc import wspolne_stacje


def concat_monthly_means_by_city(years: list[int], input_files: list[str]) -> pd.DataFrame:
    """
    Wczytuje monthly_means_by_city.csv dla każdego roku, znajduje miasta
    wspólne dla wszystkich lat i zwraca połączony DataFrame.

    Indeks wierszy jest dwupoziomowy: (Rok, Miesiąc).
    Kolumny to miasta obecne we wszystkich latach jednocześnie.

    Parameters
    ----------
    years : list of int
        Lista lat.
    input_files : list of str
        Lista ścieżek do plików monthly_means_by_city.csv (po jednym na rok).

    Returns
    -------
    pd.DataFrame
        Połączony DataFrame z przecięciem miast, indeks (Rok, Miesiąc).
    """
    dfs = {}
    for year, path in zip(years, input_files):
        df = pd.read_csv(path, index_col=[0, 1])
        df.index.names = ["Rok", "Miesiąc"]
        dfs[year] = df
        print(f"[concat_monthly_means_by_city] Wczytano {path} ({df.shape})")

    common_cities = set(dfs[years[0]].columns)
    for year in years[1:]:
        common_cities &= set(dfs[year].columns)
    common_cities = sorted(common_cities)

    print(f"[concat_monthly_means_by_city] Wspólne miasta ({len(common_cities)}): {common_cities}")

    combined = pd.concat([dfs[year][common_cities] for year in years])
    combined.index.names = ["Rok", "Miesiąc"]
    return combined


def concat_exceedance_days_by_voivodeship(years: list[int], input_files: list[str]) -> pd.DataFrame:
    """
    Wczytuje exceedance_days_by_voivodeship.csv dla każdego roku,
    znajduje województwa wspólne dla wszystkich lat
    i zwraca połączony DataFrame.

    Parameters
    ----------
    years : list of int
        Lista lat.
    input_files : list of str
        Lista ścieżek do plików exceedance_days_by_voivodeship.csv.

    Returns
    -------
    pd.DataFrame
        Połączony DataFrame z przecięciem województw.
    """
    dfs = []
    for year, path in zip(years, input_files):
        df = pd.read_csv(path, index_col=0)
        df.index = [int(year)]
        df.index.name = "Rok"
        dfs.append(df)
        print(f"[concat_exceedance_days_by_voivodeship] Wczytano {path} ({df.shape})")

    common_voiv = set(dfs[0].columns)
    for df in dfs[1:]:
        common_voiv &= set(df.columns)
    common_voiv = sorted(common_voiv)


    combined = pd.concat([df[common_voiv] for df in dfs])
    return combined


def all_years_query_summary(years: list[int], input_files: list[str]) -> pd.DataFrame:
    """
    Wczytuje summary_by_year.csv dla każdego roku i łączy je w jeden DataFrame.

    Przykład pojedynczej tabeli:
        ,PM2.5 AND Poland,air pollution AND lung,particulate matter AND health
        2015,10,433,1000

    Parameters
    ----------
    years : list of int
        Lista lat.
    input_files : list of str
        Lista ścieżek do plików summary_by_year.csv (po jednym na rok).

    Returns
    -------
    pd.DataFrame
        DataFrame gdzie indeks to lata, kolumny to zapytania,
        wartości to liczba znalezionych artykułów.
    """
    dfs = []
    for year, path in zip(years, input_files):
        df = pd.read_csv(path, index_col=0)
        df.index = [int(year)]
        df.index.name = "Rok"
        dfs.append(df)
        print(f"[all_years_query_summary] Wczytano {path} ({df.shape})")

    common_queries = set(dfs[0].columns)
    for df in dfs[1:]:
        common_queries &= set(df.columns)
    common_queries = sorted(common_queries)

    combined = pd.concat(dfs)[common_queries]
    print(f"[all_years_query_summary] Połączono {len(years)} lat, {len(common_queries)} zapytań.")
    return combined


def year_summary_by_queries_plot(all_years_query_summary_df: pd.DataFrame, output_file_name: str, show: bool = False) -> None:
    """
    Rysuje barplot liczby wyników wyszukiwania w PubMed
    pogrupowanych po latach (oś X) i zapytaniach (kolory słupków).

    Parameters
    ----------
    all_years_query_summary_df : pd.DataFrame
        DataFrame zwrócony przez all_years_query_summary().
        Indeks: lata, kolumny: zapytania, wartości: liczba artykułów.
    output_file_name : str
        Ścieżka do zapisu wykresu.
    show : bool
        True aby otworzyć wykres w oknie matplotlib.
    """
    ax = all_years_query_summary_df.plot(kind="bar", figsize=(14, 7))
    ax.set_xlabel("Rok")
    ax.set_ylabel("Liczba wyników")
    ax.set_title("Liczba publikacji w PubMed według zapytania i roku")
    ax.grid(True, axis="y")

    plt.legend(title="Zapytanie", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=0)
    plt.tight_layout()

    plt.savefig(output_file_name, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()
    print(f"[year_summary_by_queries_plot] Zapisano wykres -> {output_file_name}")