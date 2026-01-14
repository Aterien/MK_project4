import pandas as pd

def srednie_miesieczne(df_pomiary):
    miesieczne_srednie = df_pomiary.groupby([df_pomiary.index.year, df_pomiary.index.month]).mean()
    miesieczne_srednie.index.names = ['Rok','Miesiąc']
    return miesieczne_srednie

def srednie_dla_miast(miesieczne_srednie, miasta):
    wynik = {}
    for miasto in miasta:
        df_miasto = miesieczne_srednie.loc[:, miesieczne_srednie.columns.get_level_values("Miejscowość") == miasto]
        wynik[miasto] = df_miasto.mean(axis=1)
    return wynik

def srednie_po_stacjach(miesieczne_srednie):
    return miesieczne_srednie.groupby(level="Miejscowość", axis=1).mean()

def dni_przekroczenia_normy(df_pomiary, norma_dobowa, years):
    dzienne_srednie = df_pomiary.groupby([df_pomiary.index.year, df_pomiary.index.month, df_pomiary.index.day]).mean()
    dzienne_srednie.index.names = ['Rok','Miesiąc','Dzień']
    ile_dni = pd.DataFrame(index=years, columns=dzienne_srednie.columns)
    for year in years:
        df_year = dzienne_srednie.loc[year]
        ile_dni.loc[year] = (df_year > norma_dobowa).sum()
    return ile_dni

def wybierz_stacje_max_min(ile_dni_wiecej_normy, rok, ile_maxmin=3):
    max3 = ile_dni_wiecej_normy.loc[rok].sort_values(ascending=False).head(ile_maxmin)
    min3 = ile_dni_wiecej_normy.loc[rok].sort_values(ascending=False).tail(ile_maxmin)
    wybrane_stacje = max3.index.tolist() + min3.index.tolist()
    return wybrane_stacje, ile_dni_wiecej_normy[wybrane_stacje]
