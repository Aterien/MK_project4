import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def rysuj_srednie_kat_wwa(kat_srednie, wwa_srednie):
    kat_2014 = kat_srednie.xs(2014, level='Rok')
    kat_2024 = kat_srednie.xs(2024, level='Rok')
    wwa_2014 = wwa_srednie.xs(2014, level='Rok')
    wwa_2024 = wwa_srednie.xs(2024, level='Rok')

    plt.figure(figsize=(12,8))
    plt.plot(kat_2014.index, kat_2014.values, marker='*', color = "pink", label='Średnie miesięczne wartości- Katowice 2014')
    plt.plot(kat_2024.index, kat_2024.values, marker='*', color = "plum", label='Średnie miesięczne wartości- Katowice 2024')
    plt.plot(wwa_2014.index, wwa_2014.values, marker='*', color = "purple", label='Średnie miesięczne wartości- Warszawa 2014')
    plt.plot(wwa_2024.index, wwa_2024.values, marker='*', color = "navy", label='Średnie miesięczne wartości- Warszawa 2024')
    plt.xlabel('Miesiąc')
    plt.ylabel('Średnia wartość PM2.5')
    plt.title('Średnie miesięczne stężenie PM2.5 w Katowicach i Warszawie')
    plt.xticks(range(1,13))
    plt.grid(True)
    plt.legend()
    plt.show()

def rysuj_srednie_wszystkie(df_pomiary):
    miesieczne_srednie = df_pomiary.groupby([df_pomiary.index.year, df_pomiary.index.month]).mean()
    srednie_po_stacjach = miesieczne_srednie.groupby(level="Miejscowość", axis=1).mean()
    years = [2014, 2019, 2024]
    miejscowosci = srednie_po_stacjach.columns.to_list()

    fig, axes = plt.subplots((len(miejscowosci)) // 3, 3, figsize=(15, 20))
    fig.suptitle("Średnie miesięczne stężenie PM2.5 we wszystkich miejscowościach",fontsize=20)
    for nr in range(len(miejscowosci)):
        df_heat = pd.concat([srednie_po_stacjach.loc[y, miejscowosci[nr]] for y in years],axis=1).T
        df_heat.index = years
        df_heat.columns.names = ['Miesiac']
        df_heat.index.names = ['Rok']
        y = nr//3
        x = nr%3
        hm = axes[y][x].imshow(df_heat, aspect='auto',vmin=0, vmax=80)
        axes[y][x].set_title(miejscowosci[nr])
        axes[y][x].set_xticks(range(12))
        axes[y][x].set_xticklabels(range(1, 13))
        axes[y][x].set_yticks(range(len(years)))
        axes[y][x].set_yticklabels(years)
        axes[y][x].set_xlabel("Miesiąc")
        cbar = fig.colorbar(hm, ax=axes[y][x], fraction=0.046, pad=0.04)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()

def rysuj_dni_przekroczen(ile_dni_wybrane_stacje, wybrane_stacje, norma_dobowa):
    years = [2014, 2019, 2024]
    x = np.arange(6)
    width = 0.2
    plt.bar(x-width, ile_dni_wybrane_stacje.loc[2014],width,color='red')
    plt.bar(x, ile_dni_wybrane_stacje.loc[2019],width, color='green')
    plt.bar(x+width, ile_dni_wybrane_stacje.loc[2024],width,color='blue')
    plt.xticks(x, [stacja[0] for stacja in wybrane_stacje],rotation=30)
    plt.grid()
    plt.legend(years, loc='lower left')
    plt.title(f"Liczbą dni z przekroczeniem normy dobowej = {norma_dobowa} µg/m³ stężenia PM2.5  dla 6 stacji")
    plt.show()
