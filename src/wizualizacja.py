import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def wykres_porownanie_miast(srednie_miast:pd.DataFrame, lata:list[int], miasta:list[str]) -> None:

    plt.figure(figsize=(12,8))

    for miasto in miasta:
        m = srednie_miast[miasto]
        for rok in lata:
            plt.plot(m.xs(rok, level='Rok').index, m.xs(rok, level='Rok').values,
             marker='*', label=f'{miasto} {rok}')

    plt.xlabel('Miesiąc')
    plt.ylabel('Średnia wartość PM2.5')
    plt.title('Średnie miesięczne stężenie PM2.5 w Katowicach i Warszawie')
    plt.xticks(range(1,13))
    plt.grid(True)
    plt.legend()
    plt.show()

def wykres_heatmap_srednie(srednie_po_miejscach:pd.DataFrame, lata:list[int]) -> None:
    miejscowosci = srednie_po_miejscach.columns.to_list()
    fig, axes = plt.subplots((len(miejscowosci)+2)//3, 3, figsize=(15, 20))
    fig.suptitle("Średnie miesięczne stężenie PM2.5 we wszystkich miejscowościach", fontsize=20)
    for nr, miasto in enumerate(miejscowosci):
        df_heat = pd.concat([srednie_po_miejscach.loc[y, miasto] for y in lata], axis=1).T
        df_heat.index = lata
        df_heat.columns.names = ['Miesiac']
        df_heat.index.names = ['Rok']
        y = nr//3
        x = nr%3
        hm = axes[y][x].imshow(df_heat, aspect='auto', vmin=0, vmax=80)
        axes[y][x].set_title(miasto)
        axes[y][x].set_xticks(range(12))
        axes[y][x].set_xticklabels(range(1,13))
        axes[y][x].set_yticks(range(len(lata)))
        axes[y][x].set_yticklabels(lata)
        axes[y][x].set_xlabel("Miesiąc")
        axes[y][x].set_ylabel("Rok")
        fig.colorbar(hm, ax=axes[y][x], fraction=0.046, pad=0.04)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()

def wykres_przekroczenia(ile_dni_wybrane_stacje:pd.DataFrame, wybrane_stacje:list, lata:list[int], norma_dobowa:float) -> None:
    x = np.arange(len(wybrane_stacje))
    width = 0.2
    plt.figure(figsize=(10,6))
    plt.bar(x-width, ile_dni_wybrane_stacje.loc[lata[0]], width, color='red', label=lata[0])
    plt.bar(x, ile_dni_wybrane_stacje.loc[lata[1]], width, color='green', label=lata[1])
    plt.bar(x+width, ile_dni_wybrane_stacje.loc[lata[2]], width, color='blue', label=lata[2])
    plt.bar(x+2*width, ile_dni_wybrane_stacje.loc[lata[3]], width, color='purple', label=lata[3])
    plt.xticks(x, [stacja[0] for stacja in wybrane_stacje], rotation=30)
    plt.ylabel('Liczba dni z przekroczeniem normy PM2.5')
    plt.xlabel('Stacja')
    plt.title(f"Liczba dni z przekroczeniem normy dobowej = {norma_dobowa} µg/m³")
    plt.grid(True)
