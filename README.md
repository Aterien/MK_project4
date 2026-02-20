# ZTP: projekt 4
### Mykyta Khrabust

---

## Instrukcja jak uruchomić pipeline zadania 4
1) Ustawienia w configu
   - Config jest w workflows/config.yaml
   - Aby przeprowadzić analizę pomiarów PM2,5 dla określonego roku, należy dodać do configu (w years:) treść
      ```yaml
      {year}: 
        archive_id: "{archive_id}"
        filename: "{filename}"
      ```
     Gdzie:
     - {year} to rok, 
     - {archive_id} to id pobieranego archiwum na stronie 
     - {filename} to nazwa pliku .xlsx z danymi o pomiarach stężenia PM2.5 w pobranym archiwum (zwykle jest to {year}_PM25_1g.xlsx) na przykład:
     ```yaml
      2024:
        archive_id: "582"
        filename: "2024_PM25_1g.xlsx"
      ```
   - Aby dokonać mini-przegłądu literatury za wskazane lata, należy podać email (wymagany podczas wyszukiwań w bazach NCBI) oraz api key (opcjonalnie)
2) Uruchomienie Snakemake w terminalu
   Pipeline z wybranymi parametrami uruchamia się w workflows/ komendą

  ```commandline
    snakemake -s workflows/Snakefile_task4
  
  ```
### Jak weryfikuję incrementalność piplenie'a
Sprawdzam, czy wyniki nie są ponownie liczone, jeśli nie jest to konieczne, głównie na podstawie tabeli liczby wykonań reguł generowanej przez snakemake po uruchomienu pipeline'a w trybie suchym lub pełnym. 
Jeśli ani kod, ani dane w configu nie uległy zmianiom, snakemake zwróci:
```
    Building DAG of jobs...
    Nothing to be done (all requested files are present and up to date).
```
Jeśli jednak, zmienimy jeden rok na inny (na przykład 2024 -> 2019),
to snakemake wypisze:
```
    Building DAG of jobs...
    Job stats:
    job                        count
    -----------------------  -------
    all                            1
    download_year                  1
    generate_all_year_plots        1
    generate_report                1
    pm25_year                      1
    pubmed_year                    1
    total                          6

```
Zgadza się to z intuicją, po tej zmianie pipeline bedę musiał jedynie:
- pobrać dane dle tego roku (download_year)
- przeprowadzić analizę PM2.5 dla tego roku (pm25_year)
- przeprowadzić przegłąd literatury dla tego roku (pubmed_year)
- ponownie policzyć statystyki dla wielu lat (generate_all_year_plots)
- wygenerować nowy raport (generate_report )
