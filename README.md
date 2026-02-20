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

