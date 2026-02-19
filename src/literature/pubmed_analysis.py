#TODO opis skryptu
from typing import Any

from Bio import Entrez
import time
import matplotlib.pyplot as plt
import pandas as pd


def search_for_papers(query: str, year: int, entrez_email: str, entrez_api_key: str = None, retmax = 1000) -> pd.DataFrame:
    """
    Wyszukuje w PubMed artykuły pasujące do zapytania i roku,
    pobiera ich metadane i zapisuje wyniki w pliku CSV.

    Parameters
    ----------
    query : str
        Zapytanie w wyszukiwarce PubMed, np. "Lung AND Cancer" albo "PM2.5 AND pollution".
    year : int
        Rok publikacji, aby filtrować wyniki.
    entrez_email : str
        Email addres użytkownika wymagany przez NCBI Entrez API.
    entrez_api_key : str, optional
        NCBI API key (zwiększa limit szybkości z 3 do 10 żądań/sek.).
    retmax: int, optional
        Limit zwracanych wyników dla jednego zapytania. Domyślnie - 1000.

    Returns
    -------
    pd.DataFrame
        DataFrame zawierający winiki wyszukiwania z kolumnami: PMID, title, year, journal, authors, abstract.
    """
    Entrez.email = entrez_email
    if entrez_api_key:
        Entrez.api_key = entrez_api_key

    # Budowanie zapytania z filtrem roku
    full_query = f"({query}) AND {year}[PDAT]"
    print(f"[pubmed_year][{year}] Szukam zapytanie : {full_query}")

    # Wyszukiwanie - pobieramy listę PMID
    stream = Entrez.esearch(db="pubmed", term=full_query, retmax=retmax)
    search_results = Entrez.read(stream)
    stream.close()

    pmid_list = search_results["IdList"]
    print(f"[pubmed_year][{year}] Znaleziono {len(pmid_list)} artykułów.")

    if not pmid_list:
        print(f"[pubmed_year][{year}] Brak wyników dla zapytania: {full_query}")
        empty = pd.DataFrame(columns=["PMID", "title", "year", "journal", "authors", "abstract"])
        return empty

    # Pobieranie metadanych dla znalezionych PMID
    stream = Entrez.efetch(db="pubmed", id=pmid_list, rettype="xml", retmode="xml")
    records = Entrez.read(stream)
    stream.close()

    rows = []
    for record in records["PubmedArticle"]:
        article = record["MedlineCitation"]["Article"]

        pmid = str(record["MedlineCitation"]["PMID"])

        title = str(article.get("ArticleTitle", ""))

        journal = article.get("Journal", {}).get("Title", "")

        # Rok publikacji
        pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        pub_year = pub_date.get("Year", str(year))

        # Autorzy
        authors_list = article.get("AuthorList", [])
        authors = "; ".join(
            f"{a.get('LastName', '')} {a.get('ForeName', '')}".strip()
            for a in authors_list
            if "LastName" in a
        )

        # Abstrakt - może być podzielony na sekcje
        abstract_texts = article.get("Abstract", {}).get("AbstractText", [])
        if isinstance(abstract_texts, list):
            abstract = " ".join(str(a) for a in abstract_texts)
        else:
            abstract = str(abstract_texts)

        rows.append({
            "PMID": pmid,
            "title": title,
            "year": pub_year,
            "journal": journal,
            "authors": authors,
            "abstract": abstract
        })

    df = pd.DataFrame(rows, columns=["PMID", "title", "year", "journal", "authors", "abstract"])
    print(f"[pubmed_year][{year}] Pomyślnie przeanalizowano zapytanie{query}.")
    return df

def top_journals(year: int, entrez_email: str, entrez_api_key: str = None, top: int = 10, sample_size: int = 10000) -> pd.DataFrame:
    """
    Returns top journals by number of publications in PubMed for a given year.
    Based on a sample of up to sample_size records.

    Parameters
    ----------
    year : int
    entrez_email : str
    entrez_api_key : str, optional
    top : int
        Number of top journals to return.
    sample_size : int
        Max number of records to sample (default 10000).

    Returns
    -------
    pd.DataFrame with columns: journal, article_count
    """
    Entrez.email = entrez_email
    if entrez_api_key:
        Entrez.api_key = entrez_api_key

    print(f"[top_journals] Pobieranie próbki {sample_size} rekordów za rok {year}...")

    stream = Entrez.esearch(db="pubmed", term=f"{year}[PDAT]", retmax=sample_size)
    search_results = Entrez.read(stream)
    stream.close()

    pmid_list = search_results["IdList"]
    total_found = int(search_results["Count"])
    print(f"[top_journals] Znaleziono {total_found} artykułów, pobrano próbkę {len(pmid_list)}.")

    stream = Entrez.efetch(db="pubmed", id=pmid_list, rettype="xml", retmode="xml")
    records = Entrez.read(stream)
    stream.close()

    journals = []
    for record in records["PubmedArticle"]:
        article = record["MedlineCitation"]["Article"]
        journal = article.get("Journal", {}).get("Title", "Unknown")
        journals.append(journal)

    df = pd.Series(journals).value_counts().head(top).reset_index()
    df.columns = ["journal", "article_count"]
    df["year"] = year
    df["sample_size"] = len(pmid_list)
    df["total_in_pubmed"] = total_found

    print(f"[top_journals] Top {top} journals gotowe.")
    return df

def figure_top_journals_per_year(top_journals_per_year_df: pd.DataFrame, output_file_name: str) -> None:
    """
    Draws a barplot of publication counts for top journals.

    Parameters
    ----------
    top_journals_per_year_df : pd.DataFrame
        DataFrame with columns: journal, article_count (result of top_journals())
    output_file_name : str
        Path to save the figure.
    """
    fig, ax = plt.subplots(figsize=(15, 8))
    max_len = 50
    ax.bar(top_journals_per_year_df["journal"], top_journals_per_year_df["article_count"])

    ax.set_xlabel("Czasopismo")
    ax.set_ylabel("Liczba publikacji")
    ax.grid(True)
    ax.set_title(f"Top czasopisma według liczby publikacji ({top_journals_per_year_df['year'].iloc[0]})")
    ax.set_title(f"Analiza na podstawie {top_journals_per_year_df['sample_size'].iloc[0]} publikacji z {top_journals_per_year_df['total_in_pubmed'].iloc[0]}", loc="left")
    ax.set_xticks(range(len(top_journals_per_year_df)))

    labels = [j if len(j) <= max_len else j[:max_len] + "..."
              for j in top_journals_per_year_df["journal"]]
    ax.set_xticklabels(labels, rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(output_file_name)
    plt.close()
    print(f"[figure_top_journals] Zapisano wykres do {output_file_name}")