#TODO opis skryptu
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

def top_journals(year, entrez_email,entrez_api_key = None, top = 10):
    pass
def figure_papers_per_year(output_file_name):
    pass