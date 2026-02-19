#TODO opis skryptu
from Bio import Entrez
import time
import matplotlib.pyplot as plt
import pandas as pd

def search_for_papers(query: str, year: int, entrez_email: str, output_file: str, entrez_api_key: str = None, retmax = 1000) -> pd.DataFrame:
    """
    Searches PubMed for articles matching query and year,
    fetches their metadata and saves results to a CSV file.

    Parameters
    ----------
    query : str
        PubMed search query, e.g. "Lung AND Cancer" or "PM2.5 AND pollution"
    year : int
        Publication year to filter results.
    entrez_email : str
        Email address required by NCBI Entrez API.
    output_file : str
        Path to output CSV file.
    entrez_api_key : str, optional
        NCBI API key (increases rate limit from 3 to 10 requests/sec).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: PMID, title, year, journal, authors, abstract.
        Also saved to output_file as CSV.
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
        empty.to_csv(output_file, index=False)
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
    df.to_csv(output_file, index=False)
    print(f"[pubmed_year][{year}] Zapisano {len(df)} rekordów do {output_file}")
    return df

def summary_by_year(query, year, entrez_email, entrez_api_key = None):
    pass

def top_journals(year, entrez_email,entrez_api_key = None, top = 10):
    pass
def figure_papers_per_year(output_file_name):
    pass