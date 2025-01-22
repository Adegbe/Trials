import pandas as pd
from Bio import Entrez
import urllib.error

def parse_pub_date(pub_date):
    """Extract and format the publication date."""
    year = pub_date.get('Year', 'Unknown')
    month = pub_date.get('Month', '01')
    day = pub_date.get('Day', '01')
    return f"{year}-{month}-{day}"

# Set Entrez email
Entrez.email = 'adegbesamson@gmail.com'

def search_pubmed(query, db="pubmed"):
    """Search PubMed or PMC using the given query."""
    print(f"Processing query in {db}: {query}")

    try:
        # Perform the search
        handle = Entrez.esearch(db=db, term=query, retmax=50)
        record = Entrez.read(handle)
        id_list = record.get('IdList', [])

        print(f"Found {len(id_list)} results for query: {query}")
        if not id_list:
            print("No results found.")
            return pd.DataFrame()  # Return an empty DataFrame

        # Retrieve details for each result
        rows = []
        for pmid in id_list:
            try:
                handle = Entrez.efetch(db=db, id=pmid, retmode='xml')
                records = Entrez.read(handle)

                for record in records['PubmedArticle']:
                    article = record['MedlineCitation']['Article']
                    
                    # Extract details
                    title = article.get('ArticleTitle', 'Title Not Available')
                    abstract = ' '.join(article['Abstract']['AbstractText']) if 'Abstract' in article else 'Abstract Not Available'
                    authors_list = ', '.join(
                        f"{a.get('ForeName', '')} {a.get('LastName', '')}" for a in article.get('AuthorList', [])
                    ) or 'Authors Not Available'
                    journal = article['Journal'].get('Title', 'Journal Not Available')
                    pub_date = parse_pub_date(article['Journal']['JournalIssue']['PubDate'])
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    
                    # Append result
                    rows.append({
                        'PMID': pmid,
                        'Title': title,
                        'Abstract': abstract,
                        'Authors': authors_list,
                        'Journal': journal,
                        'URL': url,
                        'Publication Date': pub_date
                    })
            except Exception as e:
                print(f"Error fetching details for PMID {pmid}: {e}")
        
        return pd.DataFrame(rows)  # Return results as a DataFrame
    except Exception as e:
        print(f"Error processing query: {e}")
        return pd.DataFrame()  # Return an empty DataFrame in case of an error

# Broad queries for Alhaji Olono's example
queries = [
    '"Olono"[All Fields]',  # Broad search for last name
    '"genomic capacity"[All Fields] AND "precision health"[All Fields]',  # Broad keyword search
    '"Africa"[All Fields] AND "genomic"[All Fields]',  # Very broad keyword search
]

# Combine results for all queries from PubMed and PMC
results = pd.DataFrame()
for query in queries:
    # First search in PubMed
    df_pubmed = search_pubmed(query, db="pubmed")
    results = pd.concat([results, df_pubmed], ignore_index=True)

    # Then search in PMC if no results are found in PubMed
    if df_pubmed.empty:
        df_pmc = search_pubmed(query, db="pmc")
        results = pd.concat([results, df_pmc], ignore_index=True)

# Remove duplicates
results.drop_duplicates(subset='PMID', inplace=True)

# Save results to Excel if any results found
if not results.empty:
    output_file = 'Broad_Search_Results.xlsx'
    results.to_excel(output_file, index=False)
    print(f"Results saved to {output_file}")
else:
    print("No results found. Nothing to save.")
