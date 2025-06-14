import streamlit as st
import requests
import pandas as pd
from Bio import Entrez

# Set email for PubMed Entrez API
Entrez.email = "adegbesamson@gmail.com"

# Streamlit page config
st.set_page_config(page_title="Clinical Trials + PubMed Viewer", layout="wide")
st.title("🔬 Clinical Trials Dashboard with PubMed Summaries")

# Sidebar search filters
with st.sidebar:
    st.header("Search Filters")
    condition = st.text_input("Condition / Disease", "breast cancer")
    recruitment_status = st.selectbox("Recruitment Status", ["", "Recruiting", "Completed", "Terminated"])
    country = st.text_input("Country")

# Function to fetch PubMed abstract summary
def get_pubmed_summary(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        title = article.get('ArticleTitle', 'No Title')
        abstract = ' '.join(article['Abstract']['AbstractText']) if 'Abstract' in article else 'No Abstract'
        return f"{title}\n\n{abstract}"
    except Exception:
        return "No PubMed summary available."

# API base for ClinicalTrials.gov
API_BASE = "https://clinicaltrials.gov/api/v2/studies"
params = {
    "query.cond": condition,
    "pageSize": 20,
    "format": "json"
}

if condition:
    response = requests.get(API_BASE, params=params)

    if response.status_code == 200:
        studies = response.json().get("studies", [])
        results = []

        for study in studies:
            protocol = study.get("protocolSection", {})
            id_mod = protocol.get("identificationModule", {})
            status_mod = protocol.get("statusModule", {})
            contact_mod = protocol.get("contactsLocationsModule", {})
            reference_mod = study.get("resultsReferenceModule", {})

            trial_id = id_mod.get("nctId", "N/A")
            title = id_mod.get("briefTitle", "N/A")
            status = status_mod.get("overallStatus", "N/A")
            trial_country = contact_mod.get("locations", [{}])[0].get("country", "N/A")

            if recruitment_status and recruitment_status.lower() not in status.lower():
                continue
            if country and country.lower() not in trial_country.lower():
                continue

            pmid = reference_mod.get("referenceList", [{}])[0].get("pmid", None)
            pubmed_summary = get_pubmed_summary(pmid) if pmid else "No linked PubMed article"

            results.append({
                "NCT ID": trial_id,
                "Title": title,
                "Status": status,
                "Country": trial_country,
                "PubMed Summary": pubmed_summary
            })

        if results:
            df = pd.DataFrame(results)
            st.success(f"{len(df)} trial(s) found.")
            st.dataframe(df)

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("Download CSV", data=csv, file_name="clinical_trials_pubmed.csv", mime="text/csv")
        else:
            st.warning("No trials matched your filters.")
    else:
        st.error(f"API Error: {response.status_code}")
else:
    st.info("Please enter a condition to begin your search.")
