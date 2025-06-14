
import streamlit as st
import requests
import pandas as pd

# Streamlit page config
st.set_page_config(page_title="Clinical Trials & Publications", layout="wide")

# Page title
st.title("🔬 Clinical Trials Explorer with PubMed Integration")

# Sidebar input
with st.sidebar:
    st.header("Search Filters")
    condition = st.text_input("Condition / Disease", "breast cancer")
    recruitment_status = st.selectbox("Recruitment Status", ["", "Recruiting", "Completed", "Terminated"])
    country = st.text_input("Country (e.g., United States, China)")

# API base
API_BASE = "https://clinicaltrials.gov/api/v2/studies"
params = {
    "query.cond": condition,
    "pageSize": 20,
    "format": "json"
}

# Function to fetch PubMed links via Entrez eSearch
def get_pubmed_link(condition):
    from Bio import Entrez
    Entrez.email = "adegbesamson@gmail.com"
    try:
        handle = Entrez.esearch(db="pubmed", term=condition, retmax=1)
        record = Entrez.read(handle)
        ids = record.get("IdList", [])
        if ids:
            return f"https://pubmed.ncbi.nlm.nih.gov/{ids[0]}"
    except Exception:
        return ""
    return ""

# Run search if condition is provided
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
            design_mod = protocol.get("designModule", {})
            desc_mod = protocol.get("descriptionModule", {})
            sponsor_mod = protocol.get("sponsorCollaboratorsModule", {})
            eligibility_mod = protocol.get("eligibilityModule", {})

            trial_id = id_mod.get("nctId", "N/A")
            title = id_mod.get("briefTitle", "N/A")
            status = status_mod.get("overallStatus", "N/A")
            trial_country = contact_mod.get("locations", [{}])[0].get("country", "N/A")
            start_date = status_mod.get("startDateStruct", {}).get("date", "N/A")
            study_type = design_mod.get("studyType", "N/A")
            phase = design_mod.get("phase", "N/A")
            sponsor = sponsor_mod.get("leadSponsor", {}).get("name", "N/A")
            eligibility = eligibility_mod.get("gender", "N/A") + ", " + eligibility_mod.get("minimumAge", "N/A") + " to " + eligibility_mod.get("maximumAge", "N/A")
            summary = desc_mod.get("briefSummary", "N/A")
            pubmed_link = get_pubmed_link(condition)

            if recruitment_status and recruitment_status.lower() not in status.lower():
                continue
            if country and country.lower() not in trial_country.lower():
                continue

            results.append({
                "NCT ID": trial_id,
                "Title": title,
                "Status": status,
                "Country": trial_country,
                "Start Date": start_date,
                "Study Type": study_type,
                "Study Phase": phase,
                "Sponsor": sponsor,
                "Eligibility": eligibility,
                "Description": summary,
                "PubMed Link": pubmed_link
            })

        if results:
            df = pd.DataFrame(results)
            st.success(f"{len(df)} trial(s) found.")
            st.dataframe(df)

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download CSV", data=csv, file_name="clinical_trials_full.csv", mime="text/csv")
        else:
            st.warning("No trials matched your filters.")
    else:
        st.error(f"API Error: {response.status_code}")
else:
    st.info("Please enter a condition to begin your search.")
