import streamlit as st
import requests
import pandas as pd

# Set page config
st.set_page_config(page_title="Clinical Trials Search", layout="wide")

# App title
st.title("🔬 Clinical Trials Search Dashboard")

# Sidebar for inputs
with st.sidebar:
    st.header("Search Filters")
    condition = st.text_input("Condition / Disease", "breast cancer")
    recruitment_status = st.selectbox("Recruitment Status", ["", "Recruiting", "Completed", "Terminated"])
    country = st.text_input("Country")

# Construct API parameters
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
            sponsor_mod = protocol.get("sponsorCollaboratorsModule", {})
            design_mod = protocol.get("designModule", {})
            eligibility_mod = protocol.get("eligibilityModule", {})
            refs_mod = protocol.get("referencesModule", {})

            trial_id = id_mod.get("nctId", "N/A")
            title = id_mod.get("briefTitle", "N/A")
            status = status_mod.get("overallStatus", "N/A")
            trial_country = contact_mod.get("locations", [{}])[0].get("country", "N/A")
            start_date = status_mod.get("startDateStruct", {}).get("date", "N/A")
            completion_date = status_mod.get("completionDateStruct", {}).get("date", "N/A")
            enrollment = design_mod.get("enrollmentInfo", {}).get("count", "N/A")
            sponsor = sponsor_mod.get("leadSponsor", {}).get("name", "N/A")
            gender = eligibility_mod.get("sex", "N/A")
            min_age = eligibility_mod.get("minimumAge", "N/A")
            max_age = eligibility_mod.get("maximumAge", "N/A")
            pub_list = refs_mod.get("referenceList", [])
            pub_url = f"https://pubmed.ncbi.nlm.nih.gov/{pub_list[0]['PMID']}/" if pub_list else "N/A"

            if recruitment_status and recruitment_status.lower() not in status.lower():
                continue
            if country and country.lower() not in trial_country.lower():
                continue

            results.append({
                "NCT ID": trial_id,
                "Title": title,
                "Status": status,
                "Start Date": start_date,
                "Completion Date": completion_date,
                "Enrollment": enrollment,
                "Sponsor": sponsor,
                "Gender": gender,
                "Min Age": min_age,
                "Max Age": max_age,
                "Country": trial_country,
                "Study Publication": pub_url
            })

        if results:
            df = pd.DataFrame(results)
            st.success(f"{len(df)} trial(s) found.")
            st.dataframe(df)

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("Download CSV", data=csv, file_name="clinical_trials.csv", mime="text/csv")
        else:
            st.warning("No trials matched your filters.")
    else:
        st.error(f"API Error: {response.status_code}")
else:
    st.info("Please enter a condition to begin your search.")
