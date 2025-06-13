import streamlit as st
import requests

st.title("Clinical Trials Data Viewer")

# 👇 User input for search query
search_query = st.text_input("Enter condition (e.g., 'breast cancer')", value="breast cancer AND stroke")

API_BASE = "https://clinicaltrials.gov/api/v2/studies"
params = {
    "query.cond": search_query,
    "pageSize": 10,
    "format": "json"
}

# Fetch data only when input is entered
if search_query:
    response = requests.get(API_BASE, params=params)

    if response.status_code == 200:
        data = response.json()
        studies = data.get("studies", [])
        if studies:
            for study in studies:
                nct_id = study["protocolSection"]["identificationModule"].get("nctId")
                title = study["protocolSection"]["identificationModule"].get("briefTitle")
                status = study["protocolSection"]["statusModule"].get("overallStatus")
                st.markdown(f"**NCT ID:** {nct_id}  \n**Title:** {title}  \n**Status:** {status}")
                st.markdown("---")
        else:
            st.info("No studies found for this query.")
    else:
        st.error(f"Error: {response.status_code}, {response.text}")
