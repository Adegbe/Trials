import streamlit as st
import requests

st.title("Clinical Trials Data Viewer")

API_BASE = "https://clinicaltrials.gov/api/v2/studies"

params = {
    "query.cond": "breast cancer AND stroke",
    "pageSize": 10,
    "format": "json"
}

response = requests.get(API_BASE, params=params)

if response.status_code == 200:
    data = response.json()
    studies = data.get("studies", [])
    for study in studies:
        nct_id = study["protocolSection"]["identificationModule"].get("nctId")
        title = study["protocolSection"]["identificationModule"].get("briefTitle")
        status = study["protocolSection"]["statusModule"].get("overallStatus")
        st.markdown(f"**NCT ID:** {nct_id}  \n**Title:** {title}  \n**Status:** {status}")
        st.markdown("---")
else:
    st.error(f"Error: {response.status_code}, {response.text}")
