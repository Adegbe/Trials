import requests
import streamlit as st

# Streamlit App Title
st.title("ClinicalTrials.gov API Viewer")
st.markdown("Search for clinical trials and view results interactively.")

# Sidebar for Input
st.sidebar.header("Search Parameters")
condition = st.sidebar.text_input("Condition", "breast cancer AND stroke")
page_size = st.sidebar.slider("Number of Results", 1, 50, 10)

# API Base URL
API_BASE = "https://clinicaltrials.gov/api/v2/studies"

# API Parameters
params = {
    "query.cond": condition,
    "pageSize": page_size,
    "format": "json"
}

# Button to Trigger API Request
if st.sidebar.button("Search"):
    st.write(f"### Results for: `{condition}`")
    
    # API Request
    response = requests.get(API_BASE, params=params)

    if response.status_code == 200:
        data = response.json()
        studies = data.get("studies", [])
        
        if studies:
            # Display Results
            for study in studies:
                nct_id = study["protocolSection"]["identificationModule"].get("nctId")
                title = study["protocolSection"]["identificationModule"].get("briefTitle")
                status = study["protocolSection"]["statusModule"].get("overallStatus")
                
                # Display Study Details
                st.write(f"**NCT ID**: {nct_id}")
                st.write(f"**Title**: {title}")
                st.write(f"**Status**: {status}")
                st.write("---")
        else:
            st.write("No studies found.")
    else:
        st.error(f"Error: {response.status_code}, {response.text}")
