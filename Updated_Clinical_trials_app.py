import streamlit as st
import requests
import pandas as pd

# Streamlit page config
st.set_page_config(page_title="Clinical Trials Explorer", layout="wide")

# Page title
st.title("🔬 Clinical Trials Explorer")

# Sidebar input
with st.sidebar:
    st.header("Search Filters")
    condition = st.text_input("Condition / Disease", "breast cancer")
    recruitment_status = st.selectbox("Recruitment Status", ["All", "Recruiting", "Completed", "Terminated"])
    country = st.text_input("Country (e.g., United States, China)")

# API base
API_BASE = "https://clinicaltrials.gov/api/v2/studies"

# Run search if condition is provided
if condition:
    params = {
        "query.cond": condition,
        "pageSize": 20,
        "format": "json"
    }
    
    with st.spinner("Fetching clinical trials data..."):
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
            eligibility = f"{eligibility_mod.get('gender', 'N/A')}, {eligibility_mod.get('minimumAge', 'N/A')} to {eligibility_mod.get('maximumAge', 'N/A')}"
            summary = desc_mod.get("briefSummary", "N/A")

            # Apply filters
            if recruitment_status != "All" and recruitment_status.lower() not in status.lower():
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
                "Phase": phase,
                "Sponsor": sponsor,
                "Eligibility": eligibility,
                "Description": summary
            })

        if results:
            df = pd.DataFrame(results)
            st.success(f"✅ {len(df)} trial(s) found matching your criteria.")
            
            # Show all columns in the table
            st.dataframe(
                df,
                use_container_width=True,
                column_config={
                    "NCT ID": st.column_config.LinkColumn("NCT ID", display_text="View Trial")
                },
                hide_index=True
            )

            # Download options
            st.markdown("---")
            st.download_button(
                "📥 Download CSV",
                data=df.to_csv(index=False).encode("utf-8"),
                file_name="clinical_trials.csv",
                mime="text/csv"
            )
        else:
            st.warning("⚠️ No trials matched your filters. Try broadening your search criteria.")
    else:
        st.error(f"⚠️ API Error: {response.status_code} - {response.text}")
else:
    st.info("ℹ️ Please enter a condition to begin your search.")
