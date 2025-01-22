import requests
import streamlit as st
from Bio import Entrez

# Streamlit App Title
st.title("NCBI Database Query Tool")
st.markdown("Use this app to query NCBI databases interactively.")

# Sidebar for Input
st.sidebar.header("NCBI Query Parameters")

# Database Selection
database_options = ["nucleotide", "protein", "genome", "pubmed", "taxonomy"]
selected_database = st.sidebar.selectbox("Select a Database", database_options)

# Accession Number Input
accession_number = st.sidebar.text_input("Enter Accession Number", "NM_021803.4")

# Button to Trigger Query
if st.sidebar.button("Query"):
    if selected_database and accession_number:
        st.write(f"### Querying `{selected_database}` database for accession `{accession_number}`")
        
        try:
            # Query NCBI using Entrez
            Entrez.email = "adegbesamson@gmail.com.com"  # Replace with your email
            handle = Entrez.efetch(db=selected_database, id=accession_number, rettype="gb", retmode="text")
            result = handle.read()
            handle.close()

            # Display the Result
            st.text_area("Query Result", result, height=400)
        except Exception as e:
            st.error(f"An error occurred: {e}")
    else:
        st.warning("Please select a database and enter an accession number.")


