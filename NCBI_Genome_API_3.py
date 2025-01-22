from Bio import Entrez

Entrez.email = "adegbesamson@gmail.com"

# Fetch gene data for the accession number NM_021803.4
handle = Entrez.efetch(db="nucleotide", id="NM_021803.4", rettype="gb", retmode="text")
gene_data = handle.read()
handle.close()

# Print the data
print(gene_data)
