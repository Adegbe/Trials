import requests

url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/9606/dataset_report"
response = requests.get(url)
print(response.json())
