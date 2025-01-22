import requests

url = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/accession/NM_021803.4"
headers = {"Accept": "application/json", "api-key": "d408307bb4c1114872d63b8271f5418f9707"}

response = requests.get(url, headers=headers)

if response.status_code == 200:
    print("API Key is valid. Here's the response:")
    print(response.json())
else:
    print(f"Error: {response.status_code}, {response.text}")
