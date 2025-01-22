import requests

API_BASE = "https://clinicaltrials.gov/api/v2/studies"

params = {
    "query.cond": "breast cancer AND stroke",
    "pageSize": 10,
    "format": "json"
}

response = requests.get(API_BASE, params=params)

if response.status_code == 200:
    data = response.json()
    for study in data.get("studies", []):
        nct_id = study["protocolSection"]["identificationModule"].get("nctId")
        title = study["protocolSection"]["identificationModule"].get("briefTitle")
        status = study["protocolSection"]["statusModule"].get("overallStatus")
        print(f"NCT ID: {nct_id}, Title: {title}, Status: {status}")
else:
    print(f"Error: {response.status_code}, {response.text}")
