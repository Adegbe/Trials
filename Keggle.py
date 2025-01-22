import os
from kaggle.api.kaggle_api_extended import KaggleApi

# Authenticate using your kaggle.json file
api = KaggleApi()
api.authenticate()

# Example: Download a dataset
dataset = 'zillow/zecon'  # Replace with the dataset name
path = 'datasets/'  # Directory to download the dataset
api.dataset_download_files(dataset, path=path, unzip=True)

print(f"Dataset {dataset} downloaded to {path}")
