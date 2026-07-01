#!/usr/bin/env python3
"""
Download Harmonized Psoriasis Transcriptomics Dataset from Zenodo
DOI: 10.5281/zenodo.4009497
"""

import os
import requests
import json
from pathlib import Path

def download_zenodo_dataset():
    """
    Download the harmonized psoriasis transcriptomics dataset from Zenodo.
    Dataset: https://zenodo.org/record/4009497
    """
    
    # Zenodo API endpoint
    zenodo_id = "4009497"
    zenodo_url = f"https://zenodo.org/api/records/{zenodo_id}"
    
    # Create data directory
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    
    print(f"Downloading dataset from Zenodo (Record ID: {zenodo_id})...")
    
    try:
        # Get metadata
        response = requests.get(zenodo_url)
        response.raise_for_status()
        metadata = response.json()
        
        print("\n=== Dataset Information ===")
        print(f"Title: {metadata['metadata']['title']}")
        print(f"Description: {metadata['metadata']['description'][:200]}...")
        print(f"DOI: {metadata['metadata']['doi']}")
        print(f"Publication Date: {metadata['metadata']['publication_date']}")
        
        # Save metadata
        metadata_file = data_dir / "metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"\nMetadata saved to: {metadata_file}")
        
        # Download files
        files = metadata['files']
        print(f"\nFound {len(files)} files to download:")
        print("-" * 50)
        
        for file in files:
            filename = file['key']
            file_size = file['size'] / (1024**2)  # Convert to MB
            download_url = file['links']['self']
            
            print(f"• {filename} ({file_size:.2f} MB)")
            
            filepath = data_dir / filename
            
            # Download with progress
            print(f"  Downloading: {filename}...")
            try:
                response = requests.get(download_url, stream=True)
                response.raise_for_status()
                
                total_size = int(response.headers.get('content-length', 0))
                downloaded = 0
                
                with open(filepath, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            if total_size > 0:
                                progress = (downloaded / total_size) * 100
                                print(f"  Progress: {progress:.1f}%", end='\r')
                
                print(f"  ✓ Downloaded: {filename}")
                
            except Exception as e:
                print(f"  ✗ Error downloading {filename}: {str(e)}")
        
        print("\n" + "=" * 50)
        print("Download complete!")
        print(f"Files saved to: {data_dir.absolute()}")
        
    except requests.exceptions.RequestException as e:
        print(f"Error connecting to Zenodo: {str(e)}")
        print("Please check your internet connection and try again.")
    except Exception as e:
        print(f"Unexpected error: {str(e)}")

def list_downloaded_files():
    """List all downloaded files"""
    data_dir = Path("data")
    
    if not data_dir.exists():
        print("No data directory found. Please run download first.")
        return
    
    files = list(data_dir.glob("*"))
    
    print("\n=== Downloaded Files ===")
    for file in files:
        if file.is_file():
            size = file.stat().st_size / (1024**2)
            print(f"• {file.name} ({size:.2f} MB)")

if __name__ == "__main__":
    print("=" * 60)
    print("Psoriasis Transcriptomics Dataset Downloader")
    print("=" * 60)
    
    download_zenodo_dataset()
    list_downloaded_files()