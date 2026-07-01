#!/usr/bin/env python3
"""
Data Preprocessing and Quality Control
Load, clean, and normalize psoriasis transcriptomics data
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import os

class DataPreprocessor:
    """Preprocess psoriasis transcriptomics data"""
    
    def __init__(self, data_dir="data"):
        self.data_dir = Path(data_dir)
        self.expression_data = None
        self.metadata = None
        
    def load_data(self):
        """Load expression matrix and metadata"""
        print("Loading data from Zenodo dataset...")
        
        # Look for CSV or TSV files
        data_files = list(self.data_dir.glob("*.csv")) + list(self.data_dir.glob("*.tsv")) + list(self.data_dir.glob("*.txt"))
        
        if not data_files:
            print("No data files found. Please run download_data.py first.")
            return False
        
        print(f"Found {len(data_files)} data files")
        
        # Load the main expression file
        expr_file = None
        for file in data_files:
            if 'expression' in file.name.lower() or 'data' in file.name.lower():
                expr_file = file
                break
        
        if expr_file is None:
            expr_file = data_files[0]
        
        print(f"Loading expression data from: {expr_file.name}")
        
        try:
            if expr_file.suffix == '.csv':
                self.expression_data = pd.read_csv(expr_file, index_col=0)
            elif expr_file.suffix in ['.tsv', '.txt']:
                self.expression_data = pd.read_csv(expr_file, sep='\t', index_col=0)
            
            print(f"Expression data shape: {self.expression_data.shape}")
            print(f"Genes: {self.expression_data.shape[0]}, Samples: {self.expression_data.shape[1]}")
            
        except Exception as e:
            print(f"Error loading expression file: {str(e)}")
            return False
        
        # Load metadata if available
        metadata_file = self.data_dir / "metadata.json"
        if metadata_file.exists():
            with open(metadata_file, 'r') as f:
                self.metadata = json.load(f)
        
        return True
    
    def quality_control(self):
        """Perform quality control on expression data"""
        print("\n=== Quality Control ===")
        
        if self.expression_data is None:
            print("No expression data loaded")
            return
        
        # Remove rows with all zeros
        non_zero_genes = (self.expression_data > 0).sum(axis=1) > 0
        print(f"Genes with expression: {non_zero_genes.sum()} / {len(self.expression_data)}")
        self.expression_data = self.expression_data[non_zero_genes]
        
        # Remove columns with all zeros
        non_zero_samples = (self.expression_data > 0).sum(axis=0) > 0
        print(f"Samples with expression: {non_zero_samples.sum()} / {len(self.expression_data.columns)}")
        self.expression_data = self.expression_data[non_zero_samples]
        
        # Calculate statistics
        print("\n--- Expression Statistics ---")
        print(f"Mean expression: {self.expression_data.values.mean():.4f}")
        print(f"Median expression: {np.median(self.expression_data.values):.4f}")
        print(f"Min expression: {self.expression_data.values.min():.4f}")
        print(f"Max expression: {self.expression_data.values.max():.4f}")
        
        return self.expression_data
    
    def normalize_data(self, method='log2'):
        """Normalize expression data"""
        print(f"\n=== Normalization ({method}) ===")
        
        if self.expression_data is None:
            print("No expression data loaded")
            return
        
        if method == 'log2':
            self.expression_data = np.log2(self.expression_data + 1)
            print("Applied log2 normalization")
        elif method == 'zscore':
            self.expression_data = (self.expression_data - self.expression_data.mean()) / self.expression_data.std()
            print("Applied Z-score normalization")
        
        print(f"Normalized data shape: {self.expression_data.shape}")
        return self.expression_data
    
    def check_orm_genes(self):
        """Check for ORM1 and ORM2 genes"""
        print("\n=== Checking ORM1/ORM2 Presence ===")
        
        if self.expression_data is None:
            print("No expression data loaded")
            return
        
        genes = self.expression_data.index.tolist()
        orm_genes = [g for g in genes if 'ORM' in g.upper()]
        
        print(f"Total genes: {len(genes)}")
        print(f"ORM genes found: {len(orm_genes)}")
        
        if orm_genes:
            print("\nORM genes detected:")
            for gene in orm_genes:
                mean_expr = self.expression_data.loc[gene].mean()
                print(f"  • {gene}: mean = {mean_expr:.4f}")
        
        return orm_genes
    
    def save_processed_data(self, output_dir="results"):
        """Save processed data"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        if self.expression_data is not None:
            output_file = output_dir / "normalized_expression.csv"
            self.expression_data.to_csv(output_file)
            print(f"\nProcessed data saved to: {output_file}")
        
        return output_dir
    
    def get_summary(self):
        """Print data summary"""
        print("\n" + "="*50)
        print("DATA SUMMARY")
        print("="*50)
        
        if self.expression_data is not None:
            print(f"Genes: {self.expression_data.shape[0]}")
            print(f"Samples: {self.expression_data.shape[1]}")
            print(f"Sparsity: {(self.expression_data == 0).sum().sum() / (self.expression_data.shape[0] * self.expression_data.shape[1]) * 100:.2f}%")

def main():
    print("="*60)
    print("Psoriasis Transcriptomics Data Preprocessing")
    print("="*60)
    
    preprocessor = DataPreprocessor()
    
    if not preprocessor.load_data():
        return
    
    preprocessor.quality_control()
    preprocessor.normalize_data(method='log2')
    preprocessor.check_orm_genes()
    preprocessor.get_summary()
    preprocessor.save_processed_data()

if __name__ == "__main__":
    main()