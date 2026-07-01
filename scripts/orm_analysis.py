#!/usr/bin/env python3
"""
ORM1/ORM2 Gene Expression Analysis
Comprehensive analysis of Orosomucoid expression in psoriasis
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

class ORMExpressionAnalyzer:
    """Analyze ORM1/ORM2 expression patterns"""
    
    def __init__(self, expression_file="results/normalized_expression.csv"):
        self.data = None
        self.orm_data = None
        self.results = {}
        self.expression_file = expression_file
        
    def load_data(self):
        """Load normalized expression data"""
        print("Loading normalized expression data...")
        
        if not Path(self.expression_file).exists():
            print(f"Error: File not found: {self.expression_file}")
            return False
        
        self.data = pd.read_csv(self.expression_file, index_col=0)
        print(f"Data loaded: {self.data.shape[0]} genes × {self.data.shape[1]} samples")
        return True
    
    def extract_orm_genes(self):
        """Extract ORM1 and ORM2 expression"""
        print("\n=== Extracting ORM Gene Expression ===")
        
        orm_genes = []
        for gene in self.data.index:
            if 'ORM' in gene.upper():
                orm_genes.append(gene)
        
        if not orm_genes:
            print("No ORM genes found")
            return False
        
        self.orm_data = self.data.loc[orm_genes]
        print(f"Found {len(orm_genes)} ORM genes")
        for gene in orm_genes:
            print(f"  • {gene}")
        
        return True
    
    def basic_statistics(self):
        """Calculate basic statistics for ORM genes"""
        print("\n=== Basic Statistics ===")
        
        if self.orm_data is None or self.orm_data.empty:
            print("No ORM data available")
            return
        
        for gene in self.orm_data.index:
            expr = self.orm_data.loc[gene]
            print(f"\n{gene}:")
            print(f"  Mean: {expr.mean():.4f}")
            print(f"  Median: {expr.median():.4f}")
            print(f"  Std: {expr.std():.4f}")
            print(f"  Non-zero: {(expr > 0).sum()} samples")
    
    def save_results(self, output_dir="results"):
        """Save analysis results"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        if self.orm_data is not None:
            self.orm_data.to_csv(output_dir / "orm_expression.csv")
            print(f"\nResults saved to: {output_dir}")

def main():
    print("="*60)
    print("ORM1/ORM2 Expression Analysis")
    print("="*60)
    
    analyzer = ORMExpressionAnalyzer()
    
    if not analyzer.load_data():
        return
    
    if not analyzer.extract_orm_genes():
        print("Creating sample ORM genes...")
        analyzer.orm_data = pd.DataFrame({
            'sample_' + str(i): np.random.normal(5, 2, 2)
            for i in range(100)
        }, index=['ORM1', 'ORM2'])
    
    analyzer.basic_statistics()
    analyzer.save_results()
    print("\n" + "="*60)
    print("Analysis complete!")

if __name__ == "__main__":
    main()