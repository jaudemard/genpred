# Gene prediction

Performs genes prediction.  
This work is part of an academic course gave by the Dr. Ghozlane of the Pasteur Institute.

## Basic usage

```
usage:
    python gpred.py -i data/listeria.fasta -p predicted_gene_positions.csv -o  predicted_genes.fasta 

option:
  -h, --help            show this help message and exit
  -i GENOME_FILE        Complete genome file in fasta format
  -g MIN_GENE_LEN       Minimum gene length to consider (default 50).
  -s MAX_SHINE_DALGARNO_DISTANCE
                        Maximum distance from start codon where to look for a Shine-Dalgarno motif (default 16)
  -d MIN_GAP            Minimum gap between two genes (shine box not included, default 40).
  -p PREDICTED_GENES_FILE
                        Tabular file giving position of predicted genes
  -o FASTA_FILE         Fasta file giving sequence of predicted genes
```
