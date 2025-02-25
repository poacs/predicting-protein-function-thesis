import pandas as pd
from Bio import SeqIO

# reads a UniProt .tsv and returns a DataFrame
def load_uniprot_tsv(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    return df

# parses UniProt .fasta and returns a dictionary: key (string) = UniProt accession, val (string) = full amino acid sequence
def load_uniprot_fasta(fasta_path):
    seq_dict = {}
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # record.id has a form like "sp|P99999|PROT_HUMAN" so below splits on '|' to get the actual accession
            parts = record.id.split('|')
            if len(parts) == 3:
                accession = parts[1]
            else:
                accession = record.id
            seq_dict[accession] = str(record.seq)
    return seq_dict

# this function merges a UniProt .tsv DataFrame and a dict of .fasta seqs on the accession and 
# returns a new DataFrame with a 'sequence' column as well
def merge_tsv_and_fasta(tsv_df, seq_dict):
    accession_col = 'Entry'
    tsv_df['Sequence'] = tsv_df[accession_col].map(seq_dict)
    # drop proteins that don't match, they have NaN in the Sequence column.
    tsv_df.dropna(subset=['Sequence'], inplace=True)
    return tsv_df

def main():
    df = load_uniprot_tsv("Data/uniprotkb_data.tsv")
    seq_dict = load_uniprot_fasta("Data/uniprotkb_data.fasta")
    merged = merge_tsv_and_fasta(df, seq_dict)
    print(f"Merged: {merged.shape[0]} rows")

    # save final DataFrame
    df.to_csv("Data/dataframe.csv", index=False)
    print("Wrote Data/dataframe.csv")

if __name__ == "__main__":
    main()
