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
    # load lengths 201-400
    df_201_400 = load_uniprot_tsv("data/uniprot_201-400.tsv")
    seq_dict_201_400 = load_uniprot_fasta("data/uniprot_201-400.fasta")
    merged_201_400 = merge_tsv_and_fasta(df_201_400, seq_dict_201_400)
    print(f"Merged 201-400: {merged_201_400.shape[0]} rows")

    # load lengths 401-600
    df_401_600 = load_uniprot_tsv("data/uniprot_401-600.tsv")
    seq_dict_401_600 = load_uniprot_fasta("data/uniprot_401-600.fasta")
    merged_401_600 = merge_tsv_and_fasta(df_401_600, seq_dict_401_600)
    print(f"Merged 401-600: {merged_401_600.shape[0]} rows")

    # concatenate both
    combined_df = pd.concat([merged_201_400, merged_401_600], ignore_index=True)
    print("Total combined rows:", combined_df.shape[0])

    # save final combined DataFrame
    combined_df.to_csv("data/combined_uniprot_data.csv", index=False)
    print("Wrote data/combined_uniprot_data.csv")

if __name__ == "__main__":
    main()
