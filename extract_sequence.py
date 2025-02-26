from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd

# load the CSV file
df = pd.read_csv("data/filtered_parsed_data.csv")

records = []
for idx, row in df.iterrows():
    # use the "Entry" column as the identifier if available, otherwise use idx
    record_id = str(row["Entry"]) if "Entry" in row else f"seq{idx}"
    record = SeqRecord(Seq(row["Sequence"]), id=record_id, description="")
    records.append(record)

# write the sequences to a FASTA file
SeqIO.write(records, "data/filtered_parsed_sequences.fasta", "fasta")
print("Filtered sequences saved to data/filtered_parsed_sequences.fasta")
