import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer
from collections import Counter

# splits a GO annotation string on semicolons, strips whitespace, and returns list of terms
def parse_go_terms(go_string):
    terms = [term.strip() for term in go_string.split(';') if term.strip()]
    return terms

# retains only the GO terms from terms_list that appear in the top_terms set
def filter_top_50(terms_list, top_terms):
    return [t for t in terms_list if t in top_terms]

def main():
    # read the parsed UniProt .csv file
    df = pd.read_csv("data/dataframe.csv")
    print("Original shape:", df.shape)

    # filter out rows that have no GO (molecular function) annotation
    df_filtered = df[
        df["Gene Ontology (molecular function)"].notnull() &
        (df["Gene Ontology (molecular function)"] != "")
    ].copy()  # this avoids SettingWithCopyWarning

    print("Filtered shape:", df_filtered.shape)
    n_removed = df.shape[0] - df_filtered.shape[0]
    print(f"Removed {n_removed} rows with empty or null GO annotations.")

    # parse multi-label entries into lists of GO terms
    df_filtered["molecular_functions"] = df_filtered["Gene Ontology (molecular function)"].apply(parse_go_terms)

    # determines the top 50 most frequent GO terms overall
    all_terms = []
    for terms_list in df_filtered["molecular_functions"]:
        all_terms.extend(terms_list)

    freq = Counter(all_terms)
    most_common_50 = {term for term, count in freq.most_common(50)}
    print(f"Total unique GO terms before restricting to top 50: {len(freq)}")

    # keep only terms in that top-50 set to reduce dimensionality
    df_filtered["molecular_functions"] = df_filtered["molecular_functions"].apply(
        lambda lst: filter_top_50(lst, most_common_50)
    )

    # drop rows that have no remaining GO terms
    df_filtered = df_filtered[df_filtered["molecular_functions"].map(len) > 0]
    print("After restricting to top 50 GO terms, shape =", df_filtered.shape)

    # Multi-label binarization
    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(df_filtered["molecular_functions"])
    print("Final shape after binarization = (num_proteins, num_terms) =", Y.shape)
    print("GO terms shape = ", Y.shape, ", for ex: Y[0][:5] =", Y[0][:5])
    print("Number of unique terms used =", len(mlb.classes_))

    # save the filtered dataset
    df_filtered.to_csv("data/filtered_parsed_data.csv", index=False)
    print("Saved filtered_parsed_data.csv")

if __name__ == "__main__":
    main()
