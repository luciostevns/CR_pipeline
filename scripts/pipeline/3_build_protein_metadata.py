import re

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv, save_csv


def header_value(header: str, field: str) -> str | None:
    match = re.search(rf"(?:^| ){field}=([^ ]+)", header)
    return match.group(1) if match else None


def parse_annotation(header: str) -> str | None:
    parts = header.split(maxsplit=1)
    if len(parts) == 1:
        return None
    annotation = re.split(
        r"(?:^| )(?:OS|OX|GN|PE|SV|PROTEOME|SOURCE_DB|status)=",
        parts[1],
        maxsplit=1,
    )[0].strip()
    return annotation or None


def parse_fasta() -> pd.DataFrame:
    rows = []
    records = SeqIO.parse(DATA_DIR / "raw/all_proteomes.fasta", "fasta")
    for record in tqdm(records, desc="Parsing proteins", unit="seq"):
        header = record.description
        match = re.match(r"(?:sp|tr)\|([^|]+)\|", header)
        rows.append(
            (
                match.group(1) if match else record.id,
                header_value(header, "PROTEOME"),
                parse_annotation(header),
                header_value(header, "GN"),
                str(record.seq),
                header_value(header, "SOURCE_DB"),
            )
        )
    return pd.DataFrame(
        rows,
        columns=[
            "Protein_ID",
            "Proteome_ID",
            "Annotation",
            "Gene_name",
            "Sequence",
            "Source_database",
        ],
    )


def print_audit(proteins: pd.DataFrame, occurrences: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("RAW PATHOGEN PROTEIN DATASET AUDIT")
    print("=" * 70)
    print(f"Total raw FASTA rows: {len(proteins):,}")

    for source in ["UniProtKB", "UniParc"]:
        rows = proteins[proteins["Source_database"].eq(source)]
        print(
            f"{source} rows: {len(rows):,} "
            f"(unique Protein_IDs: {rows['Protein_ID'].nunique():,})"
        )

    known_strain = occurrences[
        ["Protein_ID", "Genus_species", "Strain"]
    ].dropna()
    known_strain = known_strain[
        known_strain["Genus_species"].astype(str).str.strip().ne("")
        & known_strain["Strain"].astype(str).str.strip().ne("")
    ]
    multi_strain = known_strain.groupby(
        ["Protein_ID", "Genus_species"]
    )["Strain"].transform("nunique").gt(1)
    print(
        "Same Protein_ID and Genus_species, different strains: "
        f"{multi_strain.sum():,} rows "
        f"({known_strain.loc[multi_strain, 'Protein_ID'].nunique():,} unique Protein_IDs)"
    )

    known_species = occurrences[["Protein_ID", "Genus_species"]].dropna()
    known_species = known_species[
        known_species["Genus_species"].astype(str).str.strip().ne("")
    ]
    cross_species = known_species.groupby("Protein_ID")[
        "Genus_species"
    ].transform("nunique").gt(1)
    print(
        "Same Protein_ID, different genus and/or species: "
        f"{cross_species.sum():,} rows "
        f"({known_species.loc[cross_species, 'Protein_ID'].nunique():,} unique Protein_IDs)"
    )
    print("=" * 70 + "\n")


def main() -> None:
    proteins = parse_fasta()
    proteomes = load_csv(DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv")
    organism_cols = ["Proteome_ID", "Scientific_name", "Genus_species", "Strain"]

    raw_occurrences = proteins[["Protein_ID", "Proteome_ID"]].merge(
        proteomes[organism_cols], on="Proteome_ID", how="left", validate="many_to_one"
    )
    print_audit(proteins, raw_occurrences)
    occurrences = raw_occurrences.drop_duplicates()

    sequences = (
        proteins.groupby(["Protein_ID", "Sequence"], as_index=False, sort=False)[
            ["Annotation", "Gene_name"]
        ]
        .first()
    )
    conflicts = sequences["Protein_ID"].duplicated(keep=False)
    if conflicts.any():
        examples = ", ".join(sequences.loc[conflicts, "Protein_ID"].unique()[:10])
        raise ValueError(f"Protein IDs with multiple sequences: {examples}")

    sequences = sequences[["Protein_ID", "Sequence", "Annotation", "Gene_name"]]
    save_csv(sequences, DATA_DIR / "intermediate/protein_sequences.csv")
    save_csv(occurrences, DATA_DIR / "intermediate/protein_occurrences.csv")

    print(f"Unique protein sequences: {len(sequences):,}")
    print(f"Protein occurrences: {len(occurrences):,}")


if __name__ == "__main__":
    main()