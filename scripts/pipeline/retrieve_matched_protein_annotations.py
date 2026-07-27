import re
import time
from collections import Counter

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

from crossreactivity.io import DATA_DIR, load_csv, save_csv


INPUT_PATH = DATA_DIR / "proccesed/iedb_match_regions_long.csv"

MISSING_IDS_PATH = (
    DATA_DIR / "intermediate/matched_proteins_missing_annotations.csv"
)

LOOKUP_PATH = (
    DATA_DIR / "intermediate/matched_protein_annotation_lookup.csv"
)

OUTPUT_PATH = (
    DATA_DIR / "proccesed/iedb_match_regions_long_annotated.csv"
)

UNIPARC_URL = "https://rest.uniprot.org/uniparc/{protein_id}/databases"
UNIPROTKB_URL = "https://rest.uniprot.org/uniprotkb/{protein_id}.json"

REQUEST_TIMEOUT = 120
REQUEST_DELAY = 0.2
CHECKPOINT_EVERY = 25


# A UPI can occur in several organisms. These columns therefore define
# the annotation target more safely than Pathogen_Protein_ID alone.
KEY_COLS = [
    "Pathogen_Protein_ID",
    "Proteome_ID",
    "Pathogen_Scientific_name",
    "Pathogen_Organism",
]

RESULT_COLS = [
    "Retrieved_Annotation",
    "Retrieved_Gene_Name",
    "Annotation_Source_Database",
    "Annotation_Source_ID",
    "Annotation_Source_Organism",
    "Annotation_Match_Level",
    "Annotation_Candidate_Count",
    "Annotation_Agreement",
    "Annotation_Retrieval_Status",
    "Annotation_Retrieval_Error",
]


def clean_text(value) -> str:
    """Return a clean string, or an empty string for missing values."""
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.lower() in {"", "na", "nan", "none", "null"}:
        return ""

    return value


def missing_text(series: pd.Series) -> pd.Series:
    """Identify missing or blank text values."""
    cleaned = series.fillna("").astype(str).str.strip()

    return (
        cleaned.eq("")
        | cleaned.str.lower().isin({"na", "nan", "none", "null"})
    )


def normalize_name(value) -> str:
    """Normalize an organism name for comparison."""
    return " ".join(clean_text(value).casefold().split())


def normalize_annotation(value) -> str:
    """
    Normalize a protein name for consensus counting.

    Sorting tokens treats descriptions such as "IS6 family transposase"
    and "Transposase, IS6 family" as equivalent.
    """
    tokens = re.findall(r"[a-z0-9]+", clean_text(value).casefold())
    return " ".join(sorted(tokens))


def row_key(row: pd.Series) -> tuple[str, ...]:
    """Create a stable key that also handles missing context values."""
    return tuple(clean_text(row[col]) for col in KEY_COLS)


def empty_result(
    status: str,
    error: str | None = None,
    match_level: str | None = None,
    candidate_count: int = 0,
    agreement: float | None = None,
) -> dict:
    """Create an unresolved lookup result."""
    return {
        "Retrieved_Annotation": None,
        "Retrieved_Gene_Name": None,
        "Annotation_Source_Database": None,
        "Annotation_Source_ID": None,
        "Annotation_Source_Organism": None,
        "Annotation_Match_Level": match_level,
        "Annotation_Candidate_Count": candidate_count,
        "Annotation_Agreement": agreement,
        "Annotation_Retrieval_Status": status,
        "Annotation_Retrieval_Error": error,
    }


def make_session() -> requests.Session:
    """Create a requests session that retries temporary API failures."""
    session = requests.Session()

    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        respect_retry_after_header=True,
        raise_on_status=False,
    )

    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    session.headers.update({
        "Accept": "application/json",
        "User-Agent": "crossreactivity-annotation-retrieval/1.0",
    })

    return session


def extract_missing_targets(match_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract unique protein/proteome/organism combinations without annotation.
    """
    required = {"Pathogen_Protein_ID", "Pathogen_Annotation"}
    missing_columns = required - set(match_df.columns)

    if missing_columns:
        raise KeyError(
            "Missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    match_df = match_df.copy()

    for col in KEY_COLS:
        if col not in match_df.columns:
            match_df[col] = pd.NA

    targets = (
        match_df.loc[
            missing_text(match_df["Pathogen_Annotation"])
            & ~missing_text(match_df["Pathogen_Protein_ID"]),
            KEY_COLS,
        ]
        .drop_duplicates()
        .sort_values(KEY_COLS, na_position="last")
        .reset_index(drop=True)
    )

    return targets


def fetch_uniparc_cross_references(
    session: requests.Session,
    protein_id: str,
) -> list[dict]:
    """
    Retrieve only UniProtKB cross-references for one UniParc ID.

    Restricting dbTypes avoids downloading thousands of unrelated
    ENA and RefSeq cross-references for common sequences.
    """
    url = UNIPARC_URL.format(protein_id=protein_id)
    params = {
        "dbTypes": "UniProtKB/Swiss-Prot,UniProtKB/TrEMBL",
        "size": 500,
    }

    cross_references = []

    while url:
        response = session.get(
            url,
            params=params,
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()

        cross_references.extend(response.json().get("results", []))

        url = response.links.get("next", {}).get("url")
        params = None

    return cross_references


def candidate_organism(candidate: dict) -> str:
    """Get the scientific name attached to a UniParc cross-reference."""
    organism = candidate.get("organism") or {}
    return clean_text(organism.get("scientificName"))


def candidate_proteomes(candidate: dict) -> set[str]:
    """Get the proteome IDs attached to a UniParc cross-reference."""
    return {
        clean_text(proteome.get("id"))
        for proteome in candidate.get("proteomes", [])
        if clean_text(proteome.get("id"))
    }


def same_species(candidate: dict, species: str) -> bool:
    """
    Match a genus-species name to a possibly strain-specific organism name.
    """
    candidate_name = normalize_name(candidate_organism(candidate))
    species = normalize_name(species)

    if not candidate_name or not species:
        return False

    return (
        candidate_name == species
        or candidate_name.startswith(species + " ")
    )


def choose_consensus(candidates: list[dict], match_level: str) -> dict:
    """
    Select an annotation when one normalized protein name has a unique mode.

    Conflicting names tied for first place are left unresolved instead of
    transferring an arbitrary annotation.
    """
    candidates = [
        candidate
        for candidate in candidates
        if clean_text(candidate.get("proteinName"))
    ]

    if not candidates:
        return empty_result(
            status="annotation_unavailable",
            match_level=match_level,
        )

    counts = Counter(
        normalize_annotation(candidate["proteinName"])
        for candidate in candidates
    )

    highest_count = max(counts.values())
    winners = [
        name
        for name, count in counts.items()
        if count == highest_count
    ]
    agreement = highest_count / len(candidates)

    if len(winners) != 1:
        return empty_result(
            status="ambiguous_annotation",
            match_level=match_level,
            candidate_count=len(candidates),
            agreement=agreement,
        )

    winner = winners[0]

    matching_candidates = [
        candidate
        for candidate in candidates
        if normalize_annotation(candidate["proteinName"]) == winner
    ]

    # After consensus is established, prefer an active record,
    # then Swiss-Prot, then the most recently updated record.
    selected = max(
        matching_candidates,
        key=lambda candidate: (
            bool(candidate.get("active")),
            candidate.get("database") == "UniProtKB/Swiss-Prot",
            clean_text(candidate.get("lastUpdated")),
        ),
    )

    return {
        "Retrieved_Annotation": clean_text(
            selected.get("proteinName")
        ) or None,
        "Retrieved_Gene_Name": clean_text(
            selected.get("geneName")
        ) or None,
        "Annotation_Source_Database": clean_text(
            selected.get("database")
        ) or None,
        "Annotation_Source_ID": clean_text(
            selected.get("id")
        ) or None,
        "Annotation_Source_Organism": (
            candidate_organism(selected) or None
        ),
        "Annotation_Match_Level": match_level,
        "Annotation_Candidate_Count": len(candidates),
        "Annotation_Agreement": agreement,
        "Annotation_Retrieval_Status": "retrieved",
        "Annotation_Retrieval_Error": None,
    }


def resolve_uniparc_annotation(
    target: pd.Series,
    cross_references: list[dict],
) -> dict:
    """
    Resolve a UPI annotation using the most specific available context.

    The order is exact proteome, exact scientific name, and then a
    same-species consensus. Active records are preferred at each level.
    """
    proteome_id = clean_text(target["Proteome_ID"])
    scientific_name = normalize_name(
        target["Pathogen_Scientific_name"]
    )
    organism = clean_text(target["Pathogen_Organism"])

    annotated = [
        candidate
        for candidate in cross_references
        if clean_text(candidate.get("proteinName"))
    ]

    tiers = []

    if proteome_id:
        exact_proteome = [
            candidate
            for candidate in annotated
            if proteome_id in candidate_proteomes(candidate)
        ]

        tiers.extend([
            (
                "exact_proteome_active",
                [
                    candidate
                    for candidate in exact_proteome
                    if candidate.get("active") is True
                ],
            ),
            ("exact_proteome_archived", exact_proteome),
        ])

    if scientific_name:
        exact_name = [
            candidate
            for candidate in annotated
            if normalize_name(candidate_organism(candidate))
            == scientific_name
        ]

        tiers.extend([
            (
                "exact_scientific_name_active",
                [
                    candidate
                    for candidate in exact_name
                    if candidate.get("active") is True
                ],
            ),
            ("exact_scientific_name_archived", exact_name),
        ])

    if organism:
        species_candidates = [
            candidate
            for candidate in annotated
            if same_species(candidate, organism)
        ]

        tiers.extend([
            (
                "same_species_active",
                [
                    candidate
                    for candidate in species_candidates
                    if candidate.get("active") is True
                ],
            ),
            ("same_species_consensus", species_candidates),
        ])

    for match_level, candidates in tiers:
        if candidates:
            return choose_consensus(candidates, match_level)

    return empty_result(
        status="annotation_unavailable",
        match_level="no_matching_cross_reference",
        candidate_count=len(annotated),
    )


def first_nested_value(
    data: dict,
    path: list[str],
) -> str | None:
    """Read a nested UniProtKB JSON value."""
    current = data

    for key in path:
        if not isinstance(current, dict):
            return None

        current = current.get(key)

    value = clean_text(current)
    return value or None


def fetch_uniprotkb_annotation(
    session: requests.Session,
    protein_id: str,
) -> dict:
    """Retrieve an annotation for a normal UniProtKB accession."""
    response = session.get(
        UNIPROTKB_URL.format(protein_id=protein_id),
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()

    payload = response.json()
    description = payload.get("proteinDescription") or {}

    annotation = first_nested_value(
        description,
        ["recommendedName", "fullName", "value"],
    )

    if annotation is None:
        for name_type in ["submissionNames", "alternativeNames"]:
            for name in description.get(name_type) or []:
                annotation = first_nested_value(
                    name,
                    ["fullName", "value"],
                )

                if annotation:
                    break

            if annotation:
                break

    gene_name = None

    for gene in payload.get("genes") or []:
        gene_name = first_nested_value(
            gene,
            ["geneName", "value"],
        )

        if gene_name:
            break

    entry_type = clean_text(payload.get("entryType")).casefold()
    source_database = (
        "UniProtKB/Swiss-Prot"
        if "reviewed" in entry_type and "unreviewed" not in entry_type
        else "UniProtKB/TrEMBL"
    )

    organism = payload.get("organism") or {}

    if annotation is None:
        return empty_result(
            status="annotation_unavailable",
            match_level="direct_uniprotkb_accession",
        )

    return {
        "Retrieved_Annotation": annotation,
        "Retrieved_Gene_Name": gene_name,
        "Annotation_Source_Database": source_database,
        "Annotation_Source_ID": (
            clean_text(payload.get("primaryAccession"))
            or protein_id
        ),
        "Annotation_Source_Organism": (
            clean_text(organism.get("scientificName")) or None
        ),
        "Annotation_Match_Level": "direct_uniprotkb_accession",
        "Annotation_Candidate_Count": 1,
        "Annotation_Agreement": 1.0,
        "Annotation_Retrieval_Status": "retrieved",
        "Annotation_Retrieval_Error": None,
    }


def load_existing_lookup() -> pd.DataFrame:
    """Load an earlier checkpoint if available."""
    if not LOOKUP_PATH.exists():
        return pd.DataFrame(columns=KEY_COLS + RESULT_COLS)

    lookup = load_csv(LOOKUP_PATH)

    for col in KEY_COLS + RESULT_COLS:
        if col not in lookup.columns:
            lookup[col] = pd.NA

    return lookup[KEY_COLS + RESULT_COLS]


def save_checkpoint(
    existing: pd.DataFrame,
    new_rows: list[dict],
) -> pd.DataFrame:
    """Save completed attempts so an interrupted run can resume."""
    frames = [existing]

    if new_rows:
        frames.append(pd.DataFrame(new_rows))

    lookup = (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates(subset=KEY_COLS, keep="last")
        .sort_values(KEY_COLS, na_position="last")
        .reset_index(drop=True)
    )

    save_csv(lookup[KEY_COLS + RESULT_COLS], LOOKUP_PATH)
    return lookup


def retrieve_annotations(targets: pd.DataFrame) -> pd.DataFrame:
    """Retrieve annotations for all unique missing protein contexts."""
    existing = load_existing_lookup()

    # Successful, unavailable, and ambiguous results are reused.
    # Only genuine request failures are attempted again.
    completed = existing[
        existing["Annotation_Retrieval_Status"] != "request_failed"
    ]
    completed_keys = {
        row_key(row)
        for _, row in completed.iterrows()
    }

    pending = targets[
        targets.apply(
            lambda row: row_key(row) not in completed_keys,
            axis=1,
        )
    ].copy()

    print(f"Unique protein-context combinations: {len(targets):,}")
    print(f"Already present in lookup: {len(targets) - len(pending):,}")
    print(
        "Unique protein IDs requiring API requests: "
        f"{pending['Pathogen_Protein_ID'].nunique():,}"
    )

    if pending.empty:
        return existing

    session = make_session()
    new_rows = []

    grouped = pending.groupby(
        "Pathogen_Protein_ID",
        sort=False,
        dropna=False,
    )

    for index, (protein_id, group) in enumerate(
        tqdm(
            grouped,
            total=grouped.ngroups,
            desc="Retrieving annotations",
            unit="protein ID",
        ),
        start=1,
    ):
        protein_id = clean_text(protein_id).upper()

        try:
            if protein_id.startswith("UPI"):
                cross_references = fetch_uniparc_cross_references(
                    session,
                    protein_id,
                )

                results = [
                    resolve_uniparc_annotation(
                        target,
                        cross_references,
                    )
                    for _, target in group.iterrows()
                ]
            else:
                result = fetch_uniprotkb_annotation(
                    session,
                    protein_id,
                )
                results = [result] * len(group)

        except requests.RequestException as error:
            results = [
                empty_result(
                    status="request_failed",
                    error=str(error),
                )
            ] * len(group)

        except (ValueError, TypeError, KeyError) as error:
            results = [
                empty_result(
                    status="response_parse_failed",
                    error=f"Could not parse API response: {error}",
                )
            ] * len(group)

        for (_, target), result in zip(group.iterrows(), results):
            new_rows.append({
                **{col: target[col] for col in KEY_COLS},
                **result,
            })

        if index % CHECKPOINT_EVERY == 0:
            existing = save_checkpoint(existing, new_rows)
            new_rows = []

        time.sleep(REQUEST_DELAY)

    return save_checkpoint(existing, new_rows)


def add_annotations_to_matches(
    match_df: pd.DataFrame,
    lookup: pd.DataFrame,
) -> pd.DataFrame:
    """Fill missing annotations without overwriting existing values."""
    match_df = match_df.copy()
    lookup = lookup.copy()

    for col in KEY_COLS:
        if col not in match_df.columns:
            match_df[col] = pd.NA

    # Temporary string keys allow missing context fields to merge safely
    # without changing the original columns.
    merge_cols = []

    for index, col in enumerate(KEY_COLS):
        merge_col = f"_annotation_key_{index}"
        merge_cols.append(merge_col)
        match_df[merge_col] = match_df[col].map(clean_text)
        lookup[merge_col] = lookup[col].map(clean_text)

    before_rows = len(match_df)

    annotated = match_df.merge(
        lookup[merge_cols + RESULT_COLS].drop_duplicates(merge_cols),
        how="left",
        on=merge_cols,
        validate="many_to_one",
    )

    if len(annotated) != before_rows:
        raise RuntimeError(
            "Annotation merge changed the number of match rows."
        )

    original_missing = missing_text(
        annotated["Pathogen_Annotation"]
    )
    retrieved_available = ~missing_text(
        annotated["Retrieved_Annotation"]
    )
    fill_annotation = original_missing & retrieved_available

    annotated.loc[
        fill_annotation,
        "Pathogen_Annotation",
    ] = annotated.loc[
        fill_annotation,
        "Retrieved_Annotation",
    ]

    if "Pathogen_Gene_Name" not in annotated.columns:
        annotated["Pathogen_Gene_Name"] = pd.NA

    fill_gene = (
        missing_text(annotated["Pathogen_Gene_Name"])
        & ~missing_text(annotated["Retrieved_Gene_Name"])
    )

    annotated.loc[
        fill_gene,
        "Pathogen_Gene_Name",
    ] = annotated.loc[
        fill_gene,
        "Retrieved_Gene_Name",
    ]

    annotated["Pathogen_Annotation_Was_Retrieved"] = fill_annotation

    originally_annotated = ~original_missing

    annotated.loc[
        originally_annotated,
        "Annotation_Source_Database",
    ] = "Original FASTA header"

    annotated.loc[
        originally_annotated,
        "Annotation_Match_Level",
    ] = "original_fasta_header"

    annotated.loc[
        originally_annotated,
        "Annotation_Retrieval_Status",
    ] = "not_needed"

    annotated = annotated.drop(
        columns=(
            merge_cols
            + ["Retrieved_Annotation", "Retrieved_Gene_Name"]
        )
    )

    print(
        "Match rows receiving a retrieved annotation: "
        f"{fill_annotation.sum():,}"
    )
    print(
        "Match rows still lacking annotation: "
        f"{missing_text(annotated['Pathogen_Annotation']).sum():,}"
    )

    return annotated


def main() -> None:
    print(f"Loading match table: {INPUT_PATH}")
    match_df = load_csv(INPUT_PATH)

    targets = extract_missing_targets(match_df)
    save_csv(targets, MISSING_IDS_PATH)

    print(
        "Unique missing protein IDs: "
        f"{targets['Pathogen_Protein_ID'].nunique():,}"
    )
    print(
        "Unique missing protein-context combinations: "
        f"{len(targets):,}"
    )
    print(f"Saved missing-protein list: {MISSING_IDS_PATH}")

    lookup = retrieve_annotations(targets)

    print("\nAnnotation retrieval status:")
    print(
        lookup["Annotation_Retrieval_Status"]
        .value_counts(dropna=False)
        .to_string()
    )
    print(f"Saved annotation lookup: {LOOKUP_PATH}")

    annotated = add_annotations_to_matches(match_df, lookup)
    save_csv(annotated, OUTPUT_PATH)

    print(f"Saved annotated match table: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()