import time

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

from crossreactivity.io import DATA_DIR, load_csv, save_csv


INPUT_PATH = (
    DATA_DIR
    / "proccesed/iedb_match_regions_long_metadata.csv"
)

# This occurrence-level lookup deliberately uses a new filename. The previous
# lookup was keyed by protein + species and is not compatible with this schema.
LOOKUP_PATH = (
    DATA_DIR
    / "intermediate/matched_uniparc_to_uniprotkb_occurrence_lookup.csv"
)

OUTPUT_PATH = (
    DATA_DIR
    / "proccesed/iedb_match_regions_long_metadata_annotated.csv"
)

UNIPARC_URL = "https://rest.uniprot.org/uniparc/{protein_id}/databases"

REQUEST_TIMEOUT = 120
REQUEST_DELAY = 0.2
CHECKPOINT_EVERY = 25


KEY_COLS = [
    "Pathogen_Protein_ID",
    "Proteome_ID",
]

TARGET_CONTEXT_COLS = [
    "Pathogen_Scientific_name",
    "Pathogen_Organism",
    "Pathogen_Strain",
]

TARGET_COLS = KEY_COLS + TARGET_CONTEXT_COLS

LOOKUP_COLS = KEY_COLS + [
    "UniProtKB_ID",
    "Retrieved_Annotation",
    "Resolution_Status",
    "Candidate_Count",
    "Resolution_Method",
    "Error_Message",
]

TERMINAL_STATUSES = {
    "resolved_unique",
    "ambiguous_multiple",
    "no_active_mapping",
    "no_occurrence_match",
}

ALL_STATUSES = TERMINAL_STATUSES | {"api_error"}


def clean_text(value) -> str:
    """Return a stripped string, or an empty string for missing values."""
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.casefold() in {"", "na", "nan", "none", "null"}:
        return ""

    return value


def missing_text(series: pd.Series) -> pd.Series:
    """Return True for missing or blank text values."""
    cleaned = series.fillna("").astype(str).str.strip()

    return (
        cleaned.eq("")
        | cleaned.str.casefold().isin(
            {"na", "nan", "none", "null"}
        )
    )


def normalize_identifier(value) -> str:
    """Normalize a protein or proteome identifier for matching."""
    return clean_text(value).upper()


def normalize_name(value) -> str:
    """Normalize an organism or annotation name for comparison."""
    return " ".join(clean_text(value).casefold().split())


def candidate_organism(candidate: dict) -> str:
    """Return the scientific name from a UniParc cross-reference."""
    organism = candidate.get("organism") or {}
    return clean_text(organism.get("scientificName"))


def candidate_proteome_ids(candidate: dict) -> set[str]:
    """Return all proteome IDs exposed by a UniParc cross-reference."""
    values = []

    for key in (
        "proteomeId",
        "proteomeID",
        "proteomeIds",
        "proteomeIDs",
    ):
        value = candidate.get(key)

        if isinstance(value, (list, tuple, set)):
            values.extend(value)
        elif value is not None:
            values.append(value)

    proteome = candidate.get("proteome")

    if isinstance(proteome, dict):
        values.extend([
            proteome.get("id"),
            proteome.get("proteomeId"),
        ])
    elif proteome is not None:
        values.append(proteome)

    return {
        normalize_identifier(value)
        for value in values
        if normalize_identifier(value)
    }


def same_species(candidate: dict, species: str) -> bool:
    """Match a species name to a possibly strain-specific organism name."""
    candidate_name = normalize_name(candidate_organism(candidate))
    species_name = normalize_name(species)

    if not candidate_name or not species_name:
        return False

    return (
        candidate_name == species_name
        or candidate_name.startswith(species_name + " ")
        or candidate_name.startswith(species_name + " (")
    )


def consensus_annotation(candidates: list[dict]) -> str | None:
    """Return a protein name only when all available names agree."""
    annotations = {}

    for candidate in candidates:
        annotation = clean_text(candidate.get("proteinName"))

        if annotation:
            annotations.setdefault(normalize_name(annotation), annotation)

    if len(annotations) == 1:
        return next(iter(annotations.values()))

    return None


def make_session() -> requests.Session:
    """Create a session that retries temporary UniProt API failures."""
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
        "User-Agent": "crossreactivity-uniparc-lookup/2.0",
    })

    return session


def extract_targets(match_df: pd.DataFrame) -> pd.DataFrame:
    """Extract every unique matched UniParc protein occurrence."""
    required_cols = set(TARGET_COLS)
    missing_cols = required_cols - set(match_df.columns)

    if missing_cols:
        raise KeyError(
            "Missing required columns: "
            + ", ".join(sorted(missing_cols))
        )

    working = match_df[TARGET_COLS].copy()
    working["Pathogen_Protein_ID"] = (
        working["Pathogen_Protein_ID"].map(normalize_identifier)
    )
    working["Proteome_ID"] = working["Proteome_ID"].map(
        normalize_identifier
    )

    is_uniparc = working["Pathogen_Protein_ID"].str.startswith("UPI")

    missing_proteome = is_uniparc & missing_text(working["Proteome_ID"])

    if missing_proteome.any():
        examples = (
            working.loc[missing_proteome, "Pathogen_Protein_ID"]
            .drop_duplicates()
            .head(20)
            .tolist()
        )
        raise ValueError(
            "Matched UniParc rows are missing Proteome_ID. Examples: "
            + ", ".join(examples)
        )

    missing_names = (
        is_uniparc
        & missing_text(working["Pathogen_Scientific_name"])
        & missing_text(working["Pathogen_Organism"])
    )

    if missing_names.any():
        examples = (
            working.loc[missing_names, KEY_COLS]
            .drop_duplicates()
            .head(20)
            .to_dict("records")
        )
        raise ValueError(
            "Matched UniParc occurrences lack both scientific and species "
            f"names. Examples: {examples}"
        )

    targets = working.loc[is_uniparc].drop_duplicates().copy()

    conflict_counts = (
        targets.groupby(KEY_COLS, dropna=False)[TARGET_CONTEXT_COLS]
        .nunique(dropna=True)
    )
    conflicting = conflict_counts.gt(1).any(axis=1)

    if conflicting.any():
        examples = conflicting[conflicting].index.tolist()[:20]
        raise ValueError(
            "Conflicting occurrence metadata were found for the same "
            "Pathogen_Protein_ID + Proteome_ID key. Examples: "
            f"{examples}"
        )

    targets = (
        targets
        .sort_values(KEY_COLS)
        .drop_duplicates(subset=KEY_COLS, keep="first")
        .reset_index(drop=True)
    )

    return targets[TARGET_COLS]


def fetch_uniprotkb_cross_references(
    session: requests.Session,
    protein_id: str,
) -> list[dict]:
    """Retrieve all UniProtKB cross-references for one UniParc ID."""
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

        cross_references.extend(
            response.json().get("results", [])
        )

        url = response.links.get("next", {}).get("url")
        params = None

    return cross_references


def active_accession_candidates(
    cross_references: list[dict],
) -> list[dict]:
    """Return active cross-references with a usable UniProtKB accession."""
    return [
        candidate
        for candidate in cross_references
        if (
            candidate.get("active") is True
            and clean_text(candidate.get("id"))
        )
    ]


def resolve_for_occurrence(
    target: pd.Series | dict,
    cross_references: list[dict],
) -> dict:
    """Resolve one UPI + proteome occurrence without arbitrary selection."""
    active_candidates = active_accession_candidates(cross_references)

    if not active_candidates:
        return {
            "UniProtKB_ID": None,
            "Retrieved_Annotation": None,
            "Resolution_Status": "no_active_mapping",
            "Candidate_Count": 0,
            "Resolution_Method": None,
            "Error_Message": None,
        }

    target_proteome = normalize_identifier(target["Proteome_ID"])
    target_scientific_name = normalize_name(
        target["Pathogen_Scientific_name"]
    )
    target_species = clean_text(target["Pathogen_Organism"])

    proteome_matches = [
        candidate
        for candidate in active_candidates
        if target_proteome in candidate_proteome_ids(candidate)
    ]

    if proteome_matches:
        candidates = proteome_matches
        method = "exact_proteome"
    else:
        scientific_name_matches = [
            candidate
            for candidate in active_candidates
            if (
                target_scientific_name
                and normalize_name(candidate_organism(candidate))
                == target_scientific_name
            )
        ]

        if scientific_name_matches:
            candidates = scientific_name_matches
            method = "exact_scientific_name"
        else:
            species_matches = [
                candidate
                for candidate in active_candidates
                if same_species(candidate, target_species)
            ]

            if species_matches:
                candidates = species_matches
                method = "unique_species"
            else:
                return {
                    "UniProtKB_ID": None,
                    "Retrieved_Annotation": None,
                    "Resolution_Status": "no_occurrence_match",
                    "Candidate_Count": 0,
                    "Resolution_Method": None,
                    "Error_Message": None,
                }

    candidates_by_accession = {}

    for candidate in candidates:
        accession = clean_text(candidate.get("id")).upper()
        candidates_by_accession.setdefault(accession, []).append(candidate)

    accession_count = len(candidates_by_accession)
    annotation = consensus_annotation(candidates)

    if accession_count == 1:
        return {
            "UniProtKB_ID": next(iter(candidates_by_accession)),
            "Retrieved_Annotation": annotation,
            "Resolution_Status": "resolved_unique",
            "Candidate_Count": 1,
            "Resolution_Method": method,
            "Error_Message": None,
        }

    return {
        "UniProtKB_ID": None,
        "Retrieved_Annotation": annotation,
        "Resolution_Status": "ambiguous_multiple",
        "Candidate_Count": accession_count,
        "Resolution_Method": method,
        "Error_Message": None,
    }


def load_existing_lookup() -> pd.DataFrame:
    """Load and validate an occurrence-level checkpoint if it exists."""
    if not LOOKUP_PATH.exists():
        return pd.DataFrame(columns=LOOKUP_COLS)

    lookup = load_csv(LOOKUP_PATH)
    missing_cols = set(LOOKUP_COLS) - set(lookup.columns)

    if missing_cols:
        raise ValueError(
            f"Existing lookup has an incompatible schema: {LOOKUP_PATH}. "
            "Move or remove it before rerunning. Missing columns: "
            + ", ".join(sorted(missing_cols))
        )

    lookup = lookup[LOOKUP_COLS].copy()
    lookup["Pathogen_Protein_ID"] = (
        lookup["Pathogen_Protein_ID"].map(normalize_identifier)
    )
    lookup["Proteome_ID"] = lookup["Proteome_ID"].map(
        normalize_identifier
    )

    unknown_statuses = (
        set(lookup["Resolution_Status"].dropna().astype(str))
        - ALL_STATUSES
    )

    if unknown_statuses:
        raise ValueError(
            "Existing lookup contains unknown Resolution_Status values: "
            + ", ".join(sorted(unknown_statuses))
        )

    return (
        lookup
        .drop_duplicates(subset=KEY_COLS, keep="last")
        .sort_values(KEY_COLS)
        .reset_index(drop=True)
    )


def save_checkpoint(
    existing: pd.DataFrame,
    new_rows: list[dict],
) -> pd.DataFrame:
    """Save all outcomes so only API errors are retried later."""
    frames = [existing]

    if new_rows:
        frames.append(pd.DataFrame(new_rows))

    lookup = pd.concat(frames, ignore_index=True)

    for col in LOOKUP_COLS:
        if col not in lookup.columns:
            lookup[col] = pd.NA

    lookup["Pathogen_Protein_ID"] = (
        lookup["Pathogen_Protein_ID"].map(normalize_identifier)
    )
    lookup["Proteome_ID"] = lookup["Proteome_ID"].map(
        normalize_identifier
    )

    lookup = (
        lookup[LOOKUP_COLS]
        .drop_duplicates(subset=KEY_COLS, keep="last")
        .sort_values(KEY_COLS)
        .reset_index(drop=True)
    )

    save_csv(lookup, LOOKUP_PATH)
    return lookup


def retrieve_mappings(targets: pd.DataFrame) -> pd.DataFrame:
    """Retrieve and resolve mappings for all matched UPI occurrences."""
    existing = load_existing_lookup()

    completed_keys = set(
        existing.loc[
            existing["Resolution_Status"].isin(TERMINAL_STATUSES),
            KEY_COLS,
        ].itertuples(index=False, name=None)
    )

    pending_mask = ~targets.apply(
        lambda row: tuple(row[col] for col in KEY_COLS)
        in completed_keys,
        axis=1,
    )
    pending = targets.loc[pending_mask].copy()

    print(f"Unique UniParc occurrence targets: {len(targets):,}")
    print(f"Already completed: {len(targets) - len(pending):,}")
    print(f"Pending or retryable: {len(pending):,}")

    if pending.empty:
        save_csv(existing[LOOKUP_COLS], LOOKUP_PATH)
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
            desc="Looking up UniParc IDs",
            unit="UniParc ID",
        ),
        start=1,
    ):
        protein_id = normalize_identifier(protein_id)

        try:
            cross_references = fetch_uniprotkb_cross_references(
                session,
                protein_id,
            )
            api_error = None
        except requests.RequestException as exc:
            cross_references = []
            api_error = clean_text(exc)[:500] or exc.__class__.__name__

        for _, target in group.iterrows():
            if api_error is None:
                result = resolve_for_occurrence(
                    target,
                    cross_references,
                )
            else:
                result = {
                    "UniProtKB_ID": None,
                    "Retrieved_Annotation": None,
                    "Resolution_Status": "api_error",
                    "Candidate_Count": pd.NA,
                    "Resolution_Method": None,
                    "Error_Message": api_error,
                }

            new_rows.append({
                "Pathogen_Protein_ID": protein_id,
                "Proteome_ID": normalize_identifier(
                    target["Proteome_ID"]
                ),
                **result,
            })

        if index % CHECKPOINT_EVERY == 0:
            existing = save_checkpoint(existing, new_rows)
            new_rows = []

        time.sleep(REQUEST_DELAY)

    return save_checkpoint(existing, new_rows)


def add_mappings_to_matches(
    match_df: pd.DataFrame,
    lookup: pd.DataFrame,
) -> pd.DataFrame:
    """Add occurrence-specific IDs while preserving canonical UniParc IDs."""
    required_cols = {
        "Pathogen_Protein_ID",
        "Proteome_ID",
        "Pathogen_Annotation",
    }
    missing_cols = required_cols - set(match_df.columns)

    if missing_cols:
        raise KeyError(
            "Missing required match-table columns: "
            + ", ".join(sorted(missing_cols))
        )

    annotated = match_df.copy()
    before_rows = len(annotated)

    annotated["_Pathogen_Protein_ID_key"] = (
        annotated["Pathogen_Protein_ID"].map(normalize_identifier)
    )
    annotated["_Proteome_ID_key"] = (
        annotated["Proteome_ID"].map(normalize_identifier)
    )

    lookup_for_merge = lookup[
        KEY_COLS + ["UniProtKB_ID", "Retrieved_Annotation"]
    ].copy()
    lookup_for_merge["Pathogen_Protein_ID"] = (
        lookup_for_merge["Pathogen_Protein_ID"].map(normalize_identifier)
    )
    lookup_for_merge["Proteome_ID"] = (
        lookup_for_merge["Proteome_ID"].map(normalize_identifier)
    )
    lookup_for_merge = lookup_for_merge.rename(columns={
        "Pathogen_Protein_ID": "_Pathogen_Protein_ID_key",
        "Proteome_ID": "_Proteome_ID_key",
        "UniProtKB_ID": "_Retrieved_UniProtKB_ID",
    })

    annotated = annotated.merge(
        lookup_for_merge,
        how="left",
        on=["_Pathogen_Protein_ID_key", "_Proteome_ID_key"],
        validate="many_to_one",
    )

    if len(annotated) != before_rows:
        raise RuntimeError(
            "The occurrence lookup merge changed the number of match rows."
        )

    if "UniProtKB_ID" not in annotated.columns:
        insert_at = annotated.columns.get_loc("Pathogen_Protein_ID") + 1
        annotated.insert(insert_at, "UniProtKB_ID", pd.NA)

    conflicting_ids = (
        ~missing_text(annotated["UniProtKB_ID"])
        & ~missing_text(annotated["_Retrieved_UniProtKB_ID"])
        & (
            annotated["UniProtKB_ID"].map(normalize_identifier)
            != annotated["_Retrieved_UniProtKB_ID"].map(
                normalize_identifier
            )
        )
    )

    if conflicting_ids.any():
        examples = annotated.loc[
            conflicting_ids,
            [
                "Pathogen_Protein_ID",
                "Proteome_ID",
                "UniProtKB_ID",
                "_Retrieved_UniProtKB_ID",
            ],
        ].head(20)
        raise ValueError(
            "Existing and retrieved UniProtKB IDs disagree. Examples:\n"
            + examples.to_string(index=False)
        )

    fill_id_mask = (
        missing_text(annotated["UniProtKB_ID"])
        & ~missing_text(annotated["_Retrieved_UniProtKB_ID"])
    )
    annotated.loc[fill_id_mask, "UniProtKB_ID"] = annotated.loc[
        fill_id_mask,
        "_Retrieved_UniProtKB_ID",
    ]

    fill_annotation_mask = (
        missing_text(annotated["Pathogen_Annotation"])
        & ~missing_text(annotated["Retrieved_Annotation"])
    )
    annotated.loc[
        fill_annotation_mask,
        "Pathogen_Annotation",
    ] = annotated.loc[
        fill_annotation_mask,
        "Retrieved_Annotation",
    ]

    # The canonical Pathogen_Protein_ID remains the original UPI. A separate
    # UniParc_ID copy is therefore redundant in the normalized workflow.
    columns_to_drop = [
        "_Pathogen_Protein_ID_key",
        "_Proteome_ID_key",
        "_Retrieved_UniProtKB_ID",
        "Retrieved_Annotation",
    ]

    if "UniParc_ID" in annotated.columns:
        columns_to_drop.append("UniParc_ID")

    annotated = annotated.drop(columns=columns_to_drop)

    print(
        "Match rows receiving a unique UniProtKB ID: "
        f"{fill_id_mask.sum():,}"
    )
    print(
        "Match rows receiving a consensus annotation: "
        f"{fill_annotation_mask.sum():,}"
    )
    print(
        "Match rows still lacking annotation: "
        f"{missing_text(annotated['Pathogen_Annotation']).sum():,}"
    )

    return annotated


def print_resolution_summary(lookup: pd.DataFrame) -> None:
    """Print occurrence-level resolution counts."""
    print("\nResolution status counts:")
    print(
        lookup["Resolution_Status"]
        .value_counts(dropna=False)
        .to_string()
    )


def main() -> None:
    print(f"Loading match table: {INPUT_PATH}")
    match_df = load_csv(INPUT_PATH)

    targets = extract_targets(match_df)
    lookup = retrieve_mappings(targets)

    print(f"Saved occurrence-level lookup table: {LOOKUP_PATH}")
    print_resolution_summary(lookup)

    annotated = add_mappings_to_matches(match_df, lookup)
    save_csv(annotated, OUTPUT_PATH)

    print(f"Saved annotated match table: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()