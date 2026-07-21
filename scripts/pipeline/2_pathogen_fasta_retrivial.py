import time
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv


UNIPROTKB_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPARC_SEARCH_URL = "https://rest.uniprot.org/uniparc/search"


def load_proteome_ids() -> list[str]:
    proteome_ids_path = DATA_DIR / "intermediate/proteome_ids.csv"

    proteome_ids = (
        load_csv(proteome_ids_path)
        .iloc[:, 0]
        .dropna()
        .astype(str)
        .str.strip()
    )

    proteome_ids = proteome_ids[proteome_ids != ""]

    return list(dict.fromkeys(proteome_ids.tolist()))


def make_session() -> requests.Session:
    session = requests.Session()

    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )

    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    session.headers.update({
        "User-Agent": "crossreactivity-proteome-download/1.0"
    })

    return session


def get_next_link(headers: dict) -> str | None:
    link_header = headers.get("Link")

    if not link_header:
        return None

    for part in link_header.split(","):
        if 'rel="next"' in part:
            return part[part.find("<") + 1:part.find(">")]

    return None


def fetch_fasta_paginated(
    session: requests.Session,
    url: str,
    proteome_id: str,
    size: int = 500,
    timeout: int = 300,
) -> tuple[list[str], dict]:
    params = {
        "query": f"proteome:{proteome_id}",
        "format": "fasta",
        "size": size,
    }

    metadata = {
        "proteome_id": proteome_id,
        "status_code": None,
        "n_pages": 0,
        "n_records": 0,
        "error": "",
        "response_preview": "",
    }

    fasta_pages = []

    try:
        response = session.get(url, params=params, timeout=timeout)

        while True:
            metadata["status_code"] = response.status_code
            metadata["response_preview"] = response.text[:300].replace("\n", " ")

            if response.status_code != 200:
                metadata["error"] = f"HTTP {response.status_code}"
                return fasta_pages, metadata

            text = response.text.strip()

            if text:
                fasta_pages.append(text)
                metadata["n_pages"] += 1
                metadata["n_records"] += text.count(">")

            next_link = get_next_link(response.headers)

            if next_link is None:
                break

            response = session.get(next_link, timeout=timeout)

        if metadata["n_records"] == 0:
            metadata["error"] = "empty response"

        return fasta_pages, metadata

    except requests.exceptions.Timeout:
        metadata["error"] = "timeout"
        return fasta_pages, metadata

    except requests.exceptions.RequestException as e:
        metadata["error"] = f"request error: {e}"
        return fasta_pages, metadata


def fetch_proteome_fasta(
    session: requests.Session,
    proteome_id: str,
) -> tuple[list[str], dict]:
    """
    Try UniProtKB first. If no records are found, fall back to UniParc.
    """

    uniprotkb_pages, uniprotkb_meta = fetch_fasta_paginated(
        session=session,
        url=UNIPROTKB_SEARCH_URL,
        proteome_id=proteome_id,
    )

    if uniprotkb_meta["n_records"] > 0:
        return uniprotkb_pages, {
            **uniprotkb_meta,
            "source_database": "UniProtKB",
            "final_status": "success",
        }

    uniparc_pages, uniparc_meta = fetch_fasta_paginated(
        session=session,
        url=UNIPARC_SEARCH_URL,
        proteome_id=proteome_id,
    )

    if uniparc_meta["n_records"] > 0:
        return uniparc_pages, {
            **uniparc_meta,
            "source_database": "UniParc",
            "final_status": "success",
            "uniprotkb_error": uniprotkb_meta["error"],
            "uniprotkb_records": uniprotkb_meta["n_records"],
        }

    return [], {
        "proteome_id": proteome_id,
        "source_database": "none",
        "final_status": "failed",
        "n_pages": 0,
        "n_records": 0,
        "status_code": uniparc_meta["status_code"],
        "error": f"UniProtKB: {uniprotkb_meta['error']} | UniParc: {uniparc_meta['error']}",
        "response_preview": uniparc_meta["response_preview"],
    }


def write_fasta_pages(
    fasta_pages: list[str],
    proteome_id: str,
    source_database: str,
    output_handle,
) -> None:
    for page in fasta_pages:
        for line in page.splitlines():
            if line.startswith(">"):
                output_handle.write(
                    f"{line} PROTEOME={proteome_id} SOURCE_DB={source_database}\n"
                )
            else:
                output_handle.write(f"{line}\n")


def main() -> None:
    proteome_ids = load_proteome_ids()

    output_file = DATA_DIR / "raw/all_proteomes.fasta"
    report_file = DATA_DIR / "raw/proteome_download_report.csv"

    output_file.parent.mkdir(parents=True, exist_ok=True)

    session = make_session()

    report_rows = []
    total_records = 0

    print(f"Downloading {len(proteome_ids)} proteomes...")

    with open(output_file, "w", encoding="utf-8") as out:
        for proteome_id in tqdm(
            proteome_ids,
            desc="Downloading proteomes",
            unit="proteome",
        ):
            fasta_pages, metadata = fetch_proteome_fasta(
                session=session,
                proteome_id=proteome_id,
            )

            if metadata["final_status"] == "success":
                write_fasta_pages(
                    fasta_pages=fasta_pages,
                    proteome_id=proteome_id,
                    source_database=metadata["source_database"],
                    output_handle=out,
                )

                total_records += metadata["n_records"]

            report_rows.append(metadata)

            time.sleep(0.2)

    pd.DataFrame(report_rows).to_csv(report_file, index=False)

    report = pd.DataFrame(report_rows)

    print(f"Saved FASTA: {output_file}")
    print(f"Saved report: {report_file}")
    print(f"Total protein records: {total_records}")
    print()
    print(report["final_status"].value_counts(dropna=False))
    print()
    print(report["source_database"].value_counts(dropna=False))


if __name__ == "__main__":
    main()