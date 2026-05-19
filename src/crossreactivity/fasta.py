

def parse_fasta_text(fasta_text: str) -> tuple[str, str] | None:
    if not fasta_text:
        return None

    lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        return None

    header = lines[0][1:]
    sequence = "".join(lines[1:])

    if not sequence:
        return None

    return header, sequence


def fasta_wrap(seq: str, width: int = 60) -> str:
    seq = str(seq).strip()
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    handle.write(f"{fasta_wrap(sequence)}\n")