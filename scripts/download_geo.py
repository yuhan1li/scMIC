"""Download GEO supplementary files for MOT-MIC analyses.

This script intentionally downloads only supplementary files and metadata pages.
Large SRA FASTQ/BAM downloads should be managed separately with fasterq-dump or
institutional mirrors.
"""

from __future__ import annotations

import argparse
import re
import urllib.request
from pathlib import Path


DEFAULT_GSE = ["GSE173958", "GSE249057", "GSE178318", "GSE277783"]


def geo_supplement_url(gse: str) -> str:
    prefix = re.sub(r"\d{3}$", "nnn", gse)
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"


def list_supplementary_files(gse: str) -> list[str]:
    url = geo_supplement_url(gse)
    with urllib.request.urlopen(url, timeout=60) as handle:
        html = handle.read().decode("utf-8", errors="ignore")
    return sorted(set(re.findall(r'href="([^"]+\.(?:gz|txt|csv|tsv|h5|h5ad|rds|RDS|mtx))"', html)))


def parse_geo_filelist_text(text: str) -> list[str]:
    """Parse GEO filelist.txt into downloadable file names."""

    names = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] in {"Archive", "File"}:
            names.append(parts[1])
    return names


def parse_geo_filelist(path_or_url: Path | str) -> list[str]:
    """Parse GEO filelist.txt from a local path or URL."""

    if isinstance(path_or_url, Path):
        text = path_or_url.read_text(errors="ignore")
    else:
        with urllib.request.urlopen(path_or_url, timeout=60) as handle:
            text = handle.read().decode("utf-8", errors="ignore")
    return parse_geo_filelist_text(text)


def download_file(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.stat().st_size > 0:
        return
    urllib.request.urlretrieve(url, output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", nargs="*", default=DEFAULT_GSE)
    parser.add_argument("--outdir", default="data/raw")
    parser.add_argument("--from-filelist", action="store_true", help="Download files listed inside GEO filelist.txt.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    for gse in args.gse:
        base = geo_supplement_url(gse)
        files = list_supplementary_files(gse)
        print(f"{gse}: {len(files)} supplementary files")
        for filename in files:
            print(" ", base + filename)
            if not args.dry_run:
                download_file(base + filename, outdir / gse / filename)
        filelist_path = outdir / gse / "filelist.txt"
        if args.from_filelist:
            listed = parse_geo_filelist(filelist_path if filelist_path.exists() else base + "filelist.txt")
            print(f"{gse}: {len(listed)} files listed in filelist.txt")
            for filename in listed:
                print(" ", base + filename)
                if not args.dry_run:
                    download_file(base + filename, outdir / gse / filename)


if __name__ == "__main__":
    main()
