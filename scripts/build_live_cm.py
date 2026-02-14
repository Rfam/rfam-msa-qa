#!/usr/bin/env python3
"""
Build a live Rfam CM database from the Rfam SVN repository.

Downloads the latest CM files for all Rfam families from the SVN repo
and concatenates them into a single pressed CM database.

Usage:
    python3 scripts/build_live_cm.py
    python3 scripts/build_live_cm.py --output /path/to/Rfam_live.cm
    python3 scripts/build_live_cm.py --delay 0.2
"""

import argparse
import re
import shutil
import subprocess
import sys
import time


SVN_BASE_URL = "https://svn.rfam.org/svn/data_repos/trunk/Families"


def check_cmpress():
    """Check that cmpress is available in PATH."""
    if shutil.which("cmpress") is None:
        print("Error: cmpress not found in PATH. Install Infernal first.", file=sys.stderr)
        print("  See: http://eddylab.org/infernal/", file=sys.stderr)
        sys.exit(1)


def fetch_family_list():
    """Fetch the list of RF##### family IDs from the SVN directory listing."""
    result = subprocess.run(
        ["wget", "-q", "-O", "-", SVN_BASE_URL + "/"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("Error: Failed to fetch family listing from SVN.", file=sys.stderr)
        print(f"  URL: {SVN_BASE_URL}/", file=sys.stderr)
        if result.stderr:
            print(f"  wget: {result.stderr.strip()}", file=sys.stderr)
        sys.exit(1)

    families = sorted(set(re.findall(r'(RF\d{5})/', result.stdout)))
    if not families:
        print("Error: No RF##### families found in SVN listing.", file=sys.stderr)
        sys.exit(1)

    return families


def download_cm(family):
    """Download the CM file for a single family. Returns CM text or None."""
    url = f"{SVN_BASE_URL}/{family}/CM"
    result = subprocess.run(
        ["wget", "-q", "-O", "-", url],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None
    return result.stdout


def rename_cm(cm_text, family):
    """Rename the NAME field in a CM file to the family ID.

    SVN CM files are built with `cmbuild -F CM SEED`, so they all have
    NAME=SEED. We replace the NAME with the family ID (e.g., RF00001)
    to ensure unique primary keys when pressed into a single database.
    """
    return re.sub(
        r'^(NAME\s+)\S+',
        rf'\g<1>{family}',
        cm_text,
        count=1,
        flags=re.MULTILINE,
    )


def build_live_cm(output_path, delay=0.1):
    """Build a live Rfam CM database from SVN."""
    check_cmpress()

    print("Fetching family list from Rfam SVN...")
    families = fetch_family_list()
    total = len(families)
    print(f"Found {total} families.")

    downloaded = 0
    skipped = 0

    with open(output_path, "w") as out:
        for i, family in enumerate(families, 1):
            cm_text = download_cm(family)
            if cm_text:
                cm_text = rename_cm(cm_text, family)
                out.write(cm_text)
                if not cm_text.endswith("\n"):
                    out.write("\n")
                downloaded += 1
            else:
                skipped += 1

            if i % 100 == 0:
                print(f"Downloaded {i}/{total} families...")

            if delay > 0:
                time.sleep(delay)

    print(f"Done. Downloaded {downloaded} CMs, skipped {skipped}.")
    print(f"Output: {output_path}")

    # Press the database
    print("Running cmpress...")
    result = subprocess.run(["cmpress", output_path], capture_output=True, text=True)
    if result.returncode != 0:
        print("Error: cmpress failed.", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        sys.exit(1)

    print("cmpress completed successfully.")
    print(f"CM database ready: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build a live Rfam CM database from the SVN repository."
    )
    parser.add_argument(
        "--output", "-o",
        default="Rfam_live.cm",
        help="Output path for the CM database (default: Rfam_live.cm)"
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.1,
        help="Delay between requests in seconds (default: 0.1)"
    )
    args = parser.parse_args()

    build_live_cm(args.output, delay=args.delay)


if __name__ == "__main__":
    main()
