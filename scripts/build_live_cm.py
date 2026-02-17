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


def extract_first_model(cm_text, family):
    """Extract the first CM model and rename its NAME to the family ID.

    SVN CM files may contain multiple models (e.g., calibrated + uncalibrated),
    each terminated by '//'. We keep only the first model (the calibrated one)
    and rename its NAME from 'SEED' to the family ID (e.g., RF00001) to ensure
    unique primary keys when pressed into a single database.
    """
    # Find the first complete model: everything up to and including the first '//'
    end = cm_text.find("\n//\n")
    if end == -1:
        # Try without trailing newline (end of file)
        end = cm_text.find("\n//")
        if end == -1:
            return None
        model = cm_text[:end + 3]  # include \n//
    else:
        model = cm_text[:end + 4]  # include \n//\n

    # Rename NAME field to the family ID
    model = re.sub(
        r'^(NAME\s+)\S+',
        rf'\g<1>{family}',
        model,
        count=1,
        flags=re.MULTILINE,
    )
    return model


def build_live_cm(output_path, delay=0.1, press=True):
    """Build a live Rfam CM database from SVN."""
    if press:
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
                model = extract_first_model(cm_text, family)
                if model:
                    out.write(model)
                    if not model.endswith("\n"):
                        out.write("\n")
                    downloaded += 1
                else:
                    skipped += 1
            else:
                skipped += 1

            if i % 100 == 0:
                print(f"Downloaded {i}/{total} families...")

            if delay > 0:
                time.sleep(delay)

    print(f"Done. Downloaded {downloaded} CMs, skipped {skipped}.")
    print(f"Output: {output_path}")

    if press:
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
    parser.add_argument(
        "--no-press",
        action="store_true",
        help="Skip running cmpress (useful when only the raw CM file is needed)"
    )
    args = parser.parse_args()

    build_live_cm(args.output, delay=args.delay, press=not args.no_press)


if __name__ == "__main__":
    main()
