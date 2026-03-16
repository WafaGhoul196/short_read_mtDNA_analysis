#!/usr/bin/env python3
"""
filter_haplo_variants.py
Remove haplogroup-defining variants from filtered pileup CSVs.

For each sample, variants listed in Haplogrep3's Found_Polys column are
separated into two files:
  - <sample>_no_haplo_variants.csv   → putative pathogenic / private variants
  - <sample>_haplo_variants.csv      → haplogroup-defining variants (removed)

Usage:
    python3 filter_haplo_variants.py \
        --haplogrep  <merged_haplogrep_results.txt> \
        --input_dir  <directory with *_data_filtred_heteroplasmy.csv> \
        --output_dir <output directory>
"""

import argparse
import os
import sys
import pandas as pd


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_haplogrep_data(haplogrep_file: str) -> dict:
    """
    Parse a merged Haplogrep3 results file.

    Returns a dict keyed by SampleID:
        { sample_id: { 'haplogroup': str, 'defining_variants': [str, ...] } }

    The SampleID in Haplogrep output matches the filename of the filtered CSV,
    e.g. the row with SampleID='SEQ001_data_filtred_heteroplasmy.csv' corresponds
    to the file SEQ001_data_filtred_heteroplasmy.csv.
    """
    print(f"Loading Haplogrep results: {haplogrep_file}")
    haplogrep_df = pd.read_csv(haplogrep_file, sep="\t")

    required = {"SampleID", "Haplogroup", "Found_Polys"}
    missing = required - set(haplogrep_df.columns)
    if missing:
        sys.exit(f"ERROR: Haplogrep file is missing columns: {missing}")

    haplogrep_data = {}
    for _, row in haplogrep_df.iterrows():
        sample_id = str(row["SampleID"]).strip()
        defining = []
        if pd.notna(row["Found_Polys"]):
            for variant in str(row["Found_Polys"]).split():
                clean = variant.replace("!", "").replace("?", "").replace("@", "").strip()
                if clean:
                    defining.append(clean)
        haplogrep_data[sample_id] = {
            "haplogroup": str(row["Haplogroup"]).strip(),
            "defining_variants": defining,
        }

    print(f"  Loaded {len(haplogrep_data)} samples from Haplogrep.")
    return haplogrep_data


def parse_haplo_variants(defining_variants: list) -> dict:
    """
    Convert a list of variant strings (e.g. ['263A', '750G']) into
    a dict of {position (int): alt_allele (str)}.
    """
    variants_to_remove = {}
    for variant in defining_variants:
        try:
            pos_str = "".join(filter(str.isdigit, variant))
            alt_str = "".join(filter(str.isalpha, variant))
            if pos_str and alt_str:
                variants_to_remove[int(pos_str)] = alt_str.upper()
        except Exception:
            print(f"  Warning: invalid variant format '{variant}', skipping.")
    return variants_to_remove


def is_haplo_variant(row, variants_to_remove: dict) -> bool:
    pos = row["Pos"]
    alt = str(row["ALT"]).upper()
    return pos in variants_to_remove and alt == variants_to_remove[pos]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Remove haplogroup-defining variants from filtered pileup CSVs."
    )
    parser.add_argument(
        "--haplogrep",
        required=True,
        help="Path to merged_haplogrep_results.txt",
    )
    parser.add_argument(
        "--input_dir",
        default=".",
        help="Directory containing *_data_filtred_heteroplasmy.csv files",
    )
    parser.add_argument(
        "--output_dir",
        default="haplogroup_filtered",
        help="Output directory for filtered files",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load Haplogrep data
    haplogrep_data = load_haplogrep_data(args.haplogrep)

    # Discover sample CSV files
    sample_files = sorted(
        f for f in os.listdir(args.input_dir)
        if f.endswith("_data_filtred_heteroplasmy.csv")
    )
    if not sample_files:
        sys.exit(f"ERROR: No '*_data_filtred_heteroplasmy.csv' files found in {args.input_dir}")

    print(f"\nFound {len(sample_files)} sample file(s).")

    # Match CSV files to Haplogrep entries
    # Haplogrep SampleID = filename without _data_filtred_heteroplasmy.csv
    matched, unmatched = [], []
    for csv_file in sample_files:
        base = csv_file.replace("_data_filtred_heteroplasmy.csv", "")
        if base in haplogrep_data:
            matched.append(csv_file)
        else:
            unmatched.append(csv_file)

    if unmatched:
        print(f"\nWarning: {len(unmatched)} sample(s) not found in Haplogrep output:")
        for f in unmatched:
            print(f"  {f}")

    if not matched:
        print("\nERROR: No samples matched between CSV files and Haplogrep output.")
        print("  Check that SampleID in the Haplogrep file matches the CSV filenames exactly.")
        if sample_files and haplogrep_data:
            print(f"\n  Example CSV filename : '{sample_files[0]}'")
            print(f"  Example Haplogrep ID : '{list(haplogrep_data.keys())[0]}'")
        sys.exit(1)

    print(f"\nProcessing {len(matched)} matched sample(s)...\n")

    total_kept = total_removed = 0

    for i, csv_file in enumerate(matched, 1):
        print(f"[{i}/{len(matched)}] {csv_file}")
        full_path = os.path.join(args.input_dir, csv_file)
        base_name = csv_file.replace("_data_filtred_heteroplasmy.csv", "")

        sample_df = pd.read_csv(full_path)
        sample_info = haplogrep_data[base_name]
        variants_to_remove = parse_haplo_variants(sample_info["defining_variants"])

        mask = sample_df.apply(is_haplo_variant, axis=1, variants_to_remove=variants_to_remove)
        removed_df  = sample_df[mask].copy()
        filtered_df = sample_df[~mask].copy()

        removed_df["Haplogroup"]    = sample_info["haplogroup"]
        removed_df["Variant_Type"]  = "Haplogroup_Defining"

        out_filtered = os.path.join(args.output_dir, f"{base_name}_no_haplo_variants.csv")
        out_removed  = os.path.join(args.output_dir, f"{base_name}_haplo_variants.csv")

        filtered_df.to_csv(out_filtered, index=False)
        removed_df.to_csv(out_removed,  index=False)

        kept    = len(filtered_df)
        removed = len(removed_df)
        total_kept    += kept
        total_removed += removed

        print(f"  Haplogroup : {sample_info['haplogroup']}")
        print(f"  Total      : {len(sample_df)}  |  Kept: {kept}  |  Removed: {removed}")
        print(f"  → {os.path.basename(out_filtered)}")

    print(f"\nDone. Total variants kept: {total_kept}, removed: {total_removed}.")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
