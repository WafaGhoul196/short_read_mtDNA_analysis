#!/usr/bin/env python3
"""
apply_annotation.py
Merge per-sample filtered variant CSVs and annotate with multiple databases.
"""

import argparse
import os
import glob
import pandas as pd


def load_sample_csv(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    df.columns = df.columns.str.strip('"')
    if "ALT" in df.columns:
        df["ALT"] = df["ALT"].str.upper()
    if "RefCall" in df.columns:
        df["RefCall"] = df["RefCall"].str.upper()
    sample_name = os.path.basename(file_path).replace("_no_haplo_variants.csv", "")
    df["sample"] = sample_name
    return df


def main():
    parser = argparse.ArgumentParser(description="Annotate filtered mtDNA variants.")
    parser.add_argument("--input_dir",  default=".",
                        help="Directory containing *_no_haplo_variants.csv files")
    parser.add_argument("--annot_dir",  default="./annotation_databases",
                        help="Directory containing annotation database CSV files")
    parser.add_argument("--output",     default="merged_annotated.csv",
                        help="Path for the output annotated CSV")
    args = parser.parse_args()

    # ── Load sample files ────────────────────────────────────────────────────
    sample_files = glob.glob(os.path.join(args.input_dir, "*_no_haplo_variants.csv"))
    if not sample_files:
        print(f"No '*_no_haplo_variants.csv' files found in {args.input_dir}")
        raise SystemExit(1)

    print(f"Loading {len(sample_files)} sample file(s)...")
    merged_df = pd.concat([load_sample_csv(f) for f in sample_files], ignore_index=True)
    print(f"  Total variants loaded: {len(merged_df)}")

    # ── Annotation databases ─────────────────────────────────────────────────
    db_files = {
        "gnomAD":     os.path.join(args.annot_dir, "gnomAD_filtred.csv"),
        "HelixMTdb":  os.path.join(args.annot_dir, "HelixMTdb_modified.csv"),
        "MitImpact":  os.path.join(args.annot_dir, "MitImpact-db-3.1.3_modified.csv"),
        "MitoMAP":    os.path.join(args.annot_dir, "MitoMAP_RNA_tRNA.csv"),
        "t_APOGEE":   os.path.join(args.annot_dir, "t_APOGEE_data_modified.csv"),
    }

    db_data = {}
    for name, path in db_files.items():
        if not os.path.isfile(path):
            print(f"  ⚠ Database not found, skipping: {path}")
            continue
        df = pd.read_csv(path)
        df.columns = df.columns.str.strip('"')
        if "Ref" in df.columns:
            df["Ref"] = df["Ref"].str.upper()
        if "Alt" in df.columns:
            df["Alt"] = df["Alt"].str.upper()
        db_data[name] = df
        print(f"  Loaded {name}: {len(df)} rows")

    # ── Merge annotations ────────────────────────────────────────────────────
    final_df = merged_df.copy()
    for name, df in db_data.items():
        print(f"  Annotating with {name}...")
        if "Ref" in df.columns and "Alt" in df.columns:
            final_df = final_df.merge(
                df, how="left",
                left_on=["Pos", "RefCall", "ALT"],
                right_on=["Pos", "Ref", "Alt"]
            )
            final_df.drop(columns=["Ref", "Alt"], inplace=True, errors="ignore")
        else:
            final_df = final_df.merge(df, how="left", on="Pos")

    # ── Save ─────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    final_df.to_csv(args.output, index=False)
    print(f"\n✔ Annotation complete → {args.output}")
    print(f"  Rows: {len(final_df)}, Columns: {len(final_df.columns)}")


if __name__ == "__main__":
    main()
