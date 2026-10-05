#!/usr/bin/env python3
# Select one tissue's subjects from the sample annotation (columns are read by name).
# Without --covariates: write the subject IDs, one per line (sample list for plink2 --keep).
# With --covariates:    write subject_id, sex, age for those subjects, no header (input for covariates.R).

import argparse
import sys
import pandas as pd


def read_table(path, columns):
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in columns if c not in df.columns]
    if missing:
        sys.exit(f"Error: {path} is missing column(s): {', '.join(missing)}")
    return df


parser = argparse.ArgumentParser(description="Select one tissue's subjects (and their sex/age) from the annotation files")
parser.add_argument("--annotation", required=True, help="Sample annotation TSV (columns: subject_id, tissue_cln)")
parser.add_argument("--tissue", required=True, help="Tissue name, as in the tissue_cln column")
parser.add_argument("--covariates", help="Subject covariates TSV (columns: subject_id, sex, age)")
parser.add_argument("--output", required=True, help="Output file")
args = parser.parse_args()

anno = read_table(args.annotation, ["subject_id", "tissue_cln"])
subjects = anno.loc[anno["tissue_cln"] == args.tissue, "subject_id"]
if subjects.empty:
    sys.exit(f"Error: No samples found for tissue '{args.tissue}' in {args.annotation}.")

if args.covariates is None:
    subjects.to_csv(args.output, index=False, header=False)
else:
    cov = read_table(args.covariates, ["subject_id", "sex", "age"]).drop_duplicates("subject_id")
    out = subjects.to_frame().merge(cov[["subject_id", "sex", "age"]], on="subject_id", how="left")
    no_cov = out.loc[out["sex"].isna() | out["age"].isna(), "subject_id"].tolist()
    if no_cov:
        sys.exit(f"Error: {len(no_cov)} subject(s) of {args.tissue} have no sex/age in {args.covariates}: {', '.join(no_cov[:10])}")
    out.to_csv(args.output, sep="\t", index=False, header=False)
