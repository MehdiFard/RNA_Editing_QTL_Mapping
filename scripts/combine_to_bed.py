#!/usr/bin/env python3

import os
import glob
import subprocess
import argparse


def combine_qqnorm_files(input_dir, tissue, min_cov, min_sam):
    prefix = f"{tissue}.edMat.{min_cov}cov.{min_sam}samps.noXYM.txt"
    out_bed = os.path.join(input_dir, f"{prefix}.qqnorm.bed.gz")

    # Find chromosome-specific files
    chr_files = sorted(
        glob.glob(os.path.join(input_dir, f"{prefix}.qqnorm_chr*")),
        key=lambda x: int(x.split("qqnorm_chr")[-1])
    )

    if not chr_files:
        print(f"No files matched in: {input_dir}")
        exit(1)

    print(f"Found {len(chr_files)} .qqnorm_chr* files.")
    print("Making the BED file.")

    # Compose shell command
    cmd = f"""
    {{
        cat {chr_files[0]} | head -n 1
        cat {' '.join(chr_files)} | grep -v '^#' | sed 's/nan/NA/g' | sort -k1,1 -k2n,2 | \\
        awk '{{print "chr"$1,$2,$3,"chr"$1"_"$3,$0}}' OFS="\\t" | cut -f 1,2,3,4,9-
    }} | sed 's/^chr#/#/g' | bgzip -c > {out_bed}
    """

    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
    subprocess.run(f"tabix -f -p bed {out_bed}", shell=True, check=True)

    print(f"BED file written to: {out_bed}")
    print(f"Indexed with tabix: {out_bed}.tbi")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge QTLtools .qqnorm_chr* files into a sorted, bgzipped BED format.")
    parser.add_argument("--input_dir", required=True, help="Directory containing the .qqnorm_chr* files")
    parser.add_argument("--tissue", required=True, help="Tissue name")
    parser.add_argument("--min_coverage", required=True, type=int, help="Minimum coverage")
    parser.add_argument("--min_samples", required=True, type=int, help="Minimum number of samples")

    args = parser.parse_args()
    combine_qqnorm_files(args.input_dir, args.tissue, args.min_coverage, args.min_samples)

### END