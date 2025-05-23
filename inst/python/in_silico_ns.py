#!/usr/bin/env python
import os
import pysam
import csv
import re

def parse_bam_list(bam_list_path):
    """
    Parses the BAM list file (TSV) and returns a dictionary mapping sample IDs to
    a dictionary with keys 'chr' and/or 'nochr' for BAM paths.
    """
    bam_paths = {}
    with open(bam_list_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            sample = parts[0]
            path = os.path.expanduser(parts[2].strip())
            chrom_hint = parts[3]
            if "chr" in chrom_hint:
                bam_paths.setdefault(sample, {})['chr'] = path
            else:
                bam_paths.setdefault(sample, {})['nochr'] = path
    return bam_paths

def load_known_junctions(bed_file):
    """
    Load known junctions and splice sites from the BED file.
    Assumes the BED file has at least 4 columns, where the 4th column
    contains an annotation like ...:...:...:i1 (or i8, etc.).
    """
    known_junc = set()
    known_splice_sites = set()
    with open(bed_file) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue  # skip headers or empty lines

            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue  # skip malformed lines

            chrom, start, end, annot = parts[:4]  # safely unpack first 4 columns
            subparts = annot.split(":")
            if len(subparts) >= 4 and subparts[3].startswith("i"):
                start, end = int(start), int(end)
                known_junc.add((chrom, start, end))
                known_splice_sites.update([(chrom, start), (chrom, end)])
    return known_junc, known_splice_sites

def detect_novel_junctions(maf_file, bam_list, bed_file, output_file, min_reads=5):
    # Required columns in the input MAF file
    required_cols = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]

    with open(maf_file) as maf:
        reader = csv.DictReader(maf, delimiter='\t')
        for col in required_cols:
            if col not in reader.fieldnames:
                raise ValueError(f"Missing required MAF column: {col}")

        # Load known junctions and splice sites from BED file
        known_junc, known_splice_sites = load_known_junctions(bed_file)
        # Load sample-to-BAM path mapping
        bam_paths = parse_bam_list(bam_list)

        results = {}  # Dictionary to collect structured result per variant

        with open(output_file, 'w') as summary_out, \
             open(output_file + ".detailed.alignment.5", 'w') as detail_out:

            for row in reader:
                try:
                    # Extract variant information from MAF
                    sample_id = row["Tumor_Sample_Barcode"].rstrip("_T")
                    chrom = row["Chromosome"].lstrip("chr")
                    pos = int(row["Start_Position"])
                    ref = row["Reference_Allele"]
                    var = row["Tumor_Seq_Allele2"]
                except (KeyError, ValueError):
                    continue

                dellen = len(ref) if var == "-" else 0
                chr_pos = f"{chrom}:{pos-20}-{pos+20}"

                # Resolve BAM file path based on whether it uses chr prefix
                bam = None
                if sample_id in bam_paths:
                    if 'chr' in bam_paths[sample_id] and os.path.exists(bam_paths[sample_id]['chr']):
                        bam = pysam.AlignmentFile(bam_paths[sample_id]['chr'], "rb")
                        chrom = "chr" + chrom
                    elif 'nochr' in bam_paths[sample_id] and os.path.exists(bam_paths[sample_id]['nochr']):
                        bam = pysam.AlignmentFile(bam_paths[sample_id]['nochr'], "rb")

                if not bam:
                    continue  # Skip if BAM could not be opened

                count_read = {}     # Map read name to alignment string
                read_info_list = [] # List to store detailed read info

                # Fetch reads overlapping region around variant
                for read in bam.fetch(chrom, pos-20, pos+20):
                    cigar = read.cigarstring
                    match = re.match(r"^(\d+)M(\d+)N(\d+)M$", cigar or "")
                    if match:
                        m1, n, m2 = map(int, match.groups())
                        x = read.reference_start + m1      # intron start (0-based)
                        y = read.reference_start + m1 + n  # intron end (0-based, exclusive)
                        if ((chrom, x, y) not in known_junc and
                            (y - x + 1 != dellen) and
                            read.reference_name == chrom and
                            ((pos - 20 <= x <= pos + 20 and (chrom, x) not in known_splice_sites) or
                             (pos - 20 <= y <= pos + 20 and (chrom, y) not in known_splice_sites))):
                            read_id = read.query_name
                            count_read[read_id] = read.to_string()
                            read_info_list.append({
                                "read_name": read_id,
                                "sequence": read.query_sequence,
                                "junction_start": x,
                                "junction_end": y
                            })

                # Only process variants with at least min_reads support
                n_support = len(count_read)
                if n_support >= min_reads:
                    # Write to summary output file
                    summary_out.write("\t".join([row.get(col, "") for col in reader.fieldnames]) + f"\t{n_support}\n")
                    # Write to detailed alignment output file
                    detail_out.write(f"{sample_id}\t{chrom}\t{pos}\t{ref}\t{var}\n")
                    for r in sorted(count_read):
                        detail_out.write(f"{count_read[r]}\n")

                    # Save structured result for R/Python integration
                    key = f"{sample_id}_{chrom}_{pos}_{ref}_{var}"
                    results[key] = {
                        "sample": sample_id,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": var,
                        "num_supporting_reads": n_support,
                        "supporting_reads": read_info_list
                    }

        return results

