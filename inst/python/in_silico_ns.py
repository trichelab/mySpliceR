import os
import pysam


def parse_bam_list(bam_list_path):
    bam_paths = {}
    with open(bam_list_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            sample = parts[0]
            path = parts[2]
            chrom_hint = parts[3]
            if "chr" in chrom_hint:
                bam_paths.setdefault(sample, {})['chr'] = path
            else:
                bam_paths.setdefault(sample, {})['nochr'] = path
    return bam_paths


def load_known_junctions(bed_file):
    known_junc = set()
    known_splice_sites = set()
    with open(bed_file) as f:
        for line in f:
            chrom, start, end, annot = line.strip().split("\t")
            parts = annot.split(":")
            if parts[3].startswith("i"):
                known_junc.add((chrom, int(start), int(end)))
                known_splice_sites.update([(chrom, int(start)), (chrom, int(end))])
    return known_junc, known_splice_sites


def detect_novel_junctions(maf_file, bam_list, bed_file, output_file, min_reads=5):
    required_cols = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]

    with open(maf_file) as maf:
        reader = csv.DictReader(maf, delimiter='\t')
        for col in required_cols:
            if col not in reader.fieldnames:
                raise ValueError(f"Missing required MAF column: {col}")

        known_junc, known_splice_sites = load_known_junctions(bed_file)
        bam_paths = parse_bam_list(bam_list)

        with open(output_file, 'w') as summary_out, \
             open(output_file + ".detailed.alignment.5", 'w') as detail_out:

            for row in reader:
                try:
                    sample_id = row["Tumor_Sample_Barcode"].rstrip("_T")
                    chrom = row["Chromosome"].lstrip("chr")
                    pos = int(row["Start_Position"])
                    ref = row["Reference_Allele"]
                    var = row["Tumor_Seq_Allele2"]
                except (KeyError, ValueError):
                    continue

                dellen = len(ref) if var == "-" else 0
                chr_pos = f"{chrom}:{pos-20}-{pos+20}"

                bam = None
                if sample_id in bam_paths:
                    if 'chr' in bam_paths[sample_id] and os.path.exists(bam_paths[sample_id]['chr']):
                        bam = pysam.AlignmentFile(bam_paths[sample_id]['chr'], "rb")
                        chrom = "chr" + chrom
                    elif 'nochr' in bam_paths[sample_id] and os.path.exists(bam_paths[sample_id]['nochr']):
                        bam = pysam.AlignmentFile(bam_paths[sample_id]['nochr'], "rb")

                if not bam:
                    continue

                count_read = {}
                for read in bam.fetch(chrom, pos-20, pos+20):
                    cigar = read.cigarstring
                    match = re.match(r"^(\d+)M(\d+)N(\d+)M$", cigar or "")
                    if match:
                        m1, n, m2 = map(int, match.groups())
                        x = read.reference_start + m1 - 1
                        y = read.reference_start + m1 + n - 2
                        if ((chrom, x, y) not in known_junc and
                            (y - x + 1 != dellen) and
                            read.reference_name == chrom and
                            ((pos - 20 <= x <= pos + 20 and (chrom, x) not in known_splice_sites) or
                             (pos - 20 <= y <= pos + 20 and (chrom, y) not in known_splice_sites))):
                            count_read[read.query_name] = read.to_string()

                if len(count_read) >= min_reads:
                    summary_out.write("\t".join([row.get(col, "") for col in reader.fieldnames]) + f"\t{len(count_read)}\n")
                    detail_out.write(f"{sample_id}\t{chrom}\t{pos}\t{ref}\t{var}\n")
                    for r in sorted(count_read):
                        detail_out.write(f"{count_read[r]}\n")



