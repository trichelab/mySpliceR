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
    known_junc, known_splice_sites = load_known_junctions(bed_file)
    bam_paths = parse_bam_list(bam_list)

    with open(maf_file) as maf, \
         open(output_file, 'w') as summary_out, \
         open(output_file + ".detailed.alignment.5", 'w') as detail_out:

        for line in maf:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split("\t")
            try:
                sample_id = fields[15].rstrip("_T")
                chrom = fields[4].lstrip("chr")
                pos = int(fields[5])
                ref = fields[10]
                var = fields[12]
            except IndexError:
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
                if cigar and re.match(r"^(\d+)M(\d+)N(\d+)M$", cigar):
                    m1, n, m2 = map(int, re.findall(r"(\d+)M|N", cigar))
                    x = read.reference_start + m1 - 1
                    y = read.reference_start + m1 + n - 2
                    if ((chrom, x, y) not in known_junc and
                        (y - x + 1 != dellen) and
                        read.reference_name == chrom and
                        ((pos - 20 <= x <= pos + 20 and (chrom, x) not in known_splice_sites) or
                         (pos - 20 <= y <= pos + 20 and (chrom, y) not in known_splice_sites))):
                        count_read[read.query_name] = read.to_string()

            if len(count_read) >= min_reads:
                summary_out.write(f"{line.strip()}\t{len(count_read)}\n")
                detail_out.write(f"{sample_id}\t{chrom}\t{pos}\t{ref}\t{var}\n")
                for r in sorted(count_read):
                    detail_out.write(f"{count_read[r]}\n")



