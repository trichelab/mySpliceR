# Re-import necessary modules and redefine the corrected function after the kernel reset

import pandas as pd
from collections import defaultdict

def parse_gtf_to_junction_bed(gtf_path, output_bed_path=None):
    exon_records = defaultdict(list)

    with open(gtf_path, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            if feature != "exon":
                continue

            start_bed = int(start) - 1  # BED start
            end_bed = int(end)          # BED end (exclusive)

            # Parse attributes
            attr_dict = {}
            for attr in attributes.strip().split(";"):
                if attr.strip() == "":
                    continue
                key_value = attr.strip().split(" ", 1)
                if len(key_value) == 2:
                    key, value = key_value
                    attr_dict[key] = value.strip('"')

            gene_name = attr_dict.get("gene_name", "NA")
            transcript_id = attr_dict.get("transcript_id", "NA")
            exon_number = attr_dict.get("exon_number", "NA")

            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

            key = (gene_name, transcript_id, strand, chrom)
            exon_records[key].append((start_bed, end_bed, exon_number))

    bed_output = []

    for (gene_name, transcript_id, strand, chrom), exons in exon_records.items():
        # Sort exons by genomic position (left to right)
        exons_sorted = sorted(exons, key=lambda x: x[0])

        # Determine exon ordering for label based on strand
        exon_order = list(range(1, len(exons_sorted) + 1))
        if strand == "-":
            exon_order = list(reversed(exon_order))

        # Write exons with correct label order
        for i, (start, end, _) in enumerate(exons_sorted):
            bed_output.append([
                chrom,
                start,
                end,
                f"{gene_name}:{transcript_id}:{strand}:e{exon_order[i]}",
                ".",
                strand
            ])

        # Write introns between consecutive exons
        for i in range(len(exons_sorted) - 1):
            intron_start = exons_sorted[i][1]       # end of exon i
            intron_end = exons_sorted[i + 1][0]     # start of exon i+1
            if intron_start < intron_end:
                # Intron label should follow transcription order
                intron_index = exon_order[i] if strand == "+" else exon_order[i + 1]
                bed_output.append([
                    chrom,
                    intron_start,
                    intron_end,
                    f"{gene_name}:{transcript_id}:{strand}:i{intron_index}",
                    ".",
                    strand
                ])

    bed_df = pd.DataFrame(bed_output, columns=["chrom", "start", "end", "name", "score", "strand"])

    if output_bed_path:
        bed_df.to_csv(output_bed_path, sep="\t", header=False, index=False)
    else:
        return bed_df

