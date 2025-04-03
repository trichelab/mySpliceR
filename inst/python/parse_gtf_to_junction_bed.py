# Updating the Python function to optionally return a DataFrame if no output path is provided.

def parse_gtf_to_junction_bed(gtf_path, output_bed_path=None):
    bed_records = []
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

            # GTF is 1-based, BED is 0-based
            start_bed = int(start) - 1
            end_bed = int(end)

            # Parse attributes field
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

            # Clean chr prefix if needed
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

            bed_name = f"{gene_name}:{transcript_id}:{strand}:e{exon_number}"

            bed_records.append([chrom, start_bed, end_bed, bed_name, ".", strand])

    # Create a DataFrame
    bed_df = pd.DataFrame(bed_records, columns=["chrom", "start", "end", "name", "score", "strand"])

    if output_bed_path:
        bed_df.to_csv(output_bed_path, sep="\t", header=False, index=False)
    else:
        return bed_df

