#!usr/bin/env python3

gtf_file = "/grphome/grp_tb/star_genome/gencode.v49.primary_assembly.annotation.gtf"
output_file = "/grphome/grp_tb/processing_scripts/results/documents/genetype_lookup.txt"

def is_qualified_entry(line, fields):
    if line.startswith("#"):
        return False
    if len(fields) < 9:
        print("Warning: Line does not have 9 fields, skipping it")
        return False
    feature_type = fields[2].lower()
    if feature_type != "gene":
        return False
    return True

def process_line(line):
    fields = line.strip().split("\t")
    if not is_qualified_entry(line, fields):
        return
    attributes = fields[8]
    attributes.split(";")
    attribute_dict = {}
    for attribute in attributes.split(";"):
        key, value = attribute.strip().split(" ", 1)
        attribute_dict[key] = value.strip('"')
    gene_id = attribute_dict.get("gene_id", "NA")
    gene_type = attribute_dict.get("gene_type", "NA")
    if gene_id == "NA":
        print("Warning: gene_id not found in attributes, skipping line")
        return
    if gene_type == "NA":
        print(f"Warning: gene_type not found for gene_id {gene_id}, setting it to NA")
    write_output(gene_id, gene_type)

def write_output(gene_id, gene_type):
    with open(output_file, "a") as f:
        f.write(f"{gene_id}\t{gene_type}\n")

def main():
    with open(gtf_file) as f:
        for line in f:
            process_line(line)

if __name__ == "__main__":
    main()