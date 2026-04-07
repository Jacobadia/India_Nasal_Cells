significant_file = "../artifacts/pure_control_sex/deg_results_significant.txt"
lookup_file = "../artifacts/genetype_lookup.txt"

def main():
    with open(significant_file, "rt") as f:
        for line in f:
            data = iter(f.readlines())
        iter(data)
    gene_ids = []
    for line in data:
        gene_id = line.split("\t")[0]
        gene_ids.append(gene_id)
    # grep for the gene id in the lookup file and print the gene name
    with open(lookup_file, "rt") as f:
        for line in f:
            line = line.strip()
            for gene_id in gene_ids:
                if gene_id in line:
                    print(line)

if __name__ == "__main__":
    main()