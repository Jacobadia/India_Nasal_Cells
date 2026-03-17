from pathways.name_mapper import GeneNameMapper

input_file = "../artifacts/lpm_pure_no_control/deg_results_full.txt"

# Read the first 100 results and print the adjusted
n = 100

mapper = GeneNameMapper()

with open(input_file, "r") as f:
    header = f.readline().strip().split("\t")
    # print("\t".join(header))
    for _ in range(n):
        line = f.readline()
        if not line:
            break
        gene_id = line.split("\t")[0]
        print(gene_id)
