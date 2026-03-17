from pathlib import Path
import sys
from typing import List
try:
    from pathways.name_mapper import GeneNameMapper
except ModuleNotFoundError:
    # Allow running this script directly from deg/.
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    from pathways.name_mapper import GeneNameMapper

input_file = "../artifacts/lpm_pure_no_control/deg_results_full.txt"

# Read the first 100 results and print the adjusted
n = 100

mapper = GeneNameMapper()

with open(input_file, "r") as f:
    header = f.readline().strip().split("\t")
    # print("\t".join(header))
    names: List[str] = [] 
    for _ in range(n):
        line = f.readline()
        if not line:
            break
        gene_id = line.split("\t")[0]
        names.append(mapper.map_gene_id(gene_id))

print("\n".join(names))
